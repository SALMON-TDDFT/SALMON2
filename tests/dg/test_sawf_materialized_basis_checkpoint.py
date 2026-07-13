#!/usr/bin/env python3
from pathlib import Path
import os
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
FC = shutil.which("gfortran")
if not FC:
    raise SystemExit("gfortran is required")

driver = r'''program check_materialized_checkpoint
  use lcfo_wannier_sawf_templates, only: t_sawf_ragged_local_basis, &
    write_sawf_materialized_basis_checkpoint, read_sawf_materialized_basis_checkpoint, &
    stitch_sawf_materialized_neighbor_pair,build_sawf_shared_buffer_point_maps
  use lcfo_wannier_sawf_templates, only: build_sawf_fragment_gauge_tree, &
    validate_sawf_materialized_neighbor_closure
  implicit none
  type(t_sawf_ragged_local_basis) :: source, loaded
  logical :: ok, reusable
  character(256) :: message
  integer :: i, shared_left(2),shared_right(2)
  integer,allocatable :: map_left(:),map_right(:)
  integer :: origins(3,4),shapes(3,4),parents(4)
  type(t_sawf_ragged_local_basis) :: left, right
  complex(8) :: q(2,2), right_before(3,2)
  real(8) :: c,s,closure_residual,alignment_residual
  allocate(source%core(3,2), source%buffer(5,2), source%gauge_unitary(2,2))
  source%representative_fragment=2; source%operation_index=3
  source%generated_independently=.true.; source%gauge_residual=2d-12
  source%core=reshape([(cmplx(dble(i),-dble(i),8),i=1,6)],[3,2])
  source%buffer=reshape([(cmplx(.1d0*dble(i),.2d0*dble(i),8),i=1,10)],[5,2])
  source%gauge_unitary=0;source%gauge_unitary(1,1)=1;source%gauge_unitary(2,2)=1
  call write_sawf_materialized_basis_checkpoint('basis.chk','same-supercell-A',4,source,ok,message)
  call req(ok,'write')
  call read_sawf_materialized_basis_checkpoint('basis.chk','same-supercell-A',4,loaded,reusable,ok,message)
  call req(ok.and.reusable,'matching fingerprint')
  call req(maxval(abs(loaded%core-source%core))<1d-14,'core roundtrip')
  call req(maxval(abs(loaded%buffer-source%buffer))<1d-14,'buffer roundtrip')
  call req(loaded%representative_fragment==2.and.loaded%operation_index==3,'provenance roundtrip')
  call read_sawf_materialized_basis_checkpoint('basis.chk','other-supercell',4,loaded,reusable,ok,message)
  call req(ok.and..not.reusable,'cross-supercell rejection')
  call read_sawf_materialized_basis_checkpoint('basis.chk','same-supercell-A',5,loaded,reusable,ok,message)
  call req(ok.and..not.reusable,'fragment mismatch rejection')
  allocate(left%core(3,2),left%buffer(3,2),left%gauge_unitary(2,2), &
    right%core(3,2),right%buffer(3,2),right%gauge_unitary(2,2))
  left%core=0;left%core(1,1)=1;left%core(2,2)=1;left%buffer=left%core
  c=cos(.23d0);s=sin(.23d0)
  q=reshape([cmplx(c,0d0,8),cmplx(s,0d0,8),cmplx(-s,0d0,8),cmplx(c,0d0,8)],[2,2])
  right%core=matmul(left%core,conjg(transpose(q)));right%buffer=right%core
  left%gauge_unitary=0;left%gauge_unitary(1,1)=1;left%gauge_unitary(2,2)=1
  right%gauge_unitary=left%gauge_unitary
  right%gauge_residual=0d0
  shared_left=[1,2];shared_right=[1,2]
  call validate_sawf_materialized_neighbor_closure(left,right,shared_left,shared_right, &
    1d0,1d-12,1d-10,closure_residual,alignment_residual,ok,message)
  call req(.not.ok.and.closure_residual>1d-3.and.alignment_residual<1d-12, &
    'nontrivial loop gauge rejected despite good alignment')
  call stitch_sawf_materialized_neighbor_pair(left,right,shared_left,shared_right, &
    1d0,1d-12,1d-10,ok,message)
  call req(ok.and.right%gauge_residual<1d-12,'materialized pair stitch')
  call req(maxval(abs(right%core-left%core))<1d-12,'right core aligned')
  call validate_sawf_materialized_neighbor_closure(left,right,shared_left,shared_right, &
    1d0,1d-12,1d-10,closure_residual,alignment_residual,ok,message)
  call req(ok.and.closure_residual<1d-12,'aligned loop gauge closure')
  right_before=right%core;right%buffer(:,2)=right%buffer(:,1)
  call stitch_sawf_materialized_neighbor_pair(left,right,shared_left,shared_right, &
    1d0,1d-12,1d-10,ok,message)
  call req(.not.ok.and.maxval(abs(right%core-right_before))<1d-15,'rank failure atomic')
  call build_sawf_shared_buffer_point_maps([8,2,1],[0,0,0],[2,2,1], &
    [2,0,0],[2,2,1],[1,0,0],map_left,map_right,ok,message)
  call req(ok.and.size(map_left)==4.and.size(map_right)==4,'shared periodic buffer map')
  call req(all(map_left==[3,4,7,8]).and.all(map_right==[1,2,5,6]),'shared map ordering')
  call build_sawf_shared_buffer_point_maps([2,2,1],[0,0,0],[2,2,1], &
    [1,0,0],[1,2,1],[1,0,0],map_left,map_right,ok,message)
  call req(.not.ok,'duplicate periodic buffer images rejected')
  origins=reshape([0,0,0, 2,0,0, 0,2,0, 2,2,0],[3,4]);shapes=2
  shapes(3,:)=1
  call build_sawf_fragment_gauge_tree([4,4,1],origins,shapes,parents,ok,message)
  call req(ok.and.all(parents==[0,1,1,3]),'deterministic face-neighbor gauge tree')
  origins(:,4)=[5,5,0]
  call build_sawf_fragment_gauge_tree([8,8,1],origins,shapes,parents,ok,message)
  call req(.not.ok,'disconnected fragment gauge tree rejected')
  write(*,'(a)')'PASS materialized SAWF basis checkpoint'
contains
  subroutine req(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label
    if(.not.condition)then;write(*,'(a)')trim(label);error stop 1;end if
  end subroutine
end program'''

with tempfile.TemporaryDirectory(prefix="sawf-materialized-checkpoint-") as td:
    td = Path(td)
    (td / "driver.f90").write_text(driver)
    (td / "CMakeLists.txt").write_text(f'''cmake_minimum_required(VERSION 3.16)
project(materialized_checkpoint LANGUAGES Fortran)
find_package(LAPACK REQUIRED)
add_executable(check "{ROOT/'src/gs/dc/lcfo_wannier_sawf_templates.f90'}" driver.f90)
target_link_libraries(check PRIVATE ${{LAPACK_LIBRARIES}})
target_compile_options(check PRIVATE -fcheck=all -fbacktrace)
''')
    env = dict(os.environ)
    for cmd in (["cmake", "-S", str(td), "-B", str(td / "b"), f"-DCMAKE_Fortran_COMPILER={FC}"],
                ["cmake", "--build", str(td / "b"), "-j", "2"]):
        result = subprocess.run(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if result.returncode:
            raise SystemExit(result.stdout)
    result = subprocess.run([str(td / "b/check")], cwd=td, env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if result.returncode:
        raise SystemExit(result.stdout)
    print(result.stdout.strip())

flux = (ROOT / "src/gs/dc/lcfo_flux.f90").read_text().lower()
generator = flux.split("subroutine generate_sawf_dmn", 1)[1].split(
    "end subroutine generate_sawf_dmn", 1)[0]
materialize = generator.index("call materialize_sawf_ragged_local_basis")
publish = generator.index("call write_sawf_materialized_basis_checkpoint")
barrier = generator.index("call comm_sync_all", publish)
assert materialize < publish < barrier
assert ".sawf-local-basis" in generator
assert "sawf local basis publication failed" in generator
tree = generator.index("call build_sawf_fragment_gauge_tree")
read_parent = generator.index("call read_sawf_materialized_basis_checkpoint", tree)
shared_map = generator.index("call build_sawf_shared_buffer_point_maps", read_parent)
stitch = generator.index("call stitch_sawf_materialized_neighbor_pair", shared_map)
rewrite = generator.index("call write_sawf_materialized_basis_checkpoint", stitch)
assert tree < read_parent < shared_map < stitch < rewrite
assert "sawf neighbor gauge stitching failed" in generator
closure = generator.index("call validate_sawf_materialized_neighbor_closure", rewrite)
assert rewrite < closure
assert "sawf all-neighbor gauge closure failed" in generator
