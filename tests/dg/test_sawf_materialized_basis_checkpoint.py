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
    write_sawf_materialized_basis_checkpoint, read_sawf_materialized_basis_checkpoint
  implicit none
  type(t_sawf_ragged_local_basis) :: source, loaded
  logical :: ok, reusable
  character(256) :: message
  integer :: i
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
