#!/usr/bin/env python3
from pathlib import Path
import os, shutil, subprocess, tempfile

ROOT = Path(__file__).resolve().parents[2]
FC = shutil.which("gfortran")
if not FC:
    raise SystemExit("gfortran is required for SAWF template checkpoint test")

driver = r'''program check_template_checkpoint
  use lcfo_wannier_sawf_templates
  implicit none
  type(t_sawf_template_fingerprint) :: fp, stale
  type(t_sawf_template_checkpoint) :: written, loaded
  logical :: ok, reuse
  integer :: i
  character(512) :: msg
  fp%geometry='geom-a'; fp%pseudopotential='ps-a'; fp%grid='8x8x8'
  fp%band_window='1:4'; fp%complete_projection_shell='s+p'
  fp%symmetry='actual-group-a'; fp%buffer='3'; fp%generator='test-v1'
  written%fingerprint=fp
  allocate(written%centers(3,2),written%spreads(2),written%d_band(2,2,1), &
           written%d_wann(2,2,1),written%gauge_unitary(2,2), &
           written%basis(4,2),written%buffer_basis(8,2),written%orbitals(4,2))
  written%centers=reshape([0d0,1d0,2d0,3d0,4d0,5d0],[3,2])
  written%spreads=[0.25d0,0.5d0]
  written%basis=reshape([(dble(i),i=1,8)],[4,2])
  written%buffer_basis=reshape([(dble(i),i=1,16)],[8,2])
  written%orbitals=cmplx(written%basis,0d0,kind=8)
  written%d_band=(0d0,0d0); written%d_band(1,1,1)=1; written%d_band(2,2,1)=1
  written%d_wann=written%d_band
  written%gauge_unitary=written%d_band(:,:,1)
  written%gauge_residual=3d-13
  call write_sawf_template_checkpoint('template.chk',written,ok,msg)
  call require(ok,'write: '//trim(msg))
  call read_sawf_template_checkpoint('template.chk',fp,loaded,reuse,ok,msg)
  call require(ok.and.reuse,'matching cache must be reused: '//trim(msg))
  call require(maxval(abs(loaded%centers-written%centers))<1d-15,'centers roundtrip')
  call require(maxval(abs(loaded%d_band-written%d_band))<1d-15,'D_band roundtrip')
  call require(maxval(abs(loaded%d_wann-written%d_wann))<1d-15,'D_wann roundtrip')
  call require(abs(loaded%gauge_residual-written%gauge_residual)<1d-20,'gauge residual roundtrip')
  stale=fp; stale%buffer='4'
  call read_sawf_template_checkpoint('template.chk',stale,loaded,reuse,ok,msg)
  call require(ok.and..not.reuse,'fingerprint mismatch must force regeneration')
  call require(index(msg,'regeneration required')>0,'mismatch diagnostic')
  write(*,'(a)') 'PASS SAWF template checkpoint roundtrip and invalidation'
contains
  subroutine require(condition,text)
    logical,intent(in)::condition; character(*),intent(in)::text
    if(.not.condition)then; write(*,'(a)') trim(text); error stop 1; end if
  end subroutine
end program'''

with tempfile.TemporaryDirectory(prefix="sawf-template-checkpoint-") as td:
    td = Path(td)
    (td / "driver.f90").write_text(driver)
    cmake = f'''cmake_minimum_required(VERSION 3.16)
project(sawf_template_checkpoint LANGUAGES Fortran)
find_package(LAPACK REQUIRED)
add_executable(checkpoint "{ROOT / 'src/gs/dc/lcfo_wannier_sawf_templates.f90'}" driver.f90)
target_link_libraries(checkpoint PRIVATE ${{LAPACK_LIBRARIES}})
'''
    (td / "CMakeLists.txt").write_text(cmake)
    env = dict(os.environ)
    subprocess.run(["cmake", "-S", str(td), "-B", str(td / "build"),
                    f"-DCMAKE_Fortran_COMPILER={FC}"], check=True, env=env,
                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    result = subprocess.run(["cmake", "--build", str(td / "build"), "-j", "2"],
                            env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if result.returncode:
        raise SystemExit(result.stdout)
    result = subprocess.run([str(td / "build/checkpoint")], cwd=td, env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if result.returncode:
        raise SystemExit(result.stdout)
    print(result.stdout.strip())
