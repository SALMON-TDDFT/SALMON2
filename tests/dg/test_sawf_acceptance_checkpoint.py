#!/usr/bin/env python3
from pathlib import Path
import os, shutil, subprocess, tempfile

ROOT=Path(__file__).resolve().parents[2]
FC=shutil.which('gfortran')
if not FC: raise SystemExit('gfortran is required')
driver=r'''program check_acceptance
 use lcfo_wannier_sawf_templates, only: t_sawf_acceptance_checkpoint, &
   write_sawf_acceptance_checkpoint,read_sawf_acceptance_checkpoint,validate_sawf_acceptance_checkpoint
 use lcfo_wannier_sawf_templates, only: admit_sawf_hierarchical_basis
 implicit none
 type(t_sawf_acceptance_checkpoint)::a,b
 logical::ok,reusable
 character(256)::msg
 allocate(a%buffer_size(3),a%center_residual(2),a%projector_residual(2), &
   a%overlap_residual(2),a%ww_residual(2),a%wp_residual(2),a%face_residual(3,2), &
   a%global_local_face_residual(3))
 a%supercell_fingerprint='cell-A';a%buffer_size=[2,3,4]
 a%center_residual=1d-9;a%projector_residual=2d-9;a%overlap_residual=3d-9
 a%ww_residual=4d-9;a%wp_residual=5d-9;a%face_residual=6d-9
 a%global_local_projector_residual=2d-9;a%global_local_overlap_residual=3d-9
 a%global_local_fixed_h_residual=4d-9;a%global_local_face_residual=5d-9
 call validate_sawf_acceptance_checkpoint(a,1d-8,ok,msg);call req(ok,'valid acceptance')
 call write_sawf_acceptance_checkpoint('accept.chk',a,ok,msg);call req(ok,'write')
 call read_sawf_acceptance_checkpoint('accept.chk','cell-A',b,reusable,ok,msg)
 call req(ok.and.reusable,'matching read')
 call validate_sawf_acceptance_checkpoint(b,1d-8,ok,msg);call req(ok,'roundtrip validate')
 call read_sawf_acceptance_checkpoint('accept.chk','cell-B',b,reusable,ok,msg)
 call req(ok.and..not.reusable,'cross-supercell rejection')
 call admit_sawf_hierarchical_basis('accept.chk','cell-A',1d-8,b,ok,msg)
 call req(ok,'hierarchical admission')
 call admit_sawf_hierarchical_basis('accept.chk','cell-B',1d-8,b,ok,msg)
 call req(.not.ok,'hierarchical provenance rejection')
 a%face_residual(2,2)=2d-8
 call validate_sawf_acceptance_checkpoint(a,1d-8,ok,msg);call req(.not.ok,'face convergence rejection')
 a%face_residual(2,2)=-1d-9
 call validate_sawf_acceptance_checkpoint(a,1d-8,ok,msg);call req(.not.ok,'negative residual rejection')
 a%face_residual(2,2)=0d0
 deallocate(a%center_residual);allocate(a%center_residual(1));a%center_residual=0
 call validate_sawf_acceptance_checkpoint(a,1d-8,ok,msg);call req(.not.ok,'incomplete series rejection')
 write(*,'(a)')'PASS SAWF operator-complete acceptance checkpoint'
contains
 subroutine req(x,t);logical,intent(in)::x;character(*),intent(in)::t
  if(.not.x)then;write(*,'(a)')trim(t);error stop 1;endif
 end subroutine
end program'''
with tempfile.TemporaryDirectory(prefix='sawf-acceptance-') as td:
 td=Path(td);(td/'driver.f90').write_text(driver)
 (td/'CMakeLists.txt').write_text(f'''cmake_minimum_required(VERSION 3.16)
project(acceptance LANGUAGES Fortran)
find_package(LAPACK REQUIRED)
add_executable(check "{ROOT/'src/gs/dc/lcfo_wannier_sawf_templates.f90'}" driver.f90)
target_link_libraries(check PRIVATE ${{LAPACK_LIBRARIES}})
target_compile_options(check PRIVATE -fcheck=all -fbacktrace)
''')
 env=dict(os.environ)
 for cmd in (["cmake","-S",str(td),"-B",str(td/'b'),f"-DCMAKE_Fortran_COMPILER={FC}"],
             ["cmake","--build",str(td/'b'),"-j","2"]):
  r=subprocess.run(cmd,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
  if r.returncode:raise SystemExit(r.stdout)
 r=subprocess.run([str(td/'b/check')],cwd=td,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
 if r.returncode:raise SystemExit(r.stdout)
 print(r.stdout.strip())
