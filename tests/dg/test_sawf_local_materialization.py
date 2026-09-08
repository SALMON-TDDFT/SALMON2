#!/usr/bin/env python3
from pathlib import Path
import os, shutil, subprocess, tempfile

ROOT = Path(__file__).resolve().parents[2]
FC = shutil.which("gfortran")
if not FC:
    raise SystemExit("gfortran is required")

driver = r'''program check_materialize
  use lcfo_wannier_sawf_templates, only: materialize_sawf_local_bases, &
    build_sawf_environment_orbits, select_sawf_environment_materialization, &
    t_sawf_ragged_local_basis,materialize_sawf_ragged_local_basis
  implicit none
  complex(8) :: reps(4,2,2),dw(2,2,2),local(4,2,3)
  complex(8)::buffer_rep(8,2)
  integer :: rep_index(3),op_index(3),point_map(4,2)
  logical :: independent(3),ok
  logical :: equivalent(3,3),defect(3),regenerate(3)
  integer :: orbit(3)
  integer :: fragment_maps(3,2),representative_fragment(3),materialize_operation(3)
  logical :: generate_independently(3)
  character(256) :: msg
  type(t_sawf_ragged_local_basis)::ragged
  reps=(0d0,0d0); reps(1,1,1)=1; reps(2,2,1)=1
  buffer_rep=0d0;buffer_rep(1:4,:)=reps(:,:,1);buffer_rep(5:8,:)=reps(:,:,1)
  reps(3,1,2)=1; reps(4,2,2)=1
  dw=(0d0,0d0); dw(1,1,1)=1; dw(2,2,1)=1
  dw(1,2,2)=1; dw(2,1,2)=1
  point_map(:,1)=[1,2,3,4]; point_map(:,2)=[3,4,1,2]
  rep_index=[1,1,2]; op_index=[1,2,1]; independent=[.true.,.false.,.true.]
  call materialize_sawf_local_bases(reps,rep_index,op_index,independent,point_map,dw,local,ok,msg)
  call require(ok,'valid materialization: '//trim(msg))
  call require(abs(local(3,2,2)-1d0)<1d-14.and.abs(local(4,1,2)-1d0)<1d-14, &
    'symmetry spatial map and D_wann action')
  call require(maxval(abs(local(:,:,3)-reps(:,:,2)))<1d-14,'defect-local representative retained')
  independent(3)=.false.; rep_index(3)=1; op_index(3)=0
  call materialize_sawf_local_bases(reps,rep_index,op_index,independent,point_map,dw,local,ok,msg)
  call require(.not.ok.and.index(msg,'operation')>0,'missing defect operation rejected')
  equivalent=.false.;equivalent(1,1)=.true.;equivalent(2,2)=.true.;equivalent(3,3)=.true.
  equivalent(1,2)=.true.;equivalent(2,1)=.true.
  equivalent(2,3)=.true.;equivalent(3,2)=.true.;defect=.false.
  call build_sawf_environment_orbits(equivalent,defect,orbit,regenerate,ok,msg)
  call require(ok.and.all(orbit==orbit(1)),'environment orbit transitive closure')
  call require(count(regenerate)==1,'one representative per transitive orbit')
  orbit=[1,1,2];fragment_maps(:,1)=[1,2,3];fragment_maps(:,2)=[2,1,3]
  call select_sawf_environment_materialization(orbit,fragment_maps,representative_fragment, &
    materialize_operation,generate_independently,ok,msg)
  call require(ok,'materialization provenance')
  call require(all(representative_fragment==[1,1,3]),'orbit representatives')
  call require(all(materialize_operation==[1,2,1]),'actual-group operations')
  call require(all(generate_independently.eqv.[.true.,.false.,.true.]),'independent environments')
  call materialize_sawf_ragged_local_basis(reps(:,:,1),buffer_rep, &
    [3,4,1,2],[5,6,7,8,1,2,3,4],dw(:,:,2),1,2,.false.,ragged,ok,msg)
  call require(ok,'ragged materialization: '//trim(msg))
  call require(all(shape(ragged%core)==[4,2]).and.all(shape(ragged%buffer)==[8,2]), &
    'ragged core/buffer shapes')
  write(*,'(a)') 'PASS representative and defect-local SAWF materialization'
contains
  subroutine require(c,t)
    logical,intent(in)::c; character(*),intent(in)::t
    if(.not.c)then; write(*,'(a)')trim(t); error stop 1; end if
  end subroutine
end program'''

with tempfile.TemporaryDirectory(prefix="sawf-materialize-") as td:
    td = Path(td); (td / "driver.f90").write_text(driver)
    (td / "CMakeLists.txt").write_text(f'''cmake_minimum_required(VERSION 3.16)
project(sawf_materialize LANGUAGES Fortran)
find_package(LAPACK REQUIRED)
add_executable(check "{ROOT/'src/gs/dc/lcfo_wannier_sawf_templates.f90'}" driver.f90)
target_link_libraries(check PRIVATE ${{LAPACK_LIBRARIES}})
''')
    env=dict(os.environ)
    for cmd in (["cmake","-S",str(td),"-B",str(td/"build"),f"-DCMAKE_Fortran_COMPILER={FC}"],
                ["cmake","--build",str(td/"build"),"-j","2"]):
        result=subprocess.run(cmd,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        if result.returncode: raise SystemExit(result.stdout)
    result=subprocess.run([str(td/"build/check")],env=env,stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT,text=True)
    if result.returncode: raise SystemExit(result.stdout)
    print(result.stdout.strip())
