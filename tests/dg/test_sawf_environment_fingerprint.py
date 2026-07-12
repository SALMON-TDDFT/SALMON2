#!/usr/bin/env python3
from pathlib import Path
import os, shutil, subprocess, tempfile

ROOT=Path(__file__).resolve().parents[2]
FLUX=(ROOT/'src/gs/dc/lcfo_flux.f90').read_text().lower()
generator=FLUX.split('subroutine generate_sawf_dmn',1)[1].split('end subroutine generate_sawf_dmn',1)[0]
FC=shutil.which('gfortran')
if not FC: raise SystemExit('gfortran is required')
driver=r'''program check_environment_fingerprint
 use lcfo_wannier_sawf_templates, only: build_sawf_supercell_fingerprint, &
   build_sawf_local_environment_fingerprint, validate_sawf_structure_class
 use lcfo_wannier_sawf_templates, only: build_sawf_file_content_digest
 use lcfo_wannier_sawf_templates, only: measure_sawf_vacuum_occupancy
 implicit none
 real(8)::lattice(3,3),xyz(3,4),relative(3,3),rotated(3,3)
 integer::species(4),local_species(3),grid(3),buffer(3),permutation(3)
 logical::pbc(3),ok
 real(8)::vacuum_by_env(3)
 integer::orbit_by_env(3)
 character(256)::keys(3)
 character(256)::digest_a,digest_b
 integer::unit
 real(8)::measured_vacuum
 character(256)::super_a,super_b,bulk_a,bulk_b,defect,interface_env,surface_env,amorphous,msg
 lattice=0;lattice(1,1)=8;lattice(2,2)=8;lattice(3,3)=8
 xyz=reshape([0d0,0d0,0d0,2d0,0d0,0d0,4d0,0d0,0d0,6d0,0d0,0d0],[3,4])
 species=[14,14,14,14];grid=[16,16,16];buffer=[2,2,2];pbc=.true.
 call build_sawf_supercell_fingerprint(lattice,pbc,species,xyz,grid,buffer, &
   'sha256:pseudo-si','bands:1-8','shells:sp3','LDA','schema:2',super_a,ok,msg)
 call req(ok,'supercell A')
 lattice(1,1)=9
 call build_sawf_supercell_fingerprint(lattice,pbc,species,xyz,grid,buffer, &
   'sha256:pseudo-si','bands:1-8','shells:sp3','LDA','schema:2',super_b,ok,msg)
 call req(ok.and.super_a/=super_b,'different supercell forbidden')
 lattice(1,1)=8;lattice(1,2)=1.25d0
 relative=reshape([-1d0,0d0,0d0,0d0,0d0,0d0,1d0,0d0,0d0],[3,3]);local_species=[14,14,14]
 call local(local_species,relative,0d0,bulk_a);call local(local_species,relative,0d0,bulk_b)
 call req(bulk_a==bulk_b,'equivalent bulk')
 permutation=[3,1,2];rotated=0
 rotated(1,:)=-relative(2,permutation);rotated(2,:)=relative(1,permutation)
 rotated(3,:)=relative(3,permutation)
 call local(local_species(permutation),rotated,0d0,bulk_b)
 call req(bulk_a==bulk_b,'rotation and atom-order invariant bulk fingerprint')
 rotated=relative;rotated(1,1)=rotated(1,1)+1d0;rotated(2,3)=rotated(2,3)+1d-10
 call local(local_species,rotated,0d0,bulk_b)
 call req(bulk_a==bulk_b,'periodic image and sub-tolerance rounding invariant')
 local_species(2)=13;call local(local_species,relative,0d0,defect)
 local_species=[14,14,8];call local(local_species,relative,0d0,interface_env)
 local_species=[14,14,14];call local(local_species,relative,.5d0,surface_env)
 relative(1,3)=1.137d0;call local(local_species,relative,0d0,amorphous)
 call req(all([bulk_a/=defect,bulk_a/=interface_env,bulk_a/=surface_env,bulk_a/=amorphous]), &
   'nonbulk environments remain independent')
 keys=[bulk_a,bulk_a,bulk_a];vacuum_by_env=0;orbit_by_env=1
 call validate_sawf_structure_class('crystal',keys,vacuum_by_env,orbit_by_env,ok,msg)
 call req(ok,'crystal class')
 keys(3)=defect;orbit_by_env(3)=2
 call validate_sawf_structure_class('defect',keys,vacuum_by_env,orbit_by_env,ok,msg)
 call req(ok,'defect class')
 call validate_sawf_structure_class('crystal',keys,vacuum_by_env,orbit_by_env,ok,msg)
 call req(.not.ok,'crystal rejects inequivalent environment')
 call validate_sawf_structure_class('surface',keys,[0d0,.5d0,1d0],orbit_by_env,ok,msg)
 call req(ok,'surface vacuum class')
 call validate_sawf_structure_class('surface',keys,vacuum_by_env,orbit_by_env,ok,msg)
 call req(.not.ok,'surface requires measured vacuum')
 call validate_sawf_structure_class('interface',keys,vacuum_by_env,orbit_by_env,ok,msg)
 call req(ok,'interface multiple environments')
 open(newunit=unit,file='pseudo-a.dat',status='replace');write(unit,'(a)')'same-content';close(unit)
 open(newunit=unit,file='pseudo-b.dat',status='replace');write(unit,'(a)')'same-content';close(unit)
 call build_sawf_file_content_digest('pseudo-a.dat',digest_a,ok,msg);call req(ok,'pseudo digest A')
 call build_sawf_file_content_digest('pseudo-b.dat',digest_b,ok,msg)
 call req(ok.and.digest_a==digest_b,'same content different path')
 open(newunit=unit,file='pseudo-a.dat',status='replace');write(unit,'(a)')'changed-content';close(unit)
 call build_sawf_file_content_digest('pseudo-a.dat',digest_b,ok,msg)
 call req(ok.and.digest_a/=digest_b,'same path changed content')
 call measure_sawf_vacuum_occupancy([0d0,1d-10,1d-3],1d-8,measured_vacuum,ok,msg)
 call req(ok.and.abs(measured_vacuum-2d0/3d0)<1d-15,'density-based vacuum occupancy')
 write(*,'(a)')'PASS exact same-supercell local-environment fingerprints'
contains
 subroutine local(z,r,vacuum,key)
  integer,intent(in)::z(:);real(8),intent(in)::r(:,:),vacuum;character(256),intent(out)::key
  call build_sawf_local_environment_fingerprint(super_a,lattice,1d-8,z,r,vacuum,key,ok,msg)
  call req(ok,'local fingerprint')
 end subroutine
 subroutine req(c,t);logical,intent(in)::c;character(*),intent(in)::t
  if(.not.c)then;write(*,'(a)')trim(t);error stop 1;endif
 end subroutine
end program'''
with tempfile.TemporaryDirectory(prefix='sawf-env-fingerprint-') as td:
 td=Path(td);(td/'driver.f90').write_text(driver)
 (td/'CMakeLists.txt').write_text(f'''cmake_minimum_required(VERSION 3.16)
project(env_fingerprint LANGUAGES Fortran)
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
 r=subprocess.run([str(td/'b/check')],env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
 if r.returncode:raise SystemExit(r.stdout)
 print(r.stdout.strip())
assert 'call build_sawf_fragment_environment_fingerprints' in generator
assert 'sawf_environment_key(ifrag)==sawf_environment_key(' in generator
assert 'call select_sawf_environment_materialization' in generator
