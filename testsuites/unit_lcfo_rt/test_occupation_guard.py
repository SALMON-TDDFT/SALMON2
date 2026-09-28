"""Exercise the actual LCFO occupation guard, including nonfinite values.
The Fujitsu compiler crash itself requires compilation on that compiler.
"""
from pathlib import Path
import os
import subprocess
import tempfile

repo = Path(__file__).resolve().parents[2]
source = (repo / 'src/xc/exx_lcfo_rt.f90').read_text()
body = source.split('  subroutine refresh_master(', 1)[1]
guard = body.split('no=size(coeff,2);nsel=size(selected);ng=size(fragment_basis,1)\n', 1)[1].split('    ! On impulse', 1)[0]
probe = '''program guard_probe
 use, intrinsic :: ieee_arithmetic
 implicit none
 type system_type
  real(8), allocatable :: rocc(:,:,:)
 end type
 type(system_type) :: system
 integer :: j, mode
 real(8) :: occupation_value
 character(16) :: arg
 allocate(system%rocc(4,1,1))
 system%rocc(:,1,1)=[0d0,0.5d0,1d0,2d0]
 call get_command_argument(1,arg)
 read(arg,*) mode
 select case(mode)
 case(1)
  system%rocc(2,1,1)=-epsilon(1d0)
 case(2)
  system%rocc(4,1,1)=2d0+4*epsilon(1d0)
 case(3)
  system%rocc(1,1,1)=ieee_value(0d0,ieee_quiet_nan)
 case(4)
  system%rocc(3,1,1)=ieee_value(0d0,ieee_positive_inf)
 case(5)
  system%rocc(4,1,1)=ieee_value(0d0,ieee_negative_inf)
 end select
''' + guard + '\nend program\n'
with tempfile.TemporaryDirectory() as folder:
    p = Path(folder)
    (p/'guard.f90').write_text(probe)
    subprocess.run([os.environ.get('FC', 'gfortran'), '-O0', '-fcheck=all', str(p/'guard.f90'), '-o', str(p/'guard')], cwd=p, check=True)
    for mode in range(6):
        run = subprocess.run([str(p/'guard'), str(mode)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        if mode == 0:
            assert run.returncode == 0, run.stderr
        else:
            assert run.returncode != 0 and 'LCFO EXX: invalid occupations' in run.stderr, (mode, run.stderr)
print('Occupation guard: valid endpoints/fractional, negative, above two, NaN, +/-Inf passed')
