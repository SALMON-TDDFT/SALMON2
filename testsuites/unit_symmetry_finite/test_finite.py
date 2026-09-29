"""Exercise the actual group validator; FC/FFLAGS may select NVHPC on MIYABI."""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


class SymmetryFinite(unittest.TestCase):
    def test_group_validation(self):
        source = (ROOT / 'src/symmetry/symmetry.f90').read_text()
        start = source.index('  subroutine symmetry_validate_group()')
        end = source.index('  end subroutine symmetry_validate_group', start)
        routine = source[start:end] + '  end subroutine symmetry_validate_group\n'
        helper_start = source.index('  pure logical function finite_real_2d(')
        helper_end = source.index('  end function', helper_start) + len('  end function')
        helper = source[helper_start:helper_end] + '\n'
        # Keep the allocatable array at module scope, as in the production code.
        module = '''module probe_group
implicit none
real(8), allocatable :: SymMatA(:,:,:)
real(8) :: Amat(3,3), Ainv(3,3)
interface salmon_all_finite
  module procedure finite_real_2d
end interface
contains
''' + routine + helper + 'end module probe_group\n'
        driver = '''program probe
use probe_group
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
implicit none
character(32) :: mode
integer :: i
real(8) :: bad
call get_command_argument(1,mode)
allocate(SymMatA(3,4,1))
SymMatA=0d0
Amat=0d0
Ainv=0d0
do i=1,3
  SymMatA(i,i,1)=1d0
  Amat(i,i)=1d0
  Ainv(i,i)=1d0
enddo
if (index(mode,'nan')>0) bad=ieee_value(0d0,ieee_quiet_nan)
if (index(mode,'posinf')>0) bad=ieee_value(0d0,ieee_positive_inf)
if (index(mode,'neginf')>0) bad=ieee_value(0d0,ieee_negative_inf)
if (index(mode,'rotation')>0) SymMatA(2,3,1)=bad
if (index(mode,'translation')>0) SymMatA(3,4,1)=bad
call symmetry_validate_group()
print *, 'VALID'
end program probe
'''
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)
            (path / 'probe.f90').write_text(module + driver)
            command = shlex.split(os.environ.get('FC', 'gfortran'))
            command += shlex.split(os.environ.get('FFLAGS', '-O2'))
            subprocess.run(command + ['probe.f90', '-o', 'probe'], cwd=path, check=True)
            for mode in ['finite'] + [a + '-' + b for a in ['nan', 'posinf', 'neginf']
                                      for b in ['rotation', 'translation']]:
                with self.subTest(mode=mode):
                    run = subprocess.run([str(path / 'probe'), mode], cwd=path,
                                         capture_output=True, text=True)
                    if mode == 'finite':
                        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
                        self.assertIn('VALID', run.stdout)
                    else:
                        self.assertNotEqual(run.returncode, 0)
                        self.assertIn('Symmetry: nonfinite operation', run.stdout + run.stderr)


if __name__ == '__main__':
    unittest.main()
