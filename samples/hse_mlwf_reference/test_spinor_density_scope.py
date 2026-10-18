"""Compile the production density routine and test its functional boundary."""
from pathlib import Path
import os
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]

class DensityScopeTest(unittest.TestCase):
    def test_nonhybrid_non_diagonal_symmetry_keeps_existing_density(self):
        source = (ROOT / 'src/so/noncollinear_module.f90').read_text()
        start = source.index('  subroutine calc_dm_noncollinear(')
        end = source.index('  end subroutine calc_dm_noncollinear', start)
        routine = source[start:end] + '  end subroutine calc_dm_noncollinear\n'
        stubs = '''module structures
 type s_dft_system
 real(8)::rocc(1,1,1)=1,wtk(1)=1
 end type
 type s_parallel_info
 integer::im_s=1,im_e=1,ik_s=1,ik_e=1,io_s=1,io_e=1,icomm_ko=0
 end type
 type s_rgrid
 integer::is(3)=1,ie(3)=1,num(3)=1
 end type
 type s_orbital
 complex(8)::zwf(1,1,1,2,1,1,1)
 end type
end module
module communication
contains
 subroutine comm_summation(a,b,n,c)
 complex(8),intent(in)::a(:,:,:,:,:)
 complex(8),intent(out)::b(:,:,:,:,:)
 integer,intent(in)::n,c
 b=a
 end subroutine
end module
module salmon_global
 character(16)::xc='pbe'
end module
module sym_sub
 logical::use_symmetry=.true.
 real(8)::SymMatA(3,4,1)=0,SymMatB(3,4,1)=0
end module
module density
 complex(8),allocatable::den_mat(:,:,:,:,:)
 complex(8),parameter::zero=(0d0,0d0)
contains
'''
        program = '''end module
program probe
 use density
 use structures
 use salmon_global
 use sym_sub
 type(s_dft_system)::system
 type(s_parallel_info)::info
 type(s_rgrid)::mg
 type(s_orbital)::psi
 call get_command_argument(1,xc)
 ! A valid non-diagonal proper rotation, previously accepted for PBE.
 SymMatA(1,2,1)=-1;SymMatA(2,1,1)=1;SymMatA(3,3,1)=1
 SymMatB=SymMatA
 psi%zwf=0;psi%zwf(1,1,1,1,1,1,1)=1
 call calc_dm_noncollinear(psi,system,info,mg)
 if(abs(den_mat(1,1,1,1,1)-1)>1d-14)stop 2
 print *, 'PBE_DENSITY_COMPATIBILITY_PASS'
end program
'''
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            (path/'probe.f90').write_text(stubs + routine + program)
            compile_result = subprocess.run([os.environ.get('FC','gfortran'), '-cpp', '-fcheck=all',
                str(ROOT/'src/xc/exx_symmetry.f90'), 'probe.f90', '-o', 'probe'],
                cwd=path, capture_output=True, text=True)
            self.assertEqual(compile_result.returncode, 0, compile_result.stderr)
            pbe = subprocess.run(['./probe','pbe'], cwd=path, capture_output=True, text=True)
            self.assertEqual(pbe.returncode, 0, pbe.stdout+pbe.stderr)
            hse = subprocess.run(['./probe','hse06'], cwd=path, capture_output=True, text=True)
            self.assertNotEqual(hse.returncode, 0)
            self.assertIn('SOC symmetry requires signed diagonal rotations',hse.stderr)
