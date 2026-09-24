"""Compile the production symmetry restart guard with minimal serial stubs."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


class RestartTests(unittest.TestCase):
    def test_metadata_roundtrip_and_rejections(self):
        compiler = shutil.which('gfortran-15') or shutil.which('gfortran')
        if compiler is None:
            self.skipTest('Fortran compiler unavailable')
        text = (ROOT/'src/io/checkpoint_restart.f90').read_text()
        start = text.index('subroutine symmetry_checkpoint_metadata(')
        stop = text.index('end subroutine symmetry_checkpoint_metadata',start)+len('end subroutine symmetry_checkpoint_metadata')
        stubs = '''module structures
type s_dft_system
integer::nk=1
real(8)::vec_k(3,1)=0d0,wtk(1)=1d0
end type
type s_parallel_info
integer::id_rko=0,icomm_rko=0
end type
end module
module sym_sub
logical::use_symmetry=.true.
real(8)::SymMatA(3,4,1)=0d0,SymMatB(3,4,1)=0d0
end module
module salmon_global
character(16)::xc='pz',tdcdft='lrc'
end module
module communication
contains
logical function comm_is_root(id)
integer::id
comm_is_root=id==0
end function
subroutine comm_bcast(status,comm)
integer::status,comm
end subroutine
end module
module extracted
contains
'''
        main = '''program test
use extracted
use structures
use sym_sub
implicit none
type(s_dft_system)::system
type(s_parallel_info)::info
character(20)::mode
integer::j
call get_command_argument(1,mode)
do j=1,3
SymMatA(j,j,1)=1d0;SymMatB(j,j,1)=1d0
enddo
select case(trim(mode))
case('write')
call symmetry_checkpoint_metadata('./',system,info,.true.)
stop
case('weight')
system%wtk=.5d0
case('group')
SymMatA(1,4,1)=.5d0
case('off')
use_symmetry=.false.
end select
call symmetry_checkpoint_metadata('./',system,info,.false.)
end program
'''
        with tempfile.TemporaryDirectory() as directory:
            d=Path(directory)
            (d/'test.f90').write_text(stubs+text[start:stop]+'\nend module\n'+main)
            subprocess.run([compiler,'-J',str(d),str(d/'test.f90'),'-o',str(d/'probe')],check=True,cwd=d)
            def run(mode):
                return subprocess.run([str(d/'probe'),mode],cwd=d,capture_output=True).returncode
            self.assertEqual(run('off'),0)  # Legacy full mesh: no sidecar.
            self.assertNotEqual(run('read'),0)  # Reduced mesh requires metadata.
            self.assertEqual(run('write'),0)
            self.assertEqual(run('read'),0)
            for mode in ('weight','group','off'):
                self.assertNotEqual(run(mode),0,mode)
            (d/'symmetry_restart.bin').write_bytes(b'broken')
            self.assertNotEqual(run('read'),0)


if __name__ == '__main__':
    unittest.main()
