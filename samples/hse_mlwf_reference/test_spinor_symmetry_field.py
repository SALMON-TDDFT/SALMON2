import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

STUBS = "module communication\ncontains\nsubroutine comm_get_globalinfo(a,b,c)\ninteger,intent(out)::a,b,c\na=0;b=0;c=1\nend subroutine\nlogical function comm_is_root(a)\ninteger,intent(in)::a\ncomm_is_root=a==0\nend function\nend module\nmodule salmon_global\ncharacter(32)::xc='hse06',theory='tddft_response'\ncharacter(1)::yn_spinorbit='y'\nreal(8)::epdir_re1(3)=[1d0,0d0,0d0]\nend module\nmodule exx_functional\ncontains\nlogical function is_hybrid(x)\ncharacter(*),intent(in)::x\nis_hybrid=x=='hse06'\nend function\nend module\n"
PROBE = "program probe\nuse sym_sub\nimplicit none\nreal(8)::a(3,3),b(3,3)\ninteger::i\na=0d0;b=0d0\ndo i=1,3\na(i,i)=i+1d0\nb(i,i)=2d0*acos(-1d0)/a(i,i)\nenddo\ncall read_sw_symmetry('yyy')\ncall init_sym_sub(a,b)\nprint *,'X_FIELD_LITTLE_GROUP_PASS'\nend program\n"

class FieldTests(unittest.TestCase):
    def test_actual_x_field_and_bounded_group(self):
        fc=os.environ.get('FC','gfortran')
        if not shutil.which(fc): self.skipTest('gfortran unavailable')
        src=Path(__file__).resolve().parents[2]/'src/symmetry/symmetry.f90'
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            (d/'stub.f90').write_text(STUBS)
            (d/'probe.f90').write_text(PROBE)
            subprocess.run([fc,'-fcheck=all','stub.f90',str(src),'probe.f90','-o','probe'],cwd=d,check=True,capture_output=True)
            def run(rotation):
                (d/'sym.dat').write_text('1 0 0 0\n0 1 0 0\n0 0 1 0\n'+rotation)
                return subprocess.run([str(d/'probe')],cwd=d,text=True,capture_output=True)
            valid=run('1 0 0 0\n0 1 0 0\n0 0 -1 0\n')
            self.assertEqual(valid.returncode,0,valid.stderr)
            self.assertIn('X_FIELD_LITTLE_GROUP_PASS',valid.stdout)
            invalid=run('-1 0 0 0\n0 1 0 0\n0 0 1 0\n')
            self.assertNotEqual(invalid.returncode,0)
            self.assertIn('changes the RT field',invalid.stderr)
            # 90 degree lattice rotations are outside the bounded diagonal spin path.
            square=STUBS.replace("theory='tddft_response'","theory='dft'")
            (d/'stub.f90').write_text(square)
            (d/'probe.f90').write_text(PROBE.replace('a(i,i)=i+1d0','a(i,i)=2d0'))
            subprocess.run([fc,'stub.f90',str(src),'probe.f90','-o','probe'],cwd=d,check=True,capture_output=True)
            (d/'sym.dat').write_text('1 0 0 0\n0 1 0 0\n0 0 1 0\n0 -1 0 0\n1 0 0 0\n0 0 1 0\n-1 0 0 0\n0 -1 0 0\n0 0 1 0\n0 1 0 0\n-1 0 0 0\n0 0 1 0\n')
            bounded=subprocess.run([str(d/'probe')],cwd=d,text=True,capture_output=True)
            self.assertNotEqual(bounded.returncode,0)
            self.assertIn('diagonal',bounded.stderr)

if __name__=='__main__':unittest.main()
