from pathlib import Path
import os,re,subprocess,tempfile,unittest
ROOT=Path(__file__).resolve().parents[2]
class AutomaticSpin(unittest.TestCase):
 def test_actual_resolver(self):
  text=(ROOT/'src/io/inputoutput.f90').read_text()
  m=re.search(r'(?ims)^  subroutine resolve_ky_spin_inputs\(.*?^  end subroutine resolve_ky_spin_inputs',text)
  self.assertIsNotNone(m,'production auto-input resolver missing')
  with tempfile.TemporaryDirectory() as tmp:
   d=Path(tmp)
   (d/'legacy_rps.dat').write_text('1 1 1 1\n0 1 1\n')
   (d/'soc_rps.dat').write_text('1 1 1 1\nKY_SOC_NLCC_V1\n')
   driver="""program test
 implicit none
 character(256)::files(2),mode
 character(20)::spin
 character::soc
 integer::n
 call get_command_argument(1,mode)
 files=['legacy_rps.dat','soc_rps.dat   ']
 soc='a';spin='auto';n=2
 select case(trim(mode))
 case('scalar');files='legacy_rps.dat'
 case('explicit-y');soc='y';files='legacy_rps.dat'
 case('explicit-spin');spin='noncollinear';files='legacy_rps.dat'
 case('bad-count');n=3
 case('off');soc='n'
 case('polarized');spin='polarized'
 end select
 call resolve_ky_spin_inputs(files,n,spin,soc)
 print *, 'RESOLVED ',trim(spin),' ',soc
 contains
"""+m.group(0)+'\nend program\n'
   (d/'probe.f90').write_text(driver)
   q=subprocess.run([os.environ.get('FC','gfortran'),'-fcheck=all',str(d/'probe.f90'),'-o',str(d/'probe')],cwd=d,capture_output=True,text=True)
   self.assertEqual(q.returncode,0,q.stderr)
   for mode,expected in [('scalar','unpolarized n'),('mixed','noncollinear y'),('explicit-y','noncollinear y'),('explicit-spin','noncollinear y')]:
    q=subprocess.run([str(d/'probe'),mode],cwd=d,capture_output=True,text=True)
    self.assertEqual(q.returncode,0,q.stderr);self.assertIn(expected,' '.join(q.stdout.split()))
   for mode in ['off','polarized','bad-count']:
    q=subprocess.run([str(d/'probe'),mode],cwd=d,capture_output=True,text=True)
    self.assertNotEqual(q.returncode,0,mode)
    if mode=='bad-count':self.assertIn('invalid element count',q.stderr+q.stdout)
if __name__=='__main__':unittest.main()
