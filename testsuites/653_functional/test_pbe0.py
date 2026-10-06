"""PBE0 exchange definition and incompatible LCFO seed rejection."""
from pathlib import Path
import subprocess,tempfile,unittest
ROOT=Path(__file__).resolve().parents[2]
class Pbe0(unittest.TestCase):
 def test_definition_and_metadata(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)
   (p/'globals.f90').write_text('''module salmon_global
character(32) :: xc
real(8) :: hse_omega=.11d0,pbeh_coulomb_radius=4d0,rvv10_b=5.3d0,rvv10_c=.0093d0,exx_mlwf_radius=0d0
integer :: rvv10_nq=32
end module
''')
   (p/'probe.f90').write_text('''program probe
use salmon_global
use exx_functional
implicit none
integer :: status
xc='pbe0'
if(.not.is_global_hybrid(xc).or..not.is_hybrid(xc))error stop 'PBE0 classification'
if(is_global_hybrid('hse06').or..not.is_hybrid('hse06'))error stop 'HSE classification'
if(is_hybrid('pbe').or.is_global_hybrid('pbe'))error stop 'PBE classification'
if(.not.is_global_hybrid('pbeh40_rvv10'))error stop 'rVV10 classification'
if(abs(exchange_fraction()-.25d0)>1d-14)error stop 'PBE0 fraction'
if(exchange_screening()/=0d0)error stop 'PBE0 screening'
call lcfo_check_functional('missing','run',status)
if(status==0)error stop 'missing PBE0 metadata accepted'
call lcfo_write_functional('meta','run',status)
if(status/=0)error stop 'write'
call lcfo_check_functional('meta','run',status)
if(status/=0)error stop 'round trip'
pbeh_coulomb_radius=5d0
call lcfo_check_functional('meta','run',status)
if(status==0)error stop 'different cutoff accepted'
pbeh_coulomb_radius=4d0
xc='pbeh40'
if(abs(exchange_fraction()-.4d0)>1d-14.or.exchange_screening()/=0d0)error stop 'PBEh regression'
call lcfo_check_functional('meta','run',status)
if(status==0)error stop 'different functional accepted'
xc='hse06'
if(abs(exchange_fraction()-.25d0)>1d-14.or.exchange_screening()/=.11d0)error stop 'HSE regression'
end program
''')
   subprocess.run(['gfortran','-ffree-line-length-none',str(p/'globals.f90'),str(ROOT/'src/xc/exx_functional.f90'),str(p/'probe.f90'),'-o',str(p/'probe')],cwd=d,check=True,capture_output=True)
   run=subprocess.run([str(p/'probe')],cwd=d,capture_output=True,text=True)
   self.assertEqual(run.returncode,0,run.stdout+run.stderr)
if __name__=='__main__':unittest.main()
