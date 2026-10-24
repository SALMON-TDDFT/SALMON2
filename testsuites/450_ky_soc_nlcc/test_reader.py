from pathlib import Path
import unittest,subprocess,tempfile,os,re
ROOT=Path(__file__).resolve().parents[2]
STUB='''module structures
 implicit none
 type s_pp_info
 integer::nrmax0=2100,lmax0=3,lmax=3,nrmax=2102
 integer::mr(1),mlps(1),zps(1),lref(1),nproj(0:7,1)
 real(8)::zion,rad(2102,1),vpp(0:2100,0:7),vpp_so(0:2100,0:7),upp(0:2100,0:7)
 real(8)::anorm(0:7,1),anorm_so(0:7,1)
 logical::core23_present(1)=.false.
 real(8)::core23_tau_tbl(2102,1)=0d0,core23_current_tbl(2102,1)=0d0
 end type
end module
module salmon_global
 implicit none
 integer::nelem=1,Lmax_ps(1)=-1
 character::yn_spinorbit='y'
end module
module inputoutput
 implicit none
 real(8),parameter::au_length_aa=.529177210903d0,au_energy_ev=27.211386245988d0
end module
module actual_reader
 implicit none
 logical::flag_beta_proj_is_given=.false.
 contains
'''
DRIVER='''program test
 use actual_reader,only:read_ps_ky
 use structures,only:s_pp_info
 use salmon_global,only:yn_spinorbit
 implicit none
 type(s_pp_info)::p
 real(8)::rc(0:3),rho(0:2100,0:2)
 logical::flags(1)=.false.
 character(256)::file,mode
 call get_command_argument(1,file)
 call get_command_argument(2,mode)
 yn_spinorbit=mode(1:1)
 p%vpp=0d0;p%vpp_so=0d0;p%upp=0d0;p%anorm=0d0;p%anorm_so=0d0
 call read_ps_ky(p,rc,rho,flags,1,file)
 print *, 'READ_PASS',p%zps(1),p%mlps(1),flags
 print *, 'CORE',rho(0,:),rho(1,:)
 print *, 'NORMS',p%anorm(:,1),p%anorm_so(:,1)
end program
'''
def extract(text,name):
 pattern=r'(?ims)^subroutine '+name+r'\(.*?^end subroutine '+name+r'\s*$'
 m=re.search(pattern,text);assert m,'missing production reader '+name
 return m.group(0)
class Reader(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.tmp=tempfile.TemporaryDirectory();cls.d=Path(cls.tmp.name);text=(ROOT/'src/atom/pp/input_pp.f90').read_text()
  source=STUB+extract(text,'read_ps_ky')+'\n'+extract(text,'read_ps_ky_soc')+'\n'+extract(text,'read_ky_core23')+'\nend module\n'+DRIVER
  (cls.d/'probe.f90').write_text(source);cmd=[os.environ.get('FC','gfortran'),'-fcheck=all',str(cls.d/'probe.f90'),'-o',str(cls.d/'probe')]
  subprocess.run(cmd,cwd=cls.d,capture_output=True,text=True,check=True)
 @classmethod
 def tearDownClass(cls):cls.tmp.cleanup()
 def run_reader(self,p,mode='y'):
  return subprocess.run([str(self.d/'probe'),str(p),mode],capture_output=True,text=True,timeout=15)
 def test_standard_W_Se_and_NLCC_units(self):
  for sp,z,l in [('W',12,3),('Se',6,2)]:
   p=ROOT/f'testsuites/pseudo/{sp}_rps.dat';q=self.run_reader(p);self.assertEqual(q.returncode,0,q.stderr)
   a=next(x for x in q.stdout.splitlines() if 'READ_PASS' in x).split();self.assertEqual([int(a[1]),int(a[2])],[z,l]);self.assertEqual(a[3],'T')
   lines=p.read_text().splitlines();nc=int(lines[3].split()[2]);rows=[list(map(float,x.split())) for x in lines[4+nc:-1]]
   core=list(map(float,next(x for x in q.stdout.splitlines() if 'CORE' in x).split()[1:]));length=.529177210903
   expected=[rows[i][j]*length**(j+1) for i in [0,1] for j in [2,3,4]]
   for x,y in zip(core,expected):self.assertAlmostEqual(x,y,delta=1e-11*max(1,abs(y)))
   norms=list(map(float,next(x for x in q.stdout.splitlines() if 'NORMS' in x).split()[1:]));h=.01/length
   for k in range(nc):
    ll=k//2 if k%2==0 else (k+1)//2;beta=[(r[5+2*k]-r[1])/27.211386245988*r[6+2*k]*length**.5 for r in rows]
    if max(map(abs,beta))<1e-12:expected=0.
    else:
     val=[r[6+2*k]*length**.5*b for r,b in zip(rows,beta)];expected=1/(h*(sum(val)-.5*(val[0]+val[-1])))
    idx=ll if k%2==0 else 8+ll;self.assertAlmostEqual(norms[idx],expected,delta=1e-10*max(1,abs(expected)))
 def test_legacy_KY_scalar_and_SOC_kept(self):
  for mode in ['n','y']:
   q=self.run_reader(ROOT/'testsuites/pseudo/Si_rps.dat',mode);self.assertEqual(q.returncode,0,q.stderr)
   a=next(x for x in q.stdout.splitlines() if 'READ_PASS' in x).split();self.assertEqual(int(a[1]),4);self.assertEqual(a[3],'F')
 def test_SOCoff_rejected(self):
  q=self.run_reader(ROOT/'testsuites/pseudo/Se_rps.dat','n');self.assertNotEqual(q.returncode,0);self.assertIn('requires yn_spinorbit=y',q.stderr+q.stdout)
 def test_corrupt_units_nan_and_truncation_rejected(self):
  text=(ROOT/'testsuites/pseudo/Se_rps.dat').read_text();lines=text.splitlines();nc=int(lines[3].split()[2]);row=4+nc;bad=list(lines);a=bad[row].split();a[5]='NaN';bad[row]=' '.join(a)
  for name,s in [('units',text.replace('ANGSTROM EV','BOHR HARTREE')),('truncated','\n'.join(lines[:-2])),('nan','\n'.join(bad)+'\n')]:
   p=self.d/(name+'_rps.dat');p.write_text(s);q=self.run_reader(p);self.assertNotEqual(q.returncode,0,name)
if __name__=='__main__':unittest.main()
