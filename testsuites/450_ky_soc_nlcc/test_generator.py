from pathlib import Path
import unittest,subprocess,tempfile,os,math,re
ROOT=Path(__file__).resolve().parents[2];SOURCE=ROOT/'utility/crps_soc_nlcc'
class Generator(unittest.TestCase):
 def test_Mo_unused_empty_f_orbitals(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)
   q=subprocess.run([os.environ.get('FC','gfortran'),'-std=legacy','-O0','-fcheck=all',
    '-Wno-align-commons',str(SOURCE/'crps_output.f90'),str(SOURCE/'crps_soc_nlcc.f'),
    '-o',str(p/'crps')],cwd=p,capture_output=True,text=True)
   self.assertEqual(q.returncode,0,q.stderr)
   inp="""&crps ZA=42,
 occA=2,2,2,4,2,2,4,4,6,2,2,4,2,3,0,0,1,9*0,
 PtoA=17,11,12,13,14,cutoff_radius_bohr=5*2.6,
 core_radius=.6,c2hi=7*5.0,output_file='Mo.dat' /
"""
   (p/'Mo.nml').write_text(inp)
   q=subprocess.run([str(p/'crps'),'Mo.nml'],cwd=p,capture_output=True,text=True,timeout=120)
   self.assertEqual(q.returncode,0,q.stderr+q.stdout[-1000:])
   self.assertTrue((p/'Mo.dat').exists(),q.stdout[-1200:])
   self.assertIn('CRPS_SOC_NLCC_WRITTEN',q.stdout)
   self.assertNotIn('occupied orbitals energy is positive',q.stdout)
   channels=(p/'Mo.dat').read_text().splitlines()[4:9]
   self.assertTrue(all(float(x.split()[4])<0 for x in channels))
   lines=(p/'Mo.dat').read_text().splitlines()
   self.assertAlmostEqual(float(lines[0].split()[3]),12.)
   for line in lines[9:-1]:
    self.assertTrue(all(math.isfinite(float(x)) for x in line.split()))
 def test_even_origin_writer(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);(p/'fixture.f90').write_text('''program test
 use crps_output,only:crps_write_v1
 implicit none
 integer,parameter::n=8,nc=3
 real(8)::rho(0:n),v(0:n,nc),u(0:n,nc),vl(0:n),occ(nc),eps(nc)
 integer::l(nc),j(nc),rc(nc),i
 do i=0,n
  rho(i)=2d0+(.01d0*i)**2
 enddo
 vl=-1d0;v=-2d0;u=0d0;l=[0,1,1];j=[1,1,3];rc=2;occ=[2d0,2d0,4d0];eps=-1d0
 call crps_write_v1('fixture.dat',n,.01d0,1,8d0,nc,1,1,l,j,rc,occ,eps,rho,vl,v,u)
end program
''')
   q=subprocess.run([os.environ.get('FC','gfortran'),'-fcheck=all',str(SOURCE/'crps_output.f90'),str(p/'fixture.f90'),'-o',str(p/'probe')],cwd=p,capture_output=True,text=True);self.assertEqual(q.returncode,0,q.stderr)
   subprocess.run([str(p/'probe')],cwd=p,check=True);rows=[list(map(float,l.split())) for l in (p/'fixture.dat').read_text().splitlines()[7:-1]]
   self.assertAlmostEqual(rows[1][3],.02,delta=1e-11);self.assertEqual(rows[0][3],0.)
   for r in rows[:7]:self.assertAlmostEqual(r[4],2.,delta=1e-9)
 def test_selected_W_Se_reproducible(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);q=subprocess.run([os.environ.get('FC','gfortran'),'-std=legacy','-O0','-fcheck=all','-fbacktrace','-Wno-align-commons',str(SOURCE/'crps_output.f90'),str(SOURCE/'crps_soc_nlcc.f'),'-o',str(p/'crps')],cwd=p,capture_output=True,text=True);self.assertEqual(q.returncode,0,q.stderr)
   for sp in ['W','Se']:
    explicit=p/f'explicit-{sp}.nml'
    explicit.write_text((SOURCE/'inputs'/f'{sp}.nml').read_text().replace('&crps',
      f"&crps NA={26 if sp=='W' else 14},Lmax={3 if sp=='W' else 2},lloc={0 if sp=='W' else 2},output_file='{sp}_soc_nlcc_rps.dat',"))
    d=p/sp;d.mkdir();q=subprocess.run([str(p/'crps'),str(explicit)],cwd=d,capture_output=True,text=True,timeout=120);self.assertEqual(q.returncode,0,q.stderr)
    actual=(d/f'{sp}_soc_nlcc_rps.dat').read_text().splitlines()
    expected=(ROOT/f'testsuites/pseudo/{sp}_rps.dat').read_text().splitlines()
    self.assertEqual(len(actual),len(expected))
    for lineno,(a,b) in enumerate(zip(actual,expected),1):
     try:
      av=list(map(float,a.split()));bv=list(map(float,b.split()))
     except ValueError:
      self.assertEqual(a,b,(sp,lineno));continue
     self.assertEqual(len(av),len(bv),(sp,lineno))
     for x,y in zip(av,bv):
      self.assertTrue(math.isfinite(x) and math.isfinite(y),(sp,lineno))
      self.assertLessEqual(abs(x-y),1e-8*max(1.,abs(y)),(sp,lineno,x,y))
    # Omitted controls must resolve to exactly the explicit species definition.
    minimal=(SOURCE/'inputs'/f'{sp}.nml').read_text()
    minimal=re.sub(r'NA=\d+,?', '', minimal)
    minimal=re.sub(r'Lmax=\d+,?', '', minimal)
    minimal=re.sub(r"output_file='[^']*'", '', minimal)
    if sp=='Se':
     minimal=minimal.replace('PtoA=10,11,12,13,14,15,16','PtoA=10,11,12,13,14')
     minimal=minimal.replace('cutoff_radius_bohr=2.6,2.7,2.7,2.3,2.3,2.5,2.5',
                             'cutoff_radius_bohr=2.6,2.7,2.7,2.3,2.3')
    auto=p/f'auto-{sp}.nml';auto.write_text(minimal)
    ad=p/f'auto-{sp}';ad.mkdir()
    q=subprocess.run([str(p/'crps'),str(auto)],cwd=ad,capture_output=True,text=True,timeout=120)
    self.assertEqual(q.returncode,0,q.stderr+q.stdout[-1500:])
    self.assertEqual((ad/f'auto-{sp}_soc_nlcc_rps.dat').read_text().splitlines(),actual)
    # Discrete header and channel identities are exact; libm-dependent radial
    # values use a portable tolerance. Byte identity is separate provenance.
    self.assertEqual(actual[0].split()[::2],expected[0].split()[::2])
    self.assertEqual(actual[3],expected[3])
    for a,b in zip(actual[4:4+int(expected[3].split()[2])],expected[4:4+int(expected[3].split()[2])]):
     self.assertEqual(a.split()[:2],b.split()[:2])
if __name__=='__main__':unittest.main()
