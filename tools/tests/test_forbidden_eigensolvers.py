from pathlib import Path
import importlib.util,unittest,tempfile,subprocess,sys,os
ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('guard',ROOT/'tools/check_forbidden_eigensolvers.py');g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)
class GuardTest(unittest.TestCase):
 def check(self,source,bad,suffix='.f90'):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/('probe'+suffix);p.write_text(source)
   self.assertEqual(bool(g.source_findings(d)),bad)
 def test_direct(self):self.check('call PZHEEV("V", "U")',True)
 def test_external(self):self.check('external :: pzheev',True)
 def test_continued(self):self.check('call pzhe&\n &ev()',True)
 def test_continued_call(self):self.check('call &\n pzheev()',True)
 def test_comments_allowed(self):self.check('! call pzheev()\nprint *, "PZHEEV forbidden"',False)
 def test_allowed_driver(self):self.check('call pzheevd()\nexternal :: pzheevd',False)
 def test_binding(self):self.check('subroutine x() bind(C,name="pzheev_")',True)
 def test_c_external(self):self.check('extern void pzheev_(void);',True,'.c')
 def test_c_comments(self):self.check('/* pzheev_(); */\nconst char *x="pzheev";',False,'.c')
 def test_missing_fails(self):
  with self.assertRaises(RuntimeError):g.source_findings('/definitely/missing/pzheev-guard')
 def test_object_reference(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'bad.c').write_text('extern void pzheev_(void); void f(void){pzheev_();}')
   subprocess.run(['cc','-c',str(p/'bad.c'),'-o',str(p/'bad.o')],check=True)
   self.assertEqual(g.object_findings([p/'bad.o']),[str(p/'bad.o')])
   (p/'good.c').write_text('extern void pzheevd_(void); void f(void){pzheevd_();}')
   subprocess.run(['cc','-c',str(p/'good.c'),'-o',str(p/'good.o')],check=True)
   self.assertEqual(g.object_findings([p/'good.o']),[])
 def test_rejected_binary_quarantined(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'bad.c').write_text('extern void pzheev_(void); void f(void){pzheev_();}')
   subprocess.run(['cc','-c',str(p/'bad.c'),'-o',str(p/'bad.o')],check=True)
   binary=p/'salmon';binary.write_bytes(b'fixture');binary.chmod(0o755)
   r=subprocess.run([sys.executable,str(ROOT/'tools/check_forbidden_eigensolvers.py'),'--object-tree',d,'--built-executable',str(binary)],capture_output=True,text=True)
   self.assertEqual(r.returncode,1);self.assertEqual(binary.stat().st_mode & 0o111,0)
 def test_multiline_dlsym(self):self.check('void *p=dlsym(\n handle,\n "pzheev_"\n);',True,'.c')
 def test_multiline_binding(self):self.check('subroutine x() bind(C, &\n &name="pzheev_")',True)
 def test_dynamic_documentation_not_code(self):self.check('const char *x="dlsym(handle, \\"pzheev_\\")";',False,'.c')
 def test_unreadable_directory(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'unreadable';p.mkdir();(p/'bad.f90').write_text('call pzheev()');p.chmod(0)
   try:
    if os.access(p,os.R_OK):self.skipTest('privileged runner can read mode000')
    with self.assertRaises(OSError):g.source_findings(d)
   finally:p.chmod(0o700)
if __name__=='__main__':unittest.main()
