from pathlib import Path
import tempfile,shutil,subprocess,json
ROOT=Path(__file__).resolve().parents[2]
def run(args):
 r=subprocess.run(args,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
 return r.returncode,r.stdout
results={}
with tempfile.TemporaryDirectory() as d:
 p=Path(d);s=p/'source';b=p/'build';(s/'src').mkdir(parents=True);(s/'tools').mkdir();(s/'cmakefiles').mkdir()
 shutil.copy2(ROOT/'tools/check_forbidden_eigensolvers.py',s/'tools')
 shutil.copy2(ROOT/'cmakefiles/pzheev_guard.cmake',s/'cmakefiles')
 (s/'CMakeLists.txt').write_text('cmake_minimum_required(VERSION 3.14)\nproject(guard_fixture C)\nadd_subdirectory(src)\n')
 (s/'src/main.c').write_text('int main(void){return 0;}')
 (s/'src/CMakeLists.txt').write_text('include("${PROJECT_SOURCE_DIR}/cmakefiles/pzheev_guard.cmake")\nadd_executable(probe main.c)\nsalmon_guard_eigensolver_target(probe)\n')
 r,o=run(['cmake','-S',str(s),'-B',str(b)]);assert r==0,o
 r,o=run(['cmake','--build',str(b)]);assert r==0,o;results['safe_build']='PASS'
 (s/'src/added.f90').write_text('subroutine bad\ncall pzheev()\nend\n')
 r,o=run(['cmake','--build',str(b)]);assert r!=0 and 'PZHEEV is prohibited' in o,o;results['incremental_reintroduction_blocked']='PASS'
 (s/'src/added.f90').unlink()
 (p/'bad.c').write_text('extern void pzheev_(void);void bad(void){pzheev_();}')
 subprocess.run(['cc','-c',str(p/'bad.c'),'-o',str(b/'src/legacy.o')],check=True)
 (s/'src/main.c').write_text('int main(void){return 1;}')
 r,o=run(['cmake','--build',str(b)]);assert r!=0 and 'PZHEEV is prohibited' in o,o
 assert not (b/'src/probe').exists() or (b/'src/probe').stat().st_mode & 0o111==0;results['old_object_blocked_and_new_binary_quarantined']='PASS'
 (s/'src/added.f90').write_text('external pzheev')
 r,o=run(['cmake','-S',str(s),'-B',str(p/'bad-config')]);assert r!=0 and 'PZHEEV prohibition check failed' in o,o;results['configure_blocked']='PASS'
print(json.dumps(results,indent=2))
