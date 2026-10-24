from pathlib import Path
import subprocess,sys,os,tempfile
root=Path(__file__).resolve().parents[2];b=Path(sys.argv[1]).resolve()
checker=Path('/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/tools/check_forbidden_eigensolvers.py')
subprocess.run(['python3',str(checker),'--source-root',str(root/'src')],check=True)
objects=sorted(x for x in (b/'src/CMakeFiles/salmon.dir').rglob('*.o') if x.name!='main.f90.o')
subprocess.run(['python3',str(checker),'--object-tree',str(b/'src/CMakeFiles/salmon.dir')],check=True)
with tempfile.TemporaryDirectory(prefix='hse-thermal-native-') as d:
 d=Path(d);o=d/'probe.o';exe=d/'probe'
 subprocess.run(['mpifort','-fopenmp','-fcheck=all','-ffree-line-length-none','-I'+str(b),'-c',str(Path(__file__).with_name('factor_orbits_probe.f90')),'-o',str(o)],cwd=d,check=True)
 subprocess.run(['python3',str(checker),'--object-tree',str(d)],check=True)
 subprocess.run(['mpifort','-fopenmp',str(o),*map(str,objects),'-L/opt/homebrew/opt/openblas/lib','-lopenblas','-L/opt/homebrew/opt/fftw/lib','-lfftw3','-L/opt/homebrew/opt/libxc/lib','-lxc','-o',str(exe)],check=True)
 env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
 subprocess.run([str(exe),'startup'],check=True,env=env)
 rejected=subprocess.run([str(exe),'changed'],env=env,text=True,capture_output=True)
 if rejected.returncode==0 or 'Thermal SOC factor history requires fixed occupations' not in rejected.stderr+rejected.stdout:
  raise RuntimeError('changed-occupation guard failed: '+rejected.stdout+rejected.stderr)
 print('PASS changed-occupation collective rejection')
