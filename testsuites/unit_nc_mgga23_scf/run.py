from pathlib import Path
import argparse,subprocess,tempfile,os,shlex
p=argparse.ArgumentParser();p.add_argument('--build-dir',type=Path,required=True);p.add_argument('--fc',default='gfortran');p.add_argument('--libs',default='-lopenblas');a=p.parse_args()
root=Path(__file__).resolve().parents[2];b=a.build_dir.resolve();checker=root/'tools/check_forbidden_eigensolvers.py'
subprocess.run(['python3',str(checker),'--source-root',str(root/'src')],check=True)
objects=sorted(x for x in (b/'src/CMakeFiles/salmon.dir').rglob('*.o') if x.name!='main.f90.o')
subprocess.run(['python3',str(checker),'--object-tree',str(b/'src/CMakeFiles/salmon.dir')],check=True)
with tempfile.TemporaryDirectory(prefix='nc23-scf-') as d:
 d=Path(d);o=d/'test.o';exe=d/'test'
 subprocess.run([a.fc,'-cpp','-fopenmp','-fcheck=all','-ffree-line-length-none','-I'+str(b),'-c',str(Path(__file__).with_name('test_scf.f90')),'-o',str(o)],check=True)
 subprocess.run(['python3',str(checker),'--object-tree',str(d)],check=True)
 subprocess.run([a.fc,'-fopenmp',str(o),*map(str,objects),*shlex.split(a.libs),'-o',str(exe)],check=True)
 subprocess.run([str(exe)],check=True,env=dict(os.environ,OMP_NUM_THREADS='1'))

 for mode in ['--nan-local','--nan-tau','--nan-current','--nan-history','--nan-density','--inf-tau','--inf-density']:
  q=subprocess.run([str(exe),mode],capture_output=True,text=True,env=dict(os.environ,OMP_NUM_THREADS='1'))
  if q.returncode==0 or 'nonfinite' not in q.stderr:raise SystemExit('nonfinite residual accepted: '+mode+'\n'+q.stderr)
  print('PASS nonfinite residual rejection:',mode)
