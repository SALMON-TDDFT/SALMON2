from pathlib import Path
import subprocess,os,tempfile,sys
mpi="--mpi" in sys.argv
fc="mpifort" if mpi else "gfortran"
root=Path(__file__).resolve().parents[2]
checker=Path('/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/tools/check_forbidden_eigensolvers.py')
subprocess.run(['python3',str(checker),'--source-root',str(root/'src')],check=True)
with tempfile.TemporaryDirectory(prefix='hse-thermal-probe-') as d:
 d=Path(d); objects=[]
 for name in ['exx_collinear_fft','exx_k_backend','exx_k_exchange','exx_spinor','exx_ace','probe']:
  src=root/'src/xc'/f'{name}.f90' if name!='probe' else Path(__file__).with_name('probe_mpi.f90' if mpi else 'probe.f90')
  obj=d/f'{name}.o';objects.append(obj)
  subprocess.run([fc,'-O1','-g','-fopenmp','-fcheck=all','-ffree-line-length-none','-I/opt/homebrew/opt/fftw/include','-c',str(src),'-o',str(obj)],cwd=d,check=True)
 subprocess.run(['python3',str(checker),'--object-tree',str(d)],check=True)
 exe=d/'probe'
 subprocess.run([fc,'-fopenmp',*map(str,objects),'-L/opt/homebrew/opt/fftw/lib','-lfftw3','-L/opt/homebrew/opt/openblas/lib','-lopenblas','-o',str(exe)],check=True)
 for mode in (['','nan'] if mpi else ['', 'nan','inf','negative','overfull']):
  subprocess.run((['mpiexec','-n','2'] if mpi else [])+[str(exe),mode],check=True,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
