import os,subprocess,tempfile
from pathlib import Path
root=Path(__file__).resolve().parents[3]
with tempfile.TemporaryDirectory(prefix='ace-k-') as tmp:
    for flags in [[],['-fopenmp']]:
        subprocess.run([os.environ.get('FC','gfortran'),'-O0','-fcheck=all',*flags,
            str(root/'src/xc/exx_ace.f90'),str(Path(__file__).with_name('probe.f90')),
            '-L'+os.environ.get('OPENBLAS_ROOT','/opt/homebrew/opt/openblas')+'/lib',
            '-lopenblas','-o',tmp+'/probe'],cwd=tmp,check=True)
        subprocess.run([tmp+'/probe'],env=dict(os.environ,OPENBLAS_NUM_THREADS='1'),check=True)

    for enabled in (False,True):
        (Path(tmp)/'config.h').write_text('#define HAVE_EXX_OPENBLAS_THREADS\n' if enabled else '')
        subprocess.run([os.environ.get('FC','gfortran'),'-O0','-fcheck=all','-fopenmp','-cpp','-I'+tmp,
            str(root/'src/xc/exx_ace.f90'),str(root/'src/xc/exx_blas_threads.f90'),
            str(Path(__file__).with_name('backend_probe.f90')),
            '-L'+os.environ.get('OPENBLAS_ROOT','/opt/homebrew/opt/openblas')+'/lib',
            '-lopenblas','-o',tmp+'/probe'],cwd=tmp,check=True)
        subprocess.run([tmp+'/probe'],env=dict(os.environ,OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='3'),check=True)

    for macro,probe in [('FUJITSU','scope_probe.f90'),('NVPL','nvpl_probe.f90')]:
        (Path(tmp)/'config.h').write_text('#define HAVE_EXX_'+macro+'_THREADS\n')
        subprocess.run([os.environ.get('FC','gfortran'),'-O0','-fcheck=all','-fopenmp','-cpp','-I'+tmp,
            str(root/'src/xc/exx_ace.f90'),str(root/'src/xc/exx_blas_threads.f90'),
            str(Path(__file__).with_name(probe)),
            '-L'+os.environ.get('OPENBLAS_ROOT','/opt/homebrew/opt/openblas')+'/lib',
            '-lopenblas','-o',tmp+'/probe'],cwd=tmp,check=True)
        subprocess.run([tmp+'/probe'],env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='3'),check=True)
