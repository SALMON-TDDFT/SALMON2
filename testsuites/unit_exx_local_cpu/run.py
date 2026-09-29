"""CPU pair workspaces with and without OpenMP, including dynamic teams."""
import os,subprocess,tempfile
from pathlib import Path
root=Path(__file__).resolve().parents[2]
fftw=Path(os.environ.get('FFTW_ROOT','/opt/homebrew/opt/fftw'))
with tempfile.TemporaryDirectory(prefix='local-cpu-') as tmp:
    exe=Path(tmp)/'probe'
    for flags in [[],['-fopenmp']]:
        subprocess.run([os.environ.get('FC','gfortran'),'-O2','-fcheck=all',*flags,
            '-I'+str(fftw/'include'),str(root/'src/xc/exx_local_fft.f90'),
            str(Path(__file__).with_name('probe.f90')),'-L'+str(fftw/'lib'),'-lfftw3','-o',str(exe)],
            cwd=tmp,check=True)
        for dynamic in ['FALSE','TRUE']:
            subprocess.run([str(exe)],env=dict(os.environ,OMP_DYNAMIC=dynamic),check=True)
