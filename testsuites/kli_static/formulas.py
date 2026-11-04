"""Compile numerical probes in new temporary directories; audit all inputs and objects."""
from pathlib import Path
from tempfile import TemporaryDirectory
import argparse, os, shlex, subprocess, sys

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--checker',required=True,type=Path)
    ap.add_argument('--lapack-flags',default='-llapack -lblas')
    ap.add_argument('--fftw-includes',default='')
    ap.add_argument('--fftw-libs',default='-lfftw3')
    args=ap.parse_args()
    checker=[sys.executable,str(args.checker.resolve())]
    suite=Path(__file__).resolve().parent
    src=suite.parents[1]/'src'
    for tree in [src,suite]:
        subprocess.run(checker+['--source-root',str(tree)],check=True)
    fc=shlex.split(os.environ.get('FC','gfortran'))
    ff=shlex.split(os.environ.get('FFLAGS','-O2'))
    for module,probe,libs in [('kli_diagonal','diagonal','lapack'),('kli_self_hartree','gaussian','fftw3')]:
        with TemporaryDirectory() as out:
            def run(cmd): subprocess.run(cmd,cwd=out,check=True)
            flags=shlex.split(args.fftw_includes) if libs=='fftw3' else []
            link=shlex.split(args.fftw_libs if libs=='fftw3' else args.lapack_flags)
            for f in [src/'xc'/(module+'.f90'),suite/(probe+'.f90')]:
                run(fc+ff+flags+['-c',str(f)])
            run(checker+['--object-tree',out])
            run(fc+[module+'.o',probe+'.o']+link+['-o','probe'])
            run(['./probe'])
if __name__=='__main__': main()
