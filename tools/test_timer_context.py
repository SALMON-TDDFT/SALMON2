"""Build and run deterministic native timer-ownership checks; no GS data."""
import argparse,json,os,subprocess
from pathlib import Path
import sys


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path)
    p.add_argument('--compiler',default='gfortran')
    p.add_argument('--checker',type=Path)
    args=p.parse_args();root=Path(__file__).resolve().parents[1]
    checker=args.checker or root/'tools/check_forbidden_eigensolvers.py'
    source=root/'runtime/sources/solver/src'
    subprocess.run([sys.executable,str(checker),'--source-root',str(source)],check=True)
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    fixture=root/'tests/native/test_timer_context.f90'
    flags=['-cpp','-fopenmp','-I'+str(out),'-J'+str(out)]
    for name,file,defines in [('stubs',fixture,['-DSALMONPY_TIMER_STUBS']),
                              ('timer',source/'misc/timer.f90',[]),('test',fixture,[])]:
        subprocess.run([args.compiler,*flags,*defines,'-c',str(file),'-o',str(out/(name+'.o'))],check=True)
    subprocess.run([sys.executable,str(checker),'--object-tree',str(out)],check=True)
    exe=out/'test-timer'
    subprocess.run([args.compiler,'-fopenmp',*[str(out/(n+'.o')) for n in ('stubs','timer','test')],'-o',str(exe)],check=True)
    subprocess.run([str(exe)],env={**os.environ,'OMP_NUM_THREADS':'2'},check=True)
    (out/'result.json').write_text(json.dumps(dict(status='pass',deterministic_clock=True,
        host_and_sibling_isolated=True,running_and_thread_timers_pause=True,subrange_preserved=True,
        GS_data_collected=False))+'\n')


if __name__=='__main__':main()
