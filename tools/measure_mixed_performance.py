"""Run four sequential nine-rank mixed-mode measurements (manual, not CI).

BUILD EXISTING_SYNTHETIC_GS NEW_OUTPUT [--launcher mpiexec]. One thread per rank.
Requires the embedded host-MPI library and LD/SBE wheel kernels already built.
No downloads, source builds, SCF or GS collection. Output path must be new.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build');parser.add_argument('fixture');parser.add_argument('output',type=Path)
    parser.add_argument('--launcher',default='mpiexec')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    env=os.environ.copy()
    for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):
        env[name]='1'
    records=[];pairs=[]
    for points,steps,stride in ((32,512,1),(64,2048,16)):
        runs=[]
        # Reverse order in the larger case to expose some order/thermal bias.
        modes=('checked','performance') if points==32 else ('performance','checked')
        for mode in modes:
            name=f'{points}points-{steps}steps-{mode}';out=args.output/name
            command=[args.launcher,'-n','9',sys.executable,str(Path(__file__).with_name('test_mixed_performance_mpi.py')),
                     args.build,args.fixture,str(out),mode,str(points),str(steps),str(stride)]
            start=time.perf_counter()
            with (args.output/(name+'.log')).open('wb') as log:
                result=subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=3600)
            wall=time.perf_counter()-start
            if result.returncode:raise RuntimeError(f'calculation failed: {name}; see its log')
            report=json.loads((out/'report.json').read_text())
            report['launcher_seconds']=wall
            report['launcher_log_bytes']=(args.output/(name+'.log')).stat().st_size
            (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
            records.append(dict(name=name,mode=mode,points=points,steps=steps,stride=stride,
                                loop_seconds=report['loop_seconds_max'],launcher_seconds=wall,
                                summed_rank_peak_RSS_bytes=report['summed_rank_peak_RSS_bytes'],
                                io=report['io'],launcher_log_bytes=report['launcher_log_bytes']))
            runs.append((mode,out,report))
            print(name,round(report['loop_seconds_max'],3),'seconds',flush=True)
        a=np.load(runs[0][1]/'traces.npz')['values'];b=np.load(runs[1][1]/'traces.npz')['values']
        exact=np.array_equal(a,b);difference=float(np.max(np.abs(a-b)))
        if not np.allclose(a,b,rtol=2e-10,atol=2e-12):raise AssertionError('mode results differ')
        observation_checks=[]
        for name in ('obs1_at_point_rt.data','obs2_at_point_rt.data'):
            x=np.loadtxt(runs[0][1]/'field'/name);y=np.loadtxt(runs[1][1]/'field'/name)
            if not np.array_equal(x,y):raise AssertionError('LD-field observation differs: '+name)
            observation_checks.append(dict(file=name,bitwise_equal=True,rows=len(x)))
        checked=next(r for mode,_,r in runs if mode=='checked')
        performance=next(r for mode,_,r in runs if mode=='performance')
        pairs.append(dict(points=points,steps=steps,bitwise_equal=exact,max_abs_difference=difference,
                          LD_field_observations=observation_checks,
                          checked_over_performance=checked['loop_seconds_max']/performance['loop_seconds_max']))
        (args.output/'summary.json').write_text(json.dumps(dict(runs=records,comparisons=pairs,
            scope='single Mac; one run per mode/size; instrumented protected updates, not a speedup guarantee'),indent=2)+'\n')
    print('mixed checked/performance comparisons PASS',flush=True)

if __name__=='__main__':main()
