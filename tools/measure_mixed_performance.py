"""Sequential nine-rank mixed-mode measurements (manual, not CI).

BUILD EXISTING_SYNTHETIC_GS NEW_OUTPUT [--study scaling|observation].
Scaling: 32 points/512 steps and64 points/2048 steps. Observation: fixed64
points/1024 steps, strides1/16, both modes, two counterbalanced repetitions.
Requires prebuilt host-MPI and LD/SBE libraries. Starts compute jobs immediately.
No source build, download, SCF or GS collection. One thread per rank.
"""
import argparse
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time
import numpy as np


def compare_runs(root, records):
    comparisons=[]
    for points,steps in sorted({(r['points'],r['steps']) for r in records}):
        cases=[r for r in records if (r['points'],r['steps'])==(points,steps)]
        reference=min(cases,key=lambda r:r['stride'])
        baseline=np.load(root/reference['name']/'traces.npz')['values']
        for case in cases:
            actual=np.load(root/case['name']/'traces.npz')['values']
            if not np.array_equal(baseline,actual):
                raise AssertionError('every-step A/E/current differs: '+case['name'])
            observations=[]
            for name in ('obs1_at_point_rt.data','obs2_at_point_rt.data'):
                full=np.loadtxt(root/reference['name']/'field'/name)
                sparse=np.loadtxt(root/case['name']/'field'/name)
                full_steps=np.rint(full[:,0]/.003).astype(np.int64)
                sparse_steps=np.rint(sparse[:,0]/.003).astype(np.int64)
                indices=np.searchsorted(full_steps,sparse_steps)
                if (np.any(indices>=len(full_steps)) or
                        not np.array_equal(full_steps[indices],sparse_steps) or
                        not np.array_equal(full[indices],sparse)):
                    raise AssertionError('LD-region E/H differs at common steps: '+case['name'])
                observations.append(dict(file=name,rows=len(sparse),common_steps_bitwise_equal=True))
            comparisons.append(dict(reference=reference['name'],case=case['name'],
                every_step_bitwise_equal=True,max_abs_difference=0.,observations=observations))
    return comparisons


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build')
    parser.add_argument('fixture')
    parser.add_argument('output',type=Path)
    parser.add_argument('--launcher',default='mpiexec')
    parser.add_argument('--study',choices=('scaling','observation'),default='scaling')
    parser.add_argument('--repeats',type=int)
    args=parser.parse_args()
    repeats=args.repeats if args.repeats is not None else (2 if args.study=='observation' else 1)
    if repeats<1:parser.error('--repeats must be positive')
    args.output.mkdir(parents=True,exist_ok=False)
    env=os.environ.copy()
    for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):
        env[name]='1'
    base=([(32,512,1,'checked'),(32,512,1,'performance'),
           (64,2048,16,'performance'),(64,2048,16,'checked')]
          if args.study=='scaling' else
          [(64,1024,1,'checked'),(64,1024,16,'performance'),
           (64,1024,16,'checked'),(64,1024,1,'performance')])
    records=[]
    for repeat in range(1,repeats+1):
        cases=base if repeat%2 else list(reversed(base))
        for points,steps,stride,mode in cases:
            name=f'{points}points-{steps}steps-stride{stride}-{mode}-repeat{repeat}'
            out=args.output/name
            command=[args.launcher,'-n','9',sys.executable,
                     str(Path(__file__).with_name('test_mixed_performance_mpi.py')),
                     args.build,args.fixture,str(out),mode,str(points),str(steps),str(stride)]
            start=time.perf_counter()
            log_path=args.output/(name+'.log')
            with log_path.open('wb') as log:
                result=subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=3600)
            wall=time.perf_counter()-start
            if result.returncode:raise RuntimeError(f'calculation failed: {name}; see its log')
            report=json.loads((out/'report.json').read_text())
            if report['schema']!='salmonpy.mixed-performance-measurement/2':
                raise RuntimeError('this controller requires the matching phase measurement tool')
            report['launcher_seconds']=wall
            try:
                report['launcher_log_bytes']=log_path.stat().st_size
            except OSError as exc:
                report['launcher_log_bytes']=None
                report['resource_measurement_status']='partial'
                report['resource_measurement_errors'].append('launcher log size: '+str(exc))
            (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
            records.append(dict(name=name,mode=mode,points=points,steps=steps,stride=stride,repeat=repeat,
                loop_seconds=report['loop_seconds_max'],launcher_seconds=wall,
                summed_rank_peak_RSS_bytes=report['summed_rank_peak_RSS_bytes'],
                resource_measurement_status=report['resource_measurement_status'],
                phase_accounting_complete=report['phase_accounting_complete'],
                io=report['io'],launcher_log_bytes=report['launcher_log_bytes']))
            # Keep completed records even if a later job fails. No loop-time
            # collectives are introduced by this offline controller.
            (args.output/'progress.json').write_text(json.dumps(records,indent=2)+'\n')
            print(name,round(report['loop_seconds_max'],3),'seconds',flush=True)
    comparisons=compare_runs(args.output,records)
    timing=[]
    keys=sorted({(r['points'],r['steps'],r['stride'],r['mode']) for r in records})
    for points,steps,stride,mode in keys:
        seconds=[r['loop_seconds'] for r in records if
                 (r['points'],r['steps'],r['stride'],r['mode'])==(points,steps,stride,mode)]
        timing.append(dict(points=points,steps=steps,stride=stride,mode=mode,
            repeats=len(seconds),median_seconds=statistics.median(seconds),
            minimum_seconds=min(seconds),maximum_seconds=max(seconds)))
    summary=dict(schema='salmonpy.mixed-performance-study/2',study=args.study,runs=records,
        comparisons=comparisons,timing=timing,
        scope='single Mac, local ranks; two repeats do not establish statistical significance',
        communication_scope='reachable Python calls, including native backup votes; excludes original backend-private/Fortran MPI; blocking times include waits')
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('mixed mode/stride/repeat comparisons PASS',flush=True)


if __name__=='__main__':main()
