"""Explicit local mixed-solver scaling campaign; installed wheel required.

BUILD RUNTIME GS NEW_OUTPUT --rank-budget N --memory-budget-gib G
Uses at most15 ranks and one thread/rank. No compilation/download/remote CI.
Every run is checked against the same-shape first run. Timing repeats are
regression measurements, not a production scaling or physical accuracy claim.
"""
import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('build','runtime','bundle','output'):p.add_argument(key,type=Path)
    p.add_argument('--python',default=sys.executable)
    p.add_argument('--launcher',default=shutil.which('mpiexec'))
    p.add_argument('--launcher-arg',action='append',default=[])
    p.add_argument('--field-shape',type=int,nargs=3,default=(96,96,96))
    p.add_argument('--steps',type=int,default=128)
    p.add_argument('--repeats',type=int,default=2)
    p.add_argument('--rank-budget',type=int,required=True)
    p.add_argument('--memory-budget-gib',type=float,required=True)
    a=p.parse_args()
    if a.launcher is None:p.error('explicit MPI launcher required')
    if not math.isfinite(a.memory_budget_gib) or a.memory_budget_gib<=0:p.error('positive finite advisory memory budget required')
    if a.rank_budget<4 or a.repeats<1 or a.steps<1:p.error('rank-budget>=4, positive repeats/steps required')
    if any(n<40 or n%2 for n in a.field_shape):p.error('even field-shape >=40 required')
    plans=[shape for shape in ((1,1,1),(1,2,1),(2,2,1),(2,2,2),(3,2,2))
           if shape[0]*shape[1]*shape[2]+3<=a.rank_budget]
    if any(n%d for shape in plans for n,d in zip(a.field_shape,shape)):
        p.error('field-shape must be divisible by every selected field decomposition')
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    paths=[str(getattr(a,key).resolve()) for key in ('build','runtime','bundle')]
    env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1'}
    summary={'schema':'salmonpy.scaling-smoke/1','field_shape':a.field_shape,'steps':a.steps,
             'repeats':a.repeats,'rank_budget':a.rank_budget,'memory_budget_gib':a.memory_budget_gib,
             'timing_scope':'coupled loop excluding setup/output/serial replay; MPI wait included',
             'rss_scope':'sum of individual process peak RSS, not simultaneous resident total',
             'runs':[]}
    first=None
    for repeat in range(a.repeats):
        for shape in plans:
            count=shape[0]*shape[1]*shape[2];ranks=count+3
            case=out/f'repeat-{repeat}-field-{count}'
            command=[a.launcher,*a.launcher_arg,'-n',str(ranks),a.python,'-I','-m','salmonpy.examples.mpi_mixed_native',
                     *paths,str(case),'--field-grid',*map(str,shape),'--field-shape',*map(str,a.field_shape),
                     '--steps',str(a.steps),'--memory-budget-gib',str(a.memory_budget_gib),
                     '--rank-budget',str(a.rank_budget),'--run-id','scaling-fixture']
            if first is not None:command+=['--compare',str(first)]
            print('Running',case.name,'ranks',ranks,flush=True)
            start=time.perf_counter()
            with (out/(case.name+'.log')).open('w') as log:
                process=subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT)
            row={'name':case.name,'field_grid':shape,'field_ranks':count,'total_ranks':ranks,
                 'command':command,'returncode':process.returncode,'end_to_end_seconds':time.perf_counter()-start}
            summary['runs'].append(row)
            if process.returncode:
                summary['status']='calculation_failed'
                (out/'scaling.json').write_text(json.dumps(summary,indent=2)+'\n')
                raise SystemExit(process.returncode)
            result=json.loads((case/'result.json').read_text())
            row.update(status=result['status'],coupled_loop_seconds=result.get('coupled_loop_seconds_max'),
                       decomposition_comparison=result.get('decomposition_comparison'),communication=result['communication'])
            try:
                memory=json.loads((case/'memory-advice.json').read_text())
                end=next(r for r in reversed(memory['reports']) if r['label']=='mixed_segment_end')
                values=[r['live_bytes'] for r in end['ranks']]
                row['sum_rank_peak_rss_bytes']=None if any(v is None for v in values) else sum(values)
                row['memory_unknown']=end['unknown']
            except Exception as exc:
                row['memory_measurement_error']=str(exc)
            if first is None:first=case
            summary['status']='running'
            (out/'scaling.json').write_text(json.dumps(summary,indent=2)+'\n')
    summary['status']='passed'
    (out/'scaling.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('Scaling regression PASS',out/'scaling.json',flush=True)


if __name__=='__main__':main()
