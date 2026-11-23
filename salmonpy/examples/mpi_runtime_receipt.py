"""Caller-launched MPI runtime verification only; never execute SALMON.

Preview is default and does not import/initialize MPI. With --run, one rank
fully validates and broadcasts an ephemeral single-job context. Every rank then
loads through the receipt; executable bytes remain hashed on every load.
"""
import argparse
import json
from pathlib import Path
import secrets
from time import perf_counter


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('runtime',type=Path)
    p.add_argument('receipt',type=Path,help='NEW receipt path outside the runtime bundle')
    p.add_argument('--job-id',required=True,help='current allocation/job ID; use a fresh local label for a local test')
    p.add_argument('--report',type=Path,help='NEW optional report JSON, without the secret token')
    p.add_argument('--loads',type=int,default=1)
    p.add_argument('--require',nargs='*',default=[])
    p.add_argument('--run',action='store_true',help='perform MPI verification only, no solver or launcher')
    a=p.parse_args()
    if a.loads<1:p.error('--loads must be positive')
    if not a.run:
        print(json.dumps(dict(mode='preview_only',starts_solver=False,initializes_mpi=False,
            runtime=str(a.runtime),receipt=str(a.receipt),job_id=a.job_id,loads=a.loads,
            message='Launch explicitly with the site MPI launcher and --run inside allocated resources.'),indent=2))
        return
    from mpi4py import MPI
    from salmonpy.runtime import prepare_runtime_job, load_runtime
    comm=MPI.COMM_WORLD
    outcome=None
    if comm.rank==0:
        try:
            if a.report is not None and a.report.exists():raise FileExistsError(a.report)
            started=perf_counter()
            validation=prepare_runtime_job(a.runtime,a.receipt,job_id=a.job_id,
                job_token=secrets.token_urlsafe(32),required=a.require)
            outcome=(validation,perf_counter()-started,None)
        except Exception as exc:
            outcome=(None,None,type(exc).__name__+': '+str(exc))
    validation,initial_seconds,error=comm.bcast(outcome,root=0)
    if error:raise RuntimeError(error)
    local=[]
    try:
        for _ in range(a.loads):
            started=perf_counter();engine=load_runtime(a.runtime,required=a.require,validation=validation)
            local.append(dict(rank=comm.rank,seconds=perf_counter()-started,
                executable_sha256=engine.sha256,verification=engine.runtime_validation_report))
        result=dict(rank=comm.rank,loads=local,error=None)
    except Exception as exc:
        result=dict(rank=comm.rank,loads=local,error=type(exc).__name__+': '+str(exc))
    # This script ends after validation: only the reporting rank needs the rows.
    # Successful completion of this gather also proves caller MPI is still usable.
    results=comm.gather(result,root=0)
    if comm.rank!=0:
        if result['error'] is not None:raise RuntimeError(result['error'])
        return
    alive=len(results)==comm.size
    report=dict(schema='salmonpy.runtime-job-example/1',job_id=a.job_id,ranks=comm.size,
        initial_full_validation_seconds=initial_seconds,caller_mpi_alive=alive,
        starts_solver=False,results=results,global_collectives={'setup_bcast':1,'final_gather':1,'allgather':0,'allreduce':0,'barrier':0},scope='validation I/O and host MPI only; no physical observable or multi-node claim')
    if comm.rank==0:
        if a.report is not None:
            a.report.parent.mkdir(parents=True,exist_ok=True)
            with a.report.open('x',encoding='utf-8') as f:json.dump(report,f,indent=2);f.write('\n')
        print(json.dumps(report,indent=2))
    if any(r['error'] is not None for r in results):raise RuntimeError('one or more ranks could not validate the runtime; see report')


if __name__=='__main__':main()
