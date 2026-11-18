"""Advisory estimates for the ordinary mixed fixture; never physics admission."""
import json
import math
import warnings
from pathlib import Path
from .memory import MemoryEstimate, MemoryBudgetWarning, memory_advice, process_memory

RESERVE=128*2**20


def _record(world, local, budget, label, records, path):
    report=memory_advice(local,budget_bytes=budget,comm=world,label=label)
    record=report.to_dict()
    records.append(record)
    error=None
    if world.rank==0 and path is not None:
        try:
            Path(path).write_text(json.dumps({'schema':'salmonpy.memory-log/1','reports':records},indent=2)+'\n')
        except Exception as exc:
            error=type(exc).__name__+': '+str(exc)
    errors=world.allgather(error)
    if any(errors):
        record['report_write_errors']=[e for e in errors if e]
        if world.rank==0:
            try:
                warnings.warn('memory advice file unavailable; continuing: '+str(error),MemoryBudgetWarning,stacklevel=2)
            except Exception as exc:
                record['report_write_warning_error']=type(exc).__name__+': '+str(exc)
    return report


def initial_advice(world,context,calculation,bundle,field_grid,model,budget,records,path,*,kind,protected_steps=False):
    rss,basis=process_memory()
    unknown=('initialization workspaces and future solver allocations are not fully modeled',)
    if protected_steps:
        unknown+=('whole-step and per-component recovery backups coexist; native extra memory is not calibrated',)
    try:
        if context.group.role=='field':
            counts=calculation.geometry.grid.counts
            local=[math.ceil(n/p)+8 for n,p in zip(counts,field_grid)]
            poles=max((len(m.strength) for m in calculation.materials if hasattr(m,'strength')),default=0)
            storage=math.prod(local)*(2048+256*poles)
            basis+='; field heuristic: halo8, (2048+256*max_poles) bytes/cell'
        elif kind=='sbe-point':
            nb,nk=model.eigen.shape
            storage=8*(16*nb*nb*nk+48)
            basis+='; SBE dense state/model allowance'
        else:
            directory=Path(bundle)/'data_for_restart'
            files=list(directory.glob('*wfn*.bin'))
            if not files:
                raise OSError('wavefunction input sizes unavailable')
            storage=8*sum(p.stat().st_size for p in files)
            basis+='; electronic heuristic: 8*whole input wavefunction bytes per rank'
        if context.group.role!='field':
            storage*=len(context.group.point_ids)
            basis+='; local states='+str(len(context.group.point_ids))
        live=None if rss is None else rss+storage
        estimate=MemoryEstimate(live,4*storage,RESERVE,basis,
            unknown)
    except Exception as exc:
        estimate=MemoryEstimate(rss,None,RESERVE,basis,('initial estimate unavailable: '+str(exc),))
    return _record(world,estimate,budget,'mixed_initial_plan',records,path)


def live_advice(world,budget,records,path,*,children=False,label='mixed_initialized'):
    rss,basis=process_memory()
    return _record(world,MemoryEstimate(rss,0,RESERVE,basis,('child worker memory excluded',) if children else ()),budget,label,records,path)


def restore_advice(world,coupled,root,data):
    """Four payload copies plus Python decoding and reserve, above observed RSS.

    This is a heuristic, not a proven bound on the native heap. Actual payload
    parsing/verification belongs to the scientific checkpoint path, not this probe.
    """
    records=getattr(coupled,'memory_reports',None)
    enabled=world.allgather(records is not None)
    if not any(enabled):
        return None
    if records is None:
        records=[]
        coupled.memory_reports=records
    from .examples._direct_checkpoint import _state_filename
    r=coupled.response
    rss,basis=process_memory()
    unknown=[]
    try:
        payload=0;snapshots=0
        if coupled.field is not None:
            payload+=(root/f'field-{world.rank}.bin').stat().st_size
        for point,state in r.states.items():
            if hasattr(state,'_process'):
                unknown.append('child worker memory excluded')
            if hasattr(state,'native'):
                payload+=(root/_state_filename(r,point,world.rank)).stat().st_size
            elif getattr(state,'checkpoint_protocol',None)=='salmonpy.direct-state/1':
                payload+=(root/_state_filename(r,point,world.rank)).stat().st_size
            else:
                raw=state._raw
                if hasattr(raw,'_initial'):
                    nb,nk=raw._initial.eigen.shape
                    snapshots+=8*(16*nb*nb*nk+48)
                else:
                    snapshots+=12*len(raw._material.config.strength)*3*8
        metadata=(root/f'rank-{world.rank}.json').stat().st_size+(root/'frames.json').stat().st_size
        extra=4*payload+32*metadata+snapshots
        estimate=MemoryEstimate(rss,extra,RESERVE,basis+'; extra=4*payload('+str(payload)+')+32*JSON('+str(metadata)+')+snapshots('+str(snapshots)+') bytes',tuple(unknown))
    except Exception as exc:
        estimate=MemoryEstimate(rss,None,RESERVE,basis,tuple(unknown)+('candidate estimate unavailable: '+str(exc),))
    return _record(world,estimate,getattr(coupled,'memory_budget_bytes',None),
                   'mixed_protected_restore',records,getattr(coupled,'memory_report_path',None))
