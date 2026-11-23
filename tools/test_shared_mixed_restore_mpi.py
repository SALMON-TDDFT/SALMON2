"""Ordinary mixed sample with a handled late shared-material restore failure.

Pass normal mpi_mixed_native arguments including --direct-resume. Optional
--kernel LIB injects an explicitly built LD/SBE kernel for a source checkout.
No propagation or numerical kernels are mocked.
"""
import argparse
import ctypes as c
import json
from pathlib import Path
import sys
import numpy as np
from mpi4py import MPI


def main():
    parser=argparse.ArgumentParser(add_help=False)
    parser.add_argument('--kernel')
    parser.add_argument('--material-key',default='hse-point',help='binding with at least two points to fault late')
    args,sample_args=parser.parse_known_args()
    if args.kernel:
        from salmonpy import ld,sbe
        from salmonpy._material_handles import LDAPI as API
        from salmonpy._material_handles import SBEAPI
        a,b=API(args.kernel),SBEAPI(args.kernel)
        ld._handles_api=lambda:a;ld._native=lambda:a.legacy
        sbe._handles_api=lambda:b;sbe._native=lambda:b.legacy
    from salmonpy.examples import mpi_mixed_native as sample
    from salmonpy.examples import _direct_checkpoint as coordinator
    from salmonpy._native_checkpoint import NativeRestoreCandidate
    from salmonpy.checkpointing import AtomicRestoreError
    original=coordinator.checkpoint
    world=MPI.COMM_WORLD
    ordinal=0

    def freeze(coupled,out,tag):
        nonlocal ordinal
        ordinal+=1
        result={'coupled':(coupled.step_index,coupled.J.tobytes(),coupled.response._A.tobytes(),coupled.response.failed)}
        owners=[]
        if coupled.field is not None:owners.append(('field',coupled.field,'sp_eh_checkpoint'))
        for point,state in coupled.response.states.items():
            if hasattr(state,'native'):owners.append((str(point),state,'sp_tddft_checkpoint'))
            elif getattr(state,'checkpoint_protocol',None)=='salmonpy.direct-state/1':
                owners.append((str(point),state,'sp_fast_hse_checkpoint'))
            else:
                token=state._raw.snapshot()
                result[str(point)]=(id(state),token._integrity,token._data)
        for key,owner,symbol in owners:
            fn=getattr(owner.lib,symbol);fn.argtypes=[c.c_char_p,c.c_int];fn.restype=c.c_int
            path=out/f'preserved-{tag}-{world.rank}-{key}.bin'
            assert owner._call(fn,str(path).encode(),0)==0
            result[key]=(id(owner),path.read_bytes())
            for name in ('A','E','current','steps','step_index','energy','failed','closed'):
                if hasattr(owner,name):
                    value=getattr(owner,name)
                    result[key+name]=(id(value),value.tobytes()) if isinstance(value,np.ndarray) else value
            root=getattr(owner,'directory',None) or owner.work
            result[key+'files']={p.name:p.read_bytes() for p in Path(root).glob('*.data')}
        return result

    def checkpoint(world,coupled,root,signature,**kwargs):
        if not kwargs.get('restore'):return original(world,coupled,root,signature,**kwargs)
        output=Path(sample_args[3]).resolve()
        before=freeze(coupled,output,'before')
        hse=[p for p,b in coupled.response.bindings.items() if b.key==args.material_key]
        assert len(hse)>=2
        last=max(hse);stage=NativeRestoreCandidate.stage
        fired=False
        def fail_late(candidate):
            nonlocal fired
            value=stage(candidate)
            state=candidate.owner
            native=getattr(state,'native',None) or getattr(state,'_native',None)
            if getattr(state,'point_id',None)==last and native.context.solver_comm.rank==native.context.solver_comm.size-1:
                fired=True
                raise RuntimeError('injected after last material candidate staged')
            return value
        NativeRestoreCandidate.stage=fail_late
        try:
            try:original(world,coupled,root,signature,**kwargs)
            except AtomicRestoreError as exc:
                assert exc.original_preserved and not exc.cleanup_errors
            else:raise AssertionError('late restore failure was not reported')
        finally:NativeRestoreCandidate.stage=stage
        assert any(world.allgather(fired))
        after=freeze(coupled,output,'after')
        assert after==before,('original changed',world.rank,[k for k in before if before[k]!=after[k]])
        result=original(world,coupled,root,signature,**kwargs)
        if world.rank==0:
            (output/'shared-restore-validation.json').write_text(json.dumps({'status':'pass',
                'late_point':last,'material_key':args.material_key,'original_states_and_outputs_preserved':True,'retry_succeeded':True,
                'replayed_steps':0},indent=2)+'\n')
        return result
    coordinator.checkpoint=checkpoint
    sys.argv=[sys.argv[0]]+sample_args
    sample.main()


if __name__=='__main__':main()
