"""Exercise the ordinary mpi_mixed_native sample's protected restore.

Pass the ordinary sample arguments, with --direct-save-at N and --steps>N+3.
Optional --kernel LIB explicitly selects the built LD/SBE kernel for source tests.
The test injects a failure after one material rank stages a valid native payload.
No sample flags or native physics are mocked. Use mpiexec with a timeout.
"""
import argparse
import ctypes as c
import json
from pathlib import Path
import struct
import shutil
import sys
from unittest.mock import patch
from contextlib import nullcontext
import numpy as np
from mpi4py import MPI


def main():
    parser=argparse.ArgumentParser(add_help=False)
    parser.add_argument('--kernel')
    parser.add_argument('--simulate-memory-probe-failure',action='store_true')
    args, sample_args=parser.parse_known_args()
    if args.kernel:
        from salmonpy import ld,sbe
        from salmonpy._ld_handles import API
        from salmonpy._sbe_handles import API as SBEAPI
        api=API(args.kernel);ld._handles_api=lambda:api;ld._native=lambda:api.legacy
        sbeapi=SBEAPI(args.kernel);sbe._handles_api=lambda:sbeapi;sbe._native=lambda:sbeapi.legacy
    from salmonpy.examples import mpi_mixed_native as sample
    from salmonpy.examples import _direct_checkpoint as coordinator
    from salmonpy._native_checkpoint import NativeRestoreCandidate
    from salmonpy.checkpointing import AtomicRestoreError
    original=coordinator.checkpoint
    world=MPI.COMM_WORLD
    completed=[]

    def freeze(coupled,root,label):
        r=coupled.response;field=coupled.field
        result={'step':coupled.step_index,'response_step':r.step_index,'failed':r.failed,
                'J':coupled.J.tobytes(),'A':r._A.tobytes(),
                'states':tuple((p,id(s),id(getattr(s,'_raw',None))) for p,s in r.states.items()),
                'adapters':tuple((p,id(s),s.A.tobytes()) for p,s in r._adapters.items())}
        native=[]
        if field is not None:
            native.append((field,'sp_eh_checkpoint'))
            result.update(field_A=field.A.tobytes(),field_E=field.E.tobytes(),field_step=field.step_index,field_failed=field.failed)
        for point,state in r.states.items():
            if hasattr(state,'native'):
                native.append((state,'sp_tddft_checkpoint'))
                result[str(point)]=(state.steps,state.A.tobytes())
            elif getattr(state,'checkpoint_protocol',None)=='salmonpy.direct-state/1':
                native.append((state,'sp_fast_hse_checkpoint'))
                result[str(point)]=(state.steps,state.A.tobytes(),state.current.tobytes(),state.failed)
            else:
                token=state._raw.snapshot()
                result[str(point)]=(state._steps,token._integrity,token._data,token.time,
                                    id(token._owner),tuple(id(s) for s in state._prepared.live))
        for index,(state,symbol) in enumerate(native):
            fn=getattr(state.lib,symbol);fn.argtypes=[c.c_char_p,c.c_int];fn.restype=c.c_int
            path=root/f'{label}-{world.rank}-{index}.bin'
            assert state._call(fn,str(path).encode(),0)==0
            value=path.read_bytes()
            result[symbol]=value
            if symbol=='sp_eh_checkpoint':
                meta=struct.unpack_from('=3i',value,24+18*4+2*4)
                assert meta[0]==1 and meta[1]==2, 'test requires both field-local LD materials'
        return result

    def same_frame(left,right):
        if left is None:
            assert right is None
            return
        assert left.keys()==right.keys()
        for key in left:
            np.testing.assert_array_equal(left[key],right[key],err_msg=key)

    def exercise(world,coupled,directory,signature,*,restore=False,frames=None,**kwargs):
        assert not restore, 'launch the save path for this live-state test'
        saved=original(world,coupled,directory,signature,frames=frames)
        assert signature['steps']>saved+3
        expected=[coupled.step() for unused in range(2)]
        root=Path(directory).parent
        before=freeze(coupled,root,'before-failure')
        called=[]
        def forbidden(candidate):
            called.append(True)
            raise AssertionError('policy rejection must precede stage')
        with patch.object(NativeRestoreCandidate,'stage',forbidden):
            try:
                original(world,coupled,directory,signature,restore=True,protected=world.rank!=world.size-1)
                raise AssertionError('mixed protected policies accepted')
            except ValueError as error:
                assert 'policies differ' in str(error)
        assert not called and not coupled.response.failed
        manifest=json.loads((Path(directory)/'complete.json').read_text())
        identity_cases=[]
        if manifest['schema']=='salmonpy.native-direct/3':
            from copy import deepcopy
            assert manifest['identity']==json.loads((root/'identity.json').read_text())
            for name in ('run','point','material','definition','initial','initial_digest','ld','ownership','bool'):
                changed=deepcopy(manifest)
                record=changed['identity']['point_metadata']['point_run']
                if name=='run':record['run_id']='other-run'
                elif name=='ld':changed['identity']['constitutive']['cell_ids']['1']='other-material'
                elif name=='ownership':changed['identity']['point_metadata']['mpi_point_map'][0]['ranks']=[world.size-1]
                elif name=='bool':changed['identity']['point_metadata']['mpi_point_map'][0]['ranks']=[True]
                else:
                    key={'point':'point_id','material':'material_id','definition':'definition_digest',
                         'initial':'initial_state_id','initial_digest':'initial_state_digest'}[name]
                    record['points'][0][key]='e'*64
                mismatch=root/('identity-mismatch-'+name)
                if world.rank==0:
                    mismatch.mkdir()
                    (mismatch/'complete.json').write_text(json.dumps(changed))
                world.Barrier()
                with patch.object(NativeRestoreCandidate,'stage',forbidden):
                    try:
                        original(world,coupled,mismatch,signature,restore=True)
                        raise AssertionError('identity mismatch accepted: '+name)
                    except ValueError as error:assert 'identity' in str(error)
                assert not called and not coupled.response.failed
                identity_cases.append(name)
            old_identity=deepcopy(coupled.constitutive_identity)
            try:
                if world.rank==world.size-1:
                    coupled.constitutive_identity['materials']['drude']['definition_digest']='f'*64
                destination=root/'rank-identity-mismatch-save'
                with patch.object(NativeRestoreCandidate,'stage',forbidden):
                    try:
                        original(world,coupled,destination,signature,frames=frames)
                        raise AssertionError('different rank identity accepted')
                    except ValueError as error:assert 'identity' in str(error)
                assert not destination.exists() and not called
                identity_cases.append('rank-agreement')
            finally:
                coupled.constitutive_identity=old_identity
        bad_history=root/'bad-history-checkpoint'
        if world.rank==0:
            shutil.copytree(directory,bad_history)
            (bad_history/'frames.json').write_text('[]')
            manifest=json.loads((bad_history/'complete.json').read_text())
            manifest['files']['frames.json']=coordinator.digest(bad_history/'frames.json')
            (bad_history/'complete.json').write_text(json.dumps(manifest))
        world.Barrier()
        with patch.object(NativeRestoreCandidate,'stage',forbidden):
            try:
                original(world,coupled,bad_history,signature,restore=True)
                raise AssertionError('invalid history accepted')
            except ValueError as error:
                assert 'history admission' in str(error)
        assert not called and not coupled.response.failed
        stage=NativeRestoreCandidate.stage
        def fail_late(candidate):
            value=stage(candidate)
            if world.rank==world.size-1:
                raise OSError('injected after valid late electronic candidate load')
            return value
        try:
            with patch.object(NativeRestoreCandidate,'stage',fail_late):
                original(world,coupled,directory,signature,restore=True)
            raise AssertionError('injected failure was accepted')
        except AtomicRestoreError as error:
            assert error.original_preserved and not error.cleanup_errors
        after=freeze(coupled,root,'after-failure')
        assert before==after, 'state changed: '+str([key for key in before if before[key]!=after[key]])
        expected.append(coupled.step())  # Original session remains usable.
        count=coupled.step_index
        assert count==saved+3
        assert original(world,coupled,directory,signature,restore=True)==saved
        assert coupled.step_index==saved
        for frame in expected:
            same_frame(coupled.step(),frame)
        assert coupled.step_index==count and not coupled.response.failed
        statuses=world.allgather(True)
        if world.rank==0:
            record={'status':'passed','ranks':world.size,'saved_step':saved,
                'failure_step':saved+2,'continued_after_failure':True,
                'all_original_state_bytes_preserved':True,'field_LD_materials':2,
                'retry_future_steps_exact':3,'mixed_policy_rejected_before_stage':True,'invalid_history_rejected_before_stage':True,'restore_replayed_steps':0,
                'identity_cases_rejected_before_load':identity_cases,
                'memory_probe_failure_injected':args.simulate_memory_probe_failure,
                'memory_budget_bytes':getattr(coupled,'memory_budget_bytes',None),
                'memory_restore_reports':len([r for r in (getattr(coupled,'memory_reports',None) or []) if r['label']=='mixed_protected_restore']),
                'backend_keys':sorted(set(signature['materials'].values())) if 'materials' in signature else []}
            np.savez(root/'protected-continuation.npz',**{key:np.asarray([frame[key] for frame in expected]) for key in expected[0]})
            (root/'protected-restore-validation.json').write_text(json.dumps(record,indent=2)+'\n')
            print('mixed protected restore MPI PASS '+json.dumps(record),flush=True)
        completed.append(True)
        return saved

    sys.argv=['mpi_mixed_native',*sample_args]
    probe=(patch('salmonpy.memory._peak_rss',side_effect=OSError('injected resource probe failure'))
           if args.simulate_memory_probe_failure and world.rank==world.size-1 else nullcontext())
    with probe,patch.object(coordinator,'checkpoint',exercise):
        sample.main()
    assert completed and world.allreduce(1)==world.size


if __name__=='__main__':
    try:main()
    except BaseException:
        import traceback
        traceback.print_exc()
        MPI.COMM_WORLD.Abort(1)
