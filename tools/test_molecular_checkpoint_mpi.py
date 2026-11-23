"""Fault/retry validation for the public molecular direct checkpoint sample.

Pass the same arguments as salmonpy.examples.molecular_feedback, including
--resume EXISTING_CHECKPOINT and a NEW_OUTPUT directory. Uses the caller's
existing GS and runtime. No SCF, build, download or CI dispatch.
"""
from pathlib import Path
import json
import os
import runpy
import shutil
import sys
from unittest.mock import patch
import numpy as np
from salmonpy.dipole_source import MolecularFeedback
from salmonpy import dipole_source
from salmonpy.isolated_tddft import _EEngine
from salmonpy.checkpointing import AtomicRestoreError

original_restore=MolecularFeedback.restore_checkpoint


def checked_restore(session,directory):
    comm=session.comm
    output=Path(sys.argv[5]);source=Path(directory)
    field=session.field
    # Deliberately exercise rejection against a live nonzero target, outside restore.
    session.start()
    for _ in range(4):session.step()
    was_started=session.started
    accepted=session.step_index
    before_E=None if field is None else field.E.copy()
    before_A=None if field is None else field.A.copy()
    original_stats=None if session.state is None else session.state.worker_statistics()
    checkpoint_fn=None
    if field is not None:
        from salmonpy.api import _bind_field_operations
        _bind_field_operations(field.lib,('sp_eh_checkpoint',))
        checkpoint_fn=field.lib.sp_eh_checkpoint
        before_file=output/('before-'+str(comm.rank)+'.bin')
        assert field._call(checkpoint_fn,os.fsencode(before_file),0)==0
        original_call=field._call
        def no_evolution(fn,*args):
            if getattr(fn,'__name__','') in ('sp_eh_step','sp_eh_step_owned','sp_eh_step_checked'):
                raise AssertionError('restore must not advance a field')
            return original_call(fn,*args)
    else:no_evolution=None
    if session.state is not None:session.state.save_direct_checkpoint(output/'before-molecule')
    def unchanged(tag):
        assert session.step_index==accepted and session.started==was_started and not session.failed
        if field is not None:
            assert session._field_step_indices()==(accepted,-1)
            np.testing.assert_array_equal(field.A,before_A)
            np.testing.assert_array_equal(field.E,before_E)
            path=output/(tag+'-'+str(comm.rank)+'.bin')
            assert field._call(checkpoint_fn,os.fsencode(path),0)==0
            assert path.read_bytes()==before_file.read_bytes()
        if session.state is not None:
            assert session.state.worker_statistics()==original_stats
            session.state.save_direct_checkpoint(output/(tag+'-molecule'))
            assert (output/(tag+'-molecule/isolated-state.bin')).read_bytes()==(output/'before-molecule/isolated-state.bin').read_bytes()
    # The hook runs after native stage/query, so no owning-group collective is skipped.
    stream_init=dipole_source._MolecularStreams.__init__
    def fail_field(owner,*args):
        if comm.rank==session.reference_root:raise OSError('injected reference stream preparation failure')
        stream_init(owner,*args)
    bridge_init=dipole_source._IsolatedDipoleTrial.__init__
    def fail_molecule(owner,*args,**kwargs):
        if comm.rank==session.total_root:raise OSError('injected molecule binding preparation failure')
        bridge_init(owner,*args,**kwargs)
    reports=[]
    from contextlib import ExitStack
    with ExitStack() as stack:
        if field is not None:stack.enter_context(patch.object(field,'_call',side_effect=no_evolution))
        stack.enter_context(patch.object(_EEngine,'step',side_effect=AssertionError('restore must not replay electrons')))
        for tag,target,hook in [('field-failure',dipole_source._MolecularStreams,fail_field),
                                ('molecule-failure',dipole_source._IsolatedDipoleTrial,fail_molecule)]:
            try:
                with patch.object(target,'__init__',hook):original_restore(session,source)
            except AtomicRestoreError as exc:assert exc.original_preserved
            else:raise AssertionError('injected preparation failure was accepted')
            unchanged(tag);reports.append(tag)
        bad=output/'wrong-step-checkpoint'
        error=None
        if comm.rank==session.total_root:
            try:
                shutil.copytree(source,bad,copy_function=os.link)
                data=json.loads((bad/'complete.json').read_text())
                data['field_step']+=1;data['molecular_step']+=1;data['history']=None
                for row in data['rows']:
                    row['session_step']+=1
                    if row['field_step'] is not None:row['field_step']+=1
                    if row['molecular_step'] is not None:row['molecular_step']+=1
                (bad/'complete.json').unlink()
                (bad/'complete.json').write_text(json.dumps(data))
            except Exception as exc:error=str(exc)
        error=comm.bcast(error,root=session.total_root)
        if error:raise RuntimeError(error)
        try:original_restore(session,bad)
        except AtomicRestoreError as exc:assert exc.original_preserved
        else:raise AssertionError('manifest/native step mismatch was accepted')
        unchanged('step-failure');reports.append('step-failure')
        result=original_restore(session,source)
    local=dict(rank=comm.rank,target_before_field_step=accepted,cases=reports,all_original_numerical_bytes_preserved=True,
               no_field_advance=True,no_electron_step=True)
    evidence=comm.gather(local,root=session.total_root)
    if comm.rank==session.total_root:
        (output/'restore-validation.json').write_text(json.dumps(dict(ranks=evidence,
            verified_steps=result['verified_steps'],field_step=result['field_step'],
            molecular_step=result['molecular_step']),indent=2)+'\n')
    return result


if __name__=='__main__':
    with patch.object(MolecularFeedback,'restore_checkpoint',checked_restore):
        runpy.run_module('salmonpy.examples.molecular_feedback',run_name='__main__')
