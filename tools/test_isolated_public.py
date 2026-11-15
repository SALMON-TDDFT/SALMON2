"""Exercise public isolated loader/replay against validated E-worker history.

Args: EXECUTABLE VALIDATION_DIR NEW_OUTPUT
VALIDATION_DIR is produced by tools/test_isolated_worker.py; no new SCF/physics
convergence calculation. At most two simultaneous one-thread children.
"""
import json
from pathlib import Path
import sys
import numpy as np
from salmonpy.isolated_tddft import load_isolated_runtime, IsolatedTDDFTInitialState
from salmonpy.field_input import FieldPacket, FieldValue
from salmonpy.errors import RuntimeFailure


def main():
    executable, validation, output = sys.argv[1:]
    validation, output = Path(validation), Path(output)
    output.mkdir(parents=True,exist_ok=False)
    runtime = load_isolated_runtime(executable,validation/'report.json')
    with np.load(validation/'dipoles.npz',allow_pickle=False) as data:
        E, dt, reference = data['fields'].copy(),float(data['dt']),data['values'].copy()
    with IsolatedTDDFTInitialState.from_salmon(validation/'first','input.inp',runtime=runtime,
                                              max_steps=25,threads=1,timeout=30) as model:
        with model.new_state(dt=dt) as a:
            with model.new_state(dt=dt) as b:
                np.testing.assert_array_equal(a.observables.dipole,reference[0,3:])
                for n,v in enumerate(E[:10],1):
                    before=b.observables.dipole.copy()
                    frame=a.step_fields(FieldPacket(E=FieldValue(v,time=n*dt,position=(0,0,0),unit='E:au')))
                    np.testing.assert_array_equal(b.observables.dipole,before)
                    np.testing.assert_array_equal(frame.dipole,reference[n,3:])
                    np.testing.assert_array_equal(b.step(v).dipole,frame.dipole)
            saved=a.snapshot()
            a.save_checkpoint(output/'checkpoint')
            expected=np.array([a.step(v).dipole for v in E[10:]])
            model.close()
            a.restore(saved)
            actual=np.array([a.step(v).dipole for v in E[10:]])
            np.testing.assert_array_equal(actual,expected)
            # Public step recovers accepted history after a lost child.
            before=a.observables.dipole.copy()
            a._engine._proc.kill();a._engine._proc.wait()
            try:a.step([0,0,0])
            except RuntimeFailure:pass
            else:raise AssertionError('killed worker did not report failed step')
            np.testing.assert_array_equal(a.observables.dipole,before)
            assert a.time==20*dt
            a.step([0,0,0])
    with IsolatedTDDFTInitialState.from_salmon(validation/'first','input.inp',runtime=runtime,
                                              max_steps=25,threads=1,timeout=30) as rebuilt:
        with rebuilt.new_state(dt=dt) as resumed:
            resumed.restore_checkpoint(output/'checkpoint')
            disk_actual=np.array([resumed.step(v).dipole for v in E[10:]])
            np.testing.assert_array_equal(disk_actual,expected)
    np.savez(output/'continuation.npz',expected=expected,actual=actual,disk_actual=disk_actual)
    (output/'report.json').write_text(json.dumps(dict(public_loader=True,independent_states=True,
        disk_rebuilt_model=True,disk_continuation_max_abs_error=0.,
        field_packet=True,restore_continuation_max_abs_error=0.,worker_failure_recovery=True,
        runtime_sha256=runtime.executable.sha256,scope='interface regression; no convergence study'),indent=2)+'\n')
    print('public isolated loader, replay, independent states and recovery PASS')


if __name__=='__main__':main()
