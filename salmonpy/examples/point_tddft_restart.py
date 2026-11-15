"""Resume a point-bound TDDFT checkpoint after rebuilding its configured model.

Requires a prepared local saved-GS catalog, managed runtime and A/dt waveform NPZ.
Creates a new output directory. One thread/worker, maximum two live workers.
"""
import argparse
import json
from pathlib import Path
import numpy as np


def run(catalog_path, initial_id, runtime, waveform, directory):
    from salmonpy.config import load_materials
    from salmonpy.builtin_materials import builtin_material_registry
    from salmonpy.initial_states import tddft_from_catalog
    from salmonpy.runtime import load_runtime
    from salmonpy.point_runs import PointRun
    from salmonpy.errors import SnapshotError
    catalog_path = Path(catalog_path).resolve()
    with np.load(waveform, allow_pickle=False) as data:
        A = np.asarray(data['A'])
        raw_dt = np.asarray(data['dt'])
        if A.dtype.kind not in 'iuf' or raw_dt.dtype.kind not in 'iuf' or raw_dt.shape != ():
            raise ValueError('real numeric A and scalar dt required')
        A = A.astype(float)
        dt = float(raw_dt)
    if A.ndim != 2 or A.shape[1] != 3 or len(A) < 3 or not np.isfinite(A).all() or np.any(A[0]) or not np.isfinite(dt) or dt <= 0:
        raise ValueError('finite A (steps+1,3), A[0]=0 and positive dt required')
    engine = load_runtime(runtime, required=['tddft', 'single-rank'])
    def model():
        catalog = load_materials(catalog_path, registry=builtin_material_registry())
        return tddft_from_catalog(catalog, initial_id, base_directory=catalog_path.parent,
                                  executable=engine, max_steps=len(A)-1, threads=1)
    root = Path(directory).resolve()
    root.mkdir(parents=True, exist_ok=False)
    report = {'status':'running','scope':'TDDFT point checkpoint interface regression',
              'threads_per_worker':1,'max_simultaneous_workers':2,'memory_enforced':False}
    split = (len(A)-1)//2
    try:
        with model() as initial, PointRun({'p':initial,'other':initial},run_id='resume-example') as calculation:
            state = calculation.new_state('p',dt=dt)
            currents = [state.current.copy()]
            for field in A[1:split+1]: currents.append(state.step(field).copy())
            state.save_checkpoint(root/'checkpoint')
            with calculation.new_state('other',dt=dt) as other:
                before = other.current.copy()
                try: other.restore_checkpoint(root/'checkpoint')
                except SnapshotError: pass
                else: raise AssertionError('cross-point restore accepted')
                np.testing.assert_array_equal(before,other.current)
                assert other.time == 0
            for field in A[split+1:]: currents.append(state.step(field).copy())
        # Different initial model object; same explicit run/point/data identity.
        with model() as initial, PointRun({'p':initial},run_id='resume-example') as calculation:
            state = calculation.new_state('p',dt=dt)
            state.restore_checkpoint(root/'checkpoint')
            np.testing.assert_array_equal(state.current,currents[split])
            resumed = [state.current.copy()]
            for field in A[split+1:]: resumed.append(state.step(field).copy())
            np.testing.assert_array_equal(resumed,currents[split:])
            report.update(status='passed',identity=state.record(),split_step=split,
                          resumed_steps=len(resumed)-1,current_max_abs_difference=0.,
                          rejected_other_point_unchanged=True)
        np.savez(root/'comparison.npz',time=np.arange(len(A))*dt,current=currents,resumed=resumed)
    except BaseException as exc:
        report.update(status='failed',error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        (root/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('catalog','initial_id','runtime','waveform','directory'):p.add_argument(name)
    a=p.parse_args()
    print(json.dumps(run(a.catalog,a.initial_id,a.runtime,a.waveform,a.directory),indent=2))


if __name__ == '__main__':main()
