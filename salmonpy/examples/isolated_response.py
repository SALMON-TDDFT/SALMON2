"""Explicit validated E runtime + native GS bundle -> named dipole archive.

Arguments: EXECUTABLE VALIDATION_JSON BUNDLE INPUT_FILE E_NPZ NEW_OUTPUT
E_NPZ contains E of shape (steps,3), sampled at dt,2dt,..., and scalar dt.
"""
import argparse
import numpy as np
from salmonpy.isolated_tddft import load_isolated_runtime, IsolatedTDDFTInitialState
from salmonpy.outputs import save_arrays


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('executable','validation','bundle','input_file','E_npz','output'):
        parser.add_argument(key)
    parser.add_argument('--checkpoint',help='new replay checkpoint directory (separate from output)')
    parser.add_argument('--direct-checkpoint',help='new direct numerical checkpoint directory')
    parser.add_argument('--resume-direct',help='direct numerical checkpoint to import before continuing')
    parser.add_argument('--stop-after',type=int,help='accepted step at which to stop and save')
    a=parser.parse_args()
    with np.load(a.E_npz,allow_pickle=False) as data:
        E=data['E'].copy();raw_dt=data['dt']
        if E.dtype.kind not in 'iuf' or E.ndim!=2 or E.shape[1]!=3 or len(E)<1 or not np.isfinite(E).all():
            raise ValueError('finite real E with shape (steps,3) required')
        if raw_dt.shape!=() or raw_dt.dtype.kind not in 'iuf' or not np.isfinite(raw_dt) or raw_dt<=0:
            raise ValueError('positive scalar atomic dt required')
        dt=float(raw_dt)
    runtime=load_isolated_runtime(a.executable,a.validation)
    with IsolatedTDDFTInitialState.from_salmon(a.bundle,a.input_file,runtime=runtime,
                                              max_steps=len(E),threads=1) as model:
        with model.new_state(dt=dt) as state:
            if a.resume_direct:
                state.restore_direct_checkpoint(a.resume_direct)
            start=state.worker_statistics()['step']
            stop=len(E) if a.stop_after is None else a.stop_after
            if not start <= stop <= len(E):
                raise ValueError('stop-after must be between restored step and E sample count')
            frames=[state.observables]
            for value in E[start:stop]:frames.append(state.step(value))
            if a.direct_checkpoint:state.save_direct_checkpoint(a.direct_checkpoint)
            if a.checkpoint: state.save_checkpoint(a.checkpoint)
    save_arrays(a.output,{'time':np.array([f.time for f in frames]),
                         'E':np.vstack((np.zeros(3) if start==0 else E[start-1],E[start:stop])),
                         'dipole':np.array([f.dipole for f in frames]),
                         'electronic_change':np.array([f.electronic_change for f in frames])},
                metadata={'backend':'isolated-PZ-E3','gauge':'length',
                          'units':{'time':'atomic_time','E':'E:au','dipole':'e*bohr','electronic_change':'e*bohr'},
                          'executable_sha256':runtime.executable.sha256,
                          'convergence':'caller responsibility; not checked by this sample'})
    print('Saved isolated dipole response:',a.output)


if __name__=='__main__':main()
