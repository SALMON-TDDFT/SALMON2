"""Saved spatial FDTD E/H -> molecule, with no current returned to FDTD.

Uses an existing GS bundle and validated molecular worker. No SCF, download,
build, self-vacuum reference or radiation calculation is performed.
"""
import argparse
import numpy as np
from salmonpy.near_field import GridFieldSeries
from salmonpy.isolated_tddft import load_isolated_runtime, IsolatedTDDFTInitialState
from salmonpy.outputs import save_arrays


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('executable','validation','bundle','input_file','grid_npz','output'):
        parser.add_argument(name)
    parser.add_argument('--direct-checkpoint',help='new direct checkpoint directory')
    parser.add_argument('--resume-direct',help='direct checkpoint to restore without replay')
    parser.add_argument('--stop-after',type=int,help='last molecular step to compute')
    args=parser.parse_args()
    with np.load(args.grid_npz,allow_pickle=False) as data:
        names=('origin','spacing','times','E','H','molecule_position','center')
        fields=GridFieldSeries(**{name:data[name] for name in names},
            time_origin=float(data['time_origin']) if 'time_origin' in data else 0.)
    dt=float(fields.times[1]);steps=len(fields.times)-2
    # Validate the complete time grid before creating a worker.
    fields.window(step_index=1,dt=dt)
    stop=steps if args.stop_after is None else args.stop_after
    if not 0<=stop<=steps:parser.error('stop-after must be between zero and field frame count minus two')
    runtime=load_isolated_runtime(args.executable,args.validation)
    with IsolatedTDDFTInitialState.from_salmon(args.bundle,args.input_file,runtime=runtime,
            max_steps=steps,threads=1,near_field=True) as model:
        with model.new_state(dt=dt,position=fields.molecule_position) as state:
            if args.resume_direct:state.restore_direct_checkpoint(args.resume_direct)
            restore_statistics=state.worker_statistics()
            start=restore_statistics['step']
            if stop<start:raise ValueError('stop-after precedes restored molecular step')
            frames=[state.observables]
            for n in range(start+1,stop+1):
                frames.append(state.step_grid_field(fields.window(step_index=n,dt=dt),step_index=n))
            if args.direct_checkpoint:state.save_direct_checkpoint(args.direct_checkpoint)
            final_statistics=state.worker_statistics()
    arrays={
        'molecular_step':np.arange(start,stop+1,dtype=np.int64),
        'time':np.array([frame.time for frame in frames]),
        'dipole':np.array([frame.dipole for frame in frames]),
        'electron_number':np.array([frame.electron_number for frame in frames]),
        'electronic_change':np.array([frame.electronic_change for frame in frames])}
    energy_available=all(frame.internal_energy is not None for frame in frames)
    if energy_available:arrays['internal_energy']=np.array([frame.internal_energy for frame in frames])
    save_arrays(args.output,arrays,metadata={
        'internal_energy_available':energy_available,
        'coupling':'FDTD_to_molecule_only','spatial_EH':True,'feedback':False,
        'restore_statistics':restore_statistics,'final_statistics':final_statistics,
        'source_time_origin':fields.time_origin,'source_frames_per_molecular_step':'n,n+1',
        'executable_sha256':runtime.executable.sha256,
        'units':{'time':'atomic_time','dipole':'e*bohr','electronic_change':'e*bohr',
                 'electron_number':'electrons','internal_energy':'hartree'},
        'convergence':'not checked; weak molecular back-action assumed'})
    print('Saved one-way molecular response:',args.output)


if __name__=='__main__':main()
