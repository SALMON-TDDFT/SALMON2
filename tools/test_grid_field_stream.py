"""Validate streamed molecular E/H against an identical prepared history."""
import argparse
import json
from pathlib import Path
import shutil
import numpy as np


def run(executable,bundle,validation,output):
    from salmonpy.near_field import GridFieldSeries
    from salmonpy.isolated_tddft import load_isolated_runtime,IsolatedTDDFTInitialState
    from salmonpy.errors import InputError
    root=Path(output);root.mkdir(parents=True,exist_ok=False);bundle=Path(bundle)
    runtime=load_isolated_runtime(executable,validation)
    dt=.003;steps=4;shape=(steps+2,2,2,2,3)
    E=np.zeros(shape);H=np.zeros(shape)
    for n in range(1,steps+2):
        E[n]=np.array([.001,-.002,.0015])*n
        E[n,:,0,:,2]+=.0001*n
        H[n]=np.array([.0002,-.0001,.00015])*n
    geometry=dict(origin=(-40.,)*3,spacing=(80.,)*3,molecule_position=(.17,-.23,.31),center=(0.,0.,0.))
    prepared=GridFieldSeries(times=np.arange(steps+2)*dt,E=E,H=H,**geometry)
    blank=GridFieldSeries(times=np.arange(steps+2)*dt,E=np.zeros(shape),H=np.zeros(shape),**geometry)
    np.savez(root/'fields.npz',times=prepared.times,E=E,H=H,**geometry)
    np.savez(root/'extra-E.npz',E=np.zeros((steps,3)),dt=dt)
    def window(n,**changes):
        return GridFieldSeries(times=[0.,dt],E=E[n:n+2],H=H[n:n+2],time_origin=n*dt,**dict(geometry,**changes))
    def make(label,series):
        folder=root/label;folder.mkdir()
        for asset in ('C_rps.dat','H_rps.dat','input.inp'):shutil.copy2(bundle/asset,folder/asset)
        shutil.copytree(bundle/'data_for_restart',folder/'data_for_restart');series.write(folder/'salmonpy-near-eh.dat')
        return folder
    control_bundle=make('prepared-bundle',prepared);stream_bundle=make('stream-bundle',blank)
    def model(folder):return IsolatedTDDFTInitialState.from_salmon(folder,'input.inp',runtime=runtime,max_steps=steps,near_field=True)
    values=[];densities=[]
    with model(control_bundle) as fixed_model,model(stream_bundle) as live_model:
        with fixed_model.new_state(dt=dt) as control,live_model.new_state(dt=dt) as live:
            before=live.worker_statistics();initial=live._engine.charge_density()
            try:live._engine._stage_grid_field(window(1,origin=(-41.,-40.,-40.)),1)
            except InputError:pass
            else:raise AssertionError('wrong geometry accepted')
            assert live.worker_statistics()==before
            np.testing.assert_array_equal(live._engine.charge_density().charge_density,initial.charge_density)
            for n in range(1,steps+1):
                live._engine._stage_grid_field(window(n),n)
                if n==1:
                    # Native finite rejection must not replace the last valid pair.
                    path=Path(live._engine._work.name)/'run/salmonpy-grid-frames.bin'
                    raw=bytearray(path.read_bytes());raw[-8:]=np.array([np.inf]).tobytes();path.write_bytes(raw)
                    status,step,updates=live._engine._command('GRID_E',9)
                    assert status!=0 and step==0 and updates==0
                live.step([0,0,0]);control.step([0,0,0])
                np.testing.assert_array_equal(live._engine.current,control._engine.current)
                values.append(live._engine.current.copy());densities.append(live._engine.charge_density().charge_density)
                if n==2:live.save_direct_checkpoint(root/'checkpoint')
        with live_model.new_state(dt=dt) as restored:
            restored.restore_direct_checkpoint(root/'checkpoint');assert restored.worker_statistics()['updates']==0
            np.testing.assert_array_equal(restored._engine.charge_density().charge_density,densities[1])
            for n in range(3,steps+1):
                restored._engine._stage_grid_field(window(n),n);restored.step([0,0,0])
                np.testing.assert_array_equal(restored._engine.current,values[n-1])
                np.testing.assert_array_equal(restored._engine.charge_density().charge_density,densities[n-1])
    from salmonpy.convergence import _comparison
    electric_only=GridFieldSeries(times=prepared.times,E=E,H=np.zeros(shape),**geometry)
    electric_bundle=make('electric-only-bundle',electric_only)
    with model(electric_bundle) as electric_model:
        with electric_model.new_state(dt=dt) as electric:
            electric_values=[]
            for _ in range(steps):
                electric.step([0,0,0]);electric_values.append(electric._engine.current.copy())
    magnetic_effect=_comparison(values,electric_values,0.,0.)
    assert magnetic_effect['absolute_error']>1e-16
    thin_geometry=dict(geometry,origin=(-17.,)*3,spacing=(34.,)*3)
    thin=GridFieldSeries(times=prepared.times,E=np.zeros(shape),H=np.zeros(shape),**thin_geometry)
    thin_bundle=make('compact-stream-bundle',thin)
    thin_control=GridFieldSeries(times=prepared.times,E=E,H=H,**thin_geometry)
    thin_control_bundle=make('compact-control-bundle',thin_control)
    with model(thin_bundle) as thin_model,model(thin_control_bundle) as thin_control_model:
        with thin_model.new_state(dt=dt) as state,thin_control_model.new_state(dt=dt) as control:
            # On this single-rank fixture the compact grid covers all actual
            # halo/projector queries. It must be accepted, not over-restricted.
            for n in range(1,steps+1):
                pair=GridFieldSeries(times=[0.,dt],E=E[n:n+2],H=H[n:n+2],time_origin=n*dt,**thin_geometry)
                state._engine._stage_grid_field(pair,n);state.step([0,0,0]);control.step([0,0,0])
                np.testing.assert_array_equal(state._engine.current,control._engine.current)
    evidence=json.loads(Path(validation).read_text())
    evidence['grid_stream']=dict(schema='salmonpy.grid-stream-validation/1',prescribed_matches=True,
                                direct_restart_exact=True,rejection_preserves_state=True)
    receipt=root/'runtime-validation.json';receipt.write_text(json.dumps(evidence,indent=2)+'\n')
    runtime=load_isolated_runtime(executable,receipt)
    with model(stream_bundle) as public_model:
        with public_model.new_state(dt=dt) as public:
            for n in range(1,steps+1):
                public.step_grid_field(window(n),step_index=n)
                np.testing.assert_array_equal(public._engine.current,values[n-1])
    np.savez(root/'response.npz',observations=values,charge_density=densities)
    summary=dict(steps=4,saved_step=2,prepared_streamed_bit_exact=True,orbital_magnetic_input=True,magnetic_effect=magnetic_effect,
        invalid_geometry_and_infinity_preserve_state=True,compact_complete_magnetic_grid_accepted=True,direct_process_exit_restart_exact=True,
        restore_extra_evolutions=0,public_method_exact=True,physical_convergence=False)
    (root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('executable','bundle','validation','output'):parser.add_argument(name)
    args=parser.parse_args();run(args.executable,args.bundle,args.validation,args.output)
