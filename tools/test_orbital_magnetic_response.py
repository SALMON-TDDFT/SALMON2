"""Read-only spinless magnetic response, Hamiltonian derivative and direct restart.

Reuse a prepared molecular GS. This is a software validation, not convergence.
"""
import argparse
import json
from pathlib import Path
import shutil
import numpy as np


def run(executable,bundle,validation,output):
    from salmonpy.isolated_tddft import load_isolated_runtime,IsolatedTDDFTInitialState
    from salmonpy.near_field import GridFieldSeries
    root=Path(output);root.mkdir(parents=True,exist_ok=False);bundle=Path(bundle)
    runtime=load_isolated_runtime(executable,validation)
    dt=.003;steps=4;counts=(24,24,24);shape=(steps+2,24,24,24,3)
    spacing=np.full(3,2.);base=np.full((3,3),-23.)
    Eorigin=base+np.eye(3);Horigin=base+1.-np.eye(3)
    E=np.zeros(shape);H=np.zeros(shape)
    for n in range(1,steps+2):
        E[n]=np.array([.001,-.002,.0015])*n
        H[n]=np.array([.02,-.01,.015])*n
        H[n,:,:,:,2]+=.00001*np.arange(24)[None,None,:]*n
    fields=GridFieldSeries(np.stack((Eorigin,Horigin)),spacing,np.arange(steps+2)*dt,E,H,
                           molecule_position=(-2.54,1.62,1.94))
    folder=root/'bundle';folder.mkdir()
    for asset in ('C_rps.dat','H_rps.dat','input.inp'):shutil.copy2(bundle/asset,folder/asset)
    shutil.copytree(bundle/'data_for_restart',folder/'data_for_restart');fields.write(folder/'salmonpy-near-eh.dat')
    def model():return IsolatedTDDFTInitialState.from_salmon(folder,'input.inp',runtime=runtime,max_steps=steps,near_field=True)
    observations=[];magnetic=[];electric=[];densities=[];audit=None
    with model() as initial:
        with initial.new_state(dt=dt) as active,initial.new_state(dt=dt) as control:
            old=active._engine.charge_density();densities.append(old.charge_density)
            previous=active._engine.orbital_magnetic_response();endpoints=[]
            for n in range(1,steps+1):
                active.step([0,0,0]);control.step([0,0,0])
                before=active.worker_statistics();response=active._engine.orbital_magnetic_response()
                np.testing.assert_array_equal(response.magnetization,active._engine.orbital_magnetic_response().magnetization)
                assert before==active.worker_statistics()
                np.testing.assert_array_equal(active._engine.current,control._engine.current)
                new=active._engine.charge_density()
                Je=fields.charge_current(new.positions,old.charge_density,new.charge_density,
                                        volume_elements=new.volume_elements,dt=dt)
                Jm=response.current_density(previous=previous)
                endpoints.append(response.current_density());previous=response
                magnetic.append(Jm);electric.append(Je);observations.append(active._engine.current.copy())
                densities.append(new.charge_density);old=new
                if n==2:
                    code,step,_=active._engine._command('ORBITAL_AUDIT_E',11)
                    assert code==0 and step==n
                    audit=np.fromfile(Path(active._engine._work.name)/'run/salmonpy-orbital-audit.bin',dtype=np.float64)
                    print('Hamiltonian derivative audit:',audit.tolist(),flush=True)
                    assert np.max(audit[1:3])<2e-8,audit
                    np.testing.assert_array_equal(response.magnetization,active._engine.orbital_magnetic_response().magnetization)
                    np.testing.assert_array_equal(active._engine.current,control._engine.current)
                    active.save_direct_checkpoint(root/'checkpoint')
        with initial.new_state(dt=dt) as restored:
            restored.restore_direct_checkpoint(root/'checkpoint')
            assert restored.worker_statistics()['updates']==0
            previous=restored._engine.orbital_magnetic_response()
            restored_magnetic=list(magnetic[:2]);restored_electric=list(electric[:2]);old=restored._engine.charge_density()
            np.testing.assert_array_equal(previous.current_density(),endpoints[1])
            for n in range(3,steps+1):
                restored.step([0,0,0]);np.testing.assert_array_equal(restored._engine.current,observations[n-1])
                response=restored._engine.orbital_magnetic_response()
                Jm=response.current_density(previous=previous)
                np.testing.assert_array_equal(Jm,magnetic[n-1]);previous=response
                new=restored._engine.charge_density()
                Je=fields.charge_current(new.positions,old.charge_density,new.charge_density,
                                        volume_elements=new.volume_elements,dt=dt)
                np.testing.assert_array_equal(Je,electric[n-1]);old=new
                restored_magnetic.append(Jm);restored_electric.append(Je)
                np.testing.assert_array_equal(restored._engine.charge_density().charge_density,densities[n])
    assert np.max(np.abs(magnetic))>1e-15
    evidence=json.loads(Path(validation).read_text())
    evidence['orbital_magnetic']=dict(schema='salmonpy.orbital-magnetic-validation/1',hamiltonian_derivative=True,
                                     nonlocal_projectors=True,direct_restart_exact=True)
    receipt=root/'runtime-validation.json';receipt.write_text(json.dumps(evidence,indent=2)+'\n')
    runtime=load_isolated_runtime(executable,receipt)
    with model() as initial:
        with initial.new_state(dt=dt) as state:
            state.step([0,0,0]);np.testing.assert_array_equal(state.orbital_magnetic_response().current_density(),endpoints[0])
    # Numerical summation is done by the shared Fortran primitive.
    from salmonpy.near_field import _kernel,_ptr
    import ctypes as c
    electric=np.array(electric);magnetic=np.array(magnetic);combined=np.empty_like(electric)
    assert _kernel('sp_grid_current_sum')(electric.size,_ptr(electric),_ptr(magnetic),_ptr(combined))==0
    np.savez(root/'currents.npz',electric=electric,magnetic=magnetic,combined=combined,dt=dt,
             start=np.array([8,8,8]),observations=observations)
    re=np.array(restored_electric);rm=np.array(restored_magnetic);rc=np.empty_like(re)
    assert _kernel('sp_grid_current_sum')(re.size,_ptr(re),_ptr(rm),_ptr(rc))==0
    np.testing.assert_array_equal(rc,combined)
    np.savez(root/'restored-currents.npz',electric=re,magnetic=rm,combined=rc,dt=dt,start=np.array([8,8,8]),observations=observations)
    summary=dict(steps=steps,saved_step=2,restore_extra_evolutions=0,magnetic_return_nonzero=True,
                 readout_preserves_state=True,direct_restart_bit_exact=True,ordinary_observables_unchanged=True,
                 derivative_audit=audit.tolist(),spin_included=False,physical_convergence=False)
    (root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('executable','bundle','validation','output'):parser.add_argument(name)
    args=parser.parse_args();run(args.executable,args.bundle,args.validation,args.output)
