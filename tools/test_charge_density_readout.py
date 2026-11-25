"""Explicit C2H2 density readout/restart validation; reuses a prepared GS."""
import argparse
import json
from pathlib import Path
import numpy as np


def run(executable,bundle,validation,output):
    from salmonpy.isolated_tddft import load_isolated_runtime,IsolatedTDDFTInitialState
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    runtime=load_isolated_runtime(executable,validation)
    dt=.003
    snapshots=[];observations=[]
    with IsolatedTDDFTInitialState.from_salmon(bundle,'input.inp',runtime=runtime,max_steps=4) as model:
        with model.new_state(dt=dt) as active,model.new_state(dt=dt) as control:
            # Bootstrap optional evidence against the new executable. The public
            # getter stays disabled until this validation record is complete.
            initial=active._engine.charge_density()
            repeated=active._engine.charge_density()
            np.testing.assert_array_equal(initial.charge_density,repeated.charge_density)
            assert active.worker_statistics()['updates']==0
            snapshots.append(initial)
            for n in range(1,5):
                active.step([.001,-.002,.0015]);control.step([.001,-.002,.0015])
                np.testing.assert_array_equal(active._engine.current,control._engine.current)
                before=active.worker_statistics();sample=active._engine.charge_density();after=active.worker_statistics()
                assert before==after
                assert sample.step_index==n
                np.testing.assert_allclose(-np.sum(sample.charge_density*sample.volume_elements),active.electron_number,rtol=2e-13)
                snapshots.append(sample);observations.append(active._engine.current.copy())
                if n==2:active.save_direct_checkpoint(output/'checkpoint')
        with model.new_state(dt=dt) as restored:
            restored.restore_direct_checkpoint(output/'checkpoint')
            assert restored.worker_statistics()['updates']==0
            np.testing.assert_array_equal(restored._engine.charge_density().charge_density,snapshots[2].charge_density)
            for n in range(3,5):
                restored.step([.001,-.002,.0015])
                np.testing.assert_array_equal(restored._engine.current,observations[n-1])
                actual=restored._engine.charge_density()
                np.testing.assert_array_equal(actual.charge_density,snapshots[n].charge_density)
                np.testing.assert_array_equal(actual.positions,snapshots[n].positions)
                np.testing.assert_array_equal(actual.volume_elements,snapshots[n].volume_elements)
    evidence=json.loads(Path(validation).read_text())
    evidence['charge_density']=dict(schema='salmonpy.charge-density-validation/1',signed_charge=True,
                                   read_only=True,direct_restart_exact=True)
    receipt=output/'runtime-validation.json';receipt.write_text(json.dumps(evidence,indent=2)+'\n')
    runtime=load_isolated_runtime(executable,receipt)
    with IsolatedTDDFTInitialState.from_salmon(bundle,'input.inp',runtime=runtime,max_steps=4) as model:
        with model.new_state(dt=dt) as state:
            np.testing.assert_array_equal(state.charge_density().charge_density,initial.charge_density)
    np.savez(output/'charge-history.npz',positions=initial.positions,volume_elements=initial.volume_elements,
             charge_density=[s.charge_density for s in snapshots],observations=observations,dt=dt)
    result=dict(read_only=True,independent_control_exact=True,direct_process_exit_restart_exact=True,
                restore_extra_evolutions=0,steps=4,saved_step=2,samples=len(initial.charge_density),
                sign='electronic charge density; negative electron number density',physical_convergence=False)
    (output/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('executable','bundle','validation','output'):parser.add_argument(name)
    args=parser.parse_args();run(args.executable,args.bundle,args.validation,args.output)
