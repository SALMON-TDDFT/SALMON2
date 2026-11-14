"""Compare explicit-step-config solver with unchanged worker on synthetic GS/2.

python tools/test_rt_step_config_runtime.py NEW_CONVERTER REF_CONVERTER FIXTURES NEW_OUTPUT
Each trajectory owns an independent native worker process: this is not a full
same-process/multiple-material context test. Existing binaries are SHA-verified.
"""
import argparse,json
from dataclasses import replace
from pathlib import Path
import numpy as np
from salmonpy.gs_native import SerialRestartConverter
from salmonpy.gs_archive import GroundStateArchive
from salmonpy.fdtd import SalmonExecutable
from salmonpy.dft import Pseudopotential


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['new_converter','reference_converter','fixtures','output']:p.add_argument(name,type=Path)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    converters=[SerialRestartConverter.from_build(x) for x in [a.reference_converter,a.new_converter]]
    engines=[]
    for x in [a.reference_converter,a.new_converter]:
        r=json.loads(x.read_text());engines.append(SalmonExecutable(r['target_executable'],r['target_sha256']))
    import salmonpy
    pp=Pseudopotential.from_file(Path(salmonpy.__file__).resolve().parent/'examples/si_intro_data/Si_rps.dat',atomic_number=14,local_channel=2)
    t=np.arange(26,dtype=float)+1
    values=.0001*np.column_stack((np.sin(.3*t),np.cos(.2*t)-1,.1*t))
    results=[]
    for case in ['k2','k3','k6','k8','weighted']:
        archive=GroundStateArchive.open(a.fixtures/case/'archive')
        for dt in [.001,.002]:
            for order in [4,6]:
                models=[];states=[]
                try:
                    for converter,engine in zip(converters,engines):
                        model=archive.tddft_initial_state(converter=converter,executable=engine,pseudopotentials=[pp],max_steps=40,threads=1,timeout=60)
                        # Explicit independent test model; the caller's archive/input is not edited.
                        from salmonpy.bands import _input_groups
                        groups=_input_groups(model._input)
                        groups['propagation']=f"\n propagator='middlepoint', yn_predictor_corrector='n', n_hamil={order}\n"
                        model=replace(model,_input='\n'.join('&'+k+v+'/' for k,v in groups.items())+'\n')
                        models.append(model);states.append(model.new_state(dt=dt))
                    left=[];right=[]
                    for index,A in enumerate(values):
                        left.append(states[0].step(A));right.append(states[1].step(A))
                        if index==5:token=states[1].snapshot()
                    left=np.array(left);right=np.array(right)
                    np.testing.assert_allclose(left,right,rtol=1e-12,atol=1e-14)
                    assert states[0].time==states[1].time==26*dt
                    models[1].close();states[1].restore(token)
                    replay=np.array([states[1].step(A) for A in values[6:]])
                    assert replay.tobytes()==right[6:].tobytes()
                    result=dict(case=case,dt=dt,taylor_order=order,steps=26,max_abs_current_difference=float(np.max(abs(left-right))),
                        bitwise=left.tobytes()==right.tobytes(),allclose=True,bitwise_replay_after_model_close=True,
                        final_time=states[1].time,currents=right.tolist())
                    results.append(result);print(case,dt,order,result['max_abs_current_difference'],flush=True)
                finally:
                    for state in states:state.close()
                    for model in models:model.close()
    (a.output/'result.json').write_text(json.dumps(dict(scope='synthetic independent-worker numerical regression, not native context isolation or SCF accuracy',
        reference_sha256=engines[0].sha256,new_sha256=engines[1].sha256,rtol=1e-12,atol=1e-14,cases=results),indent=2)+'\n')


if __name__=='__main__':main()
