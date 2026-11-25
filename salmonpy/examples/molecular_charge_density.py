"""Read current-step molecular charge density from a prepared GS.

No solver build, download or SCF is run. Permanent density storage is opt-in.
Requires an executable validation receipt declaring charge_density.
"""
import argparse
import json
from pathlib import Path
import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('worker','validation','bundle','output'):parser.add_argument(name,type=Path)
    parser.add_argument('--input',default='input.inp')
    parser.add_argument('--steps',type=int,default=4)
    parser.add_argument('--dt',type=float,default=.003)
    parser.add_argument('--E',nargs=3,type=float,default=(.001,-.002,.0015))
    parser.add_argument('--save-density',action='store_true')
    args=parser.parse_args()
    from salmonpy.isolated_tddft import load_isolated_runtime,IsolatedTDDFTInitialState
    runtime=load_isolated_runtime(args.worker,args.validation)
    args.output.mkdir(parents=True,exist_ok=False)
    with IsolatedTDDFTInitialState.from_salmon(args.bundle,args.input,runtime=runtime,max_steps=args.steps) as model:
        with model.new_state(dt=args.dt) as state:
            records=[]
            for n in range(args.steps+1):
                if n:state.step(args.E)
                sample=state.charge_density()
                records.append(dict(step_index=sample.step_index,time=sample.time,
                                    samples=len(sample.charge_density),electron_number=state.electron_number))
                if args.save_density:
                    with (args.output/f'charge-{n:08d}.npz').open('xb') as stream:
                        np.savez(stream,positions=sample.positions,charge_density=sample.charge_density,
                                 volume_elements=sample.volume_elements,step_index=n,time=sample.time)
    (args.output/'summary.json').write_text(json.dumps(dict(records=records,density_saved=args.save_density),indent=2)+'\n')


if __name__=='__main__':main()
