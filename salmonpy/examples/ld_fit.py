"""Explicit passive LD fit of synthetic data; no solver, table download or GS."""
import argparse
from dataclasses import asdict
import json
import warnings
import numpy as np
from salmonpy.ld import LorentzDrude, fit_lorentz_drude

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tensor',action='store_true',help='synthetic reciprocal tensor fit; requires SciPy')
    parser.add_argument('--poles',type=int,nargs='+',default=[2])
    parser.add_argument('--optimizer',choices=['native','scipy'],default='native')
    parser.add_argument('--starts',type=int,default=4)
    parser.add_argument('--max-iterations',type=int,default=300)
    args=parser.parse_args()
    if args.tensor:
        from salmonpy.ld import TensorLorentzDrude, fit_lorentz_drude_tensor
        strength=np.array([[[.5,-.12,0.],[-.12,.3,.02],[0.,.02,.2]]])
        background=np.diag([1.2,1.3,1.4])
        truth=TensorLorentzDrude(plasma_frequency=1.,strength=strength,
            resonance=[.7],damping=[.08],epsilon_infinity=background)
        initial=TensorLorentzDrude(plasma_frequency=1.,strength=strength*.9,
            resonance=[.65],damping=[.1],epsilon_infinity=background)
        omega=np.linspace(.1,1.2,100)
        model,report=fit_lorentz_drude_tensor(omega,truth.permittivity_tensor(omega),
            initial=initial,units='atomic_gaussian',rtol=1e-6,
            metadata={'reference':'synthetic tensor, not a real TDDFT material'})
        print(json.dumps({'model':asdict(model.config),'report':report},indent=2,allow_nan=False))
        return
    omega=np.linspace(.05,1.2,180)
    truth=LorentzDrude(plasma_frequency=1.,strength=[.03,.4],resonance=[0.,.8],
                      damping=[.04,.07],epsilon_infinity=2.5)
    epsilon=truth.permittivity(omega)
    records=[]
    for poles in args.poles:
        with warnings.catch_warnings(record=True):
            warnings.simplefilter('always')
            model,report=fit_lorentz_drude(omega,epsilon,poles=poles,drude=True,
                units='atomic_gaussian',optimizer=args.optimizer,starts=args.starts,max_iterations=args.max_iterations,
                reference='synthetic Drude + Lorentz, not a real material')
        records.append({'model':None if model is None else asdict(model),'report':report.as_dict()})
    print(json.dumps(records,indent=2,allow_nan=False))

if __name__=='__main__':main()
