"""Postprocess stored energy observations; never launch a solver or MPI.

python -m salmonpy.examples.fdtd_energy INPUT.npz NEW_OUTPUT
NPZ keys (Gaussian AU, charge current j):
 time: (N,), macro_E/charge_j: (N,M,3), point_volumes: (M,),
 material_ids: (M,) Unicode strings, normal: (3,), area: scalar,
 field_energy_change/other_outgoing_energy: scalar Hartree,
 incident_E/incident_H/reflected_E/reflected_H/transmitted_E/transmitted_H: (N,3), or (N,*surface_shape,3) with
 optional area_weights: surface_shape, in bohr**2 covering one transverse cell.
All E/H and E/j must already be collocated, waves direction-separated, ports
equal-area; surface arrays can be nonuniform. Explicit field/other-energy values are required, never
inferred from missing data. Material work includes reversible storage.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from salmonpy.analysis import absorbed_energy_density, reflection_transmission
from salmonpy.fdtd import energy_balance


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('trajectory',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--rtol',type=float,default=.01)
    parser.add_argument('--atol',type=float,default=0.)
    parser.add_argument('--window',choices=('none','hann'),default='none')
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    # Hash the same bytes that are decoded; NPZ never permits Python pickles.
    import io
    payload=args.trajectory.read_bytes()
    with np.load(io.BytesIO(payload),allow_pickle=False) as data:
        work=absorbed_energy_density(data['time'],data['macro_E'],data['charge_j'],
                                    units='atomic_gaussian',current_convention='charge',field_kind='total')
        ports=reflection_transmission(data['time'],*(data[name] for name in
            ('incident_E','incident_H','reflected_E','reflected_H','transmitted_E','transmitted_H')),
            normal=data['normal'],units='atomic_gaussian',window=args.window,
            area_weights=data['area_weights'] if 'area_weights' in data else None)
        ids=data['material_ids']
        if ids.dtype.kind!='U':raise ValueError('material_ids must be Unicode strings')
        budget=energy_balance(work,data['point_volumes'],ports,material_ids=ids,
            transverse_area=data['area'],field_energy_change=data['field_energy_change'],
            other_outgoing_energy=data['other_outgoing_energy'],units='atomic_gaussian',
            rtol=args.rtol,atol=args.atol)
    args.output.mkdir(parents=True,exist_ok=False)
    np.savez(args.output/'observations.npz',absorbed_energy_density=work,
             omega=ports.omega,R=ports.R,T=ports.T,valid=ports.valid)
    report={name:getattr(budget,name) for name in
        ('material_ids','incident_energy','reflected_energy','transmitted_energy','material_work',
         'field_energy_change','other_outgoing_energy','residual','relative_residual','within_tolerance')}
    report.update(config=dict(budget.config),port_config=dict(ports.config),
                  input_sha256=hashlib.sha256(payload).hexdigest())
    (args.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('Energy residual:',budget.residual,'hartree; within tolerance:',budget.within_tolerance)


if __name__=='__main__':main()
