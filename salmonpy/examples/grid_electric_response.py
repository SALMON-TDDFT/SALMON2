"""Prepare electric-path polarization/current from saved signed charge densities.

INPUT.npz must contain GridFieldSeries arguments origin, spacing, times, E, H,
molecule_position, center, plus charge_positions[N,3], old_charge_density[N],
new_charge_density[N] (signed e/bohr^3), and volume_elements[N] (bohr^3).
This calculation neither launches FDTD/TDDFT nor creates a feedback driver.
"""
import argparse
from pathlib import Path
import numpy as np
from salmonpy.near_field import GridFieldSeries


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    parser.add_argument('output',type=Path,help='new NPZ; never overwritten')
    parser.add_argument('--dt',type=float,required=True,help='charge-density interval in atomic time')
    args=parser.parse_args()
    with np.load(args.input,allow_pickle=False) as data:
        field=GridFieldSeries(**{name:data[name] for name in
            ('origin','spacing','times','E','H','molecule_position','center')})
        positions=data['charge_positions'];volumes=data['volume_elements']
        P=field.charge_polarization(positions,data['new_charge_density'],volume_elements=volumes)
        J=field.charge_current(positions,data['old_charge_density'],data['new_charge_density'],volume_elements=volumes,dt=args.dt)
    with args.output.open('xb') as stream:
        np.savez(stream,P=P,J=J,origin=field.origin[0],spacing=field.spacing,
                 molecule_position=field.molecule_position,center=field.center,dt=args.dt,
                 P_unit='e/bohr^2',J_unit='atomic',representation='electric-component-grid')


if __name__=='__main__':main()
