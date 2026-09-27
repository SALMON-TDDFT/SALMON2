#!/usr/bin/env python3
"""Measure missing DC force response using converged displaced H4 calculations.

This evaluates diagnostics, not an MD trajectory or a production force model.
"""
import argparse
import json
import os
from pathlib import Path
from test_dc_force import DCForceTest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,required=True)
    parser.add_argument('--mpiexec',default='mpiexec')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    os.environ['SALMON_TEST_EXE']=str(args.binary.resolve())
    os.environ['SALMON_TEST_MPIEXEC']=args.mpiexec
    case=DCForceTest();rows=[]
    for buffer in (4,2,1):
        energy,force=case.run_case(ranks=2,buffer=buffer)
        row=dict(buffer_grid=buffer,energy_Ha=energy,TS_Ha=case.last_ts,
            frozen_force_Ha_bohr=force[0][0],finite_difference=[])
        for step in (.002,.001):
            ep,_=case.run_case(ranks=2,buffer=buffer,delta=step);fp=case.last_free_energy
            em,_=case.run_case(ranks=2,buffer=buffer,delta=-step);fm=case.last_free_energy
            row['finite_difference'].append(dict(step_bohr=step,
                internal_energy_force=-(ep-em)/(2*step),E_minus_TS_force=-(fp-fm)/(2*step)))
        rows.append(row)
        args.output.write_text(json.dumps(dict(electronic_temperature_K=300,mpi_ranks=2,
            density_threshold=1e-10,units='Hartree and bohr',results=rows),indent=2)+'\n')
        print(json.dumps(row),flush=True)

if __name__=='__main__':main()
