"""Analyze stored NPZ time/current/E arrays; never launch a solver.

Preview is the default. --run executes postprocessing and writes a NEW directory.
A decomposed SALMON current file is optional and its sign/dt must be declared.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from salmonpy.analysis import optical_sensitivity, k_current_contributions
from salmonpy.outputs import read_decomposed_current


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('trajectory',type=Path,help='NPZ with time, current and E in atomic units; current is charge density')
    p.add_argument('output',type=Path)
    p.add_argument('--axis',choices=['x','y','z'],default='x')
    p.add_argument('--frequency-range',nargs=3,type=float,metavar=('MIN','MAX','COUNT'),required=True,help='angular atomic units and number of common frequencies')
    p.add_argument('--sample-counts',nargs='+',type=int)
    p.add_argument('--windows',nargs='+',choices=['none','hann','cubic'],default=['none','cubic'])
    p.add_argument('--decomposed-current',type=Path)
    p.add_argument('--decomposition-dt',type=float,help='native run dt, NOT the observation sample stride')
    p.add_argument('--decomposition-convention',choices=['charge','matter'])
    p.add_argument('--run',action='store_true',help='run POSTPROCESSING only; never starts SALMON/MPI')
    a=p.parse_args()
    low,high,count=a.frequency_range
    if not all(np.isfinite([low,high,count])) or low<0 or high<=low or count!=int(count) or count<2:
        p.error('frequency range needs 0 <= MIN < MAX and integer COUNT >= 2')
    with np.load(a.trajectory,allow_pickle=False) as source:
        time=np.array(source['time']);current=np.array(source['current']);electric=np.array(source['E'])
    axis='xyz'.index(a.axis)
    if current.ndim==2:current=current[:,axis]
    if electric.ndim==2:electric=electric[:,axis]
    counts=a.sample_counts or [max(2,len(time)//2),len(time)]
    if a.decomposed_current and (a.decomposition_dt is None or a.decomposition_convention is None):
        p.error('declare --decomposition-dt and --decomposition-convention explicitly')
    if not a.run:
        print(json.dumps(dict(mode='preview_only',starts_solver=False,samples=len(time),sample_counts=counts,
            windows=a.windows,frequency_range=a.frequency_range,output=str(a.output),
            message='--run performs observation analysis only; no SALMON or MPI is launched.'),indent=2))
        return
    if a.output.exists():raise FileExistsError(a.output)
    result=optical_sensitivity(time,current,electric,frequencies=np.linspace(low,high,int(count)),sample_counts=counts,
        windows=a.windows,current_convention='charge',field_kind='total',units='atomic_gaussian')
    arrays={'omega':result.spectra[0].omega}
    for i,s in enumerate(result.spectra):
        arrays[f'sigma_{i}']=s.values;arrays[f'valid_{i}']=s.valid
    report=dict(result.report)
    if a.decomposed_current:
        k=k_current_contributions(read_decomposed_current(a.decomposed_current),dt=a.decomposition_dt,
                                 current_convention=a.decomposition_convention,units='atomic_gaussian')
        report['k_current']=dict(k.report)
        arrays.update(k_time=k.time,k_ids=k.k_ids,k_current=k.current,k_total=k.total)
    a.output.mkdir(parents=True,exist_ok=False)
    np.savez(a.output/'observations.npz',**arrays)
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2,allow_nan=False))


if __name__=='__main__':main()
