"""Prepare direct Cartesian E/H input from saved collocated FDTD samples; never start TDDFT.

NPZ contains times[T] (atomic time from0), positions[N,3] (bohr), E/H[T,N,3]
(E:au/H:au). With --yee use H_left/H_right and H_left_times/H_right_times.
Counts, origin and spacing describe the full rectangular patch. Molecular position
is in FDTD coordinates; center is the same integration point in TDDFT coordinates.
"""
import argparse
import numpy as np
from salmonpy.near_field import GridFieldSeries
from salmonpy.field_sampling import SampledField


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('samples');p.add_argument('output')
    p.add_argument('--counts',type=int,nargs=3,required=True)
    p.add_argument('--origin',type=float,nargs=3,required=True)
    p.add_argument('--spacing',type=float,nargs=3,required=True)
    p.add_argument('--molecule-position',type=float,nargs=3,required=True)
    p.add_argument('--center',type=float,nargs=3,default=(0.,0.,0.))
    p.add_argument('--time-origin',type=float,default=0.,help='original FDTD time of relative frame0, in atomic units; provenance only')
    p.add_argument('--yee',action='store_true')
    a=p.parse_args()
    options=dict(counts=a.counts,origin=a.origin,spacing=a.spacing,molecule_position=a.molecule_position,center=a.center,time_origin=a.time_origin)
    with np.load(a.samples,allow_pickle=False) as data:
        E=SampledField('E',data['times'],data['positions'],data['E'],'E:au')
        if a.yee:
            left=SampledField('H',data['H_left_times'],data['positions'],data['H_left'],'H:au')
            right=SampledField('H',data['H_right_times'],data['positions'],data['H_right'],'H:au')
            series=GridFieldSeries.from_yee_samples(E,left,right,**options)
        else:
            series=GridFieldSeries.from_samples(E,SampledField('H',data['times'],data['positions'],data['H'],'H:au'),**options)
    series.write(a.output)
    print('Prepared direct-grid E/H input:',a.output)
    print('Molecular position (FDTD bohr):',series.molecule_position.tolist())
    print('Integration center (TDDFT bohr):',series.center.tolist())


if __name__=='__main__':main()
