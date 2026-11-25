"""Prepare spatial E/H and a fixed Ylm center; never start TDDFT.

NPZ: times[T] in atomic time, positions[N,3] in bohr, E/H[T,N,3]
in E:au/H:au, collocated in molecular vacuum. With --yee, supply H_left,
H_right and their time arrays instead of H; Fortran centers the half steps.
"""
import argparse
import numpy as np
from salmonpy.field_sampling import SampledField
from salmonpy.near_field import NearFieldSeries


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('samples');parser.add_argument('output')
    parser.add_argument('--center',type=float,nargs=3,help='fixed expansion origin in bohr; default sample centroid')
    parser.add_argument('--radius',type=float,help='basis scale in bohr; default farthest sample distance')
    parser.add_argument('--lmax',type=int,default=1,help='maximum total degree l+2k, including radial terms (0..8)')
    parser.add_argument('--decay-vector',type=float,nargs=3,default=(0.,0.,0.),help='fixed directional decay vector in inverse bohr; zero keeps polynomial basis')
    parser.add_argument('--yee',action='store_true')
    args=parser.parse_args()
    with np.load(args.samples,allow_pickle=False) as data:
        E=SampledField('E',data['times'],data['positions'],data['E'],'E:au')
        options=dict(center=args.center,radius=args.radius,lmax=args.lmax,decay_vector=args.decay_vector)
        if args.yee:
            left=SampledField('H',data['H_left_times'],data['positions'],data['H_left'],'H:au')
            right=SampledField('H',data['H_right_times'],data['positions'],data['H_right'],'H:au')
            result=NearFieldSeries.from_yee_samples(E,left,right,**options)
        else:
            H=SampledField('H',data['times'],data['positions'],data['H'],'H:au')
            result=NearFieldSeries.from_samples(E,H,**options)
    result.write(args.output)
    print('Prepared spatial E/H:',args.output)
    print('Expansion center (bohr):',result.center.tolist())
    print('Radius (bohr):',result.radius)
    print('Decay vector (inverse bohr):',result.decay_vector.tolist())
    print('Largest relative E/H residual:',result.relative_residual.max(axis=0).tolist())


if __name__=='__main__':main()
