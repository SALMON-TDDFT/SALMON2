"""Prepare a quasistatic electric-potential file; never launch TDDFT.

NPZ must contain times[T] (atomic time), positions[N,3] (bohr), and
values[T,N,3] (E:au), with aligned component positions and times.
This scalar approximation does not provide general spatial E/H coupling.
"""
import argparse
import numpy as np
from salmonpy.field_sampling import SampledField
from salmonpy.near_field import ElectricPotentialSeries


def read(path):
    with np.load(path,allow_pickle=False) as data:
        return SampledField('E',data['times'],data['positions'],data['values'],'E:au')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('samples');parser.add_argument('output')
    parser.add_argument('--center',type=float,nargs=3,required=True)
    parser.add_argument('--radius',type=float,required=True)
    parser.add_argument('--lmax',type=int,default=2)
    parser.add_argument('--reference',help='explicit matching vacuum samples to subtract')
    args=parser.parse_args()
    result=ElectricPotentialSeries.from_samples(read(args.samples),center=args.center,radius=args.radius,lmax=args.lmax,
                                               reference=None if args.reference is None else read(args.reference))
    result.write(args.output)
    print('Prepared electric-potential input:',args.output)
    print('Largest relative field residual:',float(result.relative_residual.max()))


if __name__=='__main__':main()
