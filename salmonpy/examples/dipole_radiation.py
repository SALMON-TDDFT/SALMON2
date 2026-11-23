"""Prepare conservative molecular polarization currents with the bundled Fortran kernel; no SALMON executable.

Run: python -m salmonpy.examples.dipole_radiation
See DIPOLE_RADIATION.md for the optional native MPI field connection.
"""
import argparse
import numpy as np
from salmonpy.outputs import OutputArchive
from salmonpy.dipole_source import DipoleMap


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive',help='outputs/1 directory from isolated_response')
    parser.add_argument('--cell-volume',type=float,help='physical volume per FDTD cell in bohr^3')
    args=parser.parse_args()
    if args.archive:
        if args.cell_volume is None:
            parser.error('--archive requires explicit --cell-volume; do not use the TDDFT vacuum volume')
        archive=OutputArchive.open(args.archive)
        units=archive.metadata.get('units',{})
        if units.get('time')!='atomic_time' or units.get('dipole')!='e*bohr':
            raise ValueError('explicit atomic time and physical charge dipole units required')
        time=archive.array('time');dipole=archive.array('dipole')
        if time.ndim!=1 or len(time)<2 or dipole.shape!=(len(time),3) or time[0]!=0:
            raise ValueError('time[steps+1] from zero and dipole[steps+1,3] required')
        dt=float(time[1]-time[0])
        if dt<=0 or not np.allclose(np.diff(time),dt,rtol=1.e-10,atol=0):
            raise ValueError('uniform positive time spacing required; resample explicitly')
        p=dipole[:,None,:]
        volume=args.cell_volume
    else:
        dt=.003
        t=np.arange(20)*dt
        p=np.zeros((20,1,3))
        p[:,0,2]=1.e-3*np.sin(10*t)
        volume=8. if args.cell_volume is None else args.cell_volume
    mapping=DipoleMap(('molecule-1',),(10,20),[[.25],[.75]],[volume,volume])
    # Explicit one-step hold for native zero-current initialization.
    delayed = np.concatenate((p[:1],p),axis=0)
    current = np.array([mapping.current(a,b,dt=dt) for a,b in zip(delayed[:-1],delayed[1:])])
    np.testing.assert_allclose(np.sum(current*mapping.cell_volumes[None,:,None],axis=1),
                               np.diff(delayed[:,0],axis=0)/dt,rtol=1.e-13,atol=1.e-16)
    assert np.all(current[0] == 0)
    print(f'Prepared {len(current)} intervals at dt={dt} au; explicit delay={dt} au; dipole conservation passed. FDTD propagation is not run by this sample.')


if __name__ == '__main__':
    main()
