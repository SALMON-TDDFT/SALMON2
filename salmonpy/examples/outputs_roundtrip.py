"""LD配列、または既存SOC出力の電流・磁化・スピン流を単位付きで保存する。"""
import argparse
import numpy as np
from salmonpy import LorentzDrude
from salmonpy.outputs import save_arrays, OutputArchive


def save_spin_arrays(rt_data,spin_data,directory):
    """Read existing Fortran observations only; never launch a calculation."""
    from salmonpy.analysis import read_table
    from salmonpy.outputs import read_spin_dynamics
    current=read_table(rt_data);spin=read_spin_dynamics(spin_data)
    names=('Jm_x','Jm_y','Jm_z')
    indices=[current.columns.index(name) for name in names]
    units=[current.units[i] for i in indices]
    if len(set(units))!=1 or units[0] is None:
        raise ValueError('labeled consistent matter-current units required')
    # SALMON Jm is matter current. NumPy's compiled negate supplies charge j=-Jm.
    charge=np.negative(np.column_stack([current.column(n) for n in names]))
    time_unit=current.units[current.columns.index('Time')]
    save_arrays(directory,{'time_current':current.column('Time'),'charge_current':charge,
                           'time_spin':spin.time,'magnetization':spin.magnetization,
                           'spin_current':spin.spin_current,'pauli_current':spin.pauli_current},
                metadata={'units':{'time_current':time_unit,'charge_current':units[0],
                                  'time_spin':spin.time_unit,'magnetization':spin.magnetization_unit,
                                  'spin_current':spin.current_unit,'pauli_current':spin.current_unit},
                          'current_convention':'charge j=-SALMON Jm',
                          'magnetization_definition':'native SALMON mag_x/y/z, component signs preserved; not magnetic moment',
                          'spin_current_definition':'native Pauli-matrix current, no extra hbar/2 factor',
                          'spin_current_axes':['sample','spin polarization x/y/z','flow direction x/y/z'],
                          'time_alignment':'separate native sample grids; no interpolation',
                          'sources':[current.source,spin.source],'restart_checkpoint':False})
    return OutputArchive.open(directory)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory', help='新規出力ディレクトリ')
    p.add_argument('--rt-data', help='既存の電流出力。計算を実行しません')
    p.add_argument('--spin-data', help='既存のSOC磁化・スピン流出力')
    args = p.parse_args()
    if bool(args.rt_data)!=bool(args.spin_data):p.error('--rt-data と --spin-data を両方指定してください')
    if args.rt_data:
        saved=save_spin_arrays(args.rt_data,args.spin_data,args.directory)
        print(saved.metadata)
        return
    time = np.arange(101)*.1
    field = np.zeros((len(time), 3))
    field[:, 0] = .001*np.sin(.2*time)
    material = LorentzDrude(plasma_frequency=.5, strength=[1.], resonance=[.3], damping=[.01])
    result = material.respond(field, dt=.1)
    save_arrays(args.directory, {'time': result.time, 'J': result.current},
                metadata={'units': {'time': 'a.u.', 'J': 'a.u.'},
                          'backend': 'LorentzDrude', 'restart_checkpoint': False})
    saved = OutputArchive.open(args.directory)
    np.testing.assert_array_equal(saved.array('J'), result.current)
    print(saved.metadata, saved.array('J')[-1])


if __name__ == '__main__':
    main()
