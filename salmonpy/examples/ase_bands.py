"""8原子立方Siの保存済みGSをASEの経路とバンド図に接続。"""
import argparse
from pathlib import Path
from salmonpy.dft import GroundStateResult
from salmonpy.fdtd import SalmonExecutable


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable'); p.add_argument('sha256'); p.add_argument('gs')
    p.add_argument('directory'); p.add_argument('--threads', type=int, default=1)
    args = p.parse_args()
    import matplotlib
    matplotlib.use('Agg')
    from salmonpy.ase import atoms_from_initial, run_bands
    engine = SalmonExecutable(args.executable, args.sha256)
    with GroundStateResult.open(args.gs).band_initial_state(executable=engine) as initial:
        atoms = atoms_from_initial(initial)
        path = atoms.cell.bandpath('GX,MR', npoints=17)
        result = run_bands(initial, atoms, path, args.directory,
                           bands=24, reference_ev=0., threads=args.threads)
    result.plot(filename=str(Path(args.directory)/'bands.png'), show=False)


if __name__ == '__main__':
    main()
