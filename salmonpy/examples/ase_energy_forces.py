"""ASE構造からPZ/KYのenergy・forces。最適化は精度確認後に別途行います。"""
import argparse
from salmonpy.dft import Pseudopotential, SCFSettings
from salmonpy.fdtd import SalmonExecutable


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable'); p.add_argument('sha256'); p.add_argument('pseudo')
    p.add_argument('directory'); p.add_argument('--threads', type=int, default=1)
    args = p.parse_args()
    from ase.build import bulk
    from salmonpy.ase_calculator import SalmonCalculator
    atoms = bulk('Si', 'diamond', a=5.43, cubic=True)
    atoms.calc = SalmonCalculator(executable=SalmonExecutable(args.executable, args.sha256),
        directory=args.directory,
        pseudopotentials=[Pseudopotential.from_file(args.pseudo, atomic_number=14, local_channel=2)],
        electrons=32, states=32, grid=(32, 32, 32), kmesh=(2, 2, 2),
        scf=SCFSettings(max_iterations=500, tolerance=1e-12), threads=args.threads)
    print('Energy [eV]:', atoms.get_potential_energy())
    print('Forces [eV/Å]:', atoms.get_forces())  # 同じ構造なら上の計算結果を再利用
    print('Saved GS:', atoms.calc.last_ground_state)


if __name__ == '__main__':
    main()
