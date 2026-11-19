"""ASEなしで8原子立方SiのGSと力を計算。格子は接続例であり収束値ではありません。"""
import argparse
import numpy as np
from salmonpy.dft import DFT, Pseudopotential, SCFSettings
from salmonpy.fdtd import SalmonExecutable


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable'); p.add_argument('sha256'); p.add_argument('pseudo')
    p.add_argument('directory', help='新規GS保存先')
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--no-restart',action='store_true',help='波動関数・GS再開データを保存せず、エネルギーと力を得る')
    p.add_argument('--max-iterations',type=int,default=500,help='SCF上限。少数反復は接続確認用で物理収束の検証ではありません。')
    args = p.parse_args()
    fcc = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    model = DFT(numbers=[14]*8, cell=np.eye(3)*10.26,
        fractional_positions=np.vstack((fcc, fcc+.25)),
        pseudopotentials=[Pseudopotential.from_file(args.pseudo, atomic_number=14, local_channel=2)],
        electrons=32, states=32, grid=(16, 16, 16), kmesh=(2, 2, 2), xc='pz')
    gs = model.run(SalmonExecutable(args.executable, args.sha256), args.directory,
                   scf=SCFSettings(max_iterations=args.max_iterations, tolerance=1e-10), threads=args.threads,
                   save_restart=not args.no_restart)
    print('Energy [Hartree]:', gs.energy_hartree)
    print('SCF converged:',gs.converged,'Restart saved:',gs.restart_saved)
    print('Forces [Hartree/bohr]:', gs.forces_hartree_per_bohr)


if __name__ == '__main__':
    main()
