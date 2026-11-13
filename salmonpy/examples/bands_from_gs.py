"""dft_siで生成した8原子立方SiのGSからバンドを計算・描画。"""
import argparse
from pathlib import Path
from salmonpy.dft import GroundStateResult
from salmonpy.fdtd import SalmonExecutable
from salmonpy.bands import BandPath


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable'); p.add_argument('sha256'); p.add_argument('gs')
    p.add_argument('directory'); p.add_argument('--threads', type=int, default=1)
    args = p.parse_args()
    import matplotlib
    matplotlib.use('Agg')
    engine = SalmonExecutable(args.executable, args.sha256)
    path = BandPath([[0, 0, 0], [.5, 0, 0], [.5, .5, 0]], [4, 4], ['Γ', 'X', 'M'])
    with GroundStateResult.open(args.gs).band_initial_state(executable=engine) as initial:
        result = initial.run(path, args.directory, bands=24, threads=args.threads)
    # 明示的な絶対エネルギー基準。Fermi準位を推定しません。
    ax = result.plot(reference_ev=0.)
    ax.figure.savefig(Path(args.directory)/'bands.png', dpi=150)


if __name__ == '__main__':
    main()
