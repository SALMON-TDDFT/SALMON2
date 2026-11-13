"""SALMONの原子単位・単一スピンの_k/_eigen/_tm.dataからSBEを実行。"""
import argparse
import numpy as np
from salmonpy.sbe import SBEInitialState
from salmonpy.fields import AcosPulse
from salmonpy.crystal import respond


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle'); p.add_argument('sysname'); p.add_argument('directory')
    p.add_argument('--bands', type=int, required=True, help='ファイルの全バンド数')
    p.add_argument('--kpoints', type=int, required=True)
    p.add_argument('--electrons', type=int, required=True)
    p.add_argument('--lattice', required=True, help='3x3テキスト、列が格子ベクトル、bohr')
    args = p.parse_args()
    initial = SBEInitialState.from_salmon(args.bundle, args.sysname,
        bands=args.bands, kpoints=args.kpoints, electrons=args.electrons,
        lattice=np.loadtxt(args.lattice))
    drive = AcosPulse(.001, .2, 10.).sample(np.arange(501)*.02)
    result = respond(initial, drive, output_directory=args.directory)
    print(result.current[-1])


if __name__ == '__main__':
    main()
