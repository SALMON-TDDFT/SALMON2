"""合成2バンド模型でA→Jと独立状態を試す。実物質の予測ではありません。"""
import argparse
import numpy as np
from salmonpy.sbe import SBEInitialState
from salmonpy.fields import AcosPulse
from salmonpy.crystal import respond


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory', help='新規出力ディレクトリ')
    args = p.parse_args()
    momentum = np.zeros((2, 2, 3, 1), complex)
    momentum[0, 1, 0, 0] = momentum[1, 0, 0, 0] = 0.3
    initial = SBEInitialState(eigen=[[-.5], [.4]], momentum=momentum,
        nonlocal_velocity=np.zeros_like(momentum), weights=[1.],
        lattice=np.eye(3)*4., electrons=2)
    drive = AcosPulse(.001, .2, 10., polarization=(1., 0., 0.)).sample(np.arange(501)*.02)
    result = respond(initial, drive, output_directory=args.directory)
    # 同じ初期状態から作る各状態は独立。AはA/c、出力は電荷電流密度 j（SALMONの物質電流は -j）。
    a = initial.new_state(dt=.02)
    b = initial.new_state(dt=.02)
    print('別々の波形の電流:', a.step([.001, 0., 0.]), b.step([.002, 0., 0.]))
    print('バッチ計算の最終電流:', result.current[-1])


if __name__ == '__main__':
    main()
