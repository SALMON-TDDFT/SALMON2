"""合成2バンドSBEを場の部品（j→A）と結合する。実物質の予測ではありません。

python -m salmonpy.examples.field_coupling NEW_OUTPUT_DIRECTORY
"""
import argparse
from pathlib import Path
import numpy as np
from salmonpy.sbe import SBEInitialState
from salmonpy.fields import AcosPulse
from salmonpy.maxwell import Transverse, Longitudinal, Multiscale1D, SPEED_OF_LIGHT, couple


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory', help='新規出力ディレクトリ（存在しないこと）')
    args = p.parse_args()
    out = Path(args.directory)
    if out.exists(): raise SystemExit(f'{out} already exists')
    out.mkdir(parents=True)

    momentum = np.zeros((2, 2, 3, 1), complex)
    momentum[0, 1, 2, 0] = momentum[1, 0, 2, 0] = 0.3
    material = SBEInitialState(eigen=[[-.5], [.4]], momentum=momentum,
        nonlocal_velocity=np.zeros_like(momentum), weights=[1.],
        lattice=np.eye(3)*4., electrons=2)
    pulse = AcosPulse(.001, .2, 60., polarization=(0., 0., 1.))  # 1Dではx方向が伝播方向
    dt, steps = .02, 4000

    # 1) 外場をそのまま与える（SALMON trans_longi='tr'）
    tr = couple(material, Transverse(pulse), dt=dt, steps=steps, output_directory=out/'transverse')
    # 2) バルクの縦波配置：誘起場が電流から決まる（trans_longi='lo'）
    lo = couple(material, Longitudinal(pulse), dt=dt, steps=steps, output_directory=out/'longitudinal')
    # 3) 1次元伝播：3セルの層に左から入射（SALMON多スケールの1D）。c*dt = 刻み幅で真空分散なし
    layer = Multiscale1D(cells=3, spacing=SPEED_OF_LIGHT*dt, vacuum=(50, 50), incident=pulse)
    ms = couple(material, layer, dt=dt, steps=steps, output_directory=out/'multiscale_1d')

    # 結果の current は電荷電流密度 j（SBEの状態も j を返す。SALMONの物質電流は -j）
    print('外場の最大|A_z|        :', np.abs(tr.A[:, 0, 2]).max())
    print('縦波配置の最大|A_z|    :', np.abs(lo.A[:, 0, 2]).max())
    print('1D層・各セルの最大|A_z|:', np.abs(ms.A[:, :, 2]).max(axis=0))
    print('出力:', out)


if __name__ == '__main__':
    main()
