"""保存済みGSからA→J。二つの波形を同じGSから順番に独立実行します。"""
import argparse
from pathlib import Path
import numpy as np
from salmonpy.dft import GroundStateResult
from salmonpy.fdtd import SalmonExecutable
from salmonpy.fields import AcosPulse
from salmonpy.crystal import respond


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable'); p.add_argument('sha256'); p.add_argument('gs')
    p.add_argument('directory', help='新規親ディレクトリ')
    p.add_argument('--threads', type=int, default=1)
    args = p.parse_args()
    output = Path(args.directory)
    output.mkdir(parents=True, exist_ok=False)
    engine = SalmonExecutable(args.executable, args.sha256)
    gs = GroundStateResult.open(args.gs)
    with gs.tddft_initial_state(executable=engine, max_steps=20, threads=args.threads) as initial:
        for i, amplitude in enumerate((.001, .002)):
            drive = AcosPulse(amplitude, .2, .4).sample(np.arange(21)*.02)
            result = respond(initial, drive, output_directory=output/f'run-{i}')
            print('電子物質電流 [a.u.]:', result.current[-1])
    # 短い接続試験用波形。物理評価にはパルス・dt・伝播時間の収束確認が必要。


if __name__ == '__main__':
    main()
