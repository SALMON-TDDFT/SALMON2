"""LDの選択した配列を単位付きで保存し、検証して読み直す。"""
import argparse
import numpy as np
from salmonpy import LorentzDrude
from salmonpy.outputs import save_arrays, OutputArchive


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory', help='新規出力ディレクトリ')
    args = p.parse_args()
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
