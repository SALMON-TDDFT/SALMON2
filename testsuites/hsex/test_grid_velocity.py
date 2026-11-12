"""Check the native scalar/SOC k derivative against an independent Hamiltonian."""
from pathlib import Path
import argparse
import os
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--build', type=Path, required=True)
    args = parser.parse_args()
    build = args.build.resolve()
    command = shlex.split((build / 'src/CMakeFiles/salmon.dir/link.txt').read_text())
    objects = []
    skip = False
    for item in command[1:]:
        if skip:
            skip = False
            continue
        if item == '-o':
            skip = True
        elif item.endswith('/main.f90.o'):
            continue
        elif item.endswith('.o'):
            objects.append(str((build / 'src' / item).resolve()))
        else:
            objects.append(item)
    checker = ROOT / 'tools/check_forbidden_eigensolvers.py'
    with tempfile.TemporaryDirectory(prefix='hsex-grid-velocity-') as temporary:
        path = Path(temporary)
        subprocess.run(['python3', str(checker), '--source-root', str(ROOT / 'src')], check=True)
        subprocess.run([command[0], '-fopenmp', '-fcheck=all', '-I' + str(build), '-c',
                        str(Path(__file__).with_name('fortran_grid_velocity.f90')), '-o', str(path / 'test.o')],
                       cwd=path, check=True)
        for tree in (build, path):
            subprocess.run(['python3', str(checker), '--object-tree', str(tree)], check=True)
        subprocess.run([command[0], *objects, str(path / 'test.o'), '-o', str(path / 'test')],
                       cwd=path, check=True)
        env = dict(os.environ, OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1')
        subprocess.run([str(path / 'test')], cwd=path, env=env, check=True)


if __name__ == '__main__':
    main()
