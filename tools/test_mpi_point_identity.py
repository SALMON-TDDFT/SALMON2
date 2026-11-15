"""Bounded local 3/4-rank point identity regression; no native solver required."""
import argparse
import os
from pathlib import Path
import sys
# Locate only the launcher helper; package imports still come from the wheel under -I.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_parallel_mpi_tests import run_bounded


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    args.directory.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    for case, ranks in (('success', 3), ('success', 4), ('mismatch', 3), ('mismatch_internal', 4), ('malformed', 3), ('writefailure', 3), ('ioerror', 3)):
        output = args.directory/f'{case}-{ranks}'
        # repr is Python syntax in a direct argv, never a shell fragment.
        code = f'from salmonpy.examples.mpi_point_identity import run; run({str(output)!r}, scenario={case!r})'
        with (args.directory/f'{case}-{ranks}.log').open('w') as log:
            run_bounded(['mpiexec', '-n', str(ranks), sys.executable, '-c', code], env=env, log=log, timeout=45)
        print(case, ranks, 'PASS', flush=True)


if __name__ == '__main__':
    main()
