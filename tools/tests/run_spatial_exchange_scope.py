"""Link an audited synthetic oracle probe using an existing MPI Makefiles build."""
from pathlib import Path
import os
import shlex
import subprocess
import sys
import tempfile

root = Path(__file__).resolve().parents[2]
build = Path(sys.argv[1]).resolve()
checker = root / 'tools/check_forbidden_eigensolvers.py'


def run(args, **kwargs):
    return subprocess.run([str(x) for x in args], check=True, **kwargs)


run(['python3', checker, '--source-root', root / 'src'])
run(['python3', checker, '--source-root', root / 'tools/tests'])
with tempfile.TemporaryDirectory(prefix='spatial-exchange-scope-') as directory:
    probe = Path(directory)
    run(['mpifort', '-fopenmp', '-I' + str(build), '-J' + directory, '-c',
         root / 'tools/tests/spatial_exchange_scope.f90', '-o', probe / 'scope.o'])
    link = shlex.split((build / 'src/CMakeFiles/salmon.dir/link.txt').read_text())
    main = 'CMakeFiles/salmon.dir/main.f90.o'
    assert link.count(main) == 1, 'Expected one main object in the existing SALMON link command'
    link[link.index(main)] = str(probe / 'scope.o')
    link[link.index('-o') + 1] = str(probe / 'scope')
    # Preserve the exact library and compiler flags from the validated build.
    # Audit every application object, including our replacement main.
    objects = [Path(x) if Path(x).is_absolute() else build / 'src' / x for x in link if x.endswith('.o')]
    expected = (build / 'src/CMakeFiles/salmon.dir').resolve()
    assert all(p.resolve().is_relative_to(expected) or p.resolve() == (probe / 'scope.o').resolve()
               for p in objects), 'Unexpected application object outside audited trees'
    run(['python3', checker, '--object-tree', expected])
    run(['python3', checker, '--object-tree', probe])
    run(link, cwd=build / 'src')
    env = {**os.environ, 'OMP_NUM_THREADS': '2', 'OPENBLAS_NUM_THREADS': '1'}
    cases = [(1, ''), (2, ''), (4, ''), (8, ''), (4, 'spatial')]
    cases += [(4, mode) for mode in ['nan', 'inf', 'negative', 'source_nan', 'source_inf', 'weight_mismatch']]
    maximum_error = 0.0
    for ranks, mode in cases:
        result = run(['mpiexec', '--oversubscribe', '-n', ranks, probe / 'scope'] + ([mode] if mode else []),
                     env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60)
        if mode not in ('', 'spatial'):
            assert 'GAMMA_PAIR_REJECT_PASS' in result.stdout, result.stdout
        else:
            passes = [line for line in result.stdout.splitlines() if 'SPATIAL_EXCHANGE_PASS' in line]
            assert len(passes) == 12, result.stdout
            maximum_error = max(maximum_error, *(float(line.split()[-1]) for line in passes))
        print(f'MPI{ranks} {mode or "orbital"}: PASS', flush=True)
    print(f'11 layouts/rejections passed; maximum direct-oracle error = {maximum_error:.6e}')
