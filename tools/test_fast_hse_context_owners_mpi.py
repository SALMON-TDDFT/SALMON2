"""Test internal same-model fast-PR owners: BUILD FULL_K_FIXTURE NEW_OUTPUT.

Explicitly compiles a development probe; never runs on import. Uses two independent
MPI groups, two contexts per process, dt/capacity and waveform differences. Does
not lift the public one-native-lifetime contract or certify physical convergence.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def setting(text, name, value):
    pattern = rf"\b{re.escape(name)}\s*=\s*(?:'[^']*'|[^,\s/]+)"
    updated, count = re.subn(pattern, name + '=' + value, text, flags=re.I)
    if count != 1:
        raise ValueError(f'fixture must contain one explicit {name}')
    return updated


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('build', 'fixture', 'output'):
        parser.add_argument(name, type=Path)
    parser.add_argument('--group-sizes', nargs='+', type=int, default=[1, 2])
    parser.add_argument('--threads', nargs='+', type=int, default=[1, 2])
    args = parser.parse_args()
    build = args.build.resolve()
    if build.is_file():
        build = build.parent
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    library = build / ('libsalmonpy_fast_hse.dylib' if sys.platform == 'darwin' else 'libsalmonpy_fast_hse.so')
    receipt = json.loads((build / 'build.json').read_text())
    if digest(library) != receipt['sha256']:
        raise ValueError('native library does not match its build receipt')
    source = root / 'tests/native/test_hse_context_owners.f90'
    if digest(build / 'src/tddft_context.f90') != digest(root / 'salmonpy/native_adapters/tddft_context.f90'):
        raise ValueError('explicit rebuild required for the current context adapter')
    compiler = shutil.which('mpifort')
    launcher = shutil.which('mpiexec')
    if not compiler or not launcher:
        raise ValueError('host MPI compiler and launcher required')
    probe = out / 'probe'
    compile_command = [compiler, '-fopenmp', '-fcheck=all', '-I' + str(build),
                       '-I' + str(build / 'build'), str(source), str(library), '-o', str(probe)]
    subprocess.run(compile_command, check=True)
    commands = [compile_command]
    output_probe = out / 'output-probe'
    output_source = root / 'tests/native/test_step_output_owners.f90'
    output_command = [compiler, '-fopenmp', '-fcheck=all', '-I' + str(build),
                      '-I' + str(build / 'build'), str(output_source), str(library), '-o', str(output_probe)]
    subprocess.run(output_command, check=True)
    output_work = out / 'output-unit-reuse'
    output_work.mkdir()
    with (output_work / 'run.log').open('w') as log:
        subprocess.run([str(output_probe)], cwd=output_work, stdout=log,
                       stderr=subprocess.STDOUT, check=True)
    commands += [output_command, [str(output_probe)]]
    records = []
    original = (args.fixture / 'rt.inp').read_text()
    if "yn_symmetry='n'" not in original.replace(' ', ''):
        raise ValueError('full-k non-symmetric fixture required; shared symmetry ownership remains separate')
    profiles = [('taylor', 'hse_taylor4', 'off'),
                ('taylor-learned', 'hse_taylor4', 'learned'),
                ('ptcn-learned', 'hse_ptcn', 'learned')]
    for label, propagator, factor in profiles:
        for size in args.group_sizes:
            for threads in args.threads:
                if size < 1 or threads < 1:
                    raise ValueError('positive MPI group size and thread count required')
                case = out / f'{label}-k{size}-omp{threads}'
                case.mkdir()
                text = setting(original, 'propagator', repr(propagator))
                text = setting(text, 'exx_factor_mode', repr(factor))
                text = setting(text, 'exx_factor_exact_interval', '1' if propagator == 'hse_ptcn' else '2')
                text = setting(text, 'yn_predictor_corrector', "'n'" if propagator == 'hse_ptcn' else "'y'")
                env = dict(os.environ, OMP_NUM_THREADS=str(threads), OPENBLAS_NUM_THREADS='1',
                           VECLIB_MAXIMUM_THREADS='1')
                for mode in ('single-a', 'single-b', 'dual', 'protected', 'late-reopen', 'late-close', 'late-failure', 'late-retry'):
                    work = case / mode
                    work.mkdir()
                    for group in range(2):
                        target = work / f'group-{group}'
                        shutil.copytree(args.fixture, target)
                        # Separate solver groups retain independent immutable material inputs.
                        material = text
                        if group:
                            material = material.replace("xc='hse06'", "xc='hse06',hse_omega=.12")
                        (target / 'inputfile').write_text(material)
                    command = [launcher, '-n', str(2 * size), str(probe), mode, str(size)]
                    commands.append(command)
                    with (work / 'run.log').open('w') as log:
                        subprocess.run(command, cwd=work, env=env, stdout=log,
                                       stderr=subprocess.STDOUT, check=True, timeout=180)
                checks = []
                for mode in ('dual', 'protected', 'late-reopen', 'late-close', 'late-failure', 'late-retry'):
                    for group in range(2):
                        for owner in (('b',) if mode.startswith('late-') else ('a', 'b')):
                            control = case / ('single-' + owner) / f'group-{group}'
                            actual = case / mode / f'group-{group}'
                            files = [f'{owner}-trace.dat'] + [f'{owner}-final-{rank}.probe' for rank in range(size)]
                            if mode != 'protected':
                                files += [f'{owner}_rt.data', f'{owner}_rt_energy.data']
                            for name in files:
                                if (control / name).read_bytes() != (actual / name).read_bytes():
                                    raise AssertionError(f'{case.name}/{mode}/group-{group}/{name} differs from isolated execution')
                            checks.append({'mode': mode, 'group': group, 'owner': owner,
                                           'files': {name: digest(actual / name) for name in files}})
                records.append({'profile': label, 'internal_ranks': size, 'threads': threads, 'checks': checks})
                print('Independent RT owners PASS', case.name, flush=True)
    report = {'status': 'pass', 'scope': 'internal non-symmetric same-model RT contexts; public lifetime unchanged',
              'library_sha256': digest(library), 'probe_source_sha256': digest(source),
              'cases': records, 'commands': commands,
              'caveat': 'Synthetic state isolation regression, not material convergence or unrestricted reentrancy.'}
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
