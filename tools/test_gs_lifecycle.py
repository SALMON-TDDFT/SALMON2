"""Compile/run the sequential GS lifecycle probe against an existing audited serial HSE build.

BUILD_RECORD REFERENCE_RECORD NEW_OUTPUT. Explicit local compilation only; no downloads.
"""
import argparse
import hashlib
import json
from decimal import Decimal
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys


INPUT = "&calculation\n theory='dft'\n/\n&control\n sysname='gs', base_directory='./', yn_restart='n', yn_self_checkpoint='n'\n/\n&units\n unit_system='au'\n/\n&system\n yn_periodic='y', spin='unpolarized', nelem=1, natom=1, nelec=4, nstate=2, al=8,8,8\n/\n&pseudo\n file_pseudo(1)='Si_rps.dat', izatom(1)=14, lloc_ps(1)=2\n/\n&functional\n xc='hse06', exx_pre_scf_threshold=1d8, exx_pre_scf_steps=2, yn_hse_wannier='n'\n/\n&rgrid\n num_rgrid=8,8,8\n/\n&kgrid\n num_kgrid=2,1,1\n/\n&scf\n nscf=8, convergence='rho_dne', threshold=1d8, nscf_init_mix_zero=1, iseed_number_change=0\n/\n&atomic_red_coor\n'Si' .3 .4 .45 1\n/\n"

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build', type=Path)
    parser.add_argument('reference', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--convergence', choices=('rho_dne','norm_rho','norm_rho_dng','norm_pot','norm_pot_dng'), default='rho_dne')
    parser.add_argument('--mixing', choices=('simple','broyden','pulay','simple_potential'), default='simple')
    parser.add_argument('--auto-mixing', action='store_true')
    parser.add_argument('--distinct-mixing', action='store_true')
    parser.add_argument('--initial-no-mix-steps', type=int, default=1)
    parser.add_argument('--different-normalization', action='store_true')
    args = parser.parse_args()
    if args.different_normalization and (args.convergence!='rho_dne' or args.distinct_mixing):
        parser.error('different normalization requires rho_dne and no distinct mixing')
    root = Path(__file__).resolve().parents[1]
    record = json.loads(args.build.read_text())
    executable = Path(record['executable'])
    if (record.get('schema') != 'salmonpy.hse-rt-state-build/1'
            or record.get('status') != 'built-and-audited'
            or hashlib.sha256(executable.read_bytes()).hexdigest() != record['sha256']):
        raise ValueError('unchanged audited build required')
    reference = json.loads(args.reference.read_text())
    reference_executable = Path(reference['executable'])
    if (reference.get('schema') != record['schema'] or reference.get('status') != 'built-and-audited'
            or hashlib.sha256(reference_executable.read_bytes()).hexdigest() != reference['sha256']):
        raise ValueError('unchanged audited reference required')
    build = executable.parent
    output = args.output.resolve()
    if any(output.is_relative_to(p.resolve()) for p in (root, build, Path(record['source_root']), reference_executable.parent, Path(reference['source_root']))):
        raise ValueError('new external output directory required')
    config = (build / 'config.h').read_text()
    if not re.search(r'^#define USE_HSE\b', config, re.M):
        raise ValueError('HSE build required')
    if re.search(r'^#define (USE_MPI|USE_OPENACC)\b', config, re.M):
        raise ValueError('serial CPU build required')
    output.mkdir(parents=True, exist_ok=False)
    source = output / 'src'
    source.mkdir()
    shutil.copy2(root / 'tests/native/test_gs_lifecycle.f90', source / 'probe.f90')
    flags = (build / 'src/CMakeFiles/salmon.dir/flags.make').read_text()
    options = [v for v in shlex.split(re.search(r'^Fortran_FLAGS\s*=\s*(.*)$', flags, re.M)[1])
               if not v.startswith('-J')]
    includes = shlex.split(re.search(r'^Fortran_INCLUDES\s*=\s*(.*)$', flags, re.M)[1])
    compiler = re.search(r'^CMAKE_Fortran_COMPILER:[^=]*=(.*)$',
                         (build / 'CMakeCache.txt').read_text(), re.M)[1].strip()
    checker = root / 'tools/check_forbidden_eigensolvers.py'
    commands = []
    with (output / 'build.log').open('w') as log:
        def run(command, **kwargs):
            commands.append([str(v) for v in command])
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT, **kwargs)
        run([sys.executable, checker, '--source-root', Path(record['source_root']) / 'src'])
        run([sys.executable, checker, '--source-root', source])
        run([compiler, *options, *includes, '-J' + str(output), '-c', source / 'probe.f90',
             '-o', output / 'probe.o'], cwd=output)
        run([sys.executable, checker, '--object-tree', build])
        run([sys.executable, checker, '--object-tree', output])
        link = shlex.split((build / 'src/CMakeFiles/salmon.dir/link.txt').read_text())
        link[link.index('-o') + 1] = str(output / 'probe')
        link = [str(output / 'probe.o') if v == 'CMakeFiles/salmon.dir/main.f90.o' else v for v in link]
        run(link, cwd=build / 'src')
    root_pseudo = root / 'salmonpy/examples/si_intro_data/Si_rps.dat'
    if not root_pseudo.exists():
        root_pseudo = root.parent / 'examples/si_intro_data/Si_rps.dat'
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    cases = []
    for steps in ((0,) if args.convergence.startswith('norm_pot') or args.mixing=='simple_potential' else (0, 1, 2, 3)):
        directory = output / f'steps-{steps}'
        directory.mkdir()
        for name in ('first', 'second', 'reference', 'reference_second'):
            (directory / name).mkdir()
        shutil.copy2(root_pseudo, directory / 'Si_rps.dat')
        text = INPUT.replace("convergence='rho_dne'", f"convergence='{args.convergence}'").replace('exx_pre_scf_threshold=1d8', f'exx_pre_scf_threshold={1e8 if steps else 0}')
        text = text.replace('&scf\n', f"&scf\n method_mixing='{args.mixing}', yn_auto_mixing='{'y' if args.auto_mixing else 'n'}', update_mixing_ratio=0, mixrate=0.35, alpha_mb=0.65, beta_p=0.55, nmemory_p=2,\n")
        text = text.replace('exx_pre_scf_steps=2', f'exx_pre_scf_steps={max(1, steps)}')
        text=text.replace('nscf_init_mix_zero=1', f'nscf_init_mix_zero={args.initial_no_mix_steps}')
        (directory / 'gs.inp').write_text(text)
        reference_input = text.replace("base_directory='./'", "base_directory='./reference/'")
        try:
            with (directory / 'probe.log').open('w') as log:
                subprocess.run([str(output / 'probe'), *(['normalization'] if args.different_normalization else ['distinct'] if args.distinct_mixing else [])], input=text, text=True, cwd=directory,
                               env=env, stdout=log, stderr=subprocess.STDOUT, timeout=120, check=True)
            with (directory / 'reference.log').open('w') as log:
                subprocess.run([str(reference_executable)], input=reference_input, text=True, cwd=directory,
                               env=env, stdout=log, stderr=subprocess.STDOUT, timeout=120, check=True)
            second_input=reference_input.replace('./reference/', './reference_second/')
            if args.distinct_mixing:
                second_input=second_input.replace('mixrate=0.35','mixrate=0.175').replace('alpha_mb=0.65','alpha_mb=0.325').replace('beta_p=0.55','beta_p=0.275').replace('nmemory_p=2','nmemory_p=1')
            with (directory / 'reference_second.log').open('w') as log:
                subprocess.run([str(reference_executable)], input=second_input, text=True, cwd=directory,
                               env=env, stdout=log, stderr=subprocess.STDOUT, timeout=120, check=True)
            log = (directory / 'probe.log').read_text()
            if log.count('EXX_GS_OWNER explicit') != 2 or 'Sequential GS lifecycle PASS' not in log:
                raise AssertionError('both GS calls must use an explicit owner and complete')
            def residuals(text):
                return [line.split() for line in text.splitlines() if line.strip().startswith('iter and ')]
            reference_residuals=residuals((directory/'reference.log').read_text())
            observed=residuals(log)
            second_residuals=residuals((directory/'reference_second.log').read_text())
            if not reference_residuals or not second_residuals:
                raise AssertionError('both references must contain convergence residuals')
            if args.different_normalization:
                if observed[:len(reference_residuals)]!=reference_residuals or len(observed)!=len(reference_residuals)+len(second_residuals):
                    raise AssertionError('normalization changed iteration count or first state')
                for actual, expected in zip(observed[len(reference_residuals):],second_residuals):
                    actual_value,expected_value=Decimal(actual[-1]),Decimal(expected[-1])
                    # Each E-format value was independently rounded to its last printed digit.
                    rounding=Decimal(10)**actual_value.as_tuple().exponent/2 + Decimal(10)**expected_value.as_tuple().exponent/4
                    if actual[:-1]!=expected[:-1] or abs(actual_value-expected_value/2)>rounding:
                        raise AssertionError('second state density residual must be halved (printed precision)')
            elif not reference_residuals or observed!=reference_residuals+second_residuals:
                raise AssertionError('per-iteration convergence labels/residuals differ from reference')
            trace = (directory / 'a_before.bin').read_bytes()
            if not trace or any((directory / name).read_bytes() != trace
                                for name in ('a_after.bin', 'a_after_close.bin')):
                raise AssertionError('retained A changed while running/closing B')
            files = ['gs_info.data', 'gs_eigen.data', 'gs_k.data']
            restart = directory / 'reference/data_for_restart'
            payloads = sorted(str(p.relative_to(restart)) for p in restart.rglob('*') if p.is_file())
            if not {'wfn.bin', 'rho_inout.bin', 'occupation.bin', 'info.bin'}.issubset(payloads):
                raise AssertionError('reference checkpoint payload incomplete')
            files += ['data_for_restart/' + name for name in payloads]
            hashes = {}
            second_hashes = {}
            for name in files:
                expected = (directory / 'reference' / name).read_bytes()
                for label in ('first', 'second'):
                    if (directory / label / name).read_bytes() != (directory / ('reference' if label=='first' else 'reference_second') / name).read_bytes():
                        raise AssertionError(f'{steps}: {label}/{name} differs from original')
                hashes[name] = hashlib.sha256(expected).hexdigest()
                second_hashes[name] = hashlib.sha256((directory/'reference_second'/name).read_bytes()).hexdigest()
            cases.append({'different_normalization':args.different_normalization, 'observed_residuals':observed, 'initial_no_mix_steps':args.initial_no_mix_steps, 'distinct_mixing':args.distinct_mixing, 'mixing_method':args.mixing, 'auto_mixing':args.auto_mixing, 'convergence_metric':args.convergence, 'pre_scf_steps': steps, 'two_sequential_calls_equal_reference': True, 'retained_a_unchanged': True,
                          'reference_file_sha256': hashes, 'reference_residuals':reference_residuals,
                          'reference_second_file_sha256':second_hashes,
                          'reference_second_residuals':residuals((directory/'reference_second.log').read_text())})
        finally:
            for name in ('first', 'second', 'reference', 'reference_second'):
                shutil.rmtree(directory / name / 'data_for_restart', ignore_errors=True)
            (directory / 'Si_rps.dat').unlink(missing_ok=True)
            for name in ('a_before.bin', 'a_after.bin', 'a_after_close.bin'):
                (directory / name).unlink(missing_ok=True)
    invalid = output / 'invalid-wrapper'
    invalid.mkdir()
    shutil.copy2(root_pseudo, invalid / 'Si_rps.dat')
    with (invalid / 'run.log').open('w') as log:
        rejected = subprocess.run([str(executable)], input=INPUT.replace('nscf=8', 'nscf=0'),
                                  text=True, cwd=invalid, env=env, stdout=log,
                                  stderr=subprocess.STDOUT, timeout=30)
    if rejected.returncode == 0 or 'Invalid GS run configuration' not in (invalid / 'run.log').read_text():
        raise AssertionError('wrapper must fail visibly for invalid SCF configuration')
    (invalid / 'Si_rps.dat').unlink()
    (output / 'result.json').write_text(json.dumps({'status': 'passed', 'cases': cases,
        'scope': 'same-input serial GS reuse and checkpoint equivalence; not concurrency or physical convergence',
        'commands': commands, 'build': record, 'reference': reference,
        'probe_sha256': hashlib.sha256((source / 'probe.f90').read_bytes()).hexdigest()}, indent=2) + '\n')
    print(f'Sequential GS/checkpoint/residual comparison PASS ({len(cases)} cases, {args.convergence})')


if __name__ == '__main__':
    main()
