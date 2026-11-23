"""Internal opt-in TDDFT symmetry regression; not a public Python/FDTD backend.

PROBE PSEUDO NEW_OUTPUT [--ranks 1|2] [--threads 1|2]
Compile tests/native/test_hse_symmetry_tddft.f90 explicitly against matching
MPI fast-PR modules/library. The probe attaches COMM_SELF per host MPI process.
All fixture creation and solver runs are explicit; no network or SCF is used.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('probe', 'pseudo', 'output'):
        parser.add_argument(name, type=Path)
    parser.add_argument('--ranks', type=int, choices=(1, 2), default=1)
    parser.add_argument('--threads', type=int, choices=(1, 2), default=1)
    parser.add_argument('--propagator', choices=('hse_taylor4', 'hse_taylor4_full'), default='hse_taylor4')
    args = parser.parse_args()
    probe = args.probe.resolve(strict=True)
    pseudo = args.pseudo.resolve(strict=True)
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads), OPENBLAS_NUM_THREADS='1',
               MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    rows = []
    commands = []
    fixtures = out / 'fixtures'
    fixtures.mkdir()
    for group in ('full', 'c2z', 'mirror-x'):
        subprocess.run([sys.executable, str(Path(__file__).with_name('prepare_fast_hse_symmetry_fixture.py')),
                        str(pseudo), str(fixtures / group), '--group', group,
                        '--propagator', args.propagator], check=True)

    def execute(group, label, mode, split=0, source=None, corrupt=None, edit=None):
        work = out / group / label
        work.mkdir(parents=True)
        for rank in range(args.ranks):
            target = work / f'rank-{rank}'
            shutil.copytree(fixtures / group, target)
            text = (target / 'rt.inp').read_text()
            # Different values and step counts detect accidental world collectives.
            text = text.replace("xc='hse06'", f"xc='hse06',hse_omega={.11 + .01*rank:.17g}")
            text = text.replace('nt=12', f'nt={12+rank}')
            if edit:
                text = edit(text)
            (target / 'inputfile').write_text(text)
            if source:
                data = (source / f'rank-{rank}/checkpoint.bin').read_bytes()
                if corrupt == 'truncate':
                    data = data[:len(data)//2]
                elif corrupt == 'trailing':
                    data += b'x'
                (target / 'checkpoint.bin').write_bytes(data)
        cmd = ['mpiexec', '-n', str(args.ranks), str(probe), mode, str(split), 'checkpoint.bin']
        commands.append(cmd)
        with (work / 'run.log').open('w') as log:
            subprocess.run(cmd, env=env, cwd=work, stdout=log, stderr=subprocess.STDOUT,
                           check=True, timeout=180)
        return work

    baseline = execute('full', 'continuous', 'continuous')
    for group in ('c2z', 'mirror-x'):
        continuous = execute(group, 'continuous', 'continuous')
        max_current = 0.
        max_energy = 0.
        for rank in range(args.ranks):
            full = np.loadtxt(baseline / f'rank-{rank}/api.dat')
            reduced = np.loadtxt(continuous / f'rank-{rank}/api.dat')
            if not np.isfinite(full).all() or not np.isfinite(reduced).all():
                raise ValueError('nonfinite full/reduced TDDFT trace')
            np.testing.assert_allclose(full, reduced, rtol=1e-10, atol=1e-12)
            max_current = max(max_current, float(abs(full[:, 1:4]-reduced[:, 1:4]).max()))
            max_energy = max(max_energy, float(abs(full[:, 4]-reduced[:, 4]).max()))
        rows.append(dict(group=group, full_k_current_max_abs=max_current, full_k_energy_max_abs=max_energy))
        for split in (6, 7):
            saved = execute(group, f'save-{split}', 'save', split)
            resumed = execute(group, f'resume-{split}', 'resume', split, saved)
            for rank in range(args.ranks):
                expected = (continuous / f'rank-{rank}/trace.bin').read_bytes()
                joined = ((saved / f'rank-{rank}/trace.bin').read_bytes()
                          + (resumed / f'rank-{rank}/trace.bin').read_bytes())
                assert joined == expected, 'resumed internal state differs from continuous'
                expected = np.loadtxt(continuous / f'rank-{rank}/api.dat')
                joined = np.vstack((np.loadtxt(saved / f'rank-{rank}/api.dat'),
                                    np.loadtxt(resumed / f'rank-{rank}/api.dat')))
                np.testing.assert_array_equal(joined, expected)
                assert (resumed / f'rank-{rank}/final-state.bin').read_bytes() == (continuous / f'rank-{rank}/final-state.bin').read_bytes()
                np.testing.assert_array_equal(np.loadtxt(resumed / f'rank-{rank}/calls.dat', dtype=int),
                                              [split, 12+rank-split, 12+rank])
            staged = execute(group, f'staged-{split}', 'staged', split, saved)
            for rank in range(args.ranks):
                target = staged / f'rank-{rank}'
                assert (target / 'before-stage.bin').read_bytes() == (target / 'after-stage.bin').read_bytes()
                expected = (continuous / f'rank-{rank}/trace.bin').read_bytes()
                assert ((saved / f'rank-{rank}/trace.bin').read_bytes() + (target / 'trace.bin').read_bytes()) == expected
                assert (target / 'final-state.bin').read_bytes() == (continuous / f'rank-{rank}/final-state.bin').read_bytes()
                np.testing.assert_array_equal(np.loadtxt(target / 'calls.dat', dtype=int),
                                              [split, 12+rank-split, 12+rank])
            rows.append(dict(group=group, split=split, internal_trace_exact=True,
                             staged_restore_exact=True, final_checkpoint_exact=True, stage_preserves_original=True, accepted_steps_replayed=0))
        field = execute(group, 'field-retry', 'field')
        for rank in range(args.ranks):
            target = field / f'rank-{rank}'
            assert (target / 'before-invalid.bin').read_bytes() == (target / 'after-invalid.bin').read_bytes()
            assert (target / 'trace.bin').read_bytes() == (continuous / f'rank-{rank}/trace.bin').read_bytes()
        rows.append(dict(group=group, broken_field_rejected=True, all_checkpoint_bytes_preserved=True,
                         corrected_field_retry_exact=True))
        shared = execute(group, 'shared-setting-retry', 'global')
        for rank in range(args.ranks):
            target = shared / f'rank-{rank}'
            assert (target / 'before-global.bin').read_bytes() == (target / 'after-global.bin').read_bytes()
            assert (target / 'trace.bin').read_bytes() == (continuous / f'rank-{rank}/trace.bin').read_bytes()
        rows.append(dict(group=group, shared_symmetry_drift_rejected=True, original_state_preserved=True))
        for corruption in ('truncate', 'trailing'):
            rejected = execute(group, 'reject-'+corruption, 'reject', 7,
                               out / group / 'save-7', corruption)
            for rank in range(args.ranks):
                target = rejected / f'rank-{rank}'
                assert (target / 'before-invalid-load.bin').read_bytes() == (target / 'after-invalid-load.bin').read_bytes()
                np.testing.assert_array_equal(np.loadtxt(target / 'calls.dat', dtype=int), [0, 0, 0])
        lifetime = execute(group, 'one-lifetime', 'lifetime')
        for rank in range(args.ranks):
            assert (lifetime / f'rank-{rank}/trace.bin').read_bytes() == (continuous / f'rank-{rank}/trace.bin').read_bytes()
        execute(group, 'public-default-denied', 'default')
        execute(group, 'learned-denied', 'open-reject', edit=lambda t: t.replace(args.propagator, 'hse_taylor4').replace("exx_factor_mode='off'", "exx_factor_mode='learned'"))
        execute(group, 'ptcn-denied', 'open-reject', edit=lambda t: t.replace(args.propagator, 'hse_ptcn').replace("yn_predictor_corrector='y'", "yn_predictor_corrector='n'"))
        if args.ranks == 2:
            execute(group, 'inner-k2-denied', 'inner-reject', edit=lambda t: t.replace('nproc_k=1', 'nproc_k=2'))
    result = dict(status='passed', scope='explicit internal symmetry TDDFT; synthetic states, no public Python/FDTD admission',
                  threads=args.threads, processes=args.ranks, internal_k_ranks=1, propagator=args.propagator,
                  probe_sha256=digest(probe), pseudo_sha256=digest(pseudo), cases=rows, commands=commands)
    (out / 'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print('Internal symmetry TDDFT/full-k/direct restart/field retry PASS', flush=True)


if __name__ == '__main__':
    main()
