"""Validate existing multi-k replay/checkpoint APIs with synthetic runtime fixtures.

python tools/test_gs_multik_checkpoints.py BUILD FIXTURES NEW_OUTPUT
FIXTURES is an output of test_gs_multik_runtime.py, never a collected physical GS.
No compilation or solver modification. Checkpoint outputs contain history only.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

from salmonpy.dft import Pseudopotential
from salmonpy.errors import SnapshotError
from salmonpy.fdtd import SalmonExecutable
from salmonpy.gs_archive import GroundStateArchive
from salmonpy.gs_native import SerialRestartConverter

DT = .001
CAPACITY = 40
PREFIX = 6
CASES = ('k2', 'k3', 'k6', 'k8', 'weighted')


def drive():
    t = np.arange(26, dtype=float) + 1
    return .0001 * np.column_stack((np.sin(t*.3), np.cos(t*.2)-1, t*.1))


def model(build, archive):
    converter = SerialRestartConverter.from_build(build)
    record = json.loads(build.read_text())
    engine = SalmonExecutable(record['target_executable'], record['target_sha256'])
    import salmonpy
    pp = Pseudopotential.from_file(Path(salmonpy.__file__).resolve().parent /
        'examples/si_intro_data/Si_rps.dat', atomic_number=14, local_channel=2)
    return GroundStateArchive.open(archive).tddft_initial_state(
        converter=converter, executable=engine, pseudopotentials=[pp],
        max_steps=CAPACITY, timeout=60, threads=1)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def rejected_without_change(state, checkpoint):
    pid, time, current = state._proc.pid, state.time, state.current.copy()
    try:
        state.restore_checkpoint(checkpoint)
    except SnapshotError:
        pass
    else:
        raise AssertionError('incompatible checkpoint was accepted')
    assert state._proc.pid == pid and state.time == time
    assert state.current.tobytes() == current.tobytes()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build', type=Path)
    parser.add_argument('fixtures', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--resume', choices=CASES)
    args = parser.parse_args()
    values = drive()
    if args.resume:
        # A genuinely fresh Python interpreter recreates its own initial model.
        with model(args.build, args.fixtures/args.resume/'archive') as initial:
            with initial.new_state(dt=DT) as state:
                state.restore_checkpoint(args.output/args.resume/'checkpoint')
                actual = np.array([state.step(A) for A in values[PREFIX:]])
        np.save(args.output/args.resume/'fresh.npy', actual)
        return
    args.output.mkdir(parents=True, exist_ok=False)
    summaries = []
    for name in CASES:
        out = args.output/name
        out.mkdir()
        archive = args.fixtures/name/'archive'
        with model(args.build, archive) as initial:
            with initial.new_state(dt=DT) as first, initial.new_state(dt=DT) as second:
                # Different histories are interleaved, so accidental shared state
                # cannot hide behind identical inputs.
                for A in values[:PREFIX]:
                    first.step(A)
                    second.step(-A)
                token = first.snapshot()
                first.save_checkpoint(out/'checkpoint')
                before = {p.relative_to(out/'checkpoint').as_posix(): digest(p)
                          for p in (out/'checkpoint').rglob('*') if p.is_file()}
                expected, survivor = [], []
                for A in values[PREFIX:]:
                    expected.append(first.step(A))
                    survivor.append(second.step(-A))
                expected, survivor = np.array(expected), np.array(survivor)
                np.save(out/'expected.npy', expected)
                # Close the model while its states remain alive; restore uses
                # the state's retained immutable baseline.
                initial.close()
                first.restore(token)
                replay = np.array([first.step(A) for A in values[PREFIX:]])
                assert replay.tobytes() == expected.tobytes()
                first.close()
                # Closing a sibling must not prevent further evolution.
                survivor_tail = second.step([0, 0, 0])
        with model(args.build, archive) as independent:
            with independent.new_state(dt=DT) as reference:
                for A in values[:PREFIX]: reference.step(-A)
                isolated = np.array([reference.step(-A) for A in values[PREFIX:]])
                isolated_tail = reference.step([0, 0, 0])
                assert isolated.tobytes() == survivor.tobytes()
                assert isolated_tail.tobytes() == survivor_tail.tobytes()
            with independent.new_state(dt=DT*2) as wrong_dt:
                rejected_without_change(wrong_dt, out/'checkpoint')
        # Different k meshes and weights must fail before replay changes state.
        other = 'weighted' if name == 'k2' else 'k2'
        with model(args.build, args.fixtures/other/'archive') as mismatch:
            with mismatch.new_state(dt=DT) as wrong_model:
                rejected_without_change(wrong_model, out/'checkpoint')
        command = [sys.executable, str(Path(__file__).resolve()), str(args.build.resolve()),
                   str(args.fixtures.resolve()), str(args.output.resolve()), '--resume', name]
        with (out/'fresh-process.log').open('w') as log:
            subprocess.run(command, env=dict(os.environ), stdout=log,
                           stderr=subprocess.STDOUT, check=True, timeout=180)
        actual = np.load(out/'fresh.npy', allow_pickle=False)
        assert actual.tobytes() == expected.tobytes()
        after = {p.relative_to(out/'checkpoint').as_posix(): digest(p)
                 for p in (out/'checkpoint').rglob('*') if p.is_file()}
        assert before == after
        assert not any('wfn' in p or 'pseudo' in p or 'restart' in p for p in after)
        summaries.append(dict(case=name, prefix_steps=PREFIX, continuation_steps=20,
            bitwise_memory_replay=True, bitwise_fresh_process_checkpoint=True,
            bitwise_interleaved_isolation=True, bitwise_survivor_after_sibling_close=True,
            incompatible_dt_rejected=True, incompatible_model_rejected=True,
            rejection_preserves_worker=True, checkpoint_unmodified=True,
            checkpoint_bytes=sum((out/'checkpoint'/p).stat().st_size for p in after),
            max_current_difference=float(np.max(abs(expected-actual)))))
        print(name, 'passed', flush=True)
    record = json.loads(args.build.read_text())
    result = dict(scope='synthetic multi-k API regression; no SCF/physical convergence claim',
                  target_sha256=record['target_sha256'], dt=DT, capacity=CAPACITY,
                  cases=summaries)
    (args.output/'result.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
