"""Serial user-backend checkpoint admission and exact independent continuation.

Synthetic A-driven recurrence for interface testing, not a material prediction.
No native library, MPI, external input, network or compilation is needed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from salmonpy.checkpointing import (CheckpointCapability, CheckpointTarget, CheckpointState,
    execute_checkpoint, restore_checkpoint, AtomicRestoreError)

GROUPS = [{'name': 'user', 'ranks': [0], 'role': 'material', 'point_ids': ['film']}]


class UserCurrent:
    def __init__(self, decay=.97):
        self.decay = decay
        self.steps = 0
        self.A = np.zeros(3)
        self.current = np.zeros(3)
        self.load_calls = 0

    @property
    def checkpoint_capability(self):
        identity = hashlib.sha256(json.dumps({'decay': self.decay, 'rule': 1}, sort_keys=True).encode()).hexdigest()
        return CheckpointCapability('example.user-current-json/1', ('save', 'restore'), (1,), identity=identity)

    def step(self, A):
        new = np.array(A, dtype=float, copy=True)
        self.current = self.decay*self.current - (new-self.A)
        self.A = new
        self.steps += 1
        return self.current.copy()

    def save(self, path):
        Path(path).write_text(json.dumps({'steps': self.steps, 'A': self.A.tolist(),
                                        'current': self.current.tolist()}, allow_nan=False))

    def load(self, path):
        self.load_calls += 1
        value = json.loads(Path(path).read_text())
        self.steps = value['steps']
        self.A = np.array(value['A'], dtype=float)
        self.current = np.array(value['current'], dtype=float)


def run(directory):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    total, split = 32, 17
    wave = np.array([[.001*np.sin(.2*n), 0., 0.] for n in range(1, total+1)])
    reference = UserCurrent()
    expected = np.array([reference.step(A) for A in wave])
    original = UserCurrent()
    prefix = np.array([original.step(A) for A in wave[:split]])
    payload = directory/'state.json'
    saved_run = execute_checkpoint(None, GROUPS, operation='save',
        states={'film': CheckpointTarget(original.checkpoint_capability, save=lambda: original.save(payload))})
    record = saved_run['compatibility']
    sha = hashlib.sha256(payload.read_bytes()).hexdigest()
    metadata = {'compatibility': record, 'state_sha256': sha}
    (directory/'checkpoint.json').write_text(json.dumps(metadata, indent=2))

    # A fresh host verifies bytes and metadata BEFORE calling the user's loader.
    metadata = json.loads((directory/'checkpoint.json').read_text())
    if hashlib.sha256(payload.read_bytes()).hexdigest() != metadata['state_sha256']:
        raise ValueError('user payload checksum differs')
    live = CheckpointState(states={'film': UserCurrent()})
    previous = live.states['film']
    for A in wave[:4]: previous.step(A)
    before = (previous.A.tobytes(), previous.current.tobytes(), previous.steps, previous.load_calls)

    def stage(fail=False):
        candidate = UserCurrent()
        candidate.load(payload)
        if fail:
            candidate.current[:] = 12345.
            raise OSError('failure after partial candidate load')
        return candidate

    preserved = False
    try:
        restore_checkpoint(None, GROUPS, live=live,
            states={'film': CheckpointTarget(previous.checkpoint_capability, stage_restore=lambda: stage(True))},
            saved=metadata['compatibility'], required_files=('state.json',), available_files=('state.json',))
    except AtomicRestoreError as exc:
        preserved = exc.original_preserved
    assert preserved and live.states['film'] is previous
    assert before == (previous.A.tobytes(), previous.current.tobytes(), previous.steps, previous.load_calls)
    np.testing.assert_array_equal(previous.step(wave[4]), expected[4])

    restored_run = restore_checkpoint(None, GROUPS, live=live,
        states={'film': CheckpointTarget(previous.checkpoint_capability, stage_restore=stage)},
        saved=metadata['compatibility'], required_files=('state.json',), available_files=('state.json',))
    restored = live.states['film']
    assert restored_run['previous']['states']['film'] is previous
    assert previous.steps == 5 and previous.load_calls == 0 and restored.load_calls == 1
    assert restored.steps == split
    suffix = np.array([restored.step(A) for A in wave[split:]])
    actual = np.concatenate((prefix, suffix))
    np.testing.assert_array_equal(actual, expected)
    assert original.steps == split and restored.steps == total
    assert not np.shares_memory(original.current, restored.current)

    incompatible = UserCurrent(decay=.5)
    other_live = CheckpointState(states={'film': incompatible})
    stage_calls = []
    def incompatible_stage():
        stage_calls.append('load')
        candidate = UserCurrent(decay=.5)
        candidate.load(payload)
        return candidate
    rejected = False
    try:
        restore_checkpoint(None, GROUPS, live=other_live,
            states={'film': CheckpointTarget(incompatible.checkpoint_capability, stage_restore=incompatible_stage)},
            saved=record)
    except ValueError:
        rejected = True
    assert rejected and not stage_calls and incompatible.load_calls == 0
    result = {'status': 'passed', 'steps': total, 'restored_step': split,
              'new_steps_after_load': total-split, 'maximum_error': float(abs(actual-expected).max()),
              'independent_states': True, 'changed_model_rejected_before_load': rejected,
              'failed_restore_original_preserved': preserved}
    np.savez(directory/'response.npz', A=wave, current=actual, reference=expected)
    (directory/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.output)))


if __name__ == '__main__':
    main()
