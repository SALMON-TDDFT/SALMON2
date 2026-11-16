"""Serial user-backend checkpoint admission and exact independent continuation.

Synthetic A-driven recurrence for interface testing, not a material prediction.
No native library, MPI, external input, network or compilation is needed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from salmonpy.checkpointing import CheckpointCapability, check_checkpoint_compatibility

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
    record = check_checkpoint_compatibility(None, GROUPS, states={'film': original.checkpoint_capability})
    payload = directory/'state.json'
    original.save(payload)
    sha = hashlib.sha256(payload.read_bytes()).hexdigest()
    metadata = {'compatibility': record, 'state_sha256': sha}
    (directory/'checkpoint.json').write_text(json.dumps(metadata, indent=2))

    # A fresh host verifies bytes and metadata BEFORE calling the user's loader.
    metadata = json.loads((directory/'checkpoint.json').read_text())
    if hashlib.sha256(payload.read_bytes()).hexdigest() != metadata['state_sha256']:
        raise ValueError('user payload checksum differs')
    restored = UserCurrent()
    check_checkpoint_compatibility(None, GROUPS, states={'film': restored.checkpoint_capability},
        operation='restore', saved=metadata['compatibility'],
        required_files=('state.json',), available_files=('state.json',))
    assert restored.load_calls == 0
    restored.load(payload)
    assert restored.steps == split
    suffix = np.array([restored.step(A) for A in wave[split:]])
    actual = np.concatenate((prefix, suffix))
    np.testing.assert_array_equal(actual, expected)
    assert original.steps == split and restored.steps == total
    assert not np.shares_memory(original.current, restored.current)

    incompatible = UserCurrent(decay=.5)
    rejected = False
    try:
        check_checkpoint_compatibility(None, GROUPS, states={'film': incompatible.checkpoint_capability},
            operation='restore', saved=record)
    except ValueError:
        rejected = True
    assert rejected and incompatible.load_calls == 0
    result = {'status': 'passed', 'steps': total, 'restored_step': split,
              'new_steps_after_load': total-split, 'maximum_error': float(abs(actual-expected).max()),
              'independent_states': True, 'changed_model_rejected_before_load': rejected}
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
