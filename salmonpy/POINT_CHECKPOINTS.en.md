# Resume a TDDFT point with identity checks

PointState.save_checkpoint / restore_checkpoint add an identity envelope around
the existing managed TDDFT replay checkpoint. The format is
`salmonpy.point-checkpoint/1`: point.json plus a state/ directory containing the
unchanged TDDFT history archive. GS wavefunctions are not included. The caller
must supply compatible local initial data and runtime; nothing is imported or
executed from archive paths.

```python
from salmonpy.point_runs import PointRun

with PointRun({'p': configured_initial}, run_id='experiment-1') as run:
    state = run.new_state('p', dt=0.02)
    state.step([0., 0., 1.e-5])
    state.save_checkpoint(new_checkpoint_directory)

# Explicit resume: use the same IDs, definitions, initial data and numerical settings.
with PointRun({'p': rebuilt_configured_initial}, run_id='experiment-1') as run:
    state = run.new_state('p', dt=0.02)
    state.restore_checkpoint(new_checkpoint_directory)
```

The outer loader requires exact run_id, point_id, material_id/definition_digest
and optional initial_state_id/digest. It rejects different points/runs before any
backend restore. Reusing a run ID here means resuming that run, not creating a
branch. Model objects may be newly constructed; native content/runtime/dt/capacity/
thread checks still apply. This is not a model factory or arbitrary restart converter.

The envelope hashes the inner manifest; the existing loader verifies its arrays
and atomically replays accepted A/current history. Failed replay leaves the target
state as specified by the existing TDDFT contract. Cost grows with history and a
candidate worker may coexist with the target. Keep enough memory for both.
A wrong-ID rejection starts no replay. Saving never overwrites a directory and
removes only its newly created directory on failure.

The envelope rejects symlinks, extra top-level members, duplicate JSON keys and
schema/version mismatches. Identity JSON is bounded to 64 KiB; inner manifest to
16 MiB. `max_bytes` must be a positive integer when supplied, applies to the inner
checkpoint's files/decoded payload, and also bounds prehashing of its manifest.
It is not a total RSS limit. Hashes detect accidental corruption, not authenticity;
pin the outer point.json hash via a trusted source if needed. Do not modify an
archive concurrently with loading.

Only known managed TDDFT states (including ConfiguredState wrappers) are supported.
LD/SBE disk formats, custom serializers, cross-run branching and MPI whole-session
resume are not provided. Existing in-memory snapshot/restore and direct raw TDDFT
checkpoint APIs remain unchanged. PointState disk methods now use this envelope;
a raw TDDFT checkpoint without point.json is deliberately not accepted there.

```sh
python -m salmonpy.examples.point_tddft_restart \
  inputs/materials.json Si-ground /local/runtime waveform.npz /new/output
```

NPZ has real finite A (steps+1,3), A[0]=0 and positive scalar dt. The example saves
halfway, rejects restoration into another point without changing it, rebuilds the
initial model, resumes the same point/run and compares continued currents exactly.
One thread per worker, at most two live workers. For the short Si starter this is
an interface regression, not optical convergence or platform portability evidence.
