# Saved GS selection validation (2026-11-15)

Scope: M4c, explicit local SALMONpy GS/1 selection by initial_state_id, matched to
an existing configured DFT material. This is not arbitrary restart conversion,
GS/2, distributed TDDFT, or physics-convergence certification.

## Implementation and rejection tests

`salmonpy.initial_states.tddft_from_catalog` pins ground-state.json and validates
its schema/status/file digest container/energy/iterations. The existing GS loader
checks native inventory, hashes and SCF completion. Scientific signature, ordered
positions, actual gs.inp scientific groups and pseudo bytes must match the material.
The executable hash must match the GS generator. After the existing loader copies
the GS, the snapshot is checked against the pinned inventory, and the manifest is
rechecked. A failed post-load validation closes the candidate owner.

Material identity and initial-state identity survive model/state construction.
Initial-state digest is the GS manifest hash; it is not a checkpoint digest.
The sample explicitly rejects complex, boolean and string waveform values before
runtime loading, instead of silently coercing them.

Ten focused tests pass (nine loader tests and one parameterized sample test).
Full source suite: 360 tests, one skip for a case-sensitive-filesystem-only check.
Installed wheel suite: 360 tests, two skips (same filesystem condition and optional
plotting dependency). No new CI run or native runtime build was requested.

## Packaging and real computation

Built sdist then wheel with `python -m build --no-isolation`, max two build jobs.
Installed wheel in an isolated venv; tests and native smoke use `python -I` so the
source checkout is not imported. Optional test dependencies were reused locally.
Mac Apple Silicon / Python 3.14 / macOS 26 wheel; no Linux/Windows claim.

Reproduction with the installed package and a compatible local runtime:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python -I tools/test_saved_gs_selection.py \
  /path/to/local/runtime /new/output/directory
```

The script generates a fresh packaged Si starter GS, pins it in the material
catalog, explicitly collects all its native files plus structure/potential/license,
renames the transfer bundle and original GS path, verifies the relocated bundle,
then loads the GS by ID. It runs the new catalog_tddft example and separately
interleaves two independent states with opposite drives. Python socket calls are
blocked during the smoke; this is not an OS network sandbox.

Final installed run:

- GS energy: -31.74653905752293 Ha; reference difference 0.
- Selected-current difference from direct loading: 0 at all 20 steps.
- Packaged reference comparison: pass (energy atol 1e-7; current atol 1e-9,
  rtol 1e-5, even though measured differences are zero on this host).
- Interleaving does not change the other state; replay agrees exactly.
- material_id `Si-starter`, initial_state_id `Si-ground`, both digests preserved.
- End-to-end wall time about 15 s, one thread per worker, at most two workers.
- Peak memory was not measured; the 4 GiB budget is not hard-enforced.
- Runtime executable SHA256:
  `29444c0e75447a63611d90b2811207a2698c76a4a18005169543cc6c2f039bed`.

Local diagnostic outputs: `/tmp/salmonpy-saved-gs-final-native2-20261115`.
Reports/current arrays/manifests retained. Generated restart payloads are removed
after validation to avoid accumulating GS copies; rerun the script with a new
output directory to regenerate a usable bundle. Raw test/build logs are under
`/tmp/salmonpy-saved-gs-final-*` and are not distributed artifacts.

## Remaining work

The Si starter is 16³ / 2³ k mesh / 20 short steps and is only a connection test.
Point/run IDs, MPI placement, automatic disk identity, general initial-state
prepare recipes, portable GS/2 and memory/physical-convergence controls remain
outside this change. Source and runtime hashes alone do not establish provenance
trust or cross-host ABI compatibility. No commit or push in this continuation.
