# Experimental isolated E-worker validation

## Implemented boundary

Bundled standard runtime only: explicit E1 worker mode, INIT_E/STEP_E and
SALMONPY_TDDFT_E_V1 replies containing electronic dipole change and total charge
dipole. Native propagation order is unchanged. E is placed directly into rt%E_ext
and rt%E_tot at the current step. The ordinary and periodic A paths remain.
PZ, isolated fixed ions, unpolarized scalar pseudopotentials and midpoint only;
external dipole sources are rejected. Tests use initial E=0.

The Python driver under tools is validation-only. FieldAdapter supplies local E
metadata checks. No installed public isolated loader, runtime-catalog capability,
atomic recovery or snapshot API is claimed; replay here means rebuilding from GS
and re-supplying accepted E history. No current-density normalization is provided.

## Local native results

New build: GNU Fortran 15.2, GNU C 15, OpenBLAS, serial/OpenMP one thread;
Libxc/HSE/MPI disabled. Built with two jobs outside the checkout. Source inventory,
shared-kernel equality and forbidden-eigensolver checks passed. Binary SHA256:
`abb4fce2106b92c9158e6c77d4111480367bad74cfc926947086d6f649b1c62b`.

The reference is described in isolated-reference-validation-2026-11-15.md.
C2H2 PZ, 24x24x32, dt=0.001 fs, 20 steps:

- Worker/reference dipole max absolute difference: 3.8458354204148545e-12 e·bohr.
  Test tolerance: rtol=1e-6, atol=1e-10; reference E comes from converted text.
- Repeated worker trajectories: exactly equal.
- Two live workers with opposite E, alternated through FieldAdapter: each exactly
  equals its standalone trajectory; advancing one leaves the other's observable
  and accepted-step count unchanged. Responses to opposite drives differ.
- Missing E rejected by FieldAdapter before advancing the child.
- Rebuild/replay ten accepted E steps, continue ten: exactly equal to uninterrupted.
- Existing tools/test_tddft_entrypoint.py: old/new ordinary periodic output
  byte-identical; A-worker baseline comparison, interleaving, input rejection and
  failure recovery passed with the new binary.

These are wiring/isolation checks, not converged optical observables. No GitHub
Actions were run. The native two-state test uses at most two one-thread children;
RSS is not measured or limited.

## Reproduce

Prepare/build an audited source tree explicitly with
`tools/prepare_salmon_runtime.py`, using an external directory and platform-local
compiler/BLAS configuration. Then, from the source checkout with NumPy available:

```sh
python tools/test_isolated_reference.py REFERENCE_EXECUTABLE NEW_REFERENCE_DIR
PYTHONPATH=. python tools/test_isolated_worker.py NEW_EXECUTABLE NEW_REFERENCE_DIR NEW_RESULTS_DIR
SALMONPY_RUNTIME_SOURCE=PREPARED_SOURCE PYTHONPATH=. python tools/test_tddft_entrypoint.py REFERENCE_EXECUTABLE NEW_EXECUTABLE NEW_PERIODIC_RESULTS_DIR
```

Do not substitute the old executable for NEW_EXECUTABLE: it has no E1 protocol.
Raw outputs and hashes are saved in the supplied new directories.

## Packaging and review

Source suite: 399 tests, one platform skip, passed. sdist-to-wheel build passed;
fresh isolated wheel installation: 399 tests, two platform/optional-dependency
skips, passed. Guide/AI JSON synchronization and git diff whitespace check passed.
A read-only review found no important numerical/periodic-regression issue. Its
blocking-write observation was addressed by a monitored writer timeout in the
validation driver; full E-worker integration was then rerun successfully.

The optional RT step-config overlay's base/replacement hashes were updated and
fresh preparation/source checks passed. Its initialization explicitly rejects E1
rather than claiming isolated E validation for a different propagation path.
HSE source profiles were not modified. No public capability receipt is issued for
this experimental binary.
