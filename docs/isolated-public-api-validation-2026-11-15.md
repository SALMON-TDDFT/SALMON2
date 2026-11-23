# Public isolated TDDFT loader and replay validation

## Change

Added isolated_tddft.py with load_isolated_runtime, IsolatedRuntime,
IsolatedTDDFTInitialState, IsolatedTDDFTState and DipoleObservables. Local bounded
JSON evidence binds the E1 protocol/capabilities to an exact executable hash.
The loader snapshots a supplied compatible native GS bundle. Independent states
lease the shared baseline, own separate worker copies, accept E or FieldPacket,
and return named readonly electronic-change and total-charge dipoles.

Private protocol hooks in the existing TDDFT transport parameterize mode,
commands, response prefix and response width. Periodic defaults remain A/3-current;
private EEngine uses E/6-dipole. Candidate replay preserves the subclass and
validates input-history width separately from observable-history width. Public
E wrapper does not expose periodic current, disk checkpoint or finish APIs.

## Verification

- Six Python fixture tests: independent branches; loader-close lease lifetime;
  wrong dt/position/protocol snapshot rejection before launch; missing/mistimed E;
  binary drift and malformed/false/duplicate receipt rejection; killed-worker
  recovery; injected candidate replay failure leaves original live worker intact.
  These fixture workers verify contracts, not physics.
- Actual C2H2 PZ E1 runtime, SHA256
  abb4fce2106b92c9158e6c77d4111480367bad74cfc926947086d6f649b1c62b:
  public loader/FieldPacket/direct step agree exactly with validated native dipole
  history. Ten-step snapshot restored after closing the initial loader; next ten
  steps exactly match uninterrupted. Killed worker reports failure, reconstructs
  the previously accepted endpoint, and permits subsequent stepping.
- At most two native child processes coexist during the actual tests; OMP and
  OpenBLAS threads set to one. This is lifecycle/interface regression, not a new
  SCF or physical-convergence study. Peak RSS not measured.
- Periodic tools/test_tddft_entrypoint.py: old/new ordinary output byte-identical;
  A-worker baseline, independent states, input rejection and recovery pass through
  the refactored Python transport.
- Full source suite: 408 tests completed, no failures, one platform skip.
- sdist-to-wheel and fresh install: 408 tests completed, no failures, two
  platform/optional-dependency skips. Installed public native regression passes.
- Installed isolated_response example: all archive arrays exactly match source
  example, both named dipole arrays exactly match the native E1 history.
- Read-only review found no required lifetime/replay/periodic-compatibility fix;
  explicit position and protocol mismatch tests were added from review advice.
- Guide/AI JSON synchronization and whitespace checks pass. No GitHub Actions.

## Reproduce with explicit local assets

1. tools/test_isolated_worker.py creates report.json with schema
   salmonpy.isolated-worker-validation/1 and a dipoles.npz trace. It requires the
   explicitly audited E1 executable and original C2H2 reference directory.
2. tools/test_isolated_public.py EXECUTABLE VALIDATION_DIR NEW_OUTPUT exercises
   public lifecycle/replay against that trace; it runs no SCF.
3. python -m salmonpy.examples.isolated_response EXECUTABLE VALIDATION_JSON BUNDLE
   RT_INPUT E_NPZ NEW_OUTPUT writes an outputs/1 archive. E_NPZ has E[steps,3] and dt.
   Run with -I in an installed environment to exclude checkout imports.

## Limits

Only explicit local runtime loading is provided. No E-enabled prebuilt catalog
asset has been published. A receipt is caller-supplied local validation evidence,
not a signature, source audit, general physics accuracy or host portability proof.
Native initialization still rejects unsupported input modes. Compatible native
GS records must be supplied; no periodic-to-isolated conversion or ASE isolated
GS preparation is introduced. Memory replay tokens require the same initial
model object; disk checkpoints and whole-MPI restoration remain pending.
Convergence and choice of self/local-field physics remain user responsibilities.
