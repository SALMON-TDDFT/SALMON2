# Explicit A binding packet adapter validation

Local Apple Silicon/macOS, Python 3.14. No GitHub Actions or publication.

Implemented `packetize_A_binding`, retaining provider/material identity and
execution/reduction policy. Endpoint and Yee half-step mappings are explicit;
raw endpoint clocks and charge-current signs are preserved. Per-state A history,
creation cleanup, idempotent close and live-state reuse rejection are tested.
The raw backend must honor zero initial A; available A accessors are checked.
Position validation remains the enclosing session's per-point contract.

Six new unit tests cover both phases, callback failures, malformed phase/time,
nonzero initial A cleanup, shared-state refusal and real native SBE parity.
Ten SBE steps match direct calls exactly for two independent states per phase.
Source suite: 429 tests, one skip, no failures.

Three-rank MPI tests use one field rank and two independent material groups,
64 steps at dt=0.003 atomic units. Synthetic A callbacks match previous legacy
Multiscale1D and native E/H histories exactly. Managed process-isolated TDDFT
also matches legacy A/E/j histories exactly for both field implementations.
These are synthetic-GS interface regressions, not converged-material evidence.
Maximum two TDDFT workers, OMP/BLAS one thread each.

Evidence outside checkout:
- `/tmp/salmonpy-A-binding-tests.log`
- `/tmp/salmonpy-A-binding-{ms,native}-20261115`
- `/tmp/salmonpy-A-binding-td-{legacy,packet}-20261115`
- `/tmp/salmonpy-A-binding-native-td-{legacy,packet}-20261115`
- Each MPI run has a neighboring `.log` and output `history.npz`, `timing.json`.
- Native build: `/tmp/salmonpy-h-sample-build-b-20261115/build.json`.
- Managed runtime: `/tmp/salmonpy-release-mac-runtime-20261114`.
- GS fixture: `/tmp/tddft-worker-regression-pz-n/bundle`.

Distribution check builds sdist then wheel in
`/tmp/salmonpy-A-binding-dist-20261115`, installs into a fresh environment
`/tmp/salmonpy-A-binding-installed-20261115`. Dependencies are reused offline
via the existing test dependency path; this does not test dependency resolution.
Installed tests are recorded in `/tmp/salmonpy-A-binding-installed-tests.log`.

No native kernels were modified. General interpolation/subcycling, E-driven LD
packet adaptation, isolated dipole feedback and whole-session rollback remain
outside this change. Absolute-time-dependent custom models must explicitly
validate the selected physical-clock mapping.

Installed wheel result: 429 tests, two skips, no failures. A three-rank native
E/H + two-worker TDDFT run using installed Python `-I` also matches legacy A/E/j
exactly: `/tmp/salmonpy-A-binding-installed-native-td-20261115` (adjacent log).

## Follow-up review corrections

Independent review found two ownership/collective hazards. The public adapter
now rejects collective material backends before preparation and is explicitly
limited to leader_only/one-rank material groups. Field-side MPI is unaffected.
Session admission also compares wrapped raw-state identities across different
bindings; duplicate wrappers are abandoned without closing the original owner.
Regression tests reproduce both failures and verify refusal/no double close.

Review-corrected source tests and build logs:
`/tmp/salmonpy-A-binding-reviewed-final-tests.log`,
`/tmp/salmonpy-A-binding-reviewed-build.log`.
The earlier 429-test evidence predates these two additional regressions.

Final review-corrected results: source 431 tests/one skip; installed wheel 431
tests/two skips, no failures. Rebuilt from sdist, then installed with no dependency
resolution. Installed `-I` three-rank TDDFT cases for both Multiscale1D and native
E/H match earlier legacy A/E/j bit for bit (64 steps). Evidence:
`/tmp/salmonpy-A-binding-reviewed-installed-tests.log` and
`/tmp/salmonpy-A-binding-reviewed-{ms,native}-20261115` with adjacent logs.
A second read-only review found no remaining required corrections in this scope.
