# DG Fragment-WF Restart and SCDM Initial-Gauge Design

## Goal

Remove the terminal LCFO occupation failure found by the Si64 fixed-density
continuation diagnostic, avoid repeating the multi-hour fragment-Wannier stage
while debugging downstream DG code, and reduce the cost of clean fragment-WF
generation with an SCDM-derived initial gauge.

This work does not change the accepted physical model: the ordinary-DC density
and potential remain fixed during DG interface continuation, each MPI rank owns
exactly one fragment, Wannier90 remains the final localization refinement, and
the complete Hybrid LCFO Hamiltonian is diagonalized exactly once at the end.

## Observed failure and root cause

The eight-rank Si64 run in
`/tmp/si64-dg-interface-continuation-randomA-omp1-20260906-55uyBx` completed all
six continuation points and then failed in the terminal occupation solve.  The
Schwarz occupation path was called with `300d0` and interprets that value in
kelvin, converting it internally to `k_B T` in hartree.  Its persistent state
therefore records `temperature=300` in kelvin.

The terminal LCFO path passed that stored value directly to
`solve_spectrum_occupations`, whose temperature contract is already in hartree.
The ordinary SALMON input conversion gives 300 K as approximately
`9.50043e-4` hartree.  Passing `300` hartree makes the Fermi distribution nearly
flat over the 232-state terminal spectrum.  The bounded chemical-potential
search cannot reach the 256-electron target and reports
`constant-electron occupation solve did not converge`.

The failure is therefore a temperature-unit boundary bug, not evidence that the
232-state Hybrid basis lacks electron capacity and not a direct SIPG failure.

## Selected approach

Implement three coordinated changes:

1. Correct and make explicit the kelvin/hartree boundary between the Schwarz
   occupation state and the common SALMON occupation kernel.
2. Publish and automatically reuse a strict, versioned fragment-WF checkpoint
   after successful fragment localization.
3. Construct a deterministic SCDM initial gauge within each retained DC
   subspace and pass its unitary gauge to Wannier90 for final refinement.

Wannier90 is not bypassed in the production route.  A raw-SCDM-only production
mode would require separate localization-tail and interface-accuracy evidence
and is outside this change.

## Temperature contract

The common occupation kernel and terminal generalized LCFO continue to accept
electronic temperature in hartree.  The Schwarz occupation routines continue
to accept and record temperature in kelvin because they explicitly multiply by
the Boltzmann constant in hartree per kelvin.

The terminal divided-Hybrid driver must pass the SALMON global `temperature`
value, which has already been converted to atomic units by input processing.  It
must not pass `bounded_schwarz_state%temperature`.  Names and diagnostics at the
boundary must include the unit where practical.  Tests must demonstrate that
300 K and its hartree conversion produce the same occupations and that the
previous 300-hartree misuse reproduces the observed failure or electron-count
inaccessibility.

## SCDM initial gauge

For each fragment, form a column-selection localization gauge from the retained
ordinary-DC subspace on the fragment construction grid.  Use a deterministic
rank-revealing, column-pivoted factorization of the subspace projector (or its
equivalent QR formulation) with stable global-grid-ID tie breaking.  The
selected real-space columns define localized trial anchors.  Assemble their
overlap with the retained states and take the closest unitary polar factor,
reusing the existing Gamma-point `A=<retained|trial>` polar-gauge machinery.

The resulting gauge is only a unitary rotation of the retained DC subspace:

- it must not change the retained projector;
- it must not add symmetry constraints to individual WFs;
- it must remain deterministic for an identical seed and layout;
- it must respect the existing bounded-workspace contract;
- its centers are diagnostic inputs to Wannier90, not a new physical selection
  rule.

Wannier90 receives this SCDM gauge and performs the final spread minimization.
The existing spectral and deterministic-random gauges remain available for
diagnosis and A/B comparison.  SCDM becomes the recommended production default
only after the Si64 comparison confirms projector preservation and reduced
Wannier iteration count without degrading retained-tail/interface diagnostics.

## Fragment-WF checkpoint

The checkpoint serializes the internal post-Wannier fragment basis needed by
the later projection, DG continuation, and terminal LCFO stages.  It is not a
Wannier90 text-output replay.  Each rank writes the payload for its single owned
fragment; a collective manifest makes the generation visible only after every
rank has written and validated its payload.

The format is versioned and records at least:

- format/version and numeric-kind metadata;
- MPI communicator size;
- exact rank-to-fragment mapping and its fingerprint;
- ordinary-DC seed publication ID and seed fingerprint;
- grid, cell, pseudopotential, fragment core/buffer, and boundary-condition
  fingerprints;
- retained-state inventory, ordering, selection, and basis-generation
  fingerprints;
- initial-gauge mode and algorithm/version fingerprint;
- fragment basis coefficients, centers, relevant selection metadata, and
  independent payload hashes.

The formal reuse rule is unchanged: reuse is permitted only when the MPI rank
count and exact rank-fragment association are identical.  All other recorded
physics and representation fingerprints must also match.  There is no
best-effort remapping between ranks or fragments.

The default policy is `auto`:

- a complete and fully compatible generation is reused, and Wannier90 is
  skipped;
- a missing or incompatible generation is rejected with a concise reason, then
  regenerated and republished;
- corrupt, partial, or version-unknown data are never consumed;
- explicit off/read/write controls may be supplied for validation, with strict
  read mode failing rather than silently regenerating.

Publication is transactional.  Rank-local temporary files are validated first;
the manifest is committed last.  A failed run leaves no apparently complete
generation.  A cache hit emits a receipt containing the publication ID, mapping
fingerprint, fragment count, basis fingerprint, and payload fingerprint.

The already completed Si64 run did not create this internal checkpoint.  Its
`.wout` files are retained as evidence but are insufficient to reconstruct the
exact internal fragment coefficients safely.  One clean generation is therefore
still required after this feature is implemented; later terminal/debug runs can
reuse it.

## Data flow

1. Restore the compatible ordinary-DC seed under the existing exact MPI/mapping
   checks.
2. Compute the expected fragment-WF checkpoint identity before invoking
   Wannier90.
3. In `auto` mode, collectively validate the manifest and all rank-local
   payloads.
4. On a hit, restore the internal fragment basis and continue directly to the
   projected WF+PW basis construction.
5. On a miss, construct the SCDM gauge, run Wannier90 refinement, validate the
   localized basis, and transactionally publish the checkpoint.
6. Run fixed-density interface continuation.
7. At lambda one, solve the complete generalized LCFO problem once and derive
   occupations using the electronic temperature in hartree.

## Failure handling

Checkpoint compatibility decisions are collective.  Rank disagreement,
missing peer payloads, fingerprint disagreement, truncated data, or payload-hash
failure rejects the entire generation.  Strict read mode stops with the precise
reason; auto mode regenerates on every rank.  No rank may enter Wannier90 while
another rank takes the cache-hit route.

SCDM rank deficiency, nonfinite factors, excessive polar-unitarity defect, or a
projector-preservation failure rejects the gauge before Wannier90.  Production
may fall back to the existing spectral gauge only when the selected policy
explicitly permits fallback and records that fact; tests and Si64 certification
use fail-closed SCDM mode.

## Testing and acceptance

Test-driven implementation must cover:

- a focused regression reproducing the 300-K/300-hartree unit mismatch and
  proving the corrected terminal occupation path reaches 256 electrons;
- collective temperature-unit and rank-agreement checks;
- deterministic SCDM pivots, unitary gauge, retained-projector invariance, and
  bounded workspace at 1, 2, 4, and 8 MPI ranks where applicable;
- cache write/read round trip with bitwise-stable restored metadata and
  tolerance-appropriate coefficient equality;
- mandatory rejection for changed MPI size or any rank-fragment permutation;
- rejection for changed seed, grid, fragment geometry, selection, gauge version,
  truncated payload, corrupt hash, and incomplete publication;
- collective auto-miss regeneration and collective auto-hit bypass of
  Wannier90;
- source-route checks proving that the terminal LCFO uses hartree temperature,
  performs one full diagonalization, and does no post-LCFO density update.

The final Si64 acceptance run uses the same eight-rank seed and mapping.  It must
complete the terminal LCFO, reproduce 256 electrons at 300 K, publish the
occupied checkpoint, and record Wannier90 iteration counts.  A second run must
show an exact fragment-WF cache hit on all eight ranks and skip Wannier90 while
reproducing the projected-basis and continuation fingerprints exactly and the
terminal floating state to the solver tolerance.  The occupied-checkpoint state
fingerprint is an integrity hash of every floating bit, not a tolerance-aware
physical-state identifier, so its value is recorded but is not required to be
bitwise stable across independent ScaLAPACK solves.

## Task 5 production-smoke evidence (2026-09-07)

The eight-fragment Si8 smoke reused the authoritative ordinary-DC publication
`6738787162284537300` with mapping fingerprint `398352827599223578`; all three
runs reported `scf_skipped=T`.  Complete logs and the machine-readable receipt
are under `/private/tmp/dg-fragment-wf-smoke-task5-08`.

The first production attempt exposed Wannier90 3.1's fixed 50-character
`seedname` limit: the original fragment/generation namespace was truncated
before `.win`.  The fragment scratch components are now the bounded
`fNNNNNN/gNNNNNNNN/w` form, and both the implementation and MPI fixture reject
seeds longer than the library limit.  The next attempt exposed two independent
collective-receipt defects before a complete hit could be certified:

- a fragment-local reference-frame fingerprint was passed where final-state
  validation requires one rank-agreed receipt; production now reduces the
  local fingerprints to a nonzero global-frame fingerprint;
- restore compared the total-communicator checkpoint rank with the rank in a
  one-rank fragment communicator.  The collective reader already validates the
  total rank, so restore now validates only the one-rank fragment layout and
  uses its local rank for cache metadata.

The final miss, hit, and incomplete-publication recovery receipts were:

| case | WF hit | WF publication | projected basis fingerprint | W90 `.wout` | wall time (s) |
|---|---:|---:|---:|---:|---:|
| miss | F | `4636335818477276088` | `7640891576956012809` | 8 | 40.717 |
| hit | T | `4636335818477276088` | `7640891576956012809` | 0 | 23.704 |
| manifest removed / recovery | F | `4903723031505013688` | `7640891576956012809` | 8 | 42.104 |

Every run performed exactly one terminal fixed-density/non-self-consistent LCFO,
performed no density update afterward, wrote an occupied checkpoint whose
complete native stream and payload fingerprint were independently read back,
and recorded 300 K in the Schwarz receipt.  The electron-count defect was
`6.52100596e-11`; terminal residuals were between `8.15e-14` and `8.53e-14`.

The following Task 5 verification completed successfully:

```text
python3 tests/dg/check_dg_fragment_wf_restart_route.py
python3 tests/dg/check_dg_hybrid_fragment_wannier_route.py
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/check_dg_dc_seed_route.py
python3 tests/dg/check_dg_hybrid_localization_first_inputs.py
python3 tests/dg/run_dg_fragment_scdm_gauge_mpi.py
python3 tests/dg/run_dg_fragment_wf_checkpoint_mpi.py
python3 tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_wannier_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_selection_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_subspace_mpi.py
cmake --build /tmp/salmon-task5-build -j4
python3 tests/dg/run_dg_fragment_wf_production_smoke.py \
  --binary /tmp/salmon-task5-build/salmon \
  --result-dir /tmp/dg-fragment-wf-smoke-task5-08 \
  --seed-directory /tmp/dg-fragment-wf-smoke-task5-01/dc-seed
git diff --check
```

## Task 6 Si64 SCDM and exact-reuse evidence (2026-09-07)

The controlled Si64 pair is preserved under
`/private/tmp/si64-t6-scdm-20260907-b`; the complete machine-readable receipt is
`si64_scdm_reuse_evidence.json`.  Both runs used eight MPI ranks, one OpenMP
thread, ordinary-DC publication `7047888166118007469`, and rank-fragment mapping
fingerprint `254086644876463474`.  Both reported `scf_skipped=T`.

The clean run missed the fragment-WF cache and published generation
`7884354523164536186`.  The checkpoint shards authenticate `scdm` gauge
algorithm version 1, candidate and retained rank 400, the fixed rank-to-fragment
mapping, and nonzero per-fragment selection/gauge fingerprints.  The
post-Wannier90 reconstruction defects range from `7.16e-13` to `7.40e-13`,
well below the `1e-9` acceptance bound.  Deterministic
pivot selection and the actual SCDM polar-gauge fingerprint remain covered by
the focused 1/2/4/8-rank MPI contract; the version-1 production checkpoint does
not serialize the pivot-ID list.

| fragment | random iterations | SCDM iterations | SCDM W90 wall (s) |
|---:|---:|---:|---:|
| 1 | 1342 | 631 | 1690.513 |
| 2 | 1192 | 631 | 1691.490 |
| 3 | 2640 | 632 | 1693.709 |
| 4 | 3058 | 631 | 1692.089 |
| 5 | 1927 | 631 | 1691.892 |
| 6 | 1837 | 631 | 1692.011 |
| 7 | 779 | 632 | 1693.606 |
| 8 | 783 | 631 | 1690.041 |

SCDM reduced the summed iteration count from 13,558 to 5,050, a 62.75% reduction
(`0.37247` ratio).  The complete clean run took 2642.7 s.  The exact-reuse run
hit the same publication on all ranks, created zero `.wout` files, and completed
in 893.9 s.  Both produced projected-basis fingerprint
`7640891576956012809`, identical six-point continuation fingerprints, exactly
one terminal LCFO solve, and no post-LCFO density update.

The clean/reuse terminal LCFO residuals were `2.25175905e-13` and
`2.55764319e-13`; electron defects were `1.30341959e-10` and
`1.28807187e-10`, respectively.  The independently authenticated occupied
streams had identical catalog, basis, provenance, and operator fingerprints.
Their bitwise state hashes differed because the independent ScaLAPACK solves
differed by roundoff: maximum occupation, eigenvalue, diagnostic-receipt, and
coefficient-component differences were `9.90e-13`, `9.41e-14`, `3.06e-14`, and
`5.92e-12`.  The relative Frobenius difference between the occupation-weighted
density matrices `C diag(f) C^H` was `4.87e-13`.  The `1e-9` physical-identity
acceptance applies to occupations, eigenvalues, invariant solver receipts, and
this orbital-gauge-invariant density comparison within the identical certified
basis.  Raw coefficient differences are
recorded only as a diagnostic because eigenvectors may acquire arbitrary phase
or rotate within an equally occupied degenerate subspace without changing the
physical state.

Task 6 verification completed with:

```text
python3 tests/dg/check_dg_hybrid_si64_scdm_runner.py
python3 tests/dg/check_dg_fragment_w90_initial_projection.py
python3 tests/dg/check_dg_fragment_wf_restart_route.py
python3 tests/dg/check_dg_hybrid_fragment_wannier_route.py
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/check_dg_dc_seed_route.py
python3 tests/dg/check_dg_hybrid_localization_first_inputs.py
python3 tests/dg/run_dg_fragment_scdm_gauge_mpi.py
python3 tests/dg/run_dg_fragment_wf_checkpoint_mpi.py
python3 tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_wannier_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_selection_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_subspace_mpi.py
python3 tests/dg/run_dg_hybrid_si64_scdm_reuse.py \
  --binary /tmp/salmon-task5-build/salmon \
  --result-dir /tmp/si64-t6-scdm-20260907-b --analyze-existing
cmake --build /tmp/salmon-task5-build -j4
git diff --check
```
