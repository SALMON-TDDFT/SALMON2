# Fast-PR HSE symmetry owner core — 2026-11-17

Project lead direction: HSE development follows the fast PR, with symmetry central.
The standard k-point port was deferred before changing tracked solver sources.
Standard/HSE symmetry arithmetic differs; standard reduction is not the target
for modifying the fast-PR algorithm. Conventional HSE remains unsupported.

## Implementation boundary

Explicit fast-PR EXX owners can now admit periodic, unpolarized scalar CPU,
full-q, internal-k-size1 symmetry within exx_factor_mode='off' and Taylor/ACE (plus the
existing test-only full-action comparison). Learned/PT-CN symmetry and internal
k-MPI2 remain rejected before mutation. Non-reduced learned/PT-CN is unchanged.
Automatic GS binding retains its existing policy. Python TDDFT/FDTD contexts and
numeric checkpoint paths still refuse symmetry; their prohibitions were NOT
removed. No public Python/C ABI or full symmetric-TDDFT capability is introduced.

The owned full-grid admission uses system ngrid, primitive vectors and hgs rather
than global num_rgrid. Count and per-axis dimensions must agree; scaled norms and
bounded integer products avoid overflow in that check. It runs before refresh/
action mutation. Symmetry matrix layout is checked before admitting an explicit
owner. Bad atoms/operation closure/weights remain deeper native-stop conditions;
this is not a general recovery contract. Existing global-owner arithmetic and
its compatibility grid check are retained. No FFT/communicator owners are cloned.

## Native comparisons

`tests/native/test_hse_symmetry_resources.f90` constructs independent C2z and
mirror-x owners with different cell spacings, screened-exchange omega and source
orbitals. Each has4 representative points for a full2³ mesh. Independently
enumerated full8-point systems use the SAME optimized fast-PR kernel and
independently generated invariant source orbitals. Representative actions and
exchange energy are compared; it is not a comparison against conventional HSE.
The real-space grid is4³ with two occupied scalar orbitals. These cases use
zero-translation two-point stars on a shifted mesh; nontrivial little-group
projectors/fractional translations need additional coverage before wider claims.

The old explicit-owner gate rejects symmetry (status2) after valid full-k
references are generated; this expected failure is retained. A dormant global
num_rgrid check also needed replacement for independent symmetry owners.
New tests poison raw grid/k/functional/electron/temperature inputs after the
owned values are prepared. They then interleave two owners, change B's source,
verify A unchanged, close A/continue B, reopen A and close B/continue A.
Wrong dimensions with the SAME total cell count, unsupported internal k size,
learned history and PT-CN reject with untouched action target and usable original
state. Automatic GS symmetry binding remains rejected, as before.

Final16-case matrix: full-action/ACE × ordinary star expansion/fast orbit reuse ×
MPI process count1/2 × OMP threads1/2. Each MPI process uses COMM_SELF; process
settings differ and one process executes additional refresh/action calls. This
checks independent processes/schedules, NOT internal symmetry k-MPI2. All inputs,
actions and energies checked by the comparison are finite (scalar IEEE checks).
Acceptance tolerance2e-11 times max(1,reference magnitude).

All16 cases passed. Maximum absolute action discrepancy 2.790712e-16;
maximum exchange-energy discrepancy 1.110223e-16.

The unchanged old full-k comparison file matches the new one byte-for-byte in
full-action/MPI1/OMP1. Full-action here is a comparison inside the fast-PR code,
not an additional supported production engine. Source/ACE/orbit reconstruction
algorithms are unchanged; independent-owner admission/grid ownership is new.
The tests are small exchange-operator regressions, not material convergence,
external-field dynamics or end-to-end symmetric TDDFT physics.

## Existing-path regression and distribution

- MPI2 standard PBE warm-start three-iteration GS: freshly generated energy/info,
  eigenvalues,k points and forces identical to unchanged baseline.
- Existing NON-reduced PBE TDDFT and fast-PR Taylor/PT-CN five-rank cases, plus
  seven-rank PBE/HSE mixed32-step cases: state/output rejection and retry pass;
  final arrays and every ordinary stream match unchanged normal baselines exactly.
- Seven-rank NON-reduced protected direct restore: original bytes preserved after
  late staging failure, original continues, retry/three future frames exact with
  zero accepted-history replay; identity negatives still pass.
- Source suite545 tests(one skip), source-distribution→wheel/native hash and
  source/object forbidden-eigensolver audits pass. Fresh installed suite545 tests
  (one skip), packaged native MPI2-process/OMP2 ACE-orbit probe, and existing
  installed seven-rank mixed/five-rank1D fault/retry comparisons all pass.
  JP/EN guides and AI JSON distinguish core and connection scope.

Local Apple Silicon64GB/18-core host. Native matrix uses at most2 MPI processes ×
2 OpenMP threads; normal regression at most7 ranks with OMP/BLAS1.48GiB mixed
budget remains advisory; no hard cap or remote Actions. Explicit fresh external
sources/libraries only; inventories/overrides/provenance updated together. Shared
LD/SBE kernel bytes unchanged. Earlier binaries remain fixed.

## Next HSE work

Prioritize owned field-little-group validation for actual A(t), complete symmetry
map/history/cache persistence, and reduced representative/star internal-MPI
partitioning with mixed failure/restore/retry. Keep current adapter/checkpoint
prohibitions until those pass. Broader shared initialization and stable handles
remain open. Molecular work stays paused; physical convergence is user-owned.

Persistent local evidence:
`/Users/otobetoshihito/.codex/outputs/salmonpy-fast-hse-symmetry-20261117/`.
Includes expected old rejection,old/new full-k references,matrix logs/commands,
prepared source/new-library receipts,normal regression/output inventories,
direct-restore evidence,distribution/tests and deferred standard-port notes.
These paths are local provenance,not installation requirements or a portable
runtime release promise. See paired FAST_HSE_SYMMETRY guides for user scope.
