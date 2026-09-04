# Hybrid Thermal Reoccupation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Initialize divided Hybrid SCF with newly determined 300 K occupations that preserve the total electron number, not the old DC density.

**Architecture:** Keep raw DC projection/support certification separate from the orthonormal local-state initial guess. Determine Fermi–Dirac occupations from current local Hamiltonian states using the existing common chemical potential. Reconstruct and converge the resulting new density using bounded local updates, followed by one terminal LCFO.

**Tech Stack:** Fortran, MPI, existing occupation kernel, Python MPI runners, CMake.

---

## Approved design amendment

The user rejected preserving old occupations after core orthonormalization and
approved filling low-energy states to the original total electron number, at
300 K. This supersedes the C3/C5 requirement that the post-initializer density
and individual orbitals reproduce the old DC state. The saved-DC audit remains
correct as a diagnosis of the old policy; it no longer mandates occupation-
matrix transport. No new density-matrix machinery is required.

Use existing `temperature_k=300` and its Kelvin-to-Hartree conversion; do not
hardcode 300 in the solver or interpret Kelvin as thermal energy. Preserve
the conventional DC occupation API and old density-preserving initializer for
their existing callers. The revised Hybrid route must explicitly select the
new policy; do not silently weaken a legacy check or add a general bypass flag.

Maintain total, not per-fragment, electron count. Equal fragments may receive
equal populations, but heterostructures must be free to exchange charge using
one chemical potential. Retain fractional occupation near the Fermi energy,
complete degenerate shells and guard states. Occupied count and working-space
count are different: do not truncate all computed states to Ne/2. Existing
capacity/tail diagnostics must request extension or reject exhausted spectra.

Raw DC energies/occupations remain immutable reference metadata, not current
Hybrid occupations. They may seed the trial subspace but do not determine the
final population. Current energies come from the existing bounded local
Hamiltonian update/Rayleigh–Ritz path, not WF labels or unchanged DC energies.
No preliminary full-system or convergence-to-exhaustion diagonalization.

Before initial-state construction, preserve exact selection/metric/projection
receipts and checks for raw DC core and required support reconstruction.
After orthogonalization or a Hamiltonian rotation, comparing each resulting
column with its old DC column is no longer the correct support test. Validate
operator sample completeness and representation in the fixed admitted basis,
metric orthogonality and independently reconstructed new-density electron
count instead. Report old/new density differences diagnostically, not as an
acceptance threshold. Do not remove physical support checks wholesale.

One rank per fragment, exact-rank/mapping checkpoint reuse, immutable W90 cache,
user-controlled PW cutoff, maximum three local updates per density epoch,
one-shot terminal LCFO and Exp RT remain unchanged. No new worktree, DC rerun,
material run or changes to unrelated dirty files. Execute here with checkpoints.

### Task O1: Verify the existing 300 K reoccupation kernel

**Checkpoint:** O1 passes on 1/2/4/8 ranks using the existing production
occupation adapter. The test imports the same `kB_au` used by input conversion,
uses one rank per fragment with reversed ownership, and compares against an
independent mu=0 Fermi–Dirac oracle. Unequal local populations, fractional
degenerate-edge filling, global electron conservation, column-order retention
and occupied-terminal-shell extension are checked. Legacy rank-invariant
fingerprints still agree across all four runs. Release build passes and
independent review has no Critical/Important finding. This is existing-kernel
reuse evidence only: no production occupation, initialization or main path
was changed. O2/O3 remain pending; no DC or W90 run was repeated.

Files: `tests/dg/test_dc_fragment_occupation_mpi.f90`, its MPI runner, and this plan.

1. Add a one-rank-per-fragment 300 K case, reversed rank ownership and deliberately
   unequal spectra. Unit core norms represent already orthonormal local states.
2. Compare occupations with an independent Fermi–Dirac formula at the returned
   common mu; verify total electrons, fractions near mu, charge transfer and
   original coefficient-column order. Include a degenerate Fermi edge and
   insufficient thermal-tail rejection. Do not use saved occupations as input.
3. Run `python3 tests/dg/run_dc_fragment_occupation_mpi.py` on 1/2/4/8 ranks;
   retain the existing rank-invariant legacy fingerprint checks. Since this
   tests an existing kernel, success is reuse evidence, not a new feature RED.
4. If a bug appears, follow test-driven-development before production changes.
   Run `cmake --build build-hybrid-release -j 4`, review and checkpoint O1.

### Task O2: Separate projection certification from state initialization

**First numerical checkpoint:** Added `initialize_dg_hybrid_fragment_trial`
with explicit initial count, guard count and optional energy cutoff; it takes
no occupations and returns no density or current-Hamiltonian spectrum. Both
this and the unchanged legacy public initializer use a private seed core,
which validates the mutually exclusive inventory policies collectively before
accessing optional arguments. No synthetic occupation mask is introduced.
The existing normalization, rank, shell and rollback policies remain intact.

The missing-symbol RED was observed before implementation. Fresh subspace
tests pass on 1/2/4/8 ranks, including half-core-norm normalization, explicit
energy-ordered inventory, degenerate guards, invalid/disagreeing count, NaN
and dependent-seed rejection without overwriting the previous state. Selection
tests (including legacy density-preserving admission) pass on 1/2/4/8 and
release builds. Independent review has no Critical/Important issue.

This is only the low-level numerical initializer; its generic distributed
kernel tests do not change the one-rank-per-fragment production requirement.
Raw-reference/projection/support factoring and the single-owner combined
trial-state adapter are still pending O2 work. O3 reoccupation/density and
main integration remain unimplemented. Existing dirty data/logs are retained.

Files: `src/gs/dc/dg_hybrid_fragment_admission.f90`,
`src/gs/dc/dg_hybrid_fragment_subspace.f90`,
`tests/dg/test_dg_hybrid_fragment_selection_mpi.f90`,
`tests/dg/test_dg_hybrid_fragment_subspace_mpi.f90` and corresponding runners.

**Combined trial checkpoint:** Added `prepare_dg_hybrid_selected_trial`.
The new and legacy public entries share a private raw-cache/selection/basis
binding, core projection and required-support verification path. Only after
those checks does the new entry create a temporary trial on MPI_COMM_SELF.
Every fragment must succeed before state publication. Initial counts may
differ per fragment; guard/tolerance/optional-cutoff policy agrees collectively.
Selection/catalog validation retains exactly one rank per fragment.

The returned report sets `trial_prepared=true`, not legacy `valid=true`.
It does not claim occupations, current-Hamiltonian eigenstates, thermal-tail
acceptance or an accepted density. The legacy entry still runs its original
density and final unchanged-raw-column support checks. There is no public
skip-check flag, no fake occupation data and no density rescaling.

The missing-API/report-field RED preceded implementation. The half-core-norm
raw-cache integration now passes the new trial path with physical core
orthogonality while continuing to fail the old density-preserving path.
Stale support fingerprints, missing required points, a single-fragment invalid
count and differing optional-cutoff presence fail without replacing an old
state. Raw W90 call counts remain unchanged. Fresh selection/subspace/300 K
occupation tests pass on 1/2/4/8, and release builds. Independent review found
no Critical/Important issue. A minor existing diagnostic limitation remains:
the generic gate may return an empty local reason on nonfailing ranks when
another fragment's trial initializer fails.

The numerical trial adapter is connected; current-Hamiltonian rotation and
operator-action checks, reoccupation and physical new-density validation are
still O3 integration work. Production support-provider completeness and
C5/C6/main promotion are not certified by this checkpoint. No DC/material run
was repeated; unrelated dirty files and logs remain preserved.

1. RED: introduce an explicit trial-state preparation API which has no density-
   preservation claim. A half-core-norm seed must retain its pre-initialization
   projection/support certification but be allowed to form an orthonormal
   trial state. Keep legacy density-preserving rejection tests unchanged.
2. Accept an explicit initial inventory/guard policy, independent of final
   occupations. Validate count, finite data, degeneracy and metric rank;
   preserve rollback on every rank's failure. Expose trial status, not a
   thermally accepted density. Do not fabricate occupied masks as physical data.
3. Factor the existing raw/selected binding and projection/support checks so
   both entry points reuse them. The new entry must not execute the old
   post-initializer density or unchanged-column support comparison.
4. Verify support loss still fails before state construction; rotations within
   the same admitted span preserve represented operator actions. Run selection
   and subspace MPI runners, build, review and checkpoint.

### Task O3: Connect bounded current-state updates, reoccupation and density

**First small-operator integration checkpoint:** The half-core-norm raw-cache
fixture now uses actual volume/SIPG/nonlocal assembly, the verified trial
adapter and the cached rectangular preconditioner. Its state is refreshed by
the real bounded updater through `run_dc_fragment_occupation_epoch` at
`300*kB_au`. Core weights are measured from the current physical states, not
assigned as identity weights. The target is the sum of the unchanged raw
reference's core electron counts across fragments.

The initial two-state spectrum needs a guard: the existing projected G=0 PW
candidate is extended through the real subspace extension routine, producing
three working states. Refresh calls keep one shared maximum-three-update
budget. Current-state Fermi–Dirac occupations agree with an independent stable
formula, core orthogonality is preserved, and independently reconstructed
density integrates to the target (5/10/20 electrons on 2/4/8 ranks). The new
density differs from raw DC density as intended. Reconstructing with padded
old occupations instead fails the electron-number criterion. No saved DC
occupations are supplied to the new thermal solve, and W90 counts do not grow.

Selection tests pass on 1/2/4/8 ranks (thermal integration on 2/4/8), existing
subspace and occupation regressions pass on 1/2/4/8, and release builds.
Independent review found no Critical/Important issue; its callback-collective
diagnostic suggestion was applied so a local failure returns status to the
outer occupation collective. No production source was modified in this
checkpoint; this is test-harness composition of real numerical routines, not
an installed main-route adapter.

Production transactional rollback/capacity integration, full operator support
coverage and the remaining C5/C6 gates are still open. Do not claim general
material acceptance or Task 8 completion. Existing dirty changes and logs
remain untouched; no material/DC computation was repeated.

Files: selection integration test/runner, `dc_fragment_occupation.f90` only if
needed for the adapter, then the task-owned main-route hunks in parent C6.

1. RED: use the actual volume/SIPG/nonlocal self-block and the new trial state.
   Invoke current-state refresh through the existing occupation epoch at 300 K.
   The half-core-norm example must yield the specified total electron count
   with new occupations, despite a changed density relative to raw DC.
2. Verify shared three-update budget across spectrum extension, current-energy
   occupation ordering, guard/tail acceptance, capacity failure, and no W90/DC
   rerun. Never manually rescale the resulting density or force local Ne.
3. Independently integrate physical core density and compare with the global
   target and sum of current occupations times measured core norms. Compare
   occupations with direct Fermi–Dirac values. A stale old-occupation substitution
   must fail the electron test. Preserve rollback until the full stage passes.
4. Resume outstanding C5 operator-support/cutoff gates, then parent C6 route
   integration. Do not declare Task 8 complete from the occupation tests alone.
   Material/RT validation remains in its later parent tasks.
