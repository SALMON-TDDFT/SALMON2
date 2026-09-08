# DG Fragment-WF Restart and SCDM Initial-Gauge Implementation Plan

> **For Codex:** REQUIRED SUB-SKILL: Use executing-plans to implement this plan task-by-task. Use systematic-debugging for the observed occupation failure and test-driven-development for every behavior change.

**Goal:** Fix the terminal 300 K occupation-unit bug, reduce clean fragment-Wannier cost with a deterministic SCDM gauge, and make the validated fragment-WF result safely reusable under exact MPI rank/fragment compatibility.

**Architecture:** Keep kelvin only inside the Schwarz thermal-inventory contract and pass SALMON's atomic-unit temperature to the common terminal occupation kernel. Add a fragment-local SCDM gauge builder that returns a unitary rotation of the retained DC subspace and feeds the existing Wannier90 Gamma `A` path. Add a versioned rank-local fragment-WF payload plus collectively committed manifest; production `auto` mode either restores all fragments or regenerates all fragments.

**Tech Stack:** Fortran 2008, MPI, LAPACK/BLAS, SALMON DG and Wannier90 interfaces, Python source-contract tests, CMake/CTest.

**Workspace constraint:** Work only in `/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/.worktrees/wpw-s-orthogonal-complement`. Do not create another worktree. Preserve all existing dirty changes and validation logs; stage only task-owned files or hunks.

---

### Task 1: Correct the terminal electronic-temperature unit

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_hybrid_generalized_eigensystem_mpi.f90`
- Modify: `tests/dg/check_dg_hybrid_divided_lcfo_route.py`
- Modify if needed: `tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py`

**Step 1: Write failing unit and route regressions**

Add a 256-electron, 232-state finite-temperature occupation case.  Verify that
300 K converted with the Schwarz Boltzmann constant reaches the electron target,
while passing the raw kelvin value to the hartree-temperature occupation kernel
does not satisfy the same contract.  In the production route check, reject
`electronic_temperature=bounded_schwarz_state%temperature` and require the
already converted SALMON `temperature` value.

**Step 2: Run RED verification**

Run:

```bash
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py
```

Expected: the route test fails on the stored-kelvin argument and the focused
regression demonstrates the old mismatch.

**Step 3: Make the minimal unit-boundary fix**

Pass `max(0d0,temperature)` to the terminal generalized LCFO publication call.
Retain `bounded_schwarz_state%temperature` as kelvin metadata used only by the
Schwarz occupation implementation and its diagnostics.  Add unit-explicit local
names/comments at the boundary; do not introduce a second conversion.

**Step 4: Run GREEN verification**

Run the two focused checks at their supported MPI ranks, the occupation-kernel
tests, `git diff --check`, and `cmake --build build-hybrid-release -j2`.

**Step 5: Commit and checkpoint**

Stage only Task 1 hunks and commit:

```bash
git commit -m "fix(dg): pass atomic-unit temperature to terminal LCFO"
```

Review the diff and test evidence before continuing.  Do not rerun Si64 yet.

---

### Task 2: Build a deterministic fragment SCDM gauge

**Files:**
- Create: `src/gs/dc/dg_fragment_scdm_gauge.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/test_dg_fragment_scdm_gauge_mpi.f90`
- Create: `tests/dg/run_dg_fragment_scdm_gauge_mpi.py`

**Step 1: Write the failing SCDM contract test**

Construct synthetic orthonormal retained subspaces with known localized column
choices.  Require deterministic pivot IDs with global-ID tie breaking, a square
unitary gauge, projector invariance before/after rotation, stable fingerprints,
and a reported workspace peak within a supplied byte cap.  Cover complex phases,
nearly tied pivots, rank deficiency, nonfinite input, insufficient byte limit,
and rank-disagreeing metadata.

**Step 2: Run RED verification**

Run `python3 tests/dg/run_dg_fragment_scdm_gauge_mpi.py`; expect failure because
the module is absent.

**Step 3: Implement pivot selection and polar gauge**

Select real-space columns from the retained-subspace projector using a
deterministic rank-revealing pivoted factorization without materializing a dense
global-grid projector.  Resolve numerical ties by the canonical global grid ID.
Build localized trial anchors from the selected columns and use the existing
Gamma polar-factor convention for `A=<retained|trial>`.  Reject a deficient
factorization or excessive unitarity/projector defect.

**Step 4: Certify bounded memory and collectivity**

Account for all temporary arrays before allocation.  Require communicator-wide
agreement on dimensions, fragment identity, controls, selected IDs, and final
fingerprint.  Do not gather all fragments or all global-grid columns to one rank.

**Step 5: Run GREEN verification and commit**

Run the new runner at `-np 1,2,4,8`, the existing Wannier90 Gamma gauge tests,
the fragment-Wannier tests, `git diff --check`, and the release build.  Commit:

```bash
git commit -m "feat(dg): add deterministic fragment SCDM gauge"
```

---

### Task 3: Add a strict fragment-WF checkpoint format

**Files:**
- Create: `src/gs/dc/dg_fragment_wf_checkpoint.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/test_dg_fragment_wf_checkpoint_mpi.f90`
- Create: `tests/dg/run_dg_fragment_wf_checkpoint_mpi.py`

**Step 1: Write failing round-trip and rejection tests**

Publish synthetic rank-local fragment basis payloads and a collective manifest,
then restore them and compare all metadata, coefficients, centers, and
fingerprints.  Require exact rejection for changed MPI size, rank-fragment
permutation, DC publication ID/fingerprint, grid/cell/fragment geometry,
pseudopotential, retained inventory/order, basis generation, gauge algorithm
version, and payload hash.

Add incomplete publication, truncated payload, corrupt scalar, unknown version,
missing peer file, and rank-disagreeing requested mode cases.  Verify that strict
read fails and auto mode reports a collective miss without partially restoring
any rank.

**Step 2: Run RED verification**

Run `python3 tests/dg/run_dg_fragment_wf_checkpoint_mpi.py`; expect failure
because the checkpoint module is absent.

**Step 3: Implement versioned rank-local payloads**

Define the minimal internal post-Wannier payload used downstream.  Serialize
numeric-kind/version metadata, exact mapping identity, all required provenance,
dimensions, basis coefficients, centers, selection metadata, and independent
hashes.  Use explicit rank/fragment filenames within one generation directory.

**Step 4: Implement transactional collective publication**

Write and reread rank-local temporary payloads first.  Collectively confirm all
hashes and compatibility, then have rank zero commit the manifest last.  Restore
only after manifest and every payload validate.  Ensure a failed write cannot
replace a previously complete generation.

**Step 5: Run GREEN verification and commit**

Run the checkpoint runner at `-np 1,2,4,8`, repeat with a rank-fragment
permutation, run `git diff --check`, and build.  Commit:

```bash
git commit -m "feat(dg): checkpoint localized fragment WF bases"
```

---

### Task 4: Connect SCDM and automatic WF reuse to production

**Files:**
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/dg_hybrid_fragment_wannier.f90`
- Modify: `tests/dg/check_dg_hybrid_fragment_wannier_route.py`
- Create: `tests/dg/check_dg_fragment_wf_restart_route.py`
- Modify: `tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_dg_continuation.in`

**Step 1: Write failing production-route tests**

Require the route to compute the expected checkpoint identity before Wannier90,
make one collective auto hit/miss decision, restore directly on a hit, and on a
miss build the SCDM gauge before invoking Wannier90 and publish only after all
fragment basis checks pass.  Reject any path where ranks disagree or one rank
runs Wannier90 while another restores.

Require user controls for checkpoint policy (`auto`, `read`, `write`, `off`) and
initial gauge (`scdm`, `spectral`, `random`) with validated defaults.  Preserve
the formal fragment-count-equals-rank-count and exact rank-fragment reuse guards.

**Step 2: Run RED verification**

Run the two route checks and input validation tests; expect missing SCDM/cache
wiring.

**Step 3: Add controls and collective decision logic**

Add narrowly scoped DG input controls and validation.  `auto` reuses only a
complete compatible generation and otherwise regenerates collectively.  `read`
fails on any miss; `write` always regenerates then publishes; `off` neither reads
nor writes.  Emit concise hit/miss/publication receipts and exact miss reasons.

**Step 4: Feed SCDM into the existing Gamma A path**

For `scdm`, construct the gauge from the retained DC subspace and pass it through
the current polar/unitarity-checked Wannier90 seed interface.  Preserve spectral
and deterministic random modes for comparison.  Validate the post-Wannier basis
before publication.

**Step 5: Run GREEN verification**

Run all new route/input checks, existing DC seed and Wannier90 route tests, the
new SCDM/checkpoint MPI matrices, affected DG MPI runners, `git diff --check`, and
the release build.

**Step 6: Commit**

Stage only Task 4 hunks and commit:

```bash
git commit -m "feat(dg): reuse fragment WF bases with SCDM seeding"
```

---

### Task 5: Short production smoke tests and failure recovery

**Files:**
- Modify if needed: `tests/dg/run_dg_hybrid_si64_divided_lcfo.py`
- Modify: `tests/dg/check_dg_fragment_wf_restart_route.py`
- Update: `docs/plans/2026-09-07-dg-wf-restart-scdm-design.md`

**Step 1: Exercise production miss/hit without the full Si64 cost**

Use the smallest representative fragment case with fragment count equal to MPI
rank count.  Run once to create the WF generation and again in auto mode.  Prove
from receipts that the second run skips Wannier90 and reproduces the projected
basis fingerprint.  Interrupt a controlled publication before manifest commit
and verify that the next auto run rejects the partial generation and regenerates.

**Step 2: Verify the terminal occupation path**

For the small case, require one terminal LCFO, the target electron count at 300
K, no post-LCFO density update, and a valid occupied checkpoint.  Save logs under
a timestamped `/tmp` directory.

**Step 3: Run regression suite and commit evidence**

Run all affected Python/MPI tests, the release build, and `git diff --check`.
Append commands, fingerprints, cache receipts, timing, and results to the design
document.  Commit only the runner/report changes:

```bash
git commit -m "test(dg): verify fragment WF restart path"
```

---

### Task 6: One clean Si64 SCDM run and exact reuse run

**Files:**
- Create: timestamped run directories under `/tmp`
- Modify only if required for reusable automation: `tests/dg/run_dg_hybrid_si64_divided_lcfo.py`
- Update: `docs/plans/2026-09-07-dg-wf-restart-scdm-design.md`
- Update: `docs/plans/2026-09-06-dg-interface-continuation-design.md`

**Step 1: Validate the exact seed and mapping before launch**

Use the existing Si64 ordinary-DC publication only with eight MPI ranks and the
same rank-fragment mapping.  Record publication ID
`7047888166118007469`, mapping fingerprint `254086644876463474`, and
`scf_skipped=T`.  Abort before Wannier90 on any mismatch.

**Step 2: Run one clean SCDM-seeded generation**

Use OMP threads equal to one for the controlled comparison.  Record per-fragment
Wannier90 iterations and wall times, SCDM pivot/gauge receipts, post-Wannier
validation, cache publication identity, all continuation diagnostics, terminal
LCFO residuals, 256-electron result, and occupied-checkpoint publication.

Do not loosen tolerances to force completion.  Compare the Wannier90 iteration
counts to the prior deterministic-random run: `1342,1192,2640,3058,1927,1837,779,783`.

**Step 3: Run exact automatic reuse**

Repeat with identical inputs, eight ranks, and mapping.  Require a cache hit on
all ranks, no Wannier90 invocation, identical fragment/projected-basis
fingerprints, and tolerance-consistent terminal observables and state
fingerprints.

**Step 4: Interpret the interface diagnostic**

Preserve the measured fixed-density schedule evidence: lambda-zero residual
`1.7943383143e2` and monotonic growth to `2.0819639807e3` at lambda one.  Record
that the interface term is the dominant source of the observed residual growth,
while the corrected terminal occupation result is a separate unit bug.

**Step 5: Final regression and evidence commit**

Run all affected tests and the release build, save complete logs, update both
design/results documents, and commit only report/runner changes:

```bash
git commit -m "test(dg): certify Si64 SCDM and WF reuse"
```

Do not begin Wannier90 OpenMP development in this plan.  That remains a separate
post-implementation optimization step.

