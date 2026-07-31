# Wannier+PW DG-DC and Exp Validation Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** For a gapped LDA system with integer occupations, construct a self-consistent Wannier+PW DG-DC initial state and reproduce the continuous-branch Full TDDFT induced polarization `Delta_Pz(t)` within 5 percent relative RMS using `dt=2 a.u.` exponential propagation.

**Architecture:** Construct a symmetry-validated fixed Wannier basis first: use monolithic SAWF for the small reference and representative-environment SAWF, symmetry replication, symmetry-inequivalent local regeneration, and gauge stitching for large systems. Then resolve the DG trial-space and length-gauge definitions, build shared overlap-metric and mixed-density components used identically by DG-DC and RT, implement fragment-local LDA plus global Hartree DG-DC, serialize an operator-complete checkpoint, and only then implement the production midpoint Exp path. The scalable production basis is distributed `windowed_kg`; dense/global and legacy G-only routes remain small-system oracle, reference, and smoke paths only.

**Tech Stack:** Fortran 2008, MPI, existing SALMON LDA/Hartree infrastructure, LAPACK/EigenExa, Python 3 with NumPy/Matplotlib for tests and analysis, CMake.

---

### Task 0: Define and validate the scalable SAWF basis-generation contract

**Files:**
- Create: `docs/plans/2026-07-12-scalable-sawf-contract.md`
- Create: `src/gs/dc/lcfo_wannier_sawf_templates.f90`
- Modify: `src/gs/dc/lcfo_wannier_sawf.f90`
- Modify: `src/gs/dc/lcfo_wannier_sawf_band.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Create: `tests/dg/check_scalable_sawf_contract.py`
- Create: `tests/dg/test_sawf_template_replication.py`
- Create: `tests/dg/test_sawf_neighbor_gauge_stitching.py`
- Create: `tests/dg/test_sawf_buffer_convergence.py`
- Create: `tests/dg/test_sawf_global_local_equivalence.py`

**Step 1: Write the failing symmetry-contract test**

Require the contract and input route to distinguish:

```text
small-system monolithic global SAWF reference
actual supercell symmetry group versus parent-crystal symmetry
symmetry-closed band and complete projection-shell selection
representative local environments and their symmetry orbits
same-supercell bulk-template reuse versus symmetry-inequivalent local regeneration
buffer definition and convergence tolerances
neighbor gauge-stitching convention
template provenance and cache invalidation
```

The SAWF acceptance gate must validate the group laws and covariance of
`D_band` and `D_wann`; `site_symmetry=true` alone must not pass the test.

**Step 2: Write failing representative-replication tests**

Construct a deterministic symmetric toy crystal. Generate one representative
basis and require translated/rotated copies to agree with direct construction:

```text
w_(gR) = g w_R D_wann(g)
S_(gR,gR') = D_wann(g)^H S_(R,R') D_wann(g)
H_(gR,gR') = D_wann(g)^H H_(R,R') D_wann(g)
```

Reject incomplete orbital shells, non-closed band subsets, invalid atom maps,
and a parent-crystal operation that is not a symmetry of the actual supercell.

**Step 3: Implement environment orbits and fingerprinted templates**

Identify symmetry-equivalent fragment environments under the actual supercell
group. Run SAWF only for one representative of each bulk orbit and serialize
the local orbitals plus `D_band`, `D_wann`, centers, spreads, and fingerprints
for geometry, pseudopotential, grid, band window, projection shells, symmetry
operations, buffer, and code/schema version. A fingerprint mismatch must force
regeneration. Do not silently fall back to an identity-only group.

**Step 4: Implement symmetry-inequivalent local-environment regeneration**

Fingerprint every core+buffer local environment in the actual supercell.
Treat defects, material interfaces, free surfaces/vacuum boundaries, and
amorphous neighborhoods as independent whenever their exact fingerprints are
not connected by an operation of the actual supercell symmetry group. Reuse
bulk templates only inside the same exact supercell; cross-supercell reuse is
forbidden. Local symmetry breaking is retained rather than projected back to a
parent crystal. In amorphous regions, do not merge approximately similar local
geometries: independently regenerate each distinct environment.

**Step 5: Write and implement the neighbor gauge-stitching test**

For every neighboring representative pair, compute the buffered-orbital
overlap and solve the unitary Procrustes/polar-decomposition alignment. Apply
the same unitary to all basis-dependent quantities before building WP and DG
face blocks. Reject rank-deficient overlaps, ambiguous orbital counts, or a
post-alignment residual above tolerance. Record the unitary and residual in the
template checkpoint.

**Step 6: Establish buffer convergence**

For at least three increasing buffers, compare centers, occupied projector,
neighbor overlaps, WW/WP blocks, and every DG face block. The largest two must
meet documented tolerances. A basis that is localized but whose face blocks
are not converged does not pass.

**Step 7: Prove local/global equivalence on the small validation system**

Run both the monolithic global SAWF reference and hierarchical construction.
After one global unitary alignment require agreement of the occupied projector,
`S`, fixed `H_kin+DG+V_ion`, and all face blocks. This gate must pass before the
hierarchical basis is accepted as DG-DC input.

**Step 8: Run focused tests, build, and commit**

Run:

```bash
python3 tests/dg/check_scalable_sawf_contract.py
python3 tests/dg/test_sawf_template_replication.py
python3 tests/dg/test_sawf_neighbor_gauge_stitching.py
python3 tests/dg/test_sawf_buffer_convergence.py
python3 tests/dg/test_sawf_global_local_equivalence.py
cmake --build build-mpi-eigenexa -j 2
```

Expected: all focused tests pass; the small system agrees with global SAWF and
defect/interface/surface/amorphous fixtures reject parent-group symmetry restoration.

```bash
git add docs/plans/2026-07-12-scalable-sawf-contract.md src/gs/dc tests/dg src/io
git commit -m "feat: add scalable symmetry-adapted Wannier construction"
```

### Task 1: Freeze the mathematical operator contract

**Files:**
- Create: `docs/plans/2026-07-12-wannier-pw-dg-operator-contract.md`
- Create: `tests/dg/check_wpw_operator_contract.py`
- Inspect: `src/gs/dc/lcfo_flux.f90`
- Inspect: `src/rt/dg/rt_dg_plane_wave.f90`
- Inspect: `src/rt/dg/rt_dg_fragment_hamiltonian.f90`

**Step 1: Write the failing contract test**

Require the contract document to define:

```text
fragment core domains and unique grid ownership
support of fragment-restricted Wannier functions
support and normalization of the PW enrichment sector
jump/average conventions and face orientation
surface-penalty scaling
which H and S blocks exist: WW, WP, PW, face-neighbor
global ownership versus physical fragment coupling
periodic length-gauge operator and polarization branch
whether any position-interface correction follows from the chosen discretization
```

**Step 2: Derive the position operator before fixing tests**

Start from the selected DG weak form and multiplicative position operator. Check
the discrete relation between the candidate `Z`, DG Hamiltonian, and velocity:

```text
v_z versus i[H_DG,Z]_S
```

Do not require a nonzero interface position term unless this derivation and a
small numerical model demonstrate it. Record the accepted definition and the
rejected alternatives.

**Step 3: Run the contract test**

Run `python3 tests/dg/check_wpw_operator_contract.py`.
Expected: PASS only when every required decision is explicit.

**Step 4: Commit**

```bash
git add docs/plans/2026-07-12-wannier-pw-dg-operator-contract.md tests/dg/check_wpw_operator_contract.py
git commit -m "docs: define Wannier PW DG operator contract"
```

### Task 2: Extract shared generalized-eigen and S-metric Exp algebra

**Files:**
- Create: `src/common/dg_generalized_algebra.f90`
- Modify: `src/CMakeLists.txt`
- Create: `tests/dg/fixtures/wpw_generalized_driver.F90`
- Create: `tests/dg/test_wpw_generalized_algebra.py`

**Step 1: Write failing numerical tests**

Use small deterministic complex Hermitian `H`, positive-definite `S`, and
reference NumPy solutions. Test:

```text
H C = S C epsilon
C^H S C = I
U^H S U = S
C(dt) = X exp(-i X^H H X dt) X^H S C(0)
```

Include negative tests for Euclidean propagation with `S /= I`, missing DG face
blocks, and metric rank loss.

**Step 2: Implement the common module**

Provide metric filtering, `S^-1/2`, generalized diagonalization, forward/backward
metric transforms, residuals, condition number, discarded-direction count, and
an S-unitary exponential action. The module must have no dependency on GS or RT
modules/data types so DG-DC and RT use the identical implementation. Use the
project's low-level LAPACK interface directly or move the required wrapper into
the common layer; do not import `eigen_subdiag_sub` from `src/gs`.

**Step 3: Test and commit**

Run `python3 -m unittest tests.dg.test_wpw_generalized_algebra -v` and the focused
Fortran build. Commit only the common module, build wiring, and tests.

### Task 3: Extract a shared complete mixed-density builder

**Files:**
- Create: `src/common/dg_wpw_density.f90`
- Modify: `src/rt/dg/rt_dg_density_reconstruct.f90`
- Modify: `src/CMakeLists.txt`
- Create: `tests/dg/test_wpw_mixed_density.py`
- Create: `tests/dg/fixtures/wpw_density_driver.F90`

**Step 1: Write failing numerical tests**

Require direct real-space density and decomposed density to agree:

```text
n_direct = n_WW + n_WP + n_PW
integral(n) = Tr(S P) = sum_i f_i = N_e
```

Add negative tests omitting WP or PP terms and tests with a nonidentity overlap.

**Step 2: Extract and reuse the production builder**

Move the complete mixed-density calculation out of the RT-only control flow.
Keep coefficient layout adapters in RT, but use the same density kernel in
DG-DC and RT. The common API accepts plain arrays and metadata and must not
depend on `s_dg_fragment_rt`; GS and RT provide thin adapters.

**Step 3: Test and commit**

Run the numerical tests plus the existing RT density checks. Commit the common
density module and RT adapter together.

### Task 4: Validate fragment-core LDA and global Hartree plumbing

**Files:**
- Create: `src/gs/dc/dg_wpw_lda_hartree.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/check_dg_wpw_lda_hartree.py`

**Step 1: Write failing ownership and energy tests**

On a small partitioned grid, require:

```text
every global point has exactly one core owner
buffer and halo points contribute zero to energy integration
sum_fragment E_xc_LDA = E_xc_global_grid
sum_fragment integral(n V_xc) = global-grid result
E_H = 0.5 integral_global(n V_H)
```

**Step 2: Implement using existing SALMON routines**

Call the existing LDA exchange-correlation implementation on owned core data and
the existing global Poisson/Hartree infrastructure on the assembled complete
mixed density. Do not create a new XC functional or fragment Hartree
approximation.

**Step 3: Test and commit**

Run the focused test and build. Commit this layer before adding SCF control.

### Task 5: Implement fixed-basis self-consistent DG-DC

**Files:**
- Create: `src/gs/dc/dg_wannier_pw_scf.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/main_dft.f90`
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Create: `tests/dg/check_dg_wpw_scf_route.py`
- Create: `tests/dg/test_dg_wpw_scf_fixed_point.py`

**Scalability correction:**

The dense solver in `dg_wannier_pw_scf.f90` is a small-system mathematical
oracle only.  Production must never assemble a global dense `H`, `S`, or
density matrix and must never compute the full eigenspectrum.  Production
uses block/neighbor storage with rank-local fragment and PW coefficient rows
and face-neighbor halo exchange.  The Task-2 transitional layout in
`prepare_local_fragment_pw_blocks` (every rank requests every PW row) and the
optional dense `H_mat_pw(n_pw,n_pw,nspin)` storage are explicitly not a
production matrix-free implementation.

**Step 1: Freeze the matrix-free production contract**

Create `tests/dg/check_dg_wpw_scf_route.py` first.  Require:

```text
dense reference solver is not called from main_dft
global dense H/S allocation is forbidden in the production route
H and S are consumed only through batched apply callbacks
only occupied states plus a bounded buffer are retained
rank-local invalid data is diagnosed before collective reductions
complete density is distributed to the existing global Poisson solver
```

Keep `tests/dg/test_dg_wpw_scf_fixed_point.py` as the independent NumPy-backed
small-system oracle for the equations and energy expression.

**Step 2: Complete distributed PW-row ownership**

Create `src/rt/dg/rt_dg_wpw_column_layout.f90` and
`tests/dg/test_dg_wpw_column_layout.f90`.  Define
`s_dg_wpw_column_layout` with `basis_kind`, `n_global_columns`, `n_g_modes`,
`pw_fragment_ids`, `pw_g_ids`, `pw_owner`, and `owned_column_ids`.  Use
fragment-major ids `column_id=(K-1)*n_G+G_id`.  Reject nonpositive dimensions,
invalid ranks, overflow of the default integer id space, duplicate/missing
pairs, and any basis kind other than `windowed_kg`.

Run the layout fixture with one and multiple ranks.  Require a bijection
between ids and `(K,G)`, deterministic ownership, bounded owned-column count,
and an exact inverse map.  This layout is separate from legacy G-only
`n_plane_waves/k_pw`; do not add `(K,G)` fields to `s_dg_fragment_rt`.

Represent every enrichment column by a stable pair `(K,G)`, with explicit
`pw_fragment_ids` and `pw_g_ids`; do not collapse the K index by summing all
window contributions into one G-only column.  The accepted Task-1 basis is the
direct sum `P_(K,G)=chi_K exp(iG.r)/sqrt(Omega_cell)`.  Existing diagnostic
`compute_wpw_overlap/compute_wpw_kinetic_weak` routines that accumulate every
fragment window into a single `n_pw by n_pw` G-only matrix are not the
production operator builder.

Replace the Task-2 transitional all-PW request list by bounded rank-local PW
ownership plus the remote rows required by local FP/PP action.  The PP apply
must use diagonal or explicitly sparse/block-distributed storage; production
must fail closed rather than allocate or consume a dense global `H_mat_pw`.
No rank may retain `O(n_pw)` rows merely because PW enrichment is enabled.

Store each sparse WP entry on the owner of its PW column.  That owner
halo-fetches only support-coupled W coefficients, computes its owned PW output,
and sends partial W-output sums to the unique W-row owners.  Store PP by PW
output-row owner and halo-fetch only support-neighbor PW coefficients.  This
layout reproduces both Hermitian WP/PW directions without replicating all PW
columns or all fragment rows.  `tests/dg/test_dg_wpw_distributed_block_action.py`
is the ownership-independent dense-equivalence oracle for this decomposition.

**Step 3: Extract a GS/RT-neutral block operator adapter**

Expose a plain callback contract for batched `Y=HX` and `Y=SX`.  The RT
adapter must use compact WW blocks, row-local FP blocks, distributed PP action,
and a row-local mixed-overlap apply; it must not call the global-replicated
`apply_overlap_operator_batch` or copy block data into dense global arrays.
Add a small-system equivalence fixture comparing callback action against the
dense oracle for WW, WP, PW, and face-neighbor blocks.

**Step 4: Implement distributed occupied-subspace iteration**

Use block subspace iteration with S-orthonormalization and Rayleigh--Ritz only
inside the bounded occupied trial subspace.  The retained dimension is
`n_occ + dg_wpw_scf_extra_states`; it must not scale to the full basis size.
The local coefficient storage follows existing DG row ownership.  Reject a
rank-deficient metric, a missing unoccupied buffer state, or a gap below the
positive namelist threshold.

**Step 5: Write failing route and fixed-point tests**

Require gapped LDA, integer occupations, fixed Wannier+PW basis, fixed DG
kinetic/surface blocks, updated `V_H+V_xc`, and convergence of:

```text
||n_out-n_in|| / ||n_out||
||V_out-V_in|| / ||V_out||
|E_tot(k)-E_tot(k-1)|
||Q_occ(k)-Q_occ(k-1)|| / ||Q_occ(k)||
max_i ||H C_i-S C_i epsilon_i||
||C^H S C-I||
Tr(SP)-N_e
```

**Step 6: Implement the matrix-free SCF loop**

Use the conventional GS only for the initial density and fixed basis. Reuse
Task 2 algebra inside the bounded Rayleigh--Ritz subspace, Task 3 rank-local
density reconstruction, Task 4 distributed LDA/global-Poisson plumbing, and
the existing block DG matrix builders. Mix the potential using SALMON's simple
potential policy. Reject metallic/smeared occupation in this milestone.

Use the conventional eigenvalue-sum energy:

```text
E_tot = sum_i f_i epsilon_i - E_H - integral(n V_xc)
        + E_xc_LDA + E_ion_ion
```

DG kinetic/surface/penalty terms are already in the eigenvalues and must not be
added again.

**Step 7: Run field-free fixed-point and scaling tests and commit**

The converged DG density passed through one additional SCF map must reproduce
itself.  A scaling fixture must increase the fragment count while asserting
that no allocation proportional to `n_basis**2` is made and that each apply
touches only local and face-neighbor blocks. Commit the solver, input contract,
and tests.

### Task 6: Define and implement the DG-DC checkpoint contract

**Files:**
- Create: `src/common/dg_wpw_checkpoint.f90`
- Modify: `src/gs/dc/dg_wannier_pw_scf.f90`
- Modify: `src/rt/dg/rt_dg_fragment_io.f90`
- Create: `tests/dg/test_dg_wpw_checkpoint_roundtrip.py`

**Step 1: Write the failing round-trip and mismatch tests**

The checkpoint must store:

```text
basis/fragment/PW fingerprints
S and fixed H_kin+DG+V_ion
converged V_H+V_xc matrix
N_basis, N_eff, N_eig, N_occ, occupation weights
metric cutoff and discarded directions
eigenvalues and all retained eigenvectors
penalty and face conventions
matrix checksums and provenance
```

Require `N_eig=N_eff` for the initial production validation. Corrupt fingerprints,
dimensions, checksums, or operator conventions and verify fatal rejection.

**Step 2: Implement one shared reader/writer**

DG-DC and RT must use the same serialization module and versioned schema.

**Step 3: Test and commit**

Run a write/read/write round trip and require identical metadata and matrices.

### Task 7: Prove DG-DC to RT operator handoff identity

**Files:**
- Modify: `src/rt/dg/rt_dg_fragment.f90`
- Modify: `src/rt/dg/rt_dg_observables.f90`
- Create: `tests/dg/check_dg_dc_rt_handoff.py`

**Step 1: Write the failing handoff test**

At RT initialization, before propagation, require:

```text
||H_RT(0)-H_DC||_F / ||H_DC||_F
||S_RT-S_DC||_F / ||S_DC||_F
occupied generalized residual
```

Checks occur before collective reductions where rank-local corruption can be
identified.

**Step 2: Implement checkpoint-backed initialization**

Do not project conventional orbitals or reselect the metric in RT. Load the
converged DG-DC state and verify reconstructed operators against checkpoint
checksums.

**Step 3: Run the static handoff test and commit**

Do not implement a field-on route until the static operator handoff gate passes.

### Task 8: Implement one namelist-driven production midpoint Exp route

**Files:**
- Create: `src/rt/dg/rt_dg_wpw_exp_production.f90`
- Modify: `src/rt/dg/CMakeLists.txt`
- Modify: `src/rt/dg/rt_dg_integrator_expdiag.f90`
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Create: `tests/dg/check_wpw_exp_production_route.py`
- Create: `tests/dg/test_wpw_exp_midpoint.py`

**Step 1: Audit scientific environment variables**

Classify every WPW/BPW/MIXED-Z/flux environment variable as result-changing or
diagnostic-only. Result-changing controls used by the accepted route must become
namelist variables or be removed. Diagnostic-only variables may remain clearly
labeled.

**Step 2: Write failing propagation tests**

Require a single production route using Task 2 S-metric algebra. Propagate
`N_occ` orbitals in the complete `N_eig=N_eff` response space. Every midpoint
corrector candidate must restart from the saved `C_n`; a cumulative-corrector
negative test must fail.

**Step 3: Implement midpoint TDDFT**

Keep the basis and DG surface matrix fixed. Iterate midpoint density,
`V_H+V_xc`, and potential-dependent matrix blocks to tolerance. Stop on
non-convergence or S-norm failure. Do not silently fall back to an experimental
local/split route.

**Step 4: Test, build, and commit**

Run Tasks 2, 3, 7, and 8 focused tests plus
`cmake --build build-mpi-eigenexa -j 2`.

Then run a field-off propagation gate from the converged checkpoint and require
density, occupied projector, S-norm, and `Delta_Pz` to remain stationary. Do not
continue to the length-gauge field-on task until this gate passes.

### Task 9: Validate the accepted length-gauge observable

**Files:**
- Modify: `src/rt/dg/rt_dg_fragment.f90`
- Modify: `src/rt/dg/rt_dg_fragment_ops.f90`
- Modify: `src/rt/dg/rt_dg_observables.f90`
- Modify: `src/io/write.f90`
- Create: `tests/dg/check_wpw_length_gauge_observable.py`

**Step 1: Implement only the Task 1 accepted operator**

Form the accepted WW/WP/PP blocks and only those interface corrections proven
by the operator contract. Report the unsymmetrized metric-Hermiticity residual
before any numerical symmetrization.

**Step 2: Test gauge and branch consistency**

Check field-sign oddness, `v_z` versus the accepted discrete commutator relation,
continuous polarization branch tracking, and consistency of
`Jz=d Delta_Pz/dt`.

**Step 3: Commit**

Commit the operator, observable output, and tests together.

### Task 10: Create the Stage 2D input and provenance matrix

**Files:**
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_stage2d_full_tddft_i1e11_dt2`
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_stage2d_dg_dc_wpw_pw16`
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_stage2d_dg_dc_wpw_pw64`
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_stage2d_dg_dc_wpw_pw128`
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_stage2d_global_wpw_exp_pw16_i1e11_dt2`
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_stage2d_global_wpw_exp_pw64_i1e11_dt2`
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_stage2d_global_wpw_exp_pw128_i1e11_dt2`
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/stage2d_wpw_exp_manifest.tsv`
- Create: `tests/dg/check_stage2d_wpw_exp_inputs.py`

**Step 1: Write the seven-input consistency test**

Require common geometry, conventional-GS provenance, and physical laser
parameters. For each DG-DC/RT pair require identical checkpoint ID, basis
fingerprint, PW selection, effective dimension, occupation, and metric policy.

Manifest columns include:

```text
path dg_dc_path checkpoint_id basis_id WF_count PW_count PW_cutoff_or_shell
N_basis N_eff N_eig N_occ occupation_policy S_eval_min S_cutoff S_condition
discarded_metric_dirs dt propagator_kind observable_source gauge
polarization_branch volume_normalization
```

**Step 2: Add inputs and smoke tests**

Run each DG-DC to convergence before its RT input. Run field-off RT first, then a
20-step field-on smoke test.

**Step 3: Commit**

Commit only the seven inputs, manifest, and parser test.

### Task 11: Implement Delta-Pz comparison and PW convergence

**Files:**
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/compare_stage2d_wpw_exp.py`
- Create: `tests/dg/test_compare_stage2d_wpw_exp.py`
- Create: `tests/dg/fixtures/stage2d_full_pz.data`
- Create: `tests/dg/fixtures/stage2d_wpw_pz.data`

**Step 1: Write failing analysis tests**

Require identical sampling, no extrapolation, a resolved continuous branch,
`Delta_Pz=Pz-Pz(0)`, identical differentiation for both `Jz`, and

```text
rel_rms = rms(Delta_Pz_WPW-Delta_Pz_full) / rms(Delta_Pz_full)
```

Reject unresolved polarization-quantum jumps. `Jz` remains secondary and cannot
replace the `Delta_Pz <= 0.05` gate.

**Step 2: Implement and test**

Write summary, aligned waveform, and plot outputs. Run
`python3 -m unittest tests.dg.test_compare_stage2d_wpw_exp -v`.

**Step 3: Run PW16/PW64/PW128 convergence**

At fixed `dt=2`, require at least one basis to pass `Delta_Pz_rel_rms <= 0.05`.
If none passes, stop and classify basis, operator, midpoint, or provenance error.

### Task 12: Final regression, documentation, and distributed handoff

**Files:**
- Modify: `doc/NOTE_DG.md`
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/stage2d_wpw_exp_summary.tsv`

**Step 1: Run the complete focused suite**

Run every new Task 1-11 test, relevant existing DG tests, and
`cmake --build build-mpi-eigenexa -j 2`.

**Step 2: Document the accepted configuration**

Record operator contract, DG-DC convergence, checkpoint ID, metric metadata,
accepted basis, `Delta_Pz_rel_rms`, `Jz_rel_rms`, and limitations.

**Step 3: Commit the validated summary**

Complete the scalable `windowed_kg` validation using the matrix-free distributed
ownership contract established in Task 5. Dense/global and legacy G-only runs
may validate small-system physics but cannot serve as production scaling evidence.
