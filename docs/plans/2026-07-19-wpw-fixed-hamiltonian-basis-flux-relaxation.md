# WPW Fixed-Hamiltonian Basis/Flux Relaxation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Produce a qualified Si64 W+P checkpoint by projecting the converged occupied subspace and continuing the complete DG interface operator while keeping the DC density and H0 fixed.

**Architecture:** Project the converged DC+LCFO occupied subspace into W+P, qualify the raw projection at zero interface strength, then continue the complete DG interface bilinear form from zero to one. Only the occupied subspace and explicitly mapped face blocks may change; potential rebuilding and density feedback are forbidden until a later opt-in SCF stage.

**Tech Stack:** Fortran 2008, MPI, bounded rank-local WW/WP/PP operators, matrix-free generalized CG, Python source contracts, MPI dense oracles.

---

### Task 1: Add fixed-H mode and frozen-state contract

**Files:**
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Create: `tests/dg/check_dg_wpw_fixed_h_relaxation.py`

1. Write a failing source contract requiring an explicit fixed-H production mode, frozen-field snapshots, and a prohibition on `wpw_potential_step` in that mode.
2. Run the contract and confirm RED.
3. Add controls for enabling fixed-H relaxation and complete-interface continuation, without changing repository input defaults.
4. Snapshot density, Hartree, XC, local potential, and frozen operator fingerprints after converged DC+LCFO publication.
5. Add collective invariant validation and transactional failure paths.
6. Run the contract and focused publication tests; expect PASS.

### Task 2: Make the operator decomposition explicit

**Files:**
- Modify: `src/common/dg_wpw_production_context.f90`
- Modify: `src/rt/dg/rt_dg_wpw_production_builder.f90`
- Modify: `src/common/dg_wpw_bounded_operator.f90`
- Test: `tests/dg/run_dg_wpw_production_operator_mpi.py`

1. Extend the MPI oracle to require separately addressable frozen-volume/nonlocal and complete DG-interface contributions.
2. Confirm RED before changing production code.
3. Store immutable H0 blocks and separately scaled interface blocks in the candidate operator.
4. Rebuild only the combined apply cache when continuation lambdas change.
5. Validate unchanged H0 fingerprints across continuation steps and prove that support/halo transport IDs are never scaled.
6. Compare decomposed and combined H/S application against the dense oracle; expect PASS.

### Task 3: Project the converged occupied subspace

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/gs/dc/dg_wpw_lcfo_operator_adapter.f90`
- Modify: `src/gs/dc/dg_wpw_matrix_free_scf.f90`
- Test: `tests/dg/test_dg_wpw_matrix_free_scf_mpi.f90`

1. Add a failing MPI oracle that projects a known occupied subspace into a larger nonorthogonal W+P basis and compares its S-metric projector, not individual coefficients.
2. Confirm RED.
3. Export the converged DC+LCFO occupied subspace and occupations before total-state storage is released.
4. Form W/P projection coefficients using the production metric and rank-local ownership.
5. Apply rank-revealing S-orthonormalization; reject loss of occupied rank.
6. Construct deterministic extra-state complements only after S-projecting out the occupied subspace.
7. Record projection loss, occupied-projector defect, charge, and reconstructed-density baseline.
8. Add a degenerate-subspace oracle proving invariance under unitary occupied rotations.

### Task 4: Implement zero-interface representation validation

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/gs/dc/dg_wpw_matrix_free_scf.f90`
- Test: `tests/dg/test_dg_wpw_matrix_free_scf_mpi.f90`

1. Add an MPI test that starts from the projected occupied subspace and solves a fixed generalized problem without invoking a potential map.
2. Confirm RED.
3. Add a fixed-H algebra driver whose only mutable state is the eigenspace/CG history.
4. Start with `lambda_interface=0`; forbid random replacement of the projected occupied subspace.
5. Qualify the unrelaxed projection before permitting any CG step.
6. Reconstruct density as a diagnostic and compare charge/density against the frozen reference without publishing it as a potential input.
7. Report projection error and subsequent fixed-H variational change separately.
8. Require finite eigenvalues, bounded residual, S-orthogonality, occupied-projector continuity, and unchanged frozen fingerprints.

### Task 5: Implement cooled orbital CG

**Files:**
- Modify: `src/gs/dc/dg_wpw_matrix_free_scf.f90`
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Test: `tests/dg/test_dg_wpw_matrix_free_scf_mpi.f90`

1. Add a wide-spectrum, multi-state failing oracle requiring accepted decrease of an occupied-subspace merit with trial rejection.
2. Add explicit orbital mix, minimum mix, cooling factor, growth factor, and accepted-growth interval controls.
3. Form an S-tangent block-CG direction, reject non-descent directions, and restart to steepest descent.
4. Trial `Q + alpha P`, S-orthonormalize, and accept only when residual/Rayleigh safeguards pass.
5. Accept using occupied trace Rayleigh quotient, maximum occupied residual, S-orthogonality, and projector-change safeguards; do not compare degenerate states coefficient by coefficient.
6. Lock converged degenerate clusters as invariant subspaces.
7. On rejection restore the last accepted occupied subspace and reduce alpha.
8. Fail closed when alpha or rejection limits are exhausted.
9. Run narrow-, wide-spectrum, and degenerate-cluster MPI oracles; expect PASS.

### Task 6: Add DG-interface continuation

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/rt/dg/rt_dg_wpw_production_builder.f90`
- Test: `tests/dg/run_dg_wpw_production_operator_mpi.py`

1. Add a failing oracle for `lambda_interface` zero, intermediate lambda, rejection, rollback, and one.
2. First publish an explicit mapping test: `ww_face_self`, `ww_cross_value`, and `wp_h_face` are scaled; all volume/nonlocal blocks and transport IDs are unchanged.
3. Prove that the mapped interface blocks at lambda one reproduce the pre-existing combined operator exactly.
4. Reuse the last accepted occupied subspace as the next initial state.
5. Adapt the continuation step independently from orbital mix.
6. Reject and cool on subspace merit, orthogonality, Hermitian defect, charge, or frozen-field violations.
7. Record H0/interface Rayleigh contributions and frozen-state fingerprints at every accepted step.
8. Require lambda one before checkpoint publication.

### Task 7: Qualify and publish the fixed-H checkpoint

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/rt/dg/rt_dg_wpw_checkpoint.f90`
- Modify: `tests/dg/check_dg_wpw_checkpoint_publication.py`

1. Extend the checkpoint test to require fixed-H mode, frozen-state fingerprints, final interface lambda, both tolerance profiles, projection baseline, and decomposed operator provenance.
2. Confirm RED.
3. Add schema fields without removing existing sparse operator/checkpoint data.
4. Publish only after `lambda_interface=1` and all qualification checks pass.
5. Label the checkpoint as fixed-H basis/flux-relaxed, not self-consistent WPW.

### Task 8: Si64 staged production validation

**Files:**
- Modify: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_stage2d_dg_dc_wpw_smoke`
- Modify: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/stage2d_wpw_exp_manifest.tsv`
- Modify: `tests/dg/check_stage2d_wpw_exp_inputs.py`
- Create: new run directories below `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/stage2d_wpw_runs/`

1. Record fixed-H mode, projection provenance, both tolerance profiles, orbital cooling, and interface-continuation controls in the manifest and consistency test before running.
2. Run a raw-projection Si64 preflight and record rank loss, projector defect, charge, and density baseline.
3. Run zero-interface fixed-H relaxation and keep its density change separate from projection error.
4. Run DG-interface continuation to one; stop on any frozen-field mutation or nonfinite result.
5. Require finite charge/eigenvalues, converged maximum occupied residual/S-orthogonality, bounded occupied-projector change, declared density error relative to the projection baseline, and identical frozen provenance.
6. Publish the checkpoint only after qualification.
7. Do not run field-off RT until the checkpoint review passes.

### Task 9: Focused verification and findings-first review

1. Run the fixed-H source contract.
2. Run matrix-free narrow/wide-spectrum MPI tests.
3. Run production operator, projector, owner-exchange, and face-trace MPI tests.
4. Run checkpoint publication/round-trip tests.
5. Build `build-mpi-eigenexa`.
6. Run `git diff --check`.
7. Review for accidental density feedback, H0 mutation, global dense allocation, collective mismatch, rollback gaps, and checkpoint ambiguity.

No commit, push, PR, field-off RT, laser run, or HHG interpretation is authorized by this plan.
