# LCFO HSE real-time implementation plan

**Goal:** Propagate the Si128 system in a fixed complex LCFO basis, retain sixteen fragment-local exchange evaluations, and compare dielectric response versus Wannier axial support.

**Architecture:** Separate a verified finite-basis propagation kernel from density/Hamiltonian assembly and the electromagnetic response. Import the actual complex LCFO basis and Hamiltonian rather than evolving sixteen isolated periodic fragments. Preserve an explicit full-support reference. No frozen-Hamiltonian trace is a self-consistent HSE dielectric spectrum.

**Tech stack:** Fortran complex128, LAPACK/BLAS, existing SALMON MPI and HSE/MLWF machinery; Python/numpy verification and analysis.

## Constraints
One numerical job or heavy analysis at a time; MPI16/OMP1/BLAS1 for the sixteen-fragment reference. Existing current branch is reused. Source-only cutoffs cannot silently be labeled variational TDHSE. Norm conservation alone is insufficient: check stationary density, energy, Hermiticity, field/current consistency and time-step convergence.

## Tasks
1. Generate the converged Si128 LCFO basis with all available fragment states (energy_cut=100 Ha, lambda_cut=1e-7). The existing DC restart is unsupported; rerun one bounded SCF+LCFO job. Independently read its versioned files and check basis orthogonality and eigen residuals. Keep original reference intact.
2. Implement a standalone orthonormal complex LCFO propagation kernel in src/rt/lcfo_rt_core.f90. Tests in testsuites/unit_lcfo_rt: compare Cayley propagation against exact eigenphases; verify preservation of complex off-diagonal coherence, norm, reversibility, stationary fractional density, and timestep convergence; reject non-Hermitian/nonfinite inputs and negative occupations. Unit tests precede implementation. No post-step Gram-Schmidt or eigenvalue clipping.
3. Implement strict LCFO file loading and reconstruction in samples/dc_hse/lcfo_rt_reference.py. Test damaged/truncated/stale files and actual matrix assembly against SALMON output. Use the reference for frozen-H verification only, clearly labeled as such.
4. Derive and implement self-consistent projected Hamiltonian updates: reconstruct global density and fragment density matrices from global coefficient occupations; Hartree/semilocal updates on consistent global density; common local exchange actions and their Hermitian global assembly. Preserve baseline and record initial DC-to-LCFO density discrepancy diagnostically. User explicitly accepts this discrepancy: propagate the saved LCFO state directly; do NOT require reconvergence or zero-field stationarity before starting. Establish an energy functional consistent with any cutoff before claiming TDHSE.
5. Introduce transported MLWF support and ACE in the projected update; direct-full/ACE parity on occupied columns, gauge tests, fractional occupation tests, full-radius identity and cutoff Hermiticity tests. Integrate input checks, U/cache reset and restart provenance. No relaxation of DC input guards before this is verified.
6. Couple weak transverse impulse and current consistently to the projected nonlocal Hamiltonian. Verify zero-field stationarity and finite-difference field derivative. Compare dt/2 and half kick before long propagation.
7. Run one-at-a-time full, 9, 8, 7 bohr axial support cases; use identical propagation length/window/frequency mesh; compare Re/Im epsilon_xx, peaks, weights, drift and wall time. Distinguish initial-state and propagation approximations; do not claim convergence of physical bulk Si k sampling from this supercell test.

## Progress
- Design approved by user (2026-09-26).
- Task 1 complete: si128-chain/lcfo-reference, MPI16/OMP1/BLAS1, 1148 SCF iterations, 1015.55 s wall time. All 16 fragments retain 64 core basis functions; global LCFO dimension 1024, 512 saved states. Independent eigen residual 1.0866e-14 Ha, basis orthogonality 1.5543e-15.
- Task 2 core implemented and standalone probe passes analytic phases, complex-density feedback, time-dependent callback, conservation, reversal, second-order convergence and rejection checks. Grid-density/potential projection expectation parity passes.
- Task 3 reader implemented: 9 tests pass for valid complex data, corrupted payload/footer/basis/coefficients, stale runs, native final symmetrization and three-fragment periodic halo geometry. Real Si128 reader validation also passes.
- Read-only review found missing native final symmetrization and missing geometric halo checks; both reproduced by failing tests then fixed (9/9 reader tests pass). No integration/TDHSE correctness is inferred from algebra tests.
- SALMON native build and CTest lcfo_rt_core pass. Actual Si128 frozen-H probe (4 forward + 4 backward steps, dt=0.02 au) passes: Cayley eigenphase error 6.6743e-16, orthogonality 9.3259e-15, reverse error 1.0550e-15. This is not self-consistent RT. Tasks 4–7 remain pending apart from local density/potential projection algebra. No RT spectrum exists yet.

- User steering: DC-to-LCFO density mismatch is accepted. Use the same saved initial LCFO state for all support cases; record zero-field response as a control, not an admission criterion. Remaining implementation gates concern the missing self-consistent Hamiltonian and field/current connection, not initial-density mismatch.
