# Symmetry plus shared response space implementation plan

**Goal:** Combine validated spatial symmetry with fixed-target finite-k response compression and compare with the full selected-q calculation.

**Architecture:** An opt-in diagnostic uses the little group of a fixed target k. Begin with symmorphic positive fractional-axis permutations, the exact grid operations already supported by EXX, preserving the cell metric, grid, k/q meshes and species-labelled ionic positions. Check occupied projectors, energy intertwining and actual Hpsi covariance; never copy individual orbital gauges. Compute only representative q responses in the common response space and reconstruct full polarization blocks by reciprocal G permutations. Independently computed all-q CG blocks validate symmetry-only restoration before reporting combined errors. The existing production W writer remains disabled for diagnostics.

**Tech Stack:** Native Fortran, MPI owner-only diagnostic and collective failure handling, existing Hpsi, BLAS/LAPACK. No production Python, remote job changes or fixed-binary replacement.

1. Add failing native tests for little-group q/G maps, degenerate rotated occupied gauges, rejection of broken atomic/Hamiltonian symmetry, and agreement of reconstructed Fourier-model polarization with independent all-q reference.
2. Implement target-local permutation discovery and q-star mapping in `hsex_response_symmetry.f90`, with finite checks and explicit covariance diagnostics. Fractional translations, sign changes, time reversal, general rotations and global irreducible-k GS reconstruction are outside this first step.
3. Add `yn_hsex_reduced_symmetry='n'` to functional input/default/broadcast/log/checks. `y` requires fixed-target diagnostic. Pass ionic coordinates/species through the adapter.
4. Compact representative RHS, call the existing response-space solver, reconstruct target polarization blocks including G relabeling. Compare symmetry-only exact CG restoration and combined approximate restoration. Representative solution residuals remain separate from full reconstructed observable errors.
5. Source/object eigensolver audits before new executable use. Test1/4MPI×12OMP, broken-symmetry cases and the existing off-mode regression. Use Si Nr12 Nk=Nq4, qcount6, rank256/384 on scratch; report total setup/validation and no-W-publication. Preserve prior binaries/logs.
6. Update manuals and test registry; save derived records, filtered screening diagnostics, inputs and code under canonical calc2026, without GS payloads or GS logs. No commit unless requested.
