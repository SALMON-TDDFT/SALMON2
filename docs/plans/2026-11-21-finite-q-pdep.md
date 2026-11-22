# Finite-q PDEP Implementation Plan

**Goal:** Measure PDEP including construction against the existing empty-state-free static Si Sternheimer response, at one finite q. User-provided revised specification authorizes this experiment.

**Architecture:** Keep production screening unchanged. Add an opt-in equilibrium scalar diagnostic in the existing Fortran SALMON adapter, using its real-space Hpsi callback and actual occupied states. Reference forms the full finite-q Coulomb-symmetrized polarization with the existing Sternheimer band routine. PDEP constructs a Krylov/Ritz subspace through polarization actions, never by diagonalizing the reference matrix. Run reference and ranks in separate processes to obtain uncontaminated per-rank peak RSS. Preserve folding, both q orientations, spin/volume/k weights. No symmetry reduction or optical-limit change.

**Tech Stack:** Fortran, existing MPI abstraction, Hpsi OpenMP, LAPACK ZHEEV, C getrusage. Python only offline checking/reporting. Default production path unchanged, no production jobs/binaries modified.

1. Add independent Hermitian PSD callback fixture: low-rank W action, full-rank equivalence, Ritz residuals, NaN rejection; observe missing-module failure.
2. Extract existing projected CG core with optional solve counts and true residual, preserving finite-T explicit terms and tolerances. Re-run existing Sternheimer tests.
3. Add generic Krylov/Ritz PDEP builder (oversampling min(N,2r+8), double orthogonalization, deterministic restarts for degenerate eigenspaces). Count every polarization action; return Ritz residuals.
4. Add a diagnostic finite-q adapter, all occupied sources, actual Hpsi, reciprocal folding and explicit +/-q sum matching screen_assemble. Record phases, RHS/iteration/action counts, explicit empty-state sum=0. Two holdout W actions and valence/conduction intraband transition-charge BSE-direct matrix elements. The conduction state is used only for observables, not for screening or a band sum.
5. Add optional input mode/rank/q indices and docs. Only equilibrium insulating scalar full-grid/k-MPI allowed; reference=-1, PDEP positive, default0 disabled. Reference file must match grid/cell/k/q/G order. Separate jobs for reference and ranks avoid peak-RSS contamination.
6. Audit source and all linked application objects before using new executables. Validate small Si then existing Nr12³ Nk4³ Si, one nonzero q, same cutoff and MPI/OMP. Rank sweep including full rank; 1% W-action and direct-matrix-element errors as initial criterion, report absolute errors too. No full-rank matrix retained during PDEP.
7. Collect derived results, input, source/build/audit provenance and RT-free diagnostic logs in canonical calc2026, excluding GS payloads. Document whether total-time or memory gains exist. Expand to other q only if worthwhile; no Z-vector, metals, dynamics, head/wing changes in this experiment.
