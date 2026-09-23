# Direct time-dependent Wannier diagnostic

**Goal:** Transform propagated Si Bloch orbitals to time-dependent Wannier orbitals and measure centers/spreads without modifying alpha.
**Architecture:** Postprocess existing shared SALMON wavefunction checkpoints. Construct an orthonormal projected bond-centered Wannier gauge in the initial occupied subspace (not MLWF), retain it during propagation, compare no pump/weak/strong pump at the same4³ k grid. Do not identify a single-orbital spread with exciton radius or dielectric screening.
**Tech Stack:** Existing MPI SALMON, Python/NumPy/Matplotlib; no production Fortran change.

1. Read shared wfn.bin using verified Fortran ordering and Si_k.data weights; validate norm/orthogonality and occupations.
2. Write numerical tests for projected gauge covariance under random occupied rotations, singular trial rejection, supercell transform normalization/density reconstruction, periodic center/spread.
3. Build16 bond-centered localized trial functions for Si8, project and symmetrically orthonormalize at each k. Check singular values, orthogonality and initial centers. No claim of maximum localization.
4. Run no-pump/1e8/1e13 W/cm², same alpha0=1 K0=0 polarization closure, omega=.2au, width60au, dt=.08, nt1600. Checkpoint every200 steps. Enable existing gs projection to separately diagnose electron/hole excitation. Keep raw checkpoints in a task-specific /private/tmp directory and retain reproducible inputs and reduced data in repo.
5. Reconstruct4³-cell Wannier density using actual shifted k mesh; periodic circular centers and minimum-image second moments. Field-free evolution broadens fixed-gauge orbitals through band phases; explicitly retain a no-pump control. Optionally undo the common field-free occupied rotation at each snapshot for a phase-controlled diagnostic; label the gauge.
6. Verify summed Wannier density equals occupied Bloch density; compare current traces to existing runs, excitation projection completeness, finite-supercell boundary weight. Plot moments and representative densities; state that4³ limits spatial extent, no k scan.
7. Review, record findings, commit locally. No alpha feedback is added based on an unvalidated spread-to-screening relation.

## Updated scope: previous U warm starts every ten RT steps

At user request, add occupied-subspace Marzari–Vanderbilt spread minimization on the same 4³ mesh. Reuse the converged U(k) from the previous snapshot as the next initial guess, every10 steps (0.0193511 fs). Monitor gradient convergence, orbital overlap matching and center continuity. Record the gauge-invariant spread separately. Verify the analytic gradient, mesh boundary phases and occupied-gauge covariance before analyzing the full trajectory. Compare final strong-pump warm start with an independent initial-gauge start.

This stage evaluates MLWF diagnostics offline and samples the pre-existing alpha at the same cadence. It does not yet change the production alpha update interval or feed MLWF widths into alpha: the screened interaction / alpha relation has not been specified. Reusing U is a localization initial guess, not temporal averaging of screening.

Refinement: transport previous U through the polar factor of the overlap between current Bloch and previous Wannier orbitals before minimizing. This exactly cancels a pure occupied-basis rotation, preserves the current occupied projector, and avoids spending hundreds of iterations undoing field-free phase evolution. Keep the direct-U no-pump control for comparison.
