# Global-Covariant Fixed-Rank Review Remediation Design

## Purpose

Repair the reviewed global-covariant Wannier construction without weakening the full-system atomic symmetry, target rank 384, covariance tolerances, or strict Si64 SCF evidence.

## Construction

The full-system affine group remains derived only from the instantaneous atomic configuration.  The identity operation is found from the product table rather than operation order.

The mandatory occupied seed is not a direct sum of independently diagonalized fragment orbitals.  After the converged DC calculation, the existing LCFO Hamiltonian construction and its LAPACK/EigenExa diagonalization provide one coherent set of full-system coefficients.  Each fragment evaluates its own LCFO basis contribution on its buffered box.  Contributions from every fragment buffer covering a physical grid point are accumulated onto the unique owner core.  Thus the distributed rows are restrictions of the same full-system occupied states, while no full-system real-space wavefunction is gathered.  Core points alone own quadrature; buffers are halos used to evaluate tails and later derivatives.

The rank of this coherent occupied seed is the physical occupied count (128 for ideal Si64), and its full-group closure must retain that rank within tolerance.  A rank increase at this gate is an error, because it means the reconstructed LCFO occupied projector does not carry the atomic symmetry.  Optional physical projection seeds are added only after this gate.  Every completed optional seed orbit is verified by the distributed projector leakage `||(1-P) U(g) P||`.

The closed candidate space may exceed rank 384 by at most one group orbit.  Its measured identity overlap is the metric.  The occupied block is metric-orthonormalized, and a symmetry-averaged physical localizer built from projection-seed overlap and periodic position moments selects complete invariant blocks at exactly rank 384.  The selected representation is remeasured with its measured metric.

Point actions carry lattice-wrap phases.  Gamma-only Si64 has unit phase, but the API must remain correct for general commensurate Bloch phases.

## Localization and evidence

The localization line search and gradient both use the complete periodic spread.  A sparse overlap graph may be diagnostic or preconditioning information, but cannot declare convergence.  The ideal Si64 evidence uses full-system group order, inversion, leakage, unitarity, closure, and covariance; fragment site symmetry and displaced structures are not prerequisites.  Polarization remains the primary LR/HHG observable, and every even harmonic is reported as peak, dip, or slope together with its suppression ratio.

## Acceptance

Each change follows RED-GREEN TDD.  Focused MPI verification runs on 1/2/4/8 ranks.  Critical and Important findings are re-reviewed before strict genuine-Si64 GS, V3, LR/HHG, final clean-first committed-HEAD build, commit, and dual push.
