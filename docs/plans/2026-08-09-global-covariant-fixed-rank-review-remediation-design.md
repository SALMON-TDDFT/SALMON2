# Global-Covariant Fixed-Rank Review Remediation Design

## Purpose

Repair the reviewed global-covariant Wannier construction without weakening the full-system atomic symmetry, target rank 384, covariance tolerances, or strict Si64 SCF evidence.

## Construction

The full-system affine group remains derived only from the instantaneous atomic configuration.  The identity operation is found from the product table rather than operation order.  All fragment core-owned occupied functions are processed before optional localization seeds.  Every completed seed orbit is verified by the distributed projector leakage `||(1-P) U(g) P||`.

The closed candidate space may exceed rank 384 by at most one group orbit.  Its measured identity overlap is the metric.  The occupied block is metric-orthonormalized, and a symmetry-averaged physical localizer built from projection-seed overlap and periodic position moments selects complete invariant blocks at exactly rank 384.  The selected representation is remeasured with its measured metric.

Point actions carry lattice-wrap phases.  Gamma-only Si64 has unit phase, but the API must remain correct for general commensurate Bloch phases.

## Localization and evidence

The localization line search and gradient both use the complete periodic spread.  A sparse overlap graph may be diagnostic or preconditioning information, but cannot declare convergence.  The ideal Si64 evidence uses full-system group order, inversion, leakage, unitarity, closure, and covariance; fragment site symmetry and displaced structures are not prerequisites.  Polarization remains the primary LR/HHG observable, and every even harmonic is reported as peak, dip, or slope together with its suppression ratio.

## Acceptance

Each change follows RED-GREEN TDD.  Focused MPI verification runs on 1/2/4/8 ranks.  Critical and Important findings are re-reviewed before strict genuine-Si64 GS, V3, LR/HHG, final clean-first committed-HEAD build, commit, and dual push.
