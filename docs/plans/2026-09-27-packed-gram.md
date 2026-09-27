# Packed Hermitian Gram check

Goal: reduce the all-occupied work in direct-WF orthogonality validation without reducing validation frequency or changing its1e-8 tolerance.

Design: compute only the upper triangle of C^H C with BLAS ZHERK, pack n(n+1)/2 entries, reduce them over the existing spatial communicator and compute max|Gram-I| there. Hermitian symmetry makes the upper-triangle maximum sufficient. Check all packed complex values for finiteness before max. No sparse cutoff or new physical approximation. This does not remove dense U or full occupied-column propagation.

Alternatives: skipping Gram checks would weaken validation; sparse propagation requires a separate R_prop and gauge design. The exact packed reduction is a bounded improvement to the current reference.

Steps: add an MPI2/4 comparison against a dense global Gram, including unequal/empty row partitions, orthonormal and perturbed matrices and nonfinite entries. Observe missing API failure, implement, run probe. Build and run native direct-WF regressions. Compare a frozen before/after C128 RT sequentially and record provisional timing, numerical error and the exact communication reduction. Review and update the root development note, then push.
