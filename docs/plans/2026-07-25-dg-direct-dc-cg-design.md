# Direct DG term in the existing DC ground-state loop

## Decision

Keep the conventional DC ground-state algorithm and add only the real-space
DG/SIPG interface action to its Hamiltonian applications. Do not introduce a
projected local-basis coefficient solver between DC and the DG fixed point.

## Data flow

Each fragment continues to own all `nstate_frag` orbitals on its existing
core-plus-buffer grid. `solve_orbitals` continues to call the standard
orthogonalized CG and the existing `gram_schmidt`. For every CG Hamiltonian
application to the current orbital and search direction, add the scaled
six-face DG/SIPG action after the conventional volume/nonlocal `hpsi`.

After the orbital solve, retain the unmodified DC path:

1. `calc_density` on each fragment;
2. `calc_rho_total_dcdft` for the total density;
3. the existing density mixing, Hartree, XC, and local-potential update;
4. `calc_vlocal_fragment_dcdft` back to each fragment.

All occupied and empty states requested by `nstate_frag` remain present.
There is no LCFO, EigenExa, global real-space orbital, global coefficient
diagonalization, core-metric diagonalization, or projected-density rebuild.

## Continuation

The DG scale is an SCF-level control. Start from the converged or accepted DC
state, increase the scale from zero to one with rollback on a failed stage,
and run an unmixed fixed-point check at one. Within every orbital CG
iteration, the DG action is recomputed from the current orbital or search
direction; Hartree and `vlocal` retain the normal outer-SCF update cadence.

## Topology

The SIPG surface term uses the canonical six oriented faces. Communication
must validate the reciprocal periodic neighbor for each `-x/+x`,
`-y/+y`, and `-z/+z` face. A 27-neighbor representation is not introduced
unless required independently by an existing real-space stencil.

## Checkpoint and gates

Publish a DG-only checkpoint only after the scale-one unmixed fixed point.
Store the fragment orbitals, occupations, total density, Hartree/XC/local
potentials, continuation history, six-face topology, fingerprints, and
measured acceptance diagnostics. Standard DC, LCFO, WPW, and RT publication
remain disabled.

Task 5 still requires the Si64/LDA/Gamma/non-SOI matrix and the
parent-prerequisite overlay build. If no Si64 configuration passes, stop
before Wannier Tasks 6–8.

## Rejected alternatives

- Projecting the DC orbitals onto a core metric and solving a new distributed
  coefficient problem duplicates the DC eigensolver and fails for legitimate
  rank-deficient core restrictions.
- A replicated full-system diagonalization changes the scaling and is not
  needed for the real-space DG ground state.
