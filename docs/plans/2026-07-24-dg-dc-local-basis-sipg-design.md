# DG-DC Direct Fragment-Orbital SIPG Ground-State Design

## 2026-07-24 architecture correction

The production ground-state path does not construct LCFO functions, solve a
global LCFO coefficient eigenproblem, or materialize whole-system real-space
wavefunctions.  Each DC fragment already owns
`dg_dc_candidate_orbitals_per_atom * natom_fragment` orbitals, including the
required empty states.  The DG path keeps those orbitals in the existing
`spsi` representation and reuses the existing `solve_orbitals` sequence:

1. apply the existing DC volume/nonlocal Hamiltonian;
2. add the physical-face SIPG action to both the orbital and CG search
   direction Hamiltonian applications;
3. update all fragment orbitals with the existing orthogonalized CG;
4. run the existing BLAS Gram--Schmidt;
5. reuse the existing DC density, Hartree, XC, `vlocal`, and mixing update.

The only new cross-fragment data are face values and normal-derivative stencil
layers for the same orbital indices.  No local-basis coefficient matrix is
needed by the production path.  The small generalized-eigenproblem helper
remains a test/reference utility and is forbidden by the production route
contract.

## Decision

Keep the existing DC fragment calculation as the ground-state implementation.
The DG extension adds only the SIPG interface contribution between neighboring
fragments.  Density construction, Hartree, XC, local potential, mixing,
conjugate-gradient minimization, and BLAS Gram--Schmidt remain on the existing
DC path.

This replaces the experimental nodal route that concatenated fragment-local
candidate orbitals into an axis of global real-space states.

## State and basis indices

`dg_dc_candidate_orbitals_per_atom` controls the number of fragment-local
orbitals retained from the DC solve.  With 40 orbitals per atom and eight atoms
in a fragment, each fragment retains 320 local orbitals:

```text
phi(fragment, local_orbital, local_grid)
```

These orbitals are a local basis.  They are not distinct global Kohn--Sham
bands.  The global DG states are represented by coefficients in the union of
the fragment-local bases:

```text
Psi(global_band) =
  sum(fragment, local_orbital)
    phi(fragment, local_orbital) * coefficient(fragment, local_orbital, global_band)
```

The requested global band count remains the total-system input `nstate`.
Fragment-local basis size and global band count are therefore independent.

## Hamiltonian

The volume contribution is the existing DC fragment Hamiltonian.  The only new
physical contribution is the SIPG interface matrix between physical neighboring
fragments:

```text
H_DG = H_DC_volume + H_SIPG_interface
```

The interface assembly exchanges the boundary values and normal derivatives of
the existing fragment-local orbitals.  It uses canonical physical-face
ownership, periodic neighbor topology, and the already-tested SIPG consistency,
symmetry, and penalty terms.

The assembled global matrix must be Hermitian within the Task 5 tolerance.
Auxiliary DC buffer boundaries must never be published as physical DG faces.

## Self-consistent iteration

At each SCF iteration:

1. Use the existing DC local orbitals and volume Hamiltonian.
2. Assemble or refresh the SIPG interface blocks.
3. Solve the global local-basis Hamiltonian for the requested `nstate` bands.
4. Reconstruct the fragment/core density from the global coefficients and
   fragment-local orbitals.
5. Reuse the existing DC density reduction, Hartree, XC, `vlocal`, mixing, and
   convergence machinery.
6. Repeat until the fixed-lambda stage converges.

The existing DC CG and BLAS Gram--Schmidt continue to prepare and update the
fragment-local orbitals.  No independent nodal CG, nodal Cholesky
orthogonalization, or concatenated-state expansion participates in production.

## Continuation and checkpoint

The existing adaptive lambda continuation remains default-off and transactional.
Its production callback invokes the local-basis DC iteration above.  The
checkpoint stores:

- fragment-local basis metadata and fingerprints;
- global band eigenvalues, occupations, and local-basis coefficients;
- density, Hartree, XC, and local potential;
- physical-face topology and SIPG diagnostics;
- continuation history and all acceptance diagnostics.

The checkpoint must not claim that concatenated fragment candidates are global
real-space orbitals.

## Rejected alternatives

1. **Concatenate fragment candidates as global real-space states.** This
   confuses local basis size with global band count and causes an unnecessary
   large-state eigensolve.
2. **Keep only occupied global bands.** Empty bands are required for the
   subsequent orbital-coverage and Wannier gates.
3. **Promote to LCFO/WPW.** Those consumers remain separate, default-off routes
   and are not enabled by the DG ground-state implementation.

## Verification

Tests must first fail when:

- local candidate count is multiplied into the global band count;
- auxiliary buffer boundaries are treated as physical faces;
- the SIPG interface term is omitted or double counted;
- the assembled Hamiltonian is non-Hermitian;
- density/Hartree/XC/`vlocal` bypass the existing DC path;
- checkpoint metadata labels local candidates as global nodal orbitals.

Task 5 then requires a fresh parent-prerequisite overlay build and the specified
Si64 matrix.  Tasks 6--8 remain blocked unless at least one Si64 DG
ground-state configuration passes every gate.
