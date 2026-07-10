# SAWF Symmetry-Closed LCFO Basis Design

## Goal

Construct the Wannier90 SAWF seed in an LCFO Galerkin space that is closed
under every validated fragment-compatible crystal symmetry operation, without
changing DC-SCF or the normal DC-LCFO output path.

## Evidence

The aligned C64 case has exact atom, integer-grid, and whole-fragment symmetry
maps, but its finite LCFO basis is not closed under inversion. Raising
`energy_cut` from 0.1 to 4.0 au and `nstate_frag` from 400 to 600 increased the
retained basis from 96 to 586 functions per fragment, while inversion leakage
fell only from 9.31 percent to 2.36 percent. The selected 1200-band SAWF space
still had a minimum symmetry singular value of about 6.2e-4 and a unitarity
residual of about 0.98.

For one 586-function fragment basis `B`, direct analysis gives
`dim(span(B,gB))=1172` at tolerance 1e-8 for inversion. Forty-five missing
directions have residual singular value above 0.1. A larger energy window is
therefore useful but cannot provide exact closure economically.

## Alternatives

### Continue increasing energy and state cutoffs

This improves the average leakage, but convergence is slow and approaches the
real-space complete basis. It is retained as an accuracy control, not as the
SAWF construction algorithm.

### Repair `D_band` by polar decomposition or operator symmetrization

This would make the representation unitary after projection while hiding that
the underlying state space is not invariant. It changes the representation
without supplying the missing physical basis directions and is rejected.

### Close the Galerkin basis under the actual symmetry group

This is the selected approach. For every target fragment, collect the mapped
images of source-fragment basis functions under all validated operations, then
orthonormalize their union with a rank threshold. The resulting fragment bases
carry the group action by construction. Hamiltonian, overlap, AMN, MMN, and
`D_band` are all rebuilt in that same basis.

## Scope And Isolation

The normal sequence remains unchanged when `wannier_site_symmetry='off'`:

```text
DC-SCF -> LCFO basis -> Flux Hamiltonian -> EigenExa -> normal outputs
```

When SAWF is enabled, normal DC-SCF and normal LCFO output are completed first.
A separate seed sequence then runs:

```text
validated atom-derived operations
  -> symmetry orbit of fragment bases
  -> rank-revealing orthonormalization
  -> SAWF-only Flux Hamiltonian
  -> EigenExa seed eigenstates
  -> eig/amn/mmn/dmn
  -> Wannier90 site_symmetry
```

The SAWF seed uses distinct filenames and provenance. It must not overwrite
`basis_functions.bin`, `wavefunctions.bin`, or the normal LCFO eigenvalue file.

## Symmetry-Closed Fragment Basis

Let `B_f` be the core basis of fragment `f`. A validated operation `g` supplies
an integer periodic grid map and a whole-fragment permutation `p_g`. The
candidate set for target fragment `t` is

```text
C_t = { g phi | phi in B_f, p_g(f)=t, g in G }.
```

Candidates are accumulated operation by operation in deterministic normalized
operation order. A rank-revealing eigendecomposition of the candidate Gram
matrix retains singular directions above a configurable threshold. Complete
symmetry orbits are never pruned independently: if capacity or rank policy
would split an orbit, construction fails with a rank-local diagnostic.

The algorithm validates after construction that every transformed basis is in
the corresponding target span and that all target bases have compatible group
representations. The accepted tolerance is the stricter of the basis rank
threshold and `wannier_symmetry_tolerance`.

## Hamiltonian Consistency

The closed basis is not assigned an artificial symmetrized Hamiltonian. SALMON
reapplies the same real-space local/nonlocal Hamiltonian and the same DG Flux
surface assembly used by the normal LCFO route to every closed-basis function.
The seed Hamiltonian is then diagonalized with EigenExa. This preserves the
physical operator and makes any remaining covariance defect observable.

Every seed artifact must derive from this one basis. Mixing normal-LCFO
eigenvectors with closed-basis AMN/MMN or `D_band` is fatal.

## Memory And Distribution

No permanent global real-space wavefunction array is introduced. Each rank
owns its fragment basis. Symmetry images are exchanged only between mapped
source and target fragment owners. Candidate blocks are orthogonalized on the
target owner and discarded after the final basis is formed.

Basis dimensions become dynamic and may exceed `nstate_frag`. The seed path
stores explicit per-fragment dimensions and checks integer/MPI count limits
before allocation. A configurable hard capacity prevents accidental growth to
the full real-space grid.

## Failure Policy

SAWF stops before publishing seed files when:

- an operation is not an atom/grid/whole-fragment symmetry;
- a mapped basis block is missing or non-finite;
- orbit closure exceeds the configured capacity;
- rank pruning would break an orbit;
- post-construction closure exceeds tolerance;
- the seed Hamiltonian is not Hermitian or symmetry covariant;
- seed provenance does not match the current run.

Diagnostics are emitted rank-locally before collective failure reduction.

## Verification

Focused tests first cover identity, self-mapped inversion, fragment-pair
inversion, a multi-operation group, deterministic rank pruning, capacity
failure, and a nonsymmetric input basis. Integration then checks:

1. `wannier_site_symmetry='off'` remains byte-for-byte behavior-compatible.
2. The aligned C64 seed basis has closure leakage near roundoff.
3. `D_band` is unitary and passes group/covariance checks without polar repair.
4. Wannier90 accepts `.dmn` and produces symmetry-adapted functions.
5. Field-off stationarity and long-pulse HHG are evaluated only after the
   representation gates pass.
