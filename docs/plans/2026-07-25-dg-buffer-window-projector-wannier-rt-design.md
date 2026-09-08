# DG Buffer-Window Projector, Wannier GS, and RT Design

## Goal

Add a default-off Gamma-point, real-orbital DG route that starts from the
existing DC fragment SCF, constructs a self-consistent real-space DG ground
state without LCFO or EigenExa, builds buffer-periodic Wannier orbitals,
re-solves the ground state self-consistently in the Wannier representation,
and propagates only Wannier coefficients in real time.

The normal non-DG DC route remains unchanged and continues through its
existing LCFO and EigenExa stages.

## Scope and route isolation

The two routes are distinct:

```text
Normal DC, DG disabled
  DC fragment SCF
  -> LCFO
  -> EigenExa
  -> normal publication

DG route
  DC fragment SCF
  -> real-space DG self-consistent GS
  -> buffer-periodic Wannier construction
  -> Wannier-coefficient DG self-consistent GS
  -> Wannier-coefficient RT
```

When `yn_dg_dc_local_periodic='y'`, the route must not call LCFO, EigenExa,
WPW, normal checkpoint publication, or conventional RT. It must not fall
back to any of those routes after a failure.

The initial DC calculation supplies:

- the initial density;
- Hartree, XC, and local potential;
- real-space `core + buffer` fragment orbitals;
- all occupied and empty states selected by the namelist `nstate_frag`;
- the state window from which the Wannier orbitals will later be built.

`nstate_frag` is a configured total fragment-state count. It is not a fixed
states-per-atom value.

## Why direct same-state face exchange is invalid

Conventional DC solves each fragment eigenproblem independently. Therefore,
orbital `io` in one fragment is not guaranteed to be the same state as
orbital `io` in a neighboring fragment. Even at Gamma with real orbitals,
independent solves may differ by:

- a sign for each nondegenerate state;
- state ordering;
- an arbitrary real orthogonal rotation within a degenerate or nearly
  degenerate subspace.

The restriction is stronger for chemically different neighboring fragments:
there need not be any one-to-one correspondence between their local
eigenstates.

Consequently, production DG code must not use the neighboring fragment's
same `io` value as the remote trace. The observed same-index face comparison
may remain a diagnostic, but it is not a valid Hamiltonian input.

## Real-space ownership

Each fragment retains the conventional DC real-space orbital storage on

\[
\Omega_I^{\rm ext}
=
\Omega_I^{\rm core}+\Omega_I^{\rm buffer}.
\]

All `nstate_frag` occupied and empty states are retained. Existing orbital
parallelism distributes the state axis; states are not replicated on every
rank.

DC atom and potential construction may retain its existing 27-neighbor
buffer information. The DG operator acts on exactly the six signed physical
faces. Edge and corner neighbors are not DG faces.

Every face comparison and projection must use points representing the same
physical positions:

```text
fragment I buffer <-> neighboring fragment J core
```

A core-edge to neighboring-core-edge comparison at different physical grid
points is not a valid buffer-continuation diagnostic.

## Full-window buffer-to-core projector

The neighboring trace is reconstructed from the complete configured
occupied-plus-empty window. Occupations are not used as projector weights;
they are used only in density construction.

Let \(\Phi_J\) contain all neighboring-fragment core states restricted to
the physical overlap with fragment \(I\)'s buffer, and let
\(\Psi_I^{\rm buffer}\) contain the local buffer states on the same physical
points. Because restriction to an overlap slab destroys full-box
orthonormality, the reconstruction uses the overlap metric:

\[
S_J=\Phi_J^T W\Phi_J,
\]

\[
B_{JI}=\Phi_J^T W\Psi_I^{\rm buffer},
\]

\[
C_{JI}=S_J^+B_{JI},
\]

\[
\widetilde\Psi_I^{J,\rm core}=\Phi_J C_{JI}.
\]

Here \(W\) is the physical overlap-grid weight and \(S_J^+\) is a
rank-revealing pseudoinverse. This is a local boundary reconstruction, not a
global coefficient eigensolver.

The construction is invariant to signs, state permutations, and orthogonal
rotations inside the retained neighboring state window, subject to the
rank-tolerance policy.

The projector uses a configurable depth of the buffer overlap, not only one
face layer. The reconstructed value and normal-derivative traces must come
from the same projector generation.

## Projector diagnostics and rejection

Each signed face records:

- buffer-to-neighbor-core projection residual;
- overlap-metric rank;
- smallest retained singular value;
- window escape norm;
- forward/reverse projection mismatch;
- physical grid-map identity;
- projector generation and operator fingerprint.

Insufficient rank, excessive escape norm, stale generations, nonfinite
values, or inconsistent physical grid maps reject the current stage
collectively. Production code must not silently pair state indices or
truncate the configured state window.

## Real-space DG Hamiltonian and CG

The ground-state Hamiltonian is

\[
H(\lambda)=H_{\rm DC}+\lambda H_{\rm DG},
\qquad 0\le\lambda\le1.
\]

The conventional DC volume, nonlocal, density, Hartree, XC, and local
potential implementations remain the owners of those terms. The six-face
DG action is the only additional Hamiltonian contribution.

The full-window face projector is rebuilt at an outer SCF boundary and held
fixed during one orbital-CG sweep. With the projector fixed, the DG action
is a linear operator. The same complete DG Hamiltonian must act on:

- the current orbital `xk`;
- every CG search direction `pk`;
- every refreshed `xk`.

Updating the projector during an inner CG sweep is forbidden because it
would change the Hamiltonian during the quadratic minimization.

After the orbital sweep, the existing Gram--Schmidt, occupation, core-density
assembly, Hartree, XC, and local-potential updates run in their established
order. The next outer SCF iteration then constructs a new face projector.

## DC handoff

The DG route accepts a DC handoff only after all of the following hold:

- the minimum DC SCF iteration has been reached;
- the DC density residual satisfies the configured handoff tolerance;
- the electron-count error satisfies its tolerance;
- all configured fragment orbitals are finite;
- every physical face has a valid buffer-to-core map;
- the full-window overlap metrics meet their rank policy;
- LCFO, EigenExa, Wannier construction, and normal publication have not
  started.

The handoff snapshots:

- all fragment orbitals;
- occupations;
- total and spin densities;
- Hartree potential;
- XC potential;
- local potential;
- continuation state and fingerprints.

Mixing history must not be cleared to zero. It is reset by seeding every
history slot from the accepted current density and potential. This prevents
the artificial `1 - mixrate` density jump produced by mixing a new density
with zero history.

## DG continuation and rollback

The DG strength is introduced through converged SCF stages. A stage must
converge before \(\lambda\) advances.

The accepted pre-stage residual is the rollback baseline. The first
post-increment iteration is compared with that accepted baseline; it must
not redefine an abnormal first residual as the new baseline.

On stage failure, rollback atomically restores:

- fragment orbitals;
- occupations;
- density and spin density;
- Hartree, XC, and local potential;
- mixing seeds;
- projector generation;
- accepted \(\lambda\).

The lambda step is then halved. Reaching the minimum step or rollback limit
fails the DG gate. There is no automatic fallback.

## Full-scale real-space DG fixed-point gate

Completion at \(\lambda=1\) requires all of:

- an unmixed density fixed point;
- orbital residual tolerance;
- orthogonality tolerance;
- buffer-to-core projection residual tolerance;
- window escape-norm tolerance;
- DG Hermiticity tolerance;
- forward/reverse face-balance tolerance;
- electron-count tolerance;
- complete coverage of all six signed faces;
- correct physical-periodic wrapping;
- a successful DG-checkpoint round trip.

A mixed-density residual alone cannot publish a ground state.

If the Si64 full-scale gate fails, Wannier construction and all later
Wannier and RT tasks remain blocked.

## Buffer-periodic Wannier construction

Wannier orbitals are built only after the real-space DG ground state passes
its full-scale gate.

`nstate_frag` supplies the occupied-plus-empty construction window but is not
the RT propagation dimension. Wannier orbitals are constructed on each
fragment's buffer-extended periodic space. The auxiliary periodic boundary
at the buffer exterior permits tails much wider than a core-only
construction.

Each persisted Wannier orbital includes:

- core values;
- its complete retained buffer tail;
- center and spread;
- fragment ownership;
- physical-periodic image metadata;
- source-window identity;
- gauge and projector fingerprints.

## Wannier-representation self-consistent ground state

Projecting the real-space DG Hamiltonian once is insufficient. Finite
windows, retained tails, and local buffer-periodic conditions change the
representation. The Wannier state must therefore be solved
self-consistently:

\[
\Psi_n(\mathbf r)=\sum_{Iw}C_{Iw,n}W_{Iw}(\mathbf r),
\]

\[
H_W[\rho]C=S_WC\varepsilon,
\]

\[
H_W=H_W^{\rm volume}+H_W^{\rm DG}.
\]

Each outer iteration:

1. reconstructs the core-owned density from Wannier coefficients;
2. updates the existing DC Hartree, XC, and local potential;
3. updates the volume Hamiltonian on retained Wannier tails;
4. updates all six-face Wannier DG blocks;
5. solves the required occupied and empty coefficient states iteratively;
6. checks density, coefficient, metric, face, and charge diagnostics.

The DG route must not promote this solve to LCFO or EigenExa.

The Wannier ground-state gate requires:

- a coefficient and density fixed point;
- coefficient residual tolerance;
- \(C^\dagger S_WC=I\);
- Wannier-Hamiltonian Hermiticity;
- six-face DG balance;
- electron-count tolerance;
- tail-norm and window-escape tolerances;
- bounded energy and density differences from the accepted real-space DG
  ground state;
- checkpoint round-trip success.

## Wannier-coefficient real-time propagation

Real-time propagation uses the converged Wannier orbitals, not
`nstate_frag` DC orbitals:

\[
iS_W\dot C(t)=H_W[\rho(t),A(t)]C(t).
\]

The initial implementation keeps the Wannier orbitals fixed and propagates
only their coefficients. Every Runge--Kutta stage:

1. reconstructs density from the stage coefficients;
2. updates Hartree, XC, and local potential;
3. updates external-field, nonlocal, and DG Hamiltonian terms;
4. evaluates \(S_W^{-1}H_WC\);
5. advances the stage coefficients.

The route monitors tail norm and window escape norm. Leaving the verified
fixed-Wannier space is an explicit RT failure. Automatic basis update and
fallback to conventional RT are outside this design.

## Verification strategy

### Buffer-to-core projection

Before accepting a production DG action, tests cover:

- identical and chemically different neighboring fragments;
- all six signed faces;
- physical-periodic faces;
- multiple buffer widths;
- the complete occupied-plus-empty state window;
- rank deficiency and pseudoinverse tolerance;
- sign, permutation, and degenerate-subspace rotation invariance.

The primary residual is

\[
r_{\rm proj}
=
\frac{\|\Psi_I^{\rm buffer}
-\Phi_JS_J^+\Phi_J^TW\Psi_I^{\rm buffer}\|}
{\|\Psi_I^{\rm buffer}\|}.
\]

### DG operator

Analytic fixtures verify:

- exact zero at \(\lambda=0\);
- linearity for a frozen projector;
- consistency, symmetry, and penalty terms;
- Hermiticity;
- forward/reverse face balance;
- six-face and physical-periodic ownership;
- fragment relabeling invariance;
- buffer-width convergence.

### Normal DC isolation

With DG disabled, the established SCF trace, LCFO, EigenExa, checkpoint, and
publication behaviors remain unchanged.

With DG enabled, contract tests prove that LCFO, EigenExa, WPW, normal
checkpoint publication, and conventional RT are not called.

### DG SCF

Small MPI fixtures cover:

- DC handoff;
- projector freezing inside an orbital-CG sweep;
- projector rebuilding at outer SCF boundaries;
- continuation and lambda advancement;
- first-iteration rollback against the accepted baseline;
- density and potential restoration;
- nonzero mixing reseeding;
- full-scale unmixed completion;
- checkpoint corruption and rank disagreement.

### Si64 production gate

The initial production gate uses Gamma, non-SOI LDA Si64,
`nstate_frag=400`, eight MPI ranks, and one OpenMP thread. The normal DC
reference retains `threshold=1e-9`. DG runs cover multiple handoff
tolerances, at least two buffer widths, and at least two equivalent fragment
decompositions.

Every run records:

- runtime and peak memory;
- SCF and CG counts;
- projector construction time;
- metric rank and smallest retained singular value;
- projection residual and escape norm;
- continuation and rollback history;
- density and orbital residuals;
- energy and electron count;
- orthogonality and Hermiticity;
- six-face balance;
- checkpoint publication state.

## Development requirements

Every implementation task requires:

- observed RED evidence;
- focused verification;
- specification review;
- code-quality review;
- resolution of all Critical and Important findings;
- a fresh parent-prerequisite overlay build;
- `git diff --check`;
- verification-before-completion before commit.

The parent worktree is read-only. Its uncommitted implementation must not be
copied into this branch.

