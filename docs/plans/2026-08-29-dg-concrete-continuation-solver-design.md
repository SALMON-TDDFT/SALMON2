# Concrete WF+PW DG Continuation Solver Design

## Decision

The hybrid WF+PW ground state is the self-consistent ground state of the
complete discontinuous-Galerkin operator. Build that operator from its
variational pieces; do not project SALMON's ordinary real-space `hpsi` and
then add DG faces. The ordinary finite-difference action crosses fragment
boundaries and makes the separation between broken volume and SIPG interface
terms ambiguous.

For one fixed WF+PW basis, define

\[
H_\lambda[\rho]=T_{\mathrm{broken}}+V_{\mathrm{NL}}
+V_{\mathrm{local}}[\rho]+\lambda H_{\mathrm{SIPG}},
\qquad 0\leq\lambda\leq1.
\]

Only the complete SIPG kinetic-interface operator is scaled by `lambda`.
The nonlocal pseudopotential is part of the physical volume Hamiltonian and
is present in full at every stage. When a projector crosses fragment
boundaries, assign it one canonical owner and accumulate its complete matrix
contribution exactly once.

Implement one explicit physical loop in the guarded production branch, where
the existing SALMON grid, Poisson, XC, DC, and basis objects are already in
scope. Keep lambda decisions and numerical kernels in small modules. Do not
add a callback fixture, production adapter, runtime catalog, independent trace
mixer, or another controller layer.

## Fixed variational payload

Complete and freeze the retained WF+PW basis before the lambda-zero solve.
For every fragment, retain basis values and gradients on its owned interior
grid and values and canonical-normal derivatives on its faces. The fixed
payload contains:

- global basis IDs and exactly-one row ownership;
- basis values and gradients on fragment interiors;
- the DG metric `S` and its distribution;
- the broken-volume kinetic matrix;
- the exactly-once nonlocal matrix;
- the complete coefficient-independent SIPG matrix;
- the retained-basis representation of the actual symmetry group, with one
  explicit identity operation for an identity-only system;
- basis, cutoff, selection, metric, operator, topology, and symmetry
  fingerprints.

For fragments `K`, assemble

\[
T^{K}_{ij}=\frac12\sum_{\mathbf r\in K}w_{\mathbf r}
\nabla\phi_i^{K*}(\mathbf r)\cdot\nabla\phi_j^K(\mathbf r).
\]

The sum contains only points owned by `K` and has no cross-fragment block.
All cross-fragment kinetic coupling comes from the separately assembled SIPG
operator, including consistency flux, adjoint-consistency flux, and penalty.

Assemble the nonlocal term from projector overlaps,

\[
V^{\mathrm{NL}}_{ij}=\sum_a\sum_{\mu\nu}
\langle\phi_i|\beta_{a\mu}\rangle D^a_{\mu\nu}
\langle\beta_{a\nu}|\phi_j\rangle.
\]

A projector may overlap several fragments. Its contribution is neither
fragment-localized nor lambda-scaled. Canonical ownership prevents duplicate
accounting without removing cross-fragment matrix elements.

At an SCF iteration, only the local-potential matrix changes:

\[
V^K_{ij}[\rho]=\sum_{\mathbf r\in K}w_{\mathbf r}
\phi_i^{K*}(\mathbf r)V_{\mathrm{local}}[\rho](\mathbf r)
\phi_j^K(\mathbf r).
\]

The basis, metric, broken kinetic, nonlocal matrix, and SIPG matrix remain
bitwise fixed throughout one continuation attempt.

### Concrete density-to-potential path

Do not wrap the physical updates in callbacks. The guarded production loop
follows the established DC decomposition explicitly. Each fragment owns its
current core density. Assemble those exactly-once core values into the total
real-space density using the DC fragment-to-total map. Compute only the
Hartree field on the total-system FFT distribution, using the existing total
Poisson/FFT objects, and redistribute that Hartree field to fragment grids.

Compute exchange-correlation on each fragment from its density buffer. The
accepted scope contains local or semilocal functionals, so XC needs only the
fragment values and, for a semilocal functional, its existing halo/gradient
exchange. XC does not require a second full-system density copy or FFT.
Combine the redistributed Hartree field, fragment XC field, and fixed local
ionic potential on the fragment grid, then project that combined local field
with the broken-volume local-matrix assembler.

Thus one density update is

```text
fragment core density
  -> DC core-to-total assembly
  -> total Hartree FFT
  -> total-to-fragment Hartree redistribution
  -> fragment-local/semi-local XC
  -> Hartree + XC + fixed local ionic potential
  -> WF+PW local-potential rows
```

The full-system collective is confined to the density assembly, Hartree FFT,
and Hartree redistribution. It is not used for XC, basis fields, or SIPG face
traces.

## Communication boundary

Materialize a SIPG face with one-to-one exchange between its two neighboring
fragment owners. Exchange only basis values and normal derivatives needed on
that face. Do not gather face traces or basis fields over the full
communicator. Because the basis and SIPG matrix are fixed, no face exchange
is needed during the SCF loop.

Global communication remains only where the physics or distributed algebra
requires it:

- one-time validation of exactly-one basis and projector ownership;
- SALMON's total-density Hartree/local-potential update;
- the distributed generalized eigensolver;
- scalar electron-count and residual reductions;
- collective accept or reject decisions;
- optional physical-symmetry residuals.

The procedure must not assume equivalent fragments, equal local basis sizes,
a regular neighbor graph, or nontrivial symmetry. Crystals, liquids, defects,
interfaces, surfaces, and identity-only systems use the same path.

## Three state objects

Use only three conceptual objects.

1. `fixed_payload` stores the immutable basis and matrices described above.
2. `iterate` stores one consistent SCF iterate: input and output density,
   local-potential and complete Hamiltonian rows, eigenvalues, occupations,
   occupied coefficients, `Gamma`, `Q`, interface observables, residuals,
   epochs, and fingerprints.
3. `accepted` is a deep copy of the last converged `iterate`, together with
   its lambda and the next proposed lambda step.

Here

\[
\Gamma=CfC^\dagger,
\qquad Q=C_{\mathrm{occ}}C_{\mathrm{occ}}^\dagger S.
\]

Use `Gamma` for density and physical expectation values and `Q` to track the
occupied subspace. Never use raw coefficient differences as a subspace
residual and never mix eigenvector coefficients.

## Density mixing and continuation

The immutable initial density is exactly the converged DC total density. The
first local-potential construction must read that density before any WF+PW
density reconstruction. Lambda zero is nevertheless converged as its own
finite-basis fixed point.

At fixed lambda, all ranks perform this single ordered loop:

1. build the total local potential from the current input density;
2. project only the density-dependent local potential;
3. form `T_broken + V_NL + V_local[rho] + lambda * H_SIPG`;
4. solve `H C = S C epsilon`;
5. determine cluster-consistent occupations;
6. construct `Gamma` and the `S`-metric occupied map `Q`;
7. reconstruct output density and gauge-invariant interface observables;
8. evaluate residuals belonging to this one epoch;
9. accept the fixed point, or mix only the density and repeat.

For linear damping,

\[
\rho_{\mathrm{in}}^{m+1}=\rho_{\mathrm{in}}^m+
\alpha_\rho(\rho_{\mathrm{out}}^m-\rho_{\mathrm{in}}^m),
\qquad0<\alpha_\rho\leq1.
\]

The existing SALMON density mixer may provide simple, Pulay, or Broyden
updates, but it must consume and produce density only. Recompute the local
potential from the mixed density. Interface traces are observables rebuilt
from the current `Gamma`; they are not independently mixed Hamiltonian input.

After a converged stage, propose one uniform next lambda. Stable, inexpensive
stages may grow the bounded step. Failure, sustained residual growth, or an
occupied-subspace discontinuity restores the complete accepted state and
shrinks the step. A small gap alone is not a rejection condition; it informs
occupation clustering and step reduction.

## Acceptance and rollback

Evaluate inexpensive gates first:

- generalized eigen-residual `R_H`;
- density residual `R_rho`;
- interface-observable residual `R_T`;
- metric orthogonality residual `R_S`;
- occupied-projector change;
- electron-number and occupation consistency;
- Hamiltonian and metric Hermiticity;
- finite-value and epoch consistency.

Only a candidate passing those gates receives the expensive reconstructed
real-space DG residual and, for a nontrivial actual group, operator and
complete occupied-subspace covariance checks. Individual eigenvectors need
not transform as symmetry eigenstates. Identity-only systems execute the same
checks with the explicit identity action.

A rejected trial restores density, potential, Hamiltonian, coefficients,
occupations, eigenvalues, `Gamma`, `Q`, interface observables, mixing history,
epochs, and fingerprints from `accepted`. No rejected or stale field may
survive. Convert every rank-local failure to a communicator-wide decision
before entering the next collective operation.

## Final lambda-one refresh

After apparent convergence at lambda one, perform exactly one unmixed full
refresh:

1. rebuild density and interface observables from the converged occupied
   state;
2. rebuild the local potential and local-potential matrix;
3. form the full lambda-one Hamiltonian;
4. solve the generalized eigenproblem;
5. rebuild `Gamma`, `Q`, density, and interface observables;
6. evaluate every final residual and provenance gate on that same epoch.

Publish this state only when input and output densities agree at the final
tolerance and all gates pass. Do not restore pre-refresh values or run
duplicate acceptance.

## Protected routes and removed mechanisms

Put every new call behind the explicit hybrid DG-continuation flag. Leave
ordinary GS and RT, existing DC+LCFO, Wannier90, overlapping-Wannier, and
their checkpoint formats unchanged.

The continuation branch must not use:

- an ordinary `hpsi` projection as its DG volume Hamiltonian;
- the divided-SCF one-shot final LCFO solve;
- the occupied-only `overlapping_wannier_occupied.chk` publication;
- a production adapter, runtime catalog, or callback table;
- raw eigenvector or independent trace mixing;
- face-local or fragment-local lambda;
- basis reselection during continuation.

The existing LCFO-flux weak-volume code is a numerical reference, not an
integration dependency and not a protected route to edit.

## Test strategy

Develop every behavior test-first.

1. A two-fragment analytic broken-volume test verifies fragment-interior
   kinetic and local terms, exactly-once crossing nonlocal projectors, zero
   kinetic cross block in the volume matrix, and Hermiticity.
2. A composition test verifies `H_volume + lambda H_SIPG`, uniform lambda,
   zero interface contribution at lambda zero, complete contribution at
   lambda one, and SIPG-only kinetic cross-fragment blocks.
3. The production-loop contract test starts from an exact supplied DC density,
   proves gradual density mixing, requires refresh of `Gamma`, `Q`, density,
   and traces after every solve, and forbids lambda advance before convergence.
4. Continuation tests cover adaptive step growth, complete rollback, phase
   and degenerate-space gauge invariance, small-gap handling, and collective
   rank-local failure.
5. A production-format integration test proves that the concrete loop uses
   the fixed broken-volume, nonlocal, metric, and SIPG payload and never calls
   the ordinary `hpsi` projection, one-shot solve, or occupied-only checkpoint.
6. Final acceptance uses Si64 with eight MPI ranks, `OMP_NUM_THREADS=1`, and no
   time cutoff, followed by checkpoint-identical zero-field RT stationarity
   checks. Si64 dimensions or symmetry are not solver assumptions.
