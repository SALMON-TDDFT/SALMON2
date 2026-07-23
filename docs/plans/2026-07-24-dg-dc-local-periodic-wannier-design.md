# DG-DC Local-Periodic Wannier Design

**Date:** 2026-07-24

**Status:** Approved

## Goal

Replace the fragile LCFO-to-Wannier-to-DG initialization chain with a
self-consistent real-space DG-DC ground state, then construct local Wannier
functions from the same local-periodic fragment boxes while restricting
physical ownership and observables to fragment cores.

The first gate is Si64, LDA, Gamma-only, non-SOI. Wilson links, length-gauge
coupling, polarization, and real-time propagation are outside this design.

## Motivation

The current path changes representation several times:

```text
independent DC fragments
→ LCFO
→ Wannier
→ DG/WPW
→ fixed-H correction
```

Each handoff introduces a new completeness, gauge, boundary, density, or
publication gate. A DG correction applied only after Wannier construction
cannot repair directions that were removed by an earlier local-basis
truncation.

The new path uses conventional DC only to create a stable initial density and
a deliberately rich orbital pool. It transfers to a real-space DG
variational problem before LCFO truncation or Wannier localization.

## Scope

The first implementation supports:

- Si64;
- LDA/local-density XC;
- one Gamma point;
- scalar, non-SOI orbitals;
- periodic bulk geometry;
- norm-conserving pseudopotentials already supported by the nodal complete-H
  action;
- local-periodic `core + buffer` fragment boxes;
- symmetric interior-penalty DG coupling;
- local Wannier construction and overlap-based fragment gauge alignment.

It does not change normal DC, LCFO, checkpoint, publication, or RT routes
unless the new route is explicitly enabled. Failure is collective and closed.

## Alternatives

### GS-native common nodal DG refinement

Extract the reusable nodal state, complete-H action, eigensolver, density
reconstruction, and diagnostics from the RT-owned layer into a common DG
layer. Add a GS/DC driver that performs staged DC-to-DG SCF. This is selected.

### Call RT nodal preparation from ground state

This is shorter but reverses dependencies: GS would depend on RT types and
lifecycle assumptions. It also makes Wannier export depend on an RT
preflight. This is rejected.

### Add face terms inside every independent fragment DC solve

This modifies the existing DC solver deeply while destroying its independent
fragment semantics. It is difficult to distinguish auxiliary local-periodic
boundaries from physical neighbor coupling. This is deferred.

## Two boundary systems

Each fragment owns a local periodic box:

```text
local box = physical core + buffer
```

The outer box boundary is periodic only for generating local orbitals and
local Wannier data. It is an auxiliary boundary condition.

Physical DG coupling is defined on core-to-core fragment faces. It uses the
actual neighboring fragment trace, not the opposite side of the same local
periodic box. Physical density, charge, energy contributions, and retained
Wannier centers are core-owned and counted exactly once.

The local-periodic boundary is acceptable only when increasing the buffer no
longer changes core density, eigenvalues, atomic-orbital projectability,
Wannier centers, spreads, and retained Hamiltonian elements.

## DG Hamiltonian

The existing nodal halo finite-difference action is a domain-decomposed
real-space action, not by itself a DG weak form. The ground-state operator is
defined explicitly as

\[
H_{\rm DG}=H_{\rm volume}+H_{\rm SIPG\ face}+V_{\rm nonlocal}.
\]

The face action contains:

- jumps of the orbital trace;
- averages of the normal derivative;
- symmetric consistency terms;
- a positive penalty term;
- canonical single ownership of each physical face;
- periodic physical-supercell faces;
- identical left/right contributions required for Hermiticity.

The implementation should reuse the canonical face schedule, face-trace halo,
weak-form evaluator, stencil, and nonlocal-pseudopotential infrastructure
already present in the repository. The penalty and all consistency terms
form one complete Hermitian interface operator.

## Early DC handoff

Conventional DC is an initial-density generator, not the final variational
space. It stops at the first collectively valid handoff satisfying:

- a configurable minimum number of DC iterations;
- finite density and potentials;
- correct total electron count;
- a loose density residual, initially compared at `1E-2`, `1E-3`, and
  `1E-4`;
- successful fragment eigensolves;
- no LCFO, Wannier, energy-window, or locality truncation;
- retention of the full configured candidate orbital pool.

DC mixing history is discarded at handoff because its approximate Jacobian
belongs to a different operator. The density and potentials are retained.

All retained fragment orbitals are materialized on the real-space nodal
degrees of freedom. From that point onward, the solver must be able to leave
the DC orbital span; it may not continue as a coefficient-only solve in a
frozen DC basis.

## Staged DC-to-DG SCF

Introduce a continuation parameter:

\[
H(\lambda)=H_{\rm DC}+\lambda(H_{\rm DG}-H_{\rm DC}),\qquad 0\le\lambda\le1.
\]

The complete Hermitian face operator is scaled as a unit; the penalty term is
not ramped independently.

The stages are:

1. **DC preconvergence:** run with `lambda=0` to the handoff gate.
2. **DG continuation:** at each fixed `lambda`, converge the orbital and
   density residuals to an intermediate tolerance before advancing.
3. **Full DG SCF:** hold `lambda=1` and converge to the final tolerances.
4. **Unmixed fixed-point verification:** rebuild the output density and
   potentials without mixing and repeat the full-DG eigensolve.

Advance `lambda` adaptively. Increase the step after clean convergence.
Rollback and halve it when the occupied-subspace overlap collapses, the
orbital/density residual grows beyond its allowed factor, the face residual
becomes nonfinite, or the eigensolver fails. A minimum step failure stops the
route collectively.

Do not respond to a large initial full-DG residual by converging DC more
tightly. Use it to reduce the initial continuation step. If the real-space
handoff is rank deficient, increase the retained DC candidate pool and repeat
the handoff.

## Ground-state acceptance

At `lambda=1`, require simultaneously:

\[
\max_i\|H_{\rm DG}\psi_i-\epsilon_i\psi_i\|\le\tau_\psi
\]

and

\[
\frac{\|\rho_{\rm out}-\rho_{\rm in}\|}
{\max(1,\|\rho_{\rm in}\|)}\le\tau_\rho.
\]

Also require:

- electron count 256 for the Si64 gate;
- finite density, Hartree, XC, eigenvalues, and orbitals;
- global occupied-orbital orthonormality;
- scale-aware Hermiticity of the complete operator;
- left/right cancellation of internal face flux;
- stable occupied-subspace projectors across the final iterations;
- agreement of mixed convergence and unmixed fixed-point verification;
- invariance under equivalent fragment decompositions within tolerance.

No LCFO, Wannier, checkpoint publication, or downstream route may consume an
unaccepted DG ground state.

## Orbital coverage manifest

Do not hard-code four, nine, forty, or another universal orbital count per
atom. For every species, construct a versioned orbital coverage manifest from:

- pseudopotential valence configuration;
- angular-momentum and radial projector channels;
- atomic bound and low-lying excited reference orbitals;
- local coordination and symmetry;
- projectability onto the selected DG outer eigenspace.

For each shell and radial multiplicity,

\[
N_{\rm target}^{(a)}
=\sum_{n,l}n_\zeta(n,l)(2l+1).
\]

The DC calculation may retain a much larger candidate pool, initially up to
the configured forty orbitals per atom. Those candidates form the variational
outer space and are not all automatically published as Wannier functions.
Remove only numerically dependent directions identified by a collective,
scale-aware metric-rank test before the DG solve.

Use a convergence ladder. For Si:

```text
Tier 0: 3s + 3p                         4 orbitals/atom
Tier 1: 3s + 3p + polarization d        9 orbitals/atom
Tier 2: second-zeta s/p + polarization  manifest-defined
```

Select the smallest tier for which the requested occupied and conduction
observables, density, projectability, Wannier spread, and buffer convergence
no longer change materially.

## Local Wannier construction

The accepted full-DG ground state and the manifest determine the Wannier
target dimension. For Si Tier 0, the target is four orbitals per atom:

\[
N_{\rm Wann}=64\times4=256.
\]

This contains the 128-dimensional occupied bonding space and sufficient
conduction/antibonding directions to complete the sp3 manifold. The occupied
projector is frozen and must be contained in the Wannier target:

\[
\mathcal W_{\rm target}\supseteq\mathcal P_{\rm occ}.
\]

The DG outer window may contain substantially more than 256 eigenstates,
drawn from the full retained DC candidate pool. Choose the target subspace by
projectability and rank against the manifest orbitals, not merely by taking
the lowest fixed number of eigenvalues.

For each local-periodic box:

1. build manifest atomic projections for core and buffer atoms;
2. compute the projection matrix against the DG outer eigenspace;
3. verify full target rank and singular-value bounds;
4. disentangle a target space that contains the frozen occupied projector;
5. localize within the local-periodic box;
6. retain Wannier functions whose home atom/center is core-owned;
7. use buffer functions only for neighbor overlap and matrix elements.

Align neighboring local Wannier gauges by unitary Procrustes/polar
decomposition on their shared buffer. Use a deterministic spanning tree and
diagnose closure around all fragment loops. This is an overlap-based static
gauge-stitching procedure; Wilson-link RT physics is not part of this design.

## Checkpoint and provenance

Publish a new DG-GS checkpoint only after all acceptance gates pass. Store:

- grid, fragment, core, buffer, and periodic-image geometry;
- canonical physical-face ownership and penalty controls;
- continuation history and final `lambda=1`;
- core-owned nodal orbitals, eigenvalues, occupations, density, Hartree, XC;
- pseudopotential and complete-operator fingerprint;
- residual, orthogonality, Hermiticity, face-balance, and fixed-point gates;
- orbital coverage manifest and projectability spectrum;
- local Wannier ownership, centers, spreads, gauge transforms, and closure
  defects when Wannier construction has passed.

Do not reuse the existing LCFO/WPW checkpoint schema by silently changing its
meaning.

## Failure handling

Fail collectively and transactionally for:

- invalid handoff density or electron count;
- incomplete candidate materialization;
- rank disagreement or metric rank loss;
- stale geometry/operator fingerprints;
- non-Hermitian volume, face, or nonlocal action;
- continuation rollback below the minimum step;
- orbital, density, or fixed-point nonconvergence;
- buffer-dependent core results above tolerance;
- insufficient manifest projectability or projection rank;
- occupied-projector loss during disentanglement;
- ambiguous core ownership or failed neighboring gauge alignment.

No automatic fallback to LCFO/WPW is allowed after explicit selection of this
route.

## Testing and physical gates

Begin with small deterministic fixtures:

- two fragments with an analytic SIPG matrix;
- single-versus-multiple-fragment Hermiticity and spectrum;
- internal face cancellation and physical-periodic boundary ownership;
- local-periodic auxiliary wrap that does not enter physical face coupling;
- core-only density and charge without double counting;
- DC handoff at multiple loose residuals;
- continuation advance, rollback, and minimum-step failure;
- unmixed fixed-point rejection;
- metric-rank and candidate-pool expansion;
- synthetic s/p and s/p/d coverage manifests;
- projected disentanglement with exact occupied inclusion;
- neighboring gauge alignment and loop-closure diagnostics.

Then run Si64 with:

- handoff residuals `1E-2`, `1E-3`, and `1E-4`;
- at least two buffer widths;
- at least two equivalent fragment decompositions;
- Tier 0 (`s+p`) and Tier 1 (`s+p+d`) manifests.

The first success criterion is an accepted self-consistent full-DG ground
state whose core density, energy, occupied spectrum, and face diagnostics are
stable against handoff tolerance, buffer, and decomposition. Only after that
gate passes should local Wannier results be interpreted.

Wannier acceptance additionally requires:

- full manifest projection rank;
- occupied-projector inclusion;
- stable centers and spreads with buffer;
- reproduction of the selected DG bands and occupied density;
- unique core ownership;
- neighbor-overlap consistency and small gauge-loop closure defects.

No result in this design authorizes RT or length-gauge work.

## Superseded work

The global GMRES comparison remains default-off and rejected by its B=6 gate.
The committed MINRES-QLP design and implementation plan were created but no
MINRES-QLP source implementation began. They are superseded by this
ground-state-first direction and should not be executed unless explicitly
reactivated.
