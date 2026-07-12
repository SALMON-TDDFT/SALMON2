# Wannier+PW Exp Propagation Design

## Goal

Complete the global Wannier+PW real-time propagation path by reproducing a
Full TDDFT reference calculation in the linear-response regime. The immediate
goal is numerical agreement with Full TDDFT for the selected small test system;
absolute physical fidelity of that system is outside this milestone.

## Physical and Numerical Scope

- Laser intensity: `1.0e11 W/cm^2`
- Photon energy: `1.55 eV`
- Pulse duration: `10 fs`
- Envelope: `sin^2`
- Polarization direction: `z`
- Time step: `dt = 2 a.u.`
- Wannier+PW propagator: exponential propagation in the mixed basis
- Reference: Full TDDFT started from the same converged ground-state provenance

The initial occupied states are not obtained by merely projecting the Full
TDDFT orbitals into the mixed basis. They are the occupied eigenstates of the
complete zero-field DG Hamiltonian in the Wannier+PW space, including all DG
interface/flux terms. For a non-orthogonal mixed basis, initialization solves
the generalized eigenproblem

`H_DG C = S C epsilon`.

### Basis dependence of the DG operator

For a fixed effective potential, the symmetric interior-penalty DG bilinear
form can be written schematically as

`a_DG(u,v) = sum_K <grad u,grad v>_K + <u,V_fixed v>
             - sum_F <{{grad u}},[[v]]>_F
             - sum_F <[[u]],{{grad v}}>_F
             + sum_F alpha_F <[[u]],[[v]]>_F`.

Consequently, the matrix is `H_DG[Phi]_(mu,nu)=a_DG(phi_mu,phi_nu)`.
The surface/flux contribution depends on the traces and normal derivatives of
the basis functions at fragment faces. It does not depend directly on the
time-dependent coefficient vector when the basis `Phi` is fixed.

Therefore two update levels must be distinguished:

- coefficient update with fixed Wannier+PW basis: keep the zero-field DG
  volume and surface matrices fixed and update only the external-field term;
- basis update: rebuild overlap, volume, nonlocal, Wannier-PW, PW-PW, and all DG
  surface/flux blocks before solving or propagating in the new basis.

The initial state must be a fixed point of the second level. Starting from the
imported converged potential and a trial Wannier+PW basis, repeatedly rebuild
the full DG operator, diagonalize it, reconstruct/update the basis from the
occupied subspace when the basis construction requires that update, and rebuild
the DG operator. Stop only after the occupied projector, eigenvalues, DG matrix,
and eigen-residual have converged. If the chosen global Wannier+PW basis is
defined once and is not regenerated from the DG eigenvectors, this loop reduces
to one rebuild and one generalized diagonalization; an artificial coefficient
iteration must not be introduced.

The local symmetry of individual fragments or Wannier functions is not a
constraint. For defect systems, local symmetry breaking is physical. A perfect
crystal may be used as a diagnostic to confirm that the complete Wannier space
reproduces the symmetry of the full system.

## Production Route

The first production target is the global Wannier+PW path. This removes
fragment ownership, halo exchange, and domain-boundary effects from the initial
validation. The mixed Wannier+PW Hamiltonian is propagated with an exponential
operator rather than a Taylor expansion.

Once the global calculation passes the Full TDDFT comparison, the same
formulation will be transferred to the fragment/distributed implementation.
Agreement between global and distributed calculations will then be treated as
a separate implementation test.

## Observables

The primary observable is the electronic polarization along the laser axis,
`Pz(t)`. The current is a derived observable,

`Jz(t) = dPz(t) / dt`.

The same numerical derivative and time window must be applied to both the
Wannier+PW and Full TDDFT polarization data. Current agreement is secondary and
must not replace direct comparison of `Pz(t)`.

## Acceptance Criteria

The primary acceptance criterion is a relative RMS error of at most 5 percent:

`rel_rms = rms(Pz_WPW - Pz_full) / rms(Pz_full) <= 0.05`.

The comparison window includes the pulse and a defined post-pulse interval.
Both signals must use the same time origin, sampling points, polarization sign,
and volume normalization.

Additional health checks are:

- negligible zero-field polarization drift;
- practically conserved electron number;
- stable norm/occupation behavior in the mixed basis;
- convergence as the PW count or cutoff is increased at fixed `dt = 2 a.u.`;
- consistent `Jz(t)` after applying an identical derivative procedure.

## Implementation Strategy

1. Establish a trustworthy Full TDDFT reference input and record its ground-state
   provenance, excitation mode, time step, pulse parameters, polarization
   convention, and volume normalization.
2. Reduce the existing experimental Wannier+PW `expdiag` branches to one explicit
   namelist-controlled production path. Required scientific choices must not
   remain environment-variable-only controls.
3. Assemble the complete zero-field Wannier+PW DG Hamiltonian and overlap,
   diagonalize `H_DG C = S C epsilon`, update the basis if and only if the basis
   construction depends on the resulting occupied subspace, rebuild every
   basis-dependent DG block, and iterate to a fixed point. Occupy the converged
   lowest physical eigenstates and verify the eigen-residual and field-free
   stationarity before propagation.
4. Verify the Hamiltonian blocks, position operator, and exponential update
   independently before interpreting the full waveform.
5. Run a PW-count/cutoff sequence at fixed time step and compare `Pz(t)` against
   Full TDDFT using one analysis tool and one manifest schema.
6. Promote the smallest converged Wannier+PW basis satisfying the 5 percent gate
   to the reference global configuration.
7. Only after the global gate passes, validate fragment/distributed propagation
   against the global result.

## Failure Handling

Discrepancies are classified before code changes:

- a mismatch at the first step points to full-DG eigenseed or operator
  construction;
- nonzero field-free motion points first to an incorrect DG eigenseed, missing
  DG interface terms, or an inconsistent overlap metric;
- apparent convergence obtained without rebuilding face terms after a basis
  refresh is invalid, because it diagonalizes a matrix belonging to the old
  basis;
- zero-field drift points to basis non-closure, non-Hermiticity, or propagation
  inconsistency;
- PW-count dependence points to incomplete mixed-space convergence;
- agreement early in time followed by phase drift points to eigenvalue or
  propagator error;
- global/distributed disagreement points to ownership or communication logic,
  not to the physical basis definition.

Rank-local finite-value and magnitude checks should stop on the rank that first
produces invalid data, before collective reductions obscure its origin.

## Deliverables

- one documented Full TDDFT reference input;
- one namelist-driven global Wannier+PW Exp input family;
- automated `Pz(t)` and derived `Jz(t)` comparison;
- a PW convergence table including the relative RMS gate;
- regression tests for the production Exp path;
- an eigen-residual and field-free stationarity report for the occupied
  full-DG Wannier+PW initial states;
- documentation of the accepted basis and remaining limitations.
