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

- coefficient update with fixed Wannier+PW basis: keep the kinetic, ionic, and
  DG surface matrices fixed, while updating the density-dependent potential
  matrices and external-field term;
- basis update: rebuild overlap, volume, nonlocal, Wannier-PW, PW-PW, and all DG
  surface/flux blocks before solving or propagating in the new basis.

This milestone uses one fixed global Wannier+PW basis. The basis is constructed
once from the converged reference ground-state potential and is not regenerated
from the DG eigenvectors. In that fixed basis, rebuild the complete DG operator
once and solve one generalized eigenproblem. A basis-self-consistency loop is
outside this milestone because the current production path has no validated
Wannier regeneration, gauge tracking, and operator-generation protocol.

The potential is fixed only while constructing this initial eigenseed. During
real-time TDDFT, the basis and its kinetic/DG surface matrices remain fixed, but
the density-dependent Hartree and exchange-correlation potential matrices are
updated from the propagated density. Thus the field-on Hamiltonian is

`H(t) = H_kin+DG + V_ion + V_H[n(t)] + V_xc[n(t)] + E(t) Z`.

Freezing `V_H+V_xc` during real-time propagation would define an independent-
particle approximation and is not an acceptable Full TDDFT reference match.

The local symmetry of individual fragments or Wannier functions is not a
constraint. For defect systems, local symmetry breaking is physical. A perfect
crystal may be used as a diagnostic to confirm that the complete Wannier space
reproduces the symmetry of the full system.

## Production Route

The first production target is the global-ownership Wannier+PW path. It uses the
same physical fragment boundaries and the same DG interface/flux operator as
the later distributed path, but removes MPI coefficient ownership and halo
exchange effects from the initial validation. The DG surface terms must not be
removed or folded away merely because ownership is global. The mixed
Wannier+PW Hamiltonian is propagated with an exponential operator rather than a
Taylor expansion.

For a non-orthogonal mixed basis, propagation solves

`i S dC/dt = H(t) C`.

Let `X` satisfy `X^H S X = I`. One midpoint exponential step is

`C(t+dt) = X exp[-i X^H H(t+dt/2) X dt] X^H S C(t)`.

The same overlap metric and orthogonalization convention must be used for
initial diagonalization, coefficient projection, exponential propagation, and
observable evaluation. The numerical propagator must satisfy
`U^H S U = S` within tolerance.

The midpoint density-dependent Hamiltonian is obtained with a predictor-
corrector: predict the state with the input Hamiltonian, reconstruct the
midpoint density and `V_H+V_xc`, rebuild only the basis-fixed potential matrix
blocks, and repeat until the midpoint density or Hamiltonian is converged.

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
- generalized eigen-residual and `S`-orthonormality of the initial states;
- `S`-unitarity of every exponential update;
- midpoint density/potential convergence for each real-time step;
- Hermiticity of the complete length-gauge position operator in the `S` metric.

## Implementation Strategy

1. Establish a trustworthy Full TDDFT reference input and record its ground-state
   provenance, excitation mode, time step, pulse parameters, polarization
   convention, and volume normalization.
2. Reduce the existing experimental Wannier+PW `expdiag` branches to one explicit
   namelist-controlled production path. Required scientific choices must not
   remain environment-variable-only controls.
3. Construct one fixed global Wannier+PW basis from the converged reference
   potential. Rebuild the complete zero-field DG Hamiltonian and overlap in that
   basis, solve `H_DG C = S C epsilon` once, and occupy the lowest physical
   eigenstates. Verify the eigen-residual and field-free stationarity.
4. Implement the exponential update in the same overlap metric and verify
   `S`-unitarity with small-matrix numerical tests.
5. Verify the complete WW, WP/PW, PP, and DG-interface contributions to the
   length-gauge position operator before interpreting a field-on waveform.
6. During RT, keep the basis and DG surface matrix fixed but update the
   density-dependent potential through a midpoint predictor-corrector.
7. Run a PW-count/cutoff sequence at fixed time step and compare `Pz(t)` against
   Full TDDFT using one analysis tool and one manifest schema.
8. Promote the smallest converged Wannier+PW basis satisfying the 5 percent gate
   to the reference global configuration.
9. Only after the global gate passes, validate fragment/distributed propagation
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
- norm drift in coefficient space must first be evaluated in the `S` metric;
- agreement with a frozen-potential calculation is not Full TDDFT agreement;
- a field-on mismatch with a stationary field-off eigenseed points next to the
  position operator or midpoint density/potential update;
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
