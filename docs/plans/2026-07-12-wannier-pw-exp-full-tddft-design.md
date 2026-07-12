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
3. Verify the mixed-basis initial state, Hamiltonian blocks, position operator,
   and exponential update independently before interpreting the full waveform.
4. Run a PW-count/cutoff sequence at fixed time step and compare `Pz(t)` against
   Full TDDFT using one analysis tool and one manifest schema.
5. Promote the smallest converged Wannier+PW basis satisfying the 5 percent gate
   to the reference global configuration.
6. Only after the global gate passes, validate fragment/distributed propagation
   against the global result.

## Failure Handling

Discrepancies are classified before code changes:

- a mismatch at the first step points to initial projection or operator
  construction;
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
- documentation of the accepted basis and remaining limitations.
