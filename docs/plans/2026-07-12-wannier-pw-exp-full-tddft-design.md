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
- Reference: Full TDDFT started from the documented conventional ground state
- Material scope: gapped systems with integer occupations
- Exchange-correlation scope: LDA

The initial occupied states are not obtained by merely projecting conventional
Full TDDFT orbitals into the mixed basis, nor by a one-shot diagonalization in a
frozen conventional potential. They are obtained from a separate self-
consistent DG-DC ground-state calculation in the fixed Wannier+PW basis. Each
DG-DC SCF iteration solves the complete zero-field DG Hamiltonian, including all
DG interface/flux terms. For a non-orthogonal mixed basis, it solves

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
once from the converged conventional reference ground state and is not
regenerated from the DG eigenvectors. Its overlap, kinetic, ionic nonlocal, and
DG surface/penalty matrices are basis-fixed. The DG-DC SCF updates the density,
Hartree and exchange-correlation potentials, their mixed-basis matrix elements,
and the generalized eigensystem until the input and output density agree.

The DG-DC iteration is

`n(k) -> V_eff[n(k)] -> H_DG[n(k)] -> (epsilon(k),C(k)) -> n(k+1)`.

It converges the density/potential, total energy, occupied projector, and
generalized eigen-residual. The converged condition is
`H_DG[n0] P0 S = S P0 H_DG[n0]`, where
`P0=C_occ f C_occ^H`. A Wannier basis-regeneration loop remains outside this
milestone.

### DG-DC energy and potential evaluation

DG-DC uses the same Kohn-Sham LDA total-energy expression as conventional DFT;
no separate DG energy formula is introduced. DG kinetic, interface, and penalty
terms are already included in the DG Kohn-Sham eigenvalues. With integer
occupations `f_i`, the total energy is

`E_tot = sum_i f_i epsilon_i - E_H - integral(n V_xc)
         + E_xc^LDA + E_ion-ion`.

The LDA terms are local and are evaluated on disjoint fragment core grids:

`E_xc^LDA = sum_K integral_(Omega_K) n epsilon_xc(n)`,

`integral(n V_xc) = sum_K integral_(Omega_K) n V_xc`.

Buffer and halo points are excluded from these energy integrals, and every
global grid point has exactly one fragment owner. For spin-polarized extensions,
the same ownership rule applies to the spin densities.

Hartree remains a global operation. The complete mixed-basis density, including
WW, WP/PW, and PP contributions, is assembled on the global grid, the existing
SALMON Poisson/Hartree solver produces `V_H`, and `V_H` is redistributed to the
fragments. The Hartree energy is evaluated as

`E_H = 0.5 integral_Omega n V_H`.

The production energy uses the eigenvalue-sum expression above. A component-
expectation reconstruction may be retained only as a diagnostic consistency
check.

During real-time TDDFT, the same basis and kinetic/DG surface matrices remain
fixed, but the density-dependent Hartree and exchange-correlation potential
matrices are updated from the propagated density. Thus the field-on Hamiltonian is

`H(t) = H_kin+DG + V_ion + V_H[n(t)] + V_xc[n(t)] + E(t) Z`.

Freezing `V_H+V_xc` either during DG-DC convergence or real-time propagation
would define an independent-particle approximation and is not an acceptable
Full TDDFT reference match.

The local symmetry of individual fragments or Wannier functions is not a
constraint. For defect systems, local symmetry breaking is physical. A perfect
crystal may be used as a diagnostic to confirm that the complete Wannier space
reproduces the symmetry of the full system.

## Symmetry-adapted Wannier construction and large-system scaling

For a perfect crystal, the fixed Wannier sector is constructed with
symmetry-adapted Wannier functions (SAWF). The guarantee is conditional: the
selected Bloch/band subspace must be closed under the supplied crystal symmetry
group, and the supplied `D_band` and `D_wann` matrices must satisfy the group
representation and covariance tests. Merely enabling Wannier90
`site_symmetry` is not an acceptance criterion. Complete symmetry-related
orbital shells are retained; a band or projection subset rejected by the
closure test must not be passed to DG-DC.

The small validation system uses a global SAWF construction as the reference
route. A production-scale supercell must not run one monolithic global
Wannier90 localization. Instead it uses a hierarchical route:

1. partition the system into local environments with a converged real-space
   buffer and identify their orbits under the actual supercell symmetry group;
2. generate SAWFs only for one representative of each symmetry-equivalent bulk
   environment;
3. replicate them by translations and the validated `D_wann` action;
4. regenerate Wannier functions for every symmetry-inequivalent local
   environment, including defects, material interfaces, free surfaces/vacuum
   boundaries, and amorphous neighborhoods and their buffers;
5. align independently generated neighboring bases by overlap-based unitary
   Procrustes/parallel-transport transformations before constructing overlap,
   Hamiltonian, and DG face blocks;
6. assemble only sparse local WW, WP, PP, overlap, and DG face blocks globally.

All non-bulk calculations use the symmetry group of the actual supercell, not
a parent-crystal group. Thus a bulk region can reuse SAWF templates within that
same supercell while physical local symmetry breaking is retained near defects,
interfaces, surfaces, vacuum boundaries, and amorphous regions. Template reuse
across different supercells is forbidden, including geometrically similar
supercells from separate calculations.

The optional `wannier_sawf_structure_class` namelist declaration (`auto`,
`crystal`, `defect`, `interface`, `surface`, or `amorphous`) is an intent and
reuse-ceiling control, not a symmetry guarantee. It cannot merge environments
rejected by exact fingerprints or the actual supercell group. Class-specific
checks fail early on inconsistent input; amorphous mode conservatively disables
reuse unless exact actual-group equivalence is demonstrated.
Every cached template carries a complete supercell fingerprint and a local
core+buffer environment fingerprint, together with geometry, pseudopotential, grid, band-window,
projection-shell, symmetry-group, buffer, and generator-version fingerprints.
Any mismatch forces local regeneration rather than silent reuse.

Before this hierarchical basis may initialize DG-DC, it must pass four gates:

- symmetry closure and `D_band`/`D_wann` representation residuals;
- neighbor-overlap rank and gauge-stitching residuals;
- buffer convergence of local orbitals and all DG face blocks;
- on a small system, equivalence to the monolithic global SAWF reference for
  the occupied projector and assembled mixed-basis operator.

The PW enrichment is not localized by SAWF. Its global or distributed storage
is handled separately, while its WP coupling uses the final stitched Wannier
gauge. The large-system route changes basis generation and ownership only; it
must reproduce the same accepted DG operator contract as the global route.

## Production Route

The first production target is the global-ownership Wannier+PW path. Here
"global" refers to coefficient ownership, not automatically to globally
continuous basis support. Before implementation, the operator contract must
define whether Wannier and PW functions are fragment-restricted or globally
supported, and therefore which jumps and DG face blocks exist. Global ownership
must reproduce that accepted physical operator while removing MPI ownership and
halo-exchange effects from the first validation. The mixed Wannier+PW
Hamiltonian is propagated with an exponential operator rather than a Taylor
expansion.

For a non-orthogonal mixed basis, propagation solves

`i S dC/dt = H(t) C`.

Let `X` satisfy `X^H S X = I`. One midpoint exponential step is

`C(t+dt) = X exp[-i X^H H(t+dt/2) X dt] X^H S C(t)`.

The same overlap metric and orthogonalization convention must be used for
initial diagonalization, coefficient projection, exponential propagation, and
observable evaluation. The numerical propagator must satisfy
`U^H S U = S` within tolerance.

The midpoint density-dependent Hamiltonian is obtained with a predictor-
corrector. Every corrector iteration restarts from the same input state `C_n`:

`C_(n+1)^(k+1) = U[H(n_mid^(k))] C_n`,

`n_mid^(k) = (n[C_n] + n[C_(n+1)^k])/2`.

Only the trial midpoint density/potential is iterated; the exponential must not
be accumulated repeatedly within one time step. Stop if the midpoint iteration
does not converge.

Once the global calculation passes the Full TDDFT comparison, the same
formulation will be transferred to the fragment/distributed implementation.
Agreement between global and distributed calculations will then be treated as
a separate implementation test.

## Observables

The primary observable is the induced electronic polarization on one continuous
periodic branch,

`Delta_Pz(t) = Pz(t) - Pz(0)`.

The branch is tracked continuously using the same Wannier/Berry convention in
both calculations, or reconstructed in both calculations from one consistently
defined current. The current is a derived observable,

`Jz(t) = d Delta_Pz(t) / dt`.

The same numerical derivative and time window must be applied to both the
Wannier+PW and Full TDDFT polarization data. Current agreement is secondary and
must not replace direct comparison of `Delta_Pz(t)`.

## Acceptance Criteria

The primary acceptance criterion is a relative RMS error of at most 5 percent:

`rel_rms = rms(Delta_Pz_WPW - Delta_Pz_full) / rms(Delta_Pz_full) <= 0.05`.

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
- Hermiticity of the derived length-gauge position operator in the `S` metric;
- DG-DC density, potential, total-energy, occupied-projector, and eigen-residual
  convergence;
- equality of fragment-summed LDA-XC integrals with a direct global-grid
  diagnostic;
- unique core-grid ownership with no buffer/halo double counting;
- explicit separation of mixed basis dimension, retained eigenbasis dimension,
  propagated occupied-orbital count, and occupation weights;
- recorded overlap cutoff, condition number, discarded metric directions, and
  final effective dimension.

## Implementation Strategy

1. Freeze the SAWF symmetry contract and implement the scalable representative-
   environment construction, cache, replication, symmetry-inequivalent local
   regeneration,
   gauge stitching, and global-reference equivalence gates.
2. Freeze the DG trial-space, face-term, PW-support, and periodic length-gauge
   operator contract. Derive rather than assume any position-interface term.
3. Extract common generalized-eigen and `S`-metric Exp algebra usable by GS and
   RT.
4. Extract a common complete WW+WP/PW+PP density builder and verify electron
   number in both real-space and overlap metrics.
5. Validate fragment-core LDA-XC integration and global Hartree plumbing.
6. Solve the fixed-basis DG-DC SCF problem and converge density, potential,
   energy, occupied projector, and generalized residual.
7. Serialize an operator-complete, fingerprinted DG-DC checkpoint shared by GS
   and RT.
8. Prove field-off stationarity and operator identity across the DG-DC to RT
   handoff.
9. Implement one namelist-driven midpoint Exp production route, converting or
   removing all result-changing environment controls on that route.
10. Implement and validate only the length-gauge observable accepted by the
   operator contract.
11. Create the seven-input Stage 2D matrix and provenance manifest.
12. Run the PW convergence and `Delta_Pz` 5 percent Full TDDFT comparison.
13. Complete regression/documentation, then plan distributed ownership as a
   separate milestone.

## Failure Handling

Discrepancies are classified before code changes:

- a mismatch at the first step points to DG-DC/RT handoff, operator
  construction, occupation, or metric convention;
- nonzero field-free motion points first to an incorrect DG eigenseed, missing
  DG interface terms, or an inconsistent overlap metric;
- apparent convergence obtained without rebuilding face terms after a basis
  refresh is invalid, because it diagonalizes a matrix belonging to the old
  basis;
- norm drift in coefficient space must first be evaluated in the `S` metric;
- agreement with a frozen-potential calculation is not Full TDDFT agreement;
- zero-field relaxation immediately after RT start means the DG-DC state and RT
  potential update do not represent the same stationary problem;
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
- automated continuous-branch `Delta_Pz(t)` and derived `Jz(t)` comparison;
- a PW convergence table including the relative RMS gate;
- regression tests for the production Exp path;
- an eigen-residual and field-free stationarity report for the occupied
  full-DG Wannier+PW initial states;
- documentation of the accepted basis and remaining limitations.
