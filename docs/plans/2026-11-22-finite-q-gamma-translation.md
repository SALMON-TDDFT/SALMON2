# Finite-q/eta Gamma closure and grid-translation verification

User request: verify the static order of q/eta limits and advance the
nonstationary Gamma-cell connection. The identified spurious low-energy mode
comes from broken continuous translation symmetry of the spatial grid, not
an arbitrary high-energy/unoccupied-state cutoff.

## Design

1. Measure integer/sub-grid joint translations of occupied states and sampled
   local potentials. Integer-grid covariance is a numerical control; fractional
   translations expose collocation/product aliasing. Translate nonlocal frozen
   operators unitarily when isolating this local-grid defect. Do not claim this
   component test measures all pseudopotential-resampling errors.
2. At fixed nonzero q, refine eta and check the response before continuing q.
   Use independent insulating/metallic models to test the noncommuting point
   limits and distinguish them from the integrable Gamma-cell average.
3. Fit the full box matrix F(q)=|q|^2 W(q,0), not just a scalar dielectric head.
   Nested 26/52-point quadratic fits retain singular head, wings/odd terms and
   regular terms. A spherical Gamma-cell average is 3*A/qc^2+Tr(C)/3.
   In a metal A can vanish while C survives; keeping only A is incorrect.
4. Connect the candidate to the frozen-W path behind an explicit option.
   Require eta convergence, nested-q agreement, finite/Hermitian/positive
   candidates and translation controls before publication. Insufficient q data
   or unresolved translation contamination must reject, preserving the old W.
5. Test actual Si q/eta response and translation diagnostics with the existing
   GS data; never collect GS payloads. Start small and report where the data
   cannot establish a q-to-zero limit.

Main RT stays on the real-space grid with its occupied orbitals. No
instantaneous eigenstate projection, thermal replacement, arbitrary removal
of low-energy states, or forced cancellation of the pumped Drude term.

## Verification

- Exact integer translations, norm conservation, band-limited fractional
  covariance, deliberate alias-induced pinning, and improvement on refinement.
- Scalar IEEE finite/NaN/Infinity guards; Fourier Nyquist treatment of real
  potentials must be explicit and preserve Hermiticity.
- Exact polynomial matrix Gamma averages; insulating and metallic limits;
  failed eta/q fits and insufficient samples reject.
- Audited source/object builds, native serial/MPI checks, and matched Si controls.

Reference on continuous-translation grid errors (context, not a derivation of
our test): https://pubmed.ncbi.nlm.nih.gov/32456435/ .

## Implemented and verified in this stage

- `screen_shift`: periodic trigonometric translation, complex-unitary orbitals,
  explicit real-Nyquist convention, all three axes, non-power-of-two meshes,
  and scalar IEEE validation.
- The existing opt-in coherent optical diagnostic now reports local-potential
  collocation covariance and expectation-value changes for integer/half-grid
  translations of the actual occupied orbitals. Integer covariance is required;
  fractional defects are reported, never subtracted from the response.
- `screen_q2_head`: full-matrix nested quadratic fit and spherical cell moments,
  including the surviving metallic regular term. It is a numerical building
  block only: it is not yet selected by the nonstationary full-W driver.
- Hydrodynamic finite-q/eta tests reject unresolved eta sequences and recover
  the fixed-q static response before taking the small-q limit.
- Native tests and a GNU full build passed. Two matched Si Nr=8^3, Nk=2^3
  pulse/zero-field jobs passed with 8 MPI. Their half-grid local expectation
  change is about 0.02457 eV, while integer controls are at roundoff.

Remaining: translate/resample the full nonlocal pseudopotential and inspect
low-energy response modes, obtain sufficient small-q Si data, and apply both
translation and eta/q acceptance gates before publishing a new Gamma cell.
The Nr=8^3 component diagnostic is not convergence evidence for Nr=12^3,
and it does not identify which specific response pole is spurious.
