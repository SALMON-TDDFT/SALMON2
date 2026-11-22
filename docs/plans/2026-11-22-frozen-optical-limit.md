# Frozen-response optical limit implementation plan

**Goal:** Replace the minimum finite-q Gamma approximation in an opt-in,
stationary occupied-projector calculation by analytic insulating head/wings.
The user approved the staged design in the preceding exchange.

**Architecture:** Keep the legacy experimental finite-q path reproducible.
Add `yn_hsex_frozen_optical` (default `n`) for the new audited limit path.
Use only the occupied projector: diagonalize its small projected Hamiltonian,
check the out-of-subspace residual, and use the existing matrix-free optical
Sternheimer solver. No unoccupied eigenstates and no RT eigenstate projection
are introduced. Use full kinetic/pseudopotential/HSEX k derivatives.

1. Add independent spectral tests for the optical tensor, occupied-unitary
   invariance, nonstationarity rejection, nonfinite inputs/callbacks, and eta
   extrapolation including divergent/nonconvergent examples.
2. Add a native helper for the stationary occupied-space optical response and
   three-sample eta-to-zero extrapolation. Finite-q eta convergence precedes
   the Gamma-cell long-wavelength integration. Preserve both signs of q and
   the existing normalization and local-field convention.
3. Connect optional velocity requests to the existing root/owner MPI service;
   replace the complete Gamma cell, including wings, without double counting.
   Reject any failure collectively before publishing a new W.
4. Wire the opt-in input, documentation, and audited serial/MPI tests. Audit
   source before builds and all application objects before execution.
5. Build SALMON and compare an equilibrium Si snapshot with the existing
   optical path on the identical grid, cutoff, and occupied projector.

**Explicit scope boundary:** This first stage is the stationary insulating
limit. A coherent pumped projector is not made stationary by assigning
Rayleigh quotients to its orbitals. Such a snapshot must be rejected in this
mode until the nonstationary optical Ward/derivative response is validated.
Metallic static compressibility and optical Drude limits require separate
tests and are not inferred from a two-point eta comparison or an insulating
Gamma-cell fit. An unchanged legacy mode is not evidence for those limits.

**Validation:** Native independent finite-matrix spectral sums and MPI owner
layouts (including empty owners and failed remote velocity actions), finite
and NaN/Infinity rejection, full HSE-enabled/disabled GNU builds, and a small
equilibrium Si comparison. Fujitsu/accelerator builds are reported as untested
when unavailable. No existing remote production job/binary is changed.

## Validation completed

- Independent spectral optical tensor, complex occupied-unitary invariance,
  nonstationary projector rejection, finite/NaN/Infinity and failed callback
  rejection, and linear/quadratic eta extrapolation/pole rejection passed.
- Thirteen 4-MPI cases passed, including empty occupied/W owners, remote
  velocity failure, and eta nonconvergence. The 1/4-MPI optical W checksums
  agree. A too-small eta trial exposed roundoff in exactly zero channels;
  its failed logs are retained and no acceptance threshold was relaxed.
- GNU HSE-on and HSE-off builds and source/all-application-object audits passed.
  Fujitsu/accelerator compilation remains untested.
- Equilibrium Si: Nr=8^3, Nk=Nq=2^3, boxes=2^3, cutoff=1 Ha, 8 MPI/1 OMP,
  eta=0.01/0.005/0.0025 Ha. This is a code comparison, not mesh convergence.
  Existing corrected optical epsilon=14.44777895, new=14.4477789. Whole-W
  relative Frobenius error=3.7109706554e-7 against the corrected static solver;
  maximum absolute error=5.7878852822e-8 Ha. The maximum equation residual
  upper bound=8.02091991e-9 and full response time=438.557 s. RT kept 4 occupied
  orbitals with no instantaneous eigenstate projection.
- Final binary SHA256:
  `346d3b51ef1cadf0eb82e42a0a05a26747c561f2f4c48087f989e64b3b38e268`.
  Scratch: `/private/tmp/si-frozen-optical-limit-20261122/equilibrium-final`.
  Existing corrected static reference is `hybrid-reference` in that root.
  The GS reference passes screen_eta=0; repeating its input eta is not an
  independent eta convergence test. The earlier finite-shell 2.92 value is
  not the corrected reference and is not used as evidence of new accuracy.

The explicitly scoped nonstationary/metallic follow-up is still outstanding;
no strong-pump/Drude capability is claimed from these stationary tests.
