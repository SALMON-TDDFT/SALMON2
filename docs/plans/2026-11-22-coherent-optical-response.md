# Coherent optical response implementation plan

**Goal:** Extend the already approved frozen-H research path with a conserving,
occupied-orbital electric-impulse optical diagnostic for a nonstationary state.
This executes the evolving-background validation in
`2026-11-22-frozen-screening-research.md`; it does not relax the stationary gate.

**Architecture:** Let L X=[H,X] and R(z)=(z-L)^(-1). For an electric impulse,
the velocity-gauge vector potential steps. Its derivative of the evolving
background density has Laplace transform i R [v,R rho]. The conductivity is

    sigma(z) = i/volume * Tr(v R [v,R rho] + a R rho),
    v=dH/dA, a=d2H/dA2, Im(z)>0.

The second term is the contact/diamagnetic term. Both resolvents are necessary
when [H,rho] differs from zero. Replacing the inner R rho by rho/z is only a
stationary approximation. Do not Hermitize either finite-frequency response.
Use the existing low-rank Fortran solver and operator actions, without empty
states or occupation projection. Report residuals of both solves separately;
they alone do not bound the final observable error at small Im(z).

1. Add independent finite-matrix spectral and orbital-time references. Test
   coherent and stationary states, density/velocity-gauge Ward identity with
   a bounded position operator, the exact zero-velocity Drude/contact case,
   norm/charge conservation, and rejection of nonfinite/undefined actions.
2. Implement the native two-resolvent core and general-operator traces.
3. Add exact kinetic/nonlocal-pseudopotential/HSEX second A derivatives and
   compare them with finite-A operator differences before SALMON coupling.
4. Add an opt-in RT snapshot diagnostic, preserving the main RT orbitals and
   active W. Test the previously rejected Si E=1e-6 snapshot. Keep the
   stationary optical W-update mode unchanged.
5. Decide Gamma static-W continuation only after finite-frequency charge
   conservation and the order of eta/q limits have been checked. Optical
   conductivity is not automatically a static metallic compressibility.

**Validation:** Audit sources before builds and all application objects before
new executable use. Test GNU serial/MPI, empty owners and remote failures.
Keep failed inputs/logs and do not change existing production jobs/binaries.
Collect only response results/records in calc2026, without GS payloads.

**References:** The elementary derivation above is independently checked by
orbital time propagation; it is not claimed as a formula from these sources.
Octopus developer linear-response derivation:
https://octopus-code.org/documentation/15/developers/code_documentation/linear_response/
Adiabatic/nonadiabatic distinction: Perfetto et al., PRB 92, 205304 (2015),
https://arxiv.org/abs/1507.01786 . A single-frequency snapshot inverse is not
a solution of the general two-time nonequilibrium Dyson equation.

## Implemented and verified

- Native `hsex_coherent_optical` evaluates both resolvents, general operator
  traces, and the contact contribution without empty states. Independent
  spectral sums and orbital RK4 agree to 1.3e-15 and 1.1e-12 in small models.
  A bounded-position Ward identity, fractional-occupation Drude/contact case,
  charge conservation, capped solves and invalid/undefined actions are tested.
- `hsex_projected_derivative` supplies HSEX first/second A derivatives. The
  former auxiliary callbacks lacked 1/hvol, unlike `hsex_add`; an independent
  finite-A regression failed before correction and passes afterward. The
  shared helper is now used in both W-update callbacks and the new probe.
  Pseudopotential second derivatives retain periodic images and spin components.
- The opt-in coherent RT adapter partitions k work directly over MPI and
  writes conductivity, IP epsilon and selected microscopic density response.
  It does not replace W. Native tests retain NaN/Infinity rejection.
- Zero-source G-shifted occupied columns seed the Krylov spaces. This is an
  exact source refactor, checked against spectral/time references and the
  unseeded Si response. It changes neither rho nor the residual threshold.
- Si Nr8^3, Nk2^3, 8/4 MPI, 4 propagated occupied states, E=1e-6 au: all
  diagnostic trajectories complete 1600 steps. The state previously rejected
  by the stationary mode (rate=1.57e-6 Ha) is accepted by the coherent probe.
  Actual finite-A first/second derivative errors are <=1.53e-9 / 1.94e-7.
  Maximum uniform-charge error is 2.51e-9; solver bound <=4.33e-7 (tol1e-6).
  The MPI4/8 IP epsilon maximum discrepancy is <9e-10.
- Observed single-run response times: seeded pulse/zero MPI8 19.13/21.05 s,
  versus unseeded 74.74/318.25 s. These are reference timings with possible
  co-resident workloads, not a controlled repeated performance benchmark.
  Main RT/current changes only at roundoff. Existing HSEX action, current,
  energy, scale and checkpoint regression passes with 2 MPI.

## Outstanding physical closure

The raw low-frequency current response has a large pole already in the
unexcited coarse-mesh reference: epsilon_IP(i*0.01)=-523.913 and
epsilon_IP(i*0.005)=-2141.240. The actual-H derivative checks pass, and the
pole is nearly identical before/after the weak pump. Do not identify it as
photoinduced Drude response or publish it as the static-W Gamma cell.
Separate the equilibrium conductivity-sum-rule residual and converge finite q
and eta in their appropriate order before connecting general nonequilibrium W.
The corrected static insulating long-wavelength W is a separate established
calculation and is not invalidated by this raw IP current diagnostic.

An additional independent stationary-current reference now measures the
contact+paramagnetic coefficient D_raw=-0.0042902245009 and static IP
epsilon=15.279026827 on this coarse Si grid. Removing only this equilibrium
reference pole in diagnostic post-processing gives epsilon(i*eta)=15.21254
and 15.26234, whose quadratic eta extrapolation is 15.27894 (relative
5.59e-6 from the independent static IP value). This is not the macroscopic
local-field-corrected 14.44778 of the established static W path.
The pole/reference is not silently removed in the native calculation.
Pump-minus-zero epsilon still changes from about 0.00140 to 0.00565 when eta
is halved; this cannot be treated as an eta-converged insulating static head.
It must not be relabeled as a measured carrier Drude weight without finite-q
and transient-response checks.

Source on this distinction and common optical-current pitfalls:
https://arxiv.org/abs/1701.06326 .

GNU HSE-on/off builds and prohibition audits pass. Fujitsu/accelerator are
untested. Original generated inputs that exceeded the 256-character reader
were soft-stopped by SALMON; they are preserved separately and never counted
as completed calculations. Final scripts require both end SALMON and PASS.
