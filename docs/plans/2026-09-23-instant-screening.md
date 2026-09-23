# Instantaneous screening implementation plan

Goal: opt-in instantaneous A,J,E,P estimator for strong laser runs, without time averaging.
Architecture: classical transverse a and E (exclude XC to avoid tautological E_xc/P),
trapezoidal P=-integral j, vector least-squares K=(a.j-E.P)/(a.a+E.E/omega^2).
Use alpha_eff=alpha0/[1+s*4*pi*max(K-K0,0)/omega^2]. This positive bounded closure is
phenomenological, not the real-frequency Drude dielectric function or a derived finite-q kernel.
K0 must be calibrated from a weak laser reference at the same frequency and geometry;
it is not inferred from the initially unexcited state. Weak-limit preservation requires this
calibration and must be tested; it does not follow automatically for arbitrary interband response.
Tech stack: Fortran module, SALMON real-time integration, Python/Fortran tests.

1. Add failing analytic tests for quadrature estimator, exact zero crossings, hold near no field,
   bounded correction, and zero strength. Implement pure estimator.
2. Add opt-in functional keywords, validation and log. Restrict first version to atomic units,
   transverse tddft_pulse and a non-impulsive pump; omega supplied explicitly.
3. Integrate P, endpoint classical E, estimator and alpha. Preserve fixed-alpha path exactly.
   Save P,K,alpha in version-2 XC restart; accept version-1 only with fixed screening.
   Append diagnostics only for the new mode, preserving legacy output layout.
4. Build and run existing tests plus active-feedback/restart/rejection integration checks.
5. Add sample and documentation with calibration workflow and limitations. No k-grid scan.

Outcome: implemented and tested. The classical-field estimator has a material pulse-tail limitation:
residual P with external E approaching zero contaminates K and the held post-pulse alpha. The 4³
trial shows this explicitly. Retain as opt-in research code; do not claim physical accuracy improvement.
