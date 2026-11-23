# Causal molecular feedback regression (development fixture)

Use the benchmark-and-mms planner's general-family plan as a checklist, not a
physical validation certificate. Its general model is intentional: neither the
Maxwell/dipole coupled system nor SALMON is a built-in family. Do not start a
production spatial-refinement campaign; the user owns physical convergence.

First verify an implicit-midpoint unit-mass/unit-charge harmonic oscillator
against its free analytic trajectory at four small dt values. Check second-order
error and the discrete identity ΔH = E_drive·Δp for constant drive within a step.
This is a test fixture, not a new molecular model or TDDFT propagator.

Native regression: a 4-rank job splits into two independent 2-rank contexts,
each containing one native field rank and one bookkeeping rank. Both native
fields use the same vacuum grid, source deposition and observation map. The
environment has prescribed illumination; the reference has zero illumination.
Only one molecule is present. Each field library initializes once per process.

Start with a declared one-dt hold and zero source current. After field time n*dt
(n>=1), subtract the own-source reference, evolve p_n,v_n to p_(n+1),v_(n+1)
under E_drive(n), form j_(n+1/2)=(p_(n+1)-p_n)/(V*dt), and advance both fields
with that same current. Thus field n+1 is never used to determine current needed
for its own update. Do not call this schedule the isolated-TDDFT midpoint method.

Run a separate illumination-only baseline with zero molecular current. Require
selected E to agree with that baseline within rtol=1e-10, atol=1e-14; require a
nonzero reference self-field so cancellation is not a trivial zero test. Check
oscillator ΔH−E_drive·Δp within 1e-12 absolute (atomic units) per step. Record
E_total/end-point-averaged source work, reference work and selected-drive work
separately. Record the temporal splitting discrepancy; do not claim full native
field-energy conservation, closed-surface flux or vacuum radiation reaction.

Save accepted E/p/j with times, input/library hashes and numerical residuals.
No public two-way native driver or isolated TDDFT support is implied. Completing
this fixture advances R4b-schedule but does not close R4b or implement scattering
in a material environment. Use at most 4 MPI ranks, OMP/BLAS one thread.
