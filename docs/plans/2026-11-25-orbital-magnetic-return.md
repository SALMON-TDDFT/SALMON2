# Spinless orbital magnetic return

Goal: return only the magnetic part of the existing spatial molecular Hamiltonian to FDTD, alongside the existing electric polarization current.

1. Test the adjoint of H -> A/c (including the existing nonlocal path quadrature), Yee curl divergence and atomic numerical rejection.
2. Implement shared Fortran magnetic adjoints; compute -d<E>/dH from occupied wavefunctions, local paramagnetic/diamagnetic terms and nonlocal projector phases. One complete isolated molecular domain/worker initially; no spin contribution. Read-only query, no propagation or added checkpoint state.
3. Expose the read-only H-grid magnetization and its E-grid c*curl(M). Require Yee origin alignment for current conversion; never add a full orbital transport current to dP/dt.
4. Compare the response against finite differences of the actual frozen-state Hamiltonian, including nonlocal phases. Verify process-exit restart and inject the returned currents through existing checked FDTD edge injection, comparing continuous/restarted output.
5. Update Japanese/English guides, AI JSON, checklist and measured validation record; no physics CI.

This implements the user-approved variational design. Magnetic input denotes vacuum B=H. The finite-difference audit is an explicit validation operation, not production current evaluation. All numerical work stays in Fortran.

### 2026-11-25: spinless orbital magnetic return

- Done: read-only Fortran M=-delta<E>/delta B with local paramagnetic/diamagnetic and nonlocal projector phase terms; Yee c*curl(M), adjacent-step interval averaging and one-time electric/magnetic current sum. Electric coupling remains length gauge.
- Verified: actual Hamiltonian derivative absolute error7.08e-12; independent adjoint/work and divergence checks; atomic invalid/overflow rejection; C2H2 direct process-exit restore2->4 reproduces M, interval currents and ordinary observables without evolution.
- Two FDTD ranks plus bookkeeping: regenerated restored currents and field restart3 reproduce continuous arrays exactly. Magnetic-only E is nonzero; electric+magnetic superposition agrees within4.24e-22. Source-distribution-to-wheel, installed99targeted tests, installed molecular and MPI restart tools passed. See docs/orbital-magnetic-return-validation-2026-11-25.md.
- Remaining: general self-consistent magnetic own-vacuum coordinator, spin magnetization, full composed charge-continuity/energy validation, physical convergence and scale tests. This closes the response-return component, not those broader items.
