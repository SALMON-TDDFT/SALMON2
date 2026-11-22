# Pump metallic screening implementation plan

Goal: support dephased, nonthermal pump populations in static screening and collisionless probe response without an insulating Coulomb head or duplicate Drude current.

Architecture: retain finite-q fractional Sternheimer response. Evaluate eta and eta/2 with a shared target reduction, compare screened matrices, then fit the full box-projected W near zero using a Cartesian quadratic polynomial. Validate nested inversion-symmetric stencils and fit residual before integrating the polynomial over the reciprocal parallelepiped cell. This is a controlled smooth-metal extrapolation, not a universal q=0 formula: unresolved screening, directional nonanalyticity or an unstable dielectric matrix must fail. No isotropic carrier mass or fitted scattering rate is introduced. Integer-occupation snapshots retain the insulating branch. Add the derivative of the same nonlocal HSEX operator to projected intraband current.

Tech stack: Python/NumPy/SciPy/mpi4py, Fortran/MPI, existing FFT transpose.

1. Add failing analytic tests: anisotropic regular screened kernels, cell volume dependence, unresolved Coulomb singularity, nonfinite input, eta sensitivity. Add an explicit partial-band intraband Lindhard test.
2. Implement finite-q integration and convergence guards in tools/hsex/metal_head.py. Connect to prepare_w.py and paired-eta pump_screen.py. Preserve equilibrium outputs.
3. Add finite-difference current test, implement projected nonlocal velocity in kli_hsex.f90 and connect write_intraband_current. Test distributed layouts.
4. Run fractional/full-band tests, serial/MPI comparison, prohibition audits, isolated standard compiler build and representative small RT checks. Do not change queued remote jobs or binaries.
5. Document approximation limits, numerical eta versus scattering, diagnostics and validation evidence. No claims of dynamic W, scattering or arbitrary coherent pump density response.

User approved the static dephased scope on 2026-11-05. Existing suitable development checkout is reused; no unrelated changes are staged or committed.
