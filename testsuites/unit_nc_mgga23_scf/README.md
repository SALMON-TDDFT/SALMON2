# Full gKS SCF damping

The material SCF must damp orbital-dependent local/tau/current coefficients;
ordinary density mixing does not change fields rebuilt directly from spinors.
The helper is used only for the SCF orbital solve. Energy band expectation and
core field response remain unmixed, and the final physical payload is restored
by XC refresh. The fixed-point residual covers all coefficient fields, including
changes invisible to total charge density, using global maxima across ranks.

`test_scf.f90` checks an otherwise unstable tau-only map, spin-current-only
nonconvergence and preservation of physical diagnostics. The runner audits
source and every linked application/test object. Use `run.py --build-dir BUILD
--libs LIBRARIES` with a current serial CPU build. A real material SCF remains
necessary to assess convergence; this algebraic regression is not a material
accuracy test.

The LSDA regression uses a transverse (1,i)/sqrt(2) spinor. After simple_dm,
its eigenchannel densities must remain1,0 rather than laboratory diagonals
1/2,1/2. A charge-invisible imaginary off-diagonal density change must fail
the full-spin criterion. Child cases test NaN/Infinity rejection for local,
tau, current, history and complex-density inputs.
