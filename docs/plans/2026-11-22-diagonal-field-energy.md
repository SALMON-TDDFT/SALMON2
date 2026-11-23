# Diagonal constitutive energy observations

Continue the accepted six-face energy diagnostic beyond unit-background scalar LD.
First complete positive diagonal electric/magnetic backgrounds and diagonal LD,
with zero Ohmic conductivity. Keep off-diagonal and subpixel observations separate:
their spatially interpolated mass operator requires additional derivation.

Use actual staggered update coefficients, not material-ID reconstruction:
c4=-4*pi*dt/cJ and epsilon=(c4+c5)/2=c4*(1+c1)/2. The implicit Lorentz-pole term
cancels. Recover the magnetic diagonal from the actual H curl coefficient. Weight
E energy by epsilon, H cross/centered energy by mu; reconstruct the prior H using
inverse mu. Poynting face work and midpoint LD/external-current work are unchanged.
Do not change the solver, field/LD checkpoint state or normal collective count.
Keep the legacy unit-background numerical evaluation path for exact regression.
A new read-only Fortran kernel takes coefficient pointers, with scalar finite and
positive-mass checks before publishing outputs. Unsupported observation returns a
warning/unavailable record; valid calculations continue.

Tests: failing manufactured diagonal-media kernel test, anisotropic grid spacing,
LD implicit coefficient cancellation, source/LD work, partition additivity,
static all-axis energy and zero/NaN/Inf coefficient rejection. Real MPI1/2 field
steps in dielectric epsilon/mu and scalar/diagonal LD: residual checks, no-observer
and previous-library state comparisons, direct checkpoint continuation. Reuse no
GS data or new SCF. Git-only build, audits, source/wheel checks, paired docs/JSON.
