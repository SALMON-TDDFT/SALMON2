# Stable small-gradient HSE evaluation

Goal: remove the reproduced phase0 full/reduced-k disagreement without changing
k weights, the propagator, fields, tolerances, or the fast-PR exchange action.

Root cause: Libxc 5.2.3 HSE06 mixes WPBEh(omega=0), -1/4 WPBEh(omega),
and PBE correlation. In WPBEh(0), separate E1/log derivatives cancel terms
proportional to 1/sigma, losing the finite small-gradient limit.

Implementation: evaluate the same unscreened WPBEh formula in u=s² for s<=1e-4,
combining E1/log analytically before differentiating. Reassemble only selected
points from that exchange, Libxc screened exchange, and untouched PBE
correlation. Use the analytic screened quadratic limit only below roundoff in u.
All other points retain the original Libxc call and output bytes. No gradient
floor, field adjustment, or change to supported engines. Spin scaling and core
contributions follow the existing wrapper. Keep the complete source override,
manifest hashes, and attribution; never modify pinned sources or old binaries.

Verification sequence:
1. Native low-gradient regression must fail on the original library (STOP 2).
2. Check analytic derivatives against energy differences and ordinary gradients,
   density/screening/spin sweeps, and both HSE wrapper entry points.
3. Fresh explicit build; original phase0 full/reduced-k zero/weak-field and
   FDTD mixed regression; normal legacy calculations compared to old binaries.
4. Source tests, source-distribution/wheel checks, documentation and AI JSON.

Only local execution; no GitHub Actions because of the usage budget.
