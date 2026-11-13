# Material-ID-independent staggered FDTD assignment

User request: fix ID-dependent interfaces, rerun the two-layer calculation,
preserve single-material numerical results, allow independent material choices.

Use a geometrical staircase assignment, not susceptibility/coefficient averaging.
For an E edge joining two non-vacuum cells, use the lower-coordinate cell's
material; when either endpoint is vacuum keep the legacy vacuum edge. For an H
site shared by up to four cells, choose the first non-vacuum cell in a fixed
coordinate ordering. This retains one-material/vacuum behavior and the existing
ADE poles and time-loop order. It is not subpixel smoothing or exact interface
homogenization. Meep's subpixel documentation explains why naive dispersive
averaging is not equivalent to a physical normal/tangential effective tensor.

The impulse initializer also currently sums every LD material into every active
pole. Restrict it to the edge's geometrical owner; preserve single-material loops.

1. Add a native regression that fails for ID permutations before the patch, checks
   LD/LD with unequal poles, dielectric/LD and different mu, and runs along x/y/z.
2. Edit an external prepared source only, regenerate the pinned compatibility
   patch, update its SHA and the runtime interface revision so old bundles cannot
   silently satisfy the new contract.
3. Build/audit a fresh runtime. Compare singleton cases including impulse byte for
   byte with the previous executable; test ID permutations, duplicate equivalent
   material splitting and shape ordering. Rerun the original two-layer fixture.
4. Add a usable two-layer example, document assignment and supported scope in
   English/Japanese and AI JSON; include native regression in runtime bundle CI.
5. Run Python, standalone sdist/wheel, TDDFT and native checks; review changes.
