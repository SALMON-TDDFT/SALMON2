# LCFO WF radius namelist and initial coverage warning

User request: specify integration radius in namelist; warn if sphere retains less
than 99.9% of an individual WF norm. Do not choose or change radius automatically.

Add &functional hse_lcfo_wf_radius, always bohr: -1 default means legacy environment
fallback, 0 full, positive fixed radius. Explicit namelist wins over legacy env;
log ignored env to avoid ambiguity. Validate finite value >=0 or sentinel -1.
Keep experimental LCFO/MLWF activation switches unchanged.

Compute initial per-WF sphere norm using same minimum-image distance and inclusive
boundary as source mask, summing only disjoint cores across lcfo_comm. Reuse the
already available unmasked initial grid, norms and centers. O(Nwf) reduction/scratch;
no additional global WF copy. Report all geometric fractions, protected flag and
radius to lcfo_mlwf_radius.dat. Warn once at initialization if any fraction <0.999;
protected WFs remain uncut as before. No renormalization, changes to U, or RT mask.
The warning describes initial support, not guaranteed current/dielectric accuracy
or future-time per-WF support (existing runtime total norm-loss diagnostic remains).

Tests: deterministic 3D periodic geometry and boundary/full support; native MPI
namelist/env equivalence, precedence, full support and warning, invalid setting.
Build, existing direct-WF/ACE/orbital-MPI regression. Generate/validate/pack Si
4^3/6^3/8^3/10^3 GS/RT inputs and pseudo; no large jobs. Document, review, push.
