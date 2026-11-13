# ASE Calculator energy and forces

Continue the user-approved Python workflow: expose an optional ASE Calculator over
the independent DFT runner. Implemented properties only energy/forces; no stress,
finite-temperature free energy, cell optimization, MPI launching or native warm start.

1. Report existing native forces at full precision after calc_force without changing
   arithmetic/order; read and hash in GroundStateResult, retaining old result support.
2. Implement separate salmonpy.ase_calculator (explicit optional ASE import), immutable
   configuration, fresh per-geometry run directories, ASE cache invalidation and failure
   handling. Never reuse a stale GS for changed geometry. Unsupported set()/properties
   fail explicitly; user constructs a new calculator for changed model settings.
3. Compare nonzero displaced-Si forces with central energy differences at two steps;
   assess grid/discretization discrepancy before exposing optimization. If native force
   is unsuitable, preserve it as diagnostic and use verified finite differences or
   narrow supported scope explicitly rather than presenting unverified forces.
4. Test cache/mutation/failed calculation/units; real minimal fixed-cell optimizer run;
   update AI JSON, human guide, source and installed-wheel verification.
