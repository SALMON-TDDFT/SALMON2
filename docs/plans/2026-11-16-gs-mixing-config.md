# GS mixing configuration isolation

Capture method, initial coefficients (mixrate/alpha_mb/beta_p), Pulay memory,
auto-mixing switch/ratio and initial no-mixing steps in a value-owned nested
s_gs_run_config%mixing. Copy into s_mixing before DC or ordinary initialization.
Keep adaptive coefficients and density/potential history as mutable per-state data.
Keep the existing 21-slot storage and all formulas, ordering and reductions.
Legacy init_mixing callers capture globals at initialization; retained GS never
reads these globals after capture. Reject invalid explicit settings before activation.

Verify poisoned legacy inputs, independent retained traces, output/restart bytes
and residual traces against the original solver for simple/Broyden/Pulay and
potential mixing. Check auto adjustment/reset separately. No molecular work,
remote CI, or physical convergence claim. Update runtime guides and AI catalog.
