# GS electron-count normalization implementation plan

Goal: own density-convergence normalization inputs per retained GS without changing
occupations, numerical formulas, reduction order or DC precedence.

1. Extend lifecycle probe first: captured normalization total/spin counts, selection
   tests for unpolarized/polarized/fallback and nonpositive divisor rejection.
2. Add a nested value config and a pure count resolver in structures. Capture legacy
   counts in gs_input_config; copy as part of owned GS config. SCF uses owned inputs
   when present, legacy capture otherwise. DC still overrides with dc%elec_num_tot.
3. Reject nonpositive density-metric divisors before activation for non-DC contexts.
   Do not reinterpret these as occupation/material electron-number setters.
4. Include configuration in retained trace. Compare original residuals and every
   restart payload for all five metrics; separately test a doubled normalization
   divisor produces half the residual with unchanged electronic payloads.
5. Local explicit build, Python/wheel tests, reviewer, Japanese/English/AI guides.
   Keep molecular work paused and no remote CI. Record spin/DC tests vs code review.
