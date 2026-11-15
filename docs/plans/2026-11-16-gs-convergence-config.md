# GS-owned convergence metric

Continue approved GS input isolation with one bounded extraction: add the
existing five convergence metric names to s_gs_run_config. Capture namelist
values once, validate the explicit selection before activation, and route SCF
metric calculations and labels through the owned copy. Keep numerical sums,
normalization, thresholds, update order and legacy fallback unchanged.

Pre-SCF semantic checks must inspect the selected metric, not the shared input.
Reject unknown metrics before native allocation. Poison legacy convergence after
capture in the lifecycle probe; caller mutations must not alter retained state.
Compare all five supported metrics with original separate-process output. This
does not isolate electron normalization, mixing, material inputs or MPI resources.
No new XC coverage, molecular work or Actions execution.
