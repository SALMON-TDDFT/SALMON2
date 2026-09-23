# Bare-K silicon pump-probe plan

**Goal:** Determine changes in the exciton-like Si peak for alpha0=1, reference=0, using 4³ k points.
**Architecture:** Existing pump-only subtraction and delayed z impulse. Keep the existing pump (omega=.2 au, full width60 au, 1e13 W/cm²) and dt=.08. First check whether the probe reactivates the external-field screening estimator; distinguish frozen pump screening from full probe-dependent screening. Use matched post-probe observation windows and signed kick amplitudes.
**Tech Stack:** Existing SALMON MPI build, NumPy and Matplotlib.

1. Short full-feedback diagnostic: probe threshold80 au, amplitudes .001 and .0005, nt1200; subtract existing identical pump-only trajectory. Record alpha jumps and normalized differential currents, pre-probe equality.
2. Resolve probe screening convention using user preference. Frozen alpha response is a partial response at the pump-prepared screening; full feedback must pass amplitude linearity before optical interpretation.
3. Production pump-probe at post-pump delays, matched no-pump reference, plus amplitude check. Use at least ~20 fs observation for comparison with earlier spectrum. If the new no-pump model has no stable linear response, explicitly distinguish the old fixed-alpha exciton reference from a same-model equilibrium reference.
4. Extract Im epsilon from differential number current, polynomial window and actual probe threshold. Compare peak location/height/area in 2–4 eV with equal windows; assess window sensitivity from common truncations. Do not infer binding energies or resolved widths below window resolution.
5. Save reproducible inputs, spectra, metrics, plots and limitations. No k-point scan or push.

Implementation decision: use optional tdcdft_screen_stop=60 au for frozen-response diagnostics,
default -1 unchanged. Checkpoint v3 adds cutoff record; v1/v2 default to disabled.
Added integration tests: delayed impulse does not unfreeze K/alpha; P continues;
pre-cutoff trajectory identical; restart and changed-cutoff rejection; v2 compatibility.
Full-feedback short test failed amplitude linearity (~50.6% normalized-current difference),
with alpha returning to 1 immediately after the probe. Production comparison is explicitly
frozen-screening at t_probe=80 au, total time960 au. Four trajectories use the same4³ grid:
pumped .001, pumped .0005, ground alpha1, ground alpha equal to pumped final value.
