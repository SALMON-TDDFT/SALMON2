# HSE–TDCDFT linear-response comparison

**Goal:** Run the HSE impulse to a useful observation time and compare with the existing stabilized TDCDFT equilibrium impulse at identical field strength and finite observation window.

**Architecture:** Extend the atomic-restart PT-CN driver with a selectable blocked *full-kernel* backend. Resume the original full-HSE trajectory from step130; do not introduce a cutoff quench. Use the verified DistanceExchange with radius=None. Finish at step2750 (880 au =21.28618fs), matching the post-probe duration of the completed TDCDFT g10_ground trace. A detached sequential job runs propagation, then a deterministic spectrum comparison. Expose live status and results; do not claim a spectrum before completion.

**Reference:** Si8, 12³ grid, shifted4³ k, impulse A_z/c=1e-4. HSE dt=.32; native TDCDFT dt=.08, alpha0=.2, beta0, gamma.001, ELF stride10. TDCDFT delayed probe begins80au after an unpumped preparation. Own self-consistent GS for each functional; identical cell and pseudopotential must be verified. No k convergence. Use full native sampling for each Fourier quadrature and the same cubic window, not interpolated spectral input. No scissors shift, peak fitting or additional damping. Show observation-window dependence before interpreting peaks.

## Tasks

1. Driver selection: test invalid method and blocked mode bypassing MLWF localization; implement backend choice without changing default, kernel, dt, restart physics or acceptance gates. Record execution method per row/status. Resume one accepted step to verify the actual path before a long run.
2. Comparison: freeze/copy the completed TDCDFT current and provenance; parse post-probe time and signed A step. Implement same-window analysis, output CSV/PNG and peak/window metrics. Reject missing/short/failed HSE, nonmatching duration, nonuniform traces or wrong impulse. Use a transient-only plot while HSE is short, with no peak claims.
3. Tests: synthetic oscillator at two time steps and shifted probe origins; amplitude normalization and input rejection; full sample suite. Review driver changes for checkpoint coherence.
4. Launch a detached sequential worker from the continued checkpoint to2750, then run comparison. Persist command/PID/code revision and a job status. Estimate remaining time from actual new steps. No recurring automation or promised notification.

## Decisions

- Full kernel is chosen because it is only slightly slower than Rc16 and permits reuse of the already converged GS and1fs response without a model change.
- Primary comparator is the latest stabilized ELF-Proca g10 equilibrium reference, not the older direct a2/a0 model. A shorter-window comparison will be reported for finite-time sensitivity.
- Launch duration is an initial finite-window comparison, not proof of timestep or spectral convergence. dt=.32 was accepted for the HSE pilot; longer-time accuracy remains to be assessed.

## Implementation ledger

- Driver selection: observed failing new-method tests, then failure/checkpoint tests passed. Real resumed step130→131 completed with same energy baseline and9.26e-12 residual.
- Input audit: TDCDFT external A first changes at80.08au, post-probe duration880au, signed amplitude1e-4. Identical pseudopotential SHA verified against HSE export source. Current and provenance copied to persistent ignored calculations storage.
- Analysis: synthetic two-dt oscillator, delayed origin, input gates and complete report generation pass. Real1.014fs transient preview rendered and visually inspected; no spectrum extracted from it.
- Job recovery: observed analysis-retry regression failure, implemented skip of already-completed propagation; added failing then passing test for changed-input fingerprint rejection.
- Review: no launch blocker; addressed analysis-only fingerprint guard. Method switch does not change physical checkpoint expectations.
- Long calculation and final spectra remain pending until the launched worker completes; do not mark the scientific comparison complete based on tests or launch alone.
