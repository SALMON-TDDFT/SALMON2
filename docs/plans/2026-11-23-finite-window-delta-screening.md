# Finite-window differential screening implementation plan

Goal: keep the existing equilibrium W0, and evaluate only pump-induced delta-chi with a finite causal delay window. User-approved design; Drude transport relaxation tau=0.5fs is the exponential reference. Apply window to auxiliary response only; propagate the original occupied RT orbitals coherently.

Definitions (atomic units): Gaussian w(s)=exp[-s^2/(2 sigma^2)]. Exponential w(s)=exp[-s/A], A=sigma sqrt(pi/2), matching integral area A. User time parameter is the common area A (physical exponential decay time); Gaussian sigma=A sqrt(2/pi). Thus baseline time0.5fs means eta=1/A≈0.04838Ha and sigma≈0.39894fs. Compare A=.25,.5,1fs, Gaussian/exponential.

Response: F(delta)=-i integral_0^infinity w(s) exp(-i delta s) ds. Exponential F=1/(i/A-delta). Gaussian F=-sqrt(2)*sigma*Dawson(delta*sigma/sqrt(2))-i*A*exp(-delta^2*sigma^2/2). Finite delay, not eta->0. Project paired response Hermitian explicitly as the existing experimental closure, with discarded fractions logged; keep reality/finite/positivity/reciprocity guards.

Equilibrium cache built at RT initialization from H0 and original occupied states: three static eta samples plus one finite-window response, sharing H-action bases. Store root-only per-q chi_eq_static and chi_eq_window; cache layout, G, parameters and baseline W. Reference static extrapolation/positivity must pass. Compare computed static W_eq against supplied W0 at <=1e-6 relative, otherwise reject (boxed W cannot uniquely reconstruct arbitrary generating chi). No silent independent LDA/TB reference substitution. Window mode currently rejects restarts and independent LDA W; old modes unaffected.

At update: chi_eff=chi_eq_static+Herm[chi_window(t)-chi_window(eq)]. Difference evaluated before addition to preserve exact zero. Invert chi_eff with existing RPA convention, then Gamma/R transform. Require reference valid/layout identical. Return W0+[W_eff-W_eq] only after W_eq/W0 validation, preserving supplied baseline and no-pump roundoff. Root chi cache is not replicated; final W remains pair-distributed. Existing k/q MPI workers and original RT orbital count retained.

Task1 (independent low-level): hsex_frozen_window module, mixed spectral filters in frozen_solve, certified finite-window error via Duhamel bound; tests with dense spectral/time integration, Gaussian/exponential area matching, zero pole, nonfinite, capacity rejection. Interface frozen_solve optional window_kind(:) int (0=resolvent,1=Gaussian,2=exponential), window_time(:) real area A, each size(z). Old calls unchanged. Only bounded solver path supports mixed/window channels, reject incompatible mode. For window response error normalization is A*||source||, not an original resolvent residual. Bound A*source_out+B*||rhs||*(H-left-defect*rnorm+H-right-defect*lnorm), B=2*A^2/pi for Gaussian and A^2 for exponential, plus roundoff. Never relabel resolvent residuals as Gaussian certification.

Task2 (root integration): optional s_frozen_window_cache in full-W driver, initialization four channels (3static+1window), update one window channel. hsex_frozen_window='none'/'gaussian'/'exponential', hsex_frozen_window_time in current time units, default0.5fs; only frozen+hermitian+pairs, no stationary optical/restart/independent reference. Global namelist defaults/broadcast/conversion/logging/validation/manual. hsex_update_frozen step0 builds reference with no W publication/mutation; initialization_rt calls after V_local creation.

Task3: audited native MPI regression, empty owners/caches/mismatch/failure paths, full application source/object audits. No changed remote production/fixed binaries. Native and Mac small Si8MPI×1OMP/BLAS, max2SALMON jobs, initial W0 from previous valid static extrapolated reference. Verify no-pump W0, Gaussian vs exponential area-time0.5fs, and width.25/.5/1fs dependence at pumped1fs. Only extend RT if physical/numerical guards pass. No GS/restart payload collection. Canonical calc2026 records, source/manual/tests. No commit requested.

Implementation status (2026-11-23): finite-window module and mixed-channel solver,
root-only differential-response cache, W0 anchoring, inputs/manual/sample and
native MPI tests implemented. Gaussian/exponential filters passed independent
spectral/quadrature and Duhamel-bound checks. Cache lifecycle/collective rejection
passed with 1 and4 MPI; quality/spec reviews passed. Full HSE-on/off builds and
source/all-object prohibition audits passed. Mac small-Si results are recorded
under canonical calc2026/si-window-delta-20261123; no GS data collected. The
0.5 fs time is the user's empirical Drude baseline, not an orbital collision term.
