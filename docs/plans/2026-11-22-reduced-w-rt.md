# Reduced finite-k W and real-time spectral validation plan

**Goal:** Test the approximate response compression by the resulting Si RT spectrum, rather than requiring every response orbital to reach1e-6 residual.

**Architecture:** Opt-in equilibrium scalar `hsex_w_solver='reduced'` uses all k/q with the existing screened-W assembly/output. Gather only occupied source projectors; for each owned target, form validated little-group q stars, batch representative q RHS (default4), solve shared projected energy-shift systems at bounded rank, reconstruct all local polarization blocks and retain the exact existing optical head/wings path. Rank-cap status2 is logged and accepted only in this explicitly approximate mode; all other failures reject collectively. Default solvers and remote jobs remain unchanged.

**Tech Stack:** Native Fortran, existing real-space Hpsi, MPI occupied-state collection/k distribution, BLAS/LAPACK, existing optical/box W routines. Python only for offline comparison and plots.

1. Add failing whole-W native pipeline calls using new optional reduced controls. Full response-space dimension on a small grid must reproduce scalar dense/matrix-free W, with symmetry off/on and4MPI12OMP.
2. Extend screen_compute with optional reduced rank/batch/tolerance/symmetry and ionic metadata. No dense real-grid H; memory bounded by q batches plus occupied source manifold and per-q small chi. Count capped batches, H actions, symmetry probes and residuals; preserve all numeric finite and MPI rejection rules.
3. Add reduced solver input validity/log/default/broadcast/manual controls; hsex_reduced_target remains0 for full W. Publish the normal W only for successful full assembly.
4. Build source/object audits and small whole-W regressions. Produce same-GS Si Nr12 Nk=Nq4 cutoff3Ha box6 reference and reduced W in fresh scratch; compare optical tensor, W norm and matrix differences. Never collect GS payloads/logs.
5. Use occupied-only TB-mBJ+HSEX RT, same GS, propagator/kick/dt/time window, fixed W. Reference and approximate kick/no-field are paired. Reuse a completed reference only after checking metadata and a short replay with the current executable; otherwise rerun it. Default target16384 steps atdt.05au, cubic SALMON window, zero-field subtraction.
6. Compare1–6eV Re/Im dielectric curves, current difference, first stable peak position/height, and75%/full time-window consistency. Im vertical axis positive. Report actual errors and resolution; do not call the approximation spectrally validated before RT completes.
7. Save figures, derived metrics, filtered diagnostics, inputs and code/build records under canonical calc2026; no GS payloads. No automatic commit, remote production replacement or cancellation.
