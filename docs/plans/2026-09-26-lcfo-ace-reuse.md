# LCFO ACE reuse implementation and validation

User-authorized extension of the native LCFO RT plan: reuse ACE between selected
physical time steps and measure speed/accuracy; append results to the current note.

Design: opt-in SALMON_LCFO_RT_ACE_INTERVAL, positive integer, default1 unchanged.
Initial exchange is always built. For impulse, force refresh at step1, then every interval from step1. For smooth
laser fields, reuse the initially constructed GS operator, refreshing at interval
multiples. At scheduled steps, refresh
both predicted and accepted endpoints using existing midpoint averaging. Other
steps keep the last accepted ACE; native local fields still update. Invalid ACE
or changed occupations force an ordinary refresh. MLWF frames still undergo
polar transport at skipped exchange refreshes, avoiding phase drift. Evaluate
0.5 trace of the retained ACE at the current orbitals for energy bookkeeping;
this is explicitly a frozen-operator diagnostic, not current-density HSE energy.
No extrapolation or adaptive trigger in this first controlled comparison.

1. Add MPI2 regression for default/interval1 identity, fewer builds at2/4,
   finite current/energy, transported MLWF reuse, invalid interval rejection.
   Observe failure against previous binary before changing production code.
2. Implement configuration, physical-step schedule, reused-ACE energy and logs.
3. Build with existing no-loop-vectorization flags; run native and transport tests.
4. Run sequential MPI16 Si128 full-support MLWF comparisons, dt0.16, intervals
   1/2/4/8 with equal steps, thread counts and continuity diagnostics disabled.
   Measure wall time, exchange rebuild count, current error and diagnostic drift.
5. Record limitations and measured results in docs/inputs/lcfo-rt-development.md
   and append the current user-facing outputs/si128-mlwf-support/report.md.
6. Resume the existing full/9/8/7 long comparison using the original per-step
   update reference; do not silently promote a short-pilot reuse setting.

User refinement: impulse step1 must rebuild before its first predictor, not merely
at the predicted endpoint. Added current-orbital forwarding to stage0 and forced
rebuild before capturing the initial ACE/MLWF frame. Laser input with a smooth
Acos2 envelope confirms step1 retains initial ACE. Native MPI2 regressions pass;
review found no scheduling, optional-argument or frame-lifecycle blockers.
