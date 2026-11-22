# Periodic laser-driven W update timing plan

**Goal:** Measure actual wall time for W updates every1fs, separating projection, W computation/IO and normal real-space propagation.

**Architecture:** Preserve the legacy one-shot step with interval0. Add an integer step interval and one pure due predicate shared by RT projection scheduling and W update guards, so every event uses the same endpoint H/A and nonthermal occupations. Add rank-max elapsed phase timers. Allow the existing dense-finite-q/real-space-optical hybrid path for scalar internal pump updates; do not route partial occupations into the GS-only reduced solver. Restart scheduling uses absolute RT step indices and the existing active-W checkpoint.

**Tech Stack:** Native Fortran for scheduling/projection/W/MPI/output. Offline Python only for timing extraction. New binaries source/object audited, existing fixed binaries untouched.

1. Failing pure schedule tests: legacy one-shot, periodic2/5/8, disabled/invalid guards, absolute restart index and large integer limits.
2. Add hsex_screen_interval default0, namelist/default/broadcast/log/validation, common schedule module, shared projection/update guards. Periodic mode requires internal update and no external probe W.
3. Add projection/W/refresh timers and internal W derivative/build/save phase timers using local system_clock differences with MPI max reductions. Preserve numerical behavior, populations and occupied RT counts.
4. Enable hybrid optical Hpsi in the existing nonthermal RT update path and validate against existing single-update results on the small fixture. Add repeated-update/restart smoke test where feasible.
5. Timing campaign: existing Si Nr12³ Nk=Nq4³ GS, box6/cutoff3Ha,4occupied propagated states,8initial projection states, Mac8MPI1OMP. Weak4fs laser at3.4eV, Epeak1e-6au, dt≈1fs/800, total2fs and updates at800/1600. Weak excitation keeps q4 in the supported insulating-head regime; do not represent this as a strongly excited/metallic benchmark. Retain nonthermal occupations without rescaling.
6. Compare same-pulse fixed W and periodic W propagation. Require2 W files, correct event times, finite currents, capture and unchanged RT orbital count. Report per-event projection/build/save/total and cost per physicalfs; derived extrapolations labelled estimates. Save only inputs, derived records and filtered diagnostics under canonical calc2026, no GS payloads/logs.
