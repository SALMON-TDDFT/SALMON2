# Hybrid unused-code cleanup implementation plan

**Goal:** Remove the approved unreachable test scaffolding and unused declarations without changing active exchange, prediction, diagnostics, GPU, or reference paths.

**Architecture:** Delete modules with no source/sample/test callers. Keep live fixed-Q trace diagnostics, but remove unused accepted-teacher APIs and storage. Preserve the diagnostic V1 output format by emitting the same false acceptance flag. Simplify the Taylor stage interface to its only used argument. No namelist or checkpoint changes.

**Tech stack:** Fortran, CMake, GNU MPI/OpenMP HSE and serial builds.

1. Confirm baseline HEAD 4647f42b8 and clean status; inspect callers and compiler unused warnings.
2. Remove exx_surrogate_controller/shadow/stage/guard/export and their CMake entries.
3. Remove trace_endpoint/trace_teacher and their unreachable accepted-entry path. Preserve trace_candidate, fitting, prediction and output. Remove unused declarations/imports and Taylor stage arguments with caller updates.
4. Build MPI HSE, serial HSE and no-HSE configurations. Re-run 21 existing regressions and bounded complex exchange action oracles. Exercise live trace capture/output separately and check the V1 false acceptance marker.
5. Review the diff, repeat static Coding Rules/unused checks, record verification. Leave uncommitted until requested.

## Verification completed

Removed 5 unreachable modules and unused teacher APIs/storage. The diagnostic V1 acceptance marker remains false without allocating an acceptance array. Current call sites use exx_taylor_stage(stage). Three build configurations and 21 regression checks passed. Bounds-checked 1/3/4/9-rank exchange comparisons passed. Pre/post-cleanup 16- and 80-step H2 runs produced byte-identical current, energy, response, trace and raw factor outputs; the 80-step shadow fit/comparison lines also matched. Unused-variable/argument warnings are zero for the audited exchange modules. Read-only review found no issues. Runtime switches and restart formats are unchanged. Changes remain uncommitted.
