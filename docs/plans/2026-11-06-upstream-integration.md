# SALMON upstream integration implementation plan

Goal: merge develop-2.0.0 376635eb into the committed local EXX baseline b3388e45 without changing remote jobs or unrelated dirty worktrees.

1. Fetch upstream and isolate the integration in codex/integrate-upstream-20261106.
2. Review conflicts and auto-merged removals. Retain local EXX names, optional modes, restart parity/metadata, custom KY SOC/core23, development evidence and guards. Incorporate upstream pseudopotential/PDOS/SCF/GPU fixes.
3. Replace upstream array IEEE inquiry with an explicit scalar loop and enforce scalar helper interface.
4. Run source prohibition audit before configure/build; full CPU HSE/MPI compilation and application object audit, focused PP/SOC/restart/EXX/finite regressions.
5. Inspect final differences and record validation limitations. Commit only the integration, with both upstream and local parents. No push or remote executable replacement.

Review: independent read-only review completed; CMake option/default ordering corrected. Inherited auto-FFT MPI layout mismatch reproduced on pre-merge baseline and documented separately. No relaxation of physics gates.
