# Si no-corrector fresh RT implementation plan

**Goal:** Execute user's authorization to restart Si Pump, Nr16/Nr20 and PT-CN from step1 with balanced FFT, Taylor predictor-corrector disabled.
**Architecture:** Existing balanced-pair code plus three admission-check changes validated by Mac experiment. Immutable new executable/source/build roots, old jobs and evidence preserved until fresh preparations pass. Reuse validated GS only; no RT checkpoint restart. Nr16 and CN use328nodes1312MPI (164x8 workers,512 pairs each); Nr20 keeps640nodes2560MPI; middle bundle keeps5472nodes21888MPI per case,10cases,total54720nodes. OMP12. CN keeps its own solver and projection; Nr16 projection preserved, Nr20/Pump projection off.
**Tech Stack:** Fortran/CMake, Python3.6-compatible harness, Fujitsu MPI/PJM through native Firefox OOD only.

1. Preserve Mac evidence and reproduce rejection of PC=n in old executable. Change only three admission predicates and corresponding messages. Audit source, build/link, audit all application objects. Commit code.
2. Independent remote source copy of balanced validated tree, exact two-file transfer and byte hashes, source audit before Fujitsu build; all application-object audit after link. New binary immutable SHA.
3. Minimal native verification for new no-PC path plus CN balanced route: numerical/MPI/OMP/capacity/child+launcher/actual streams. No threshold relaxation. Confirm actual PJM outcome independently.
4. Prepare independent step1 Nr16/Nr20/CN roots and no-PC middle bundle. Preserve GS paths/manifests, physical pulses/dt/nt, projection policy. Start-aware validation for fresh zero history. Verify all static dependencies, input line length, syntax, binary SHA, quota, locks and duplicate jobs.
5. Only after readiness recheck four old IDs QUE/unstarted, cancel each once, record scheduler replies/history, submit each replacement once. Ambiguous receipt never retry. Keep old outputs/source/locks/checkpoints immutable.
6. Independent scheduler checks; update existing96-gs monitor and latest instructions. Collect selected non-GS preparation/results records using native OOD Download into canonical calc2026 and verify every size/SHA. Cloud sync unconfirmed.
