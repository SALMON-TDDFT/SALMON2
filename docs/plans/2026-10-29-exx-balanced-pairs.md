# Nk16 Pump balanced pair implementation plan

> Execute in the existing isolated exx-worker-autolayout worktree using executing-plans.

**Goal:** Replace only queued Nk16 Pump 52911342 with a validated executable assigning all 27648 MPI ranks 24 or 25 distinct r/r-prime pairs. Nr16, Nr20, CN and WSe2 jobs are protected.

**Architecture:** Opt-in SALMON_EXX_WORKER_BALANCED_PAIRS=1 selects jagged row strips: distribute ranks over actual symmetry-reduced rows, then split each full r-prime row into contiguous quotient/remainder strips. Merge identical column strips into communication groups. Reject layouts with empty workers, pair-count spread above one or density scratch above 64MiB. Keep the old path as the default and preserve all physics and k ownership.

**Tech stack:** Fortran MPI/OpenMP/FFTW, Python independent layout oracle, native Firefox OOD.

1. Add a failing independent coverage/load/memory/invalid-geometry test and scalar/SOC action comparison.
2. Implement the bounded strip planner and opt-in worker communication, retaining the default path and scalar IEEE checks.
3. Audit full source before builds; build GNU probes and application; audit all linked application objects. Run action equality, finite rejection, OMP and native GS/RT tests.
4. Commit exact source. Transfer through authenticated native OOD, build Fujitsu independently with all source/application-object audits. Minimal new-feature test only, require all-rank evidence and actual PJM normal final.
5. Prepare new independent Nk16 production root with original physics/resources, new validated binary, static dependencies, full-rank balanced assignment gates, and existing conditional pump-kick authorization migrated to the new root.
6. Immediately before replacement confirm old job state; if still QUE/unstarted cancel once and verify CCL. Submit replacement once with durable lock/record. If started, preserve outputs and establish checkpoint continuation feasibility before interrupting. Never touch CN/Nr16/Nr20.
7. Collect non-GS preparation/results through native OOD, verify size/SHA and store directly in OneDrive calc2026. Extend only monitor96 and report new ID. No GS payload collection, no cloud-sync claim.

The human approved uniform all-rank allocation and explicitly narrowed replacement to Nk16 Pump only. Initial rectangular autolayout oracle gives 24576 active ranks, 3072 zero-work ranks; this is the regression to fix.
