# Expanded mixed validation on the 64GB / 18-core Mac

Preserve the4GiB/4-core beginner recipe. Host capacity supplied by the user:
64GB RAM,18 cores. Local campaign uses at most15 MPI ranks, OMP/BLAS1 and48GiB
advisory budget (no hard admission stop). Never parallelize campaign jobs.

Add only a sample CLI field-shape option: even counts >=40 retain the existing
synthetic geometry and source. Default40³ and centered point coordinates are
unchanged. Input signatures already bind actual grid and positions. Add explicit
shape to result metadata and time only the coupled loop (excluding setup, output
and expensive serial reference replay). Existing compare tolerances apply.

First test formerly unmeasured field2+SBE1+TDDFT2 with5 ranks and protected restore.
Strong-scaling fixtures:96³ field1/2/4/8/12 with one SBE rank and two orbital TDDFT
ranks (4/5/7/11/15 total). Compare full A/E/j histories at the same shape/input;
report end-to-end and coupled-loop timing separately, MPI communication/wait
metrics and aggregate per-rank RSS peaks. Repeat timed runs if useful; single-run
numbers are smoke evidence, not a publication-quality scalability result.
No grid-convergence or production physical-accuracy claim. TDDFT grid and k mesh
remain the small existing verified regression input. Larger electronic systems
and cross-platform tests remain separate work.
