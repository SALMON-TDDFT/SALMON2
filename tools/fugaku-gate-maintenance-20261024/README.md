# Fugaku job-gate maintenance, 2026-10-24

Archived, executed maintenance procedures for the existing Si512 Jacobi and WSe2 jobs. These are one-time operator scripts with explicit remote paths and fail-closed preconditions, not production startup scripts. Do not run them as a general setup step or automatically retry them.

The SALMON executable and physical inputs were unchanged. The source baseline is fda550b5 on codex/fix-cartesian-soc-symmetry; the Jacobi candidates remain separately versioned and are not promoted into this branch by this commit.

## Si512 Jacobi

- Accept sparse actual PLE stdout while retaining root output, every child exit and launcher completion evidence.
- Keep full-region versus fixed-27-region approximation differences as recorded physical metrics, not acceptance thresholds, per the human instruction.
- Ignore exactly one final memory sample only when both RSS and Hugetlb are null, earlier valid samples exist, and the child exited zero. Intermediate/partial/all-missing/invalid records and capacity violations remain failures.
- prepare_remote.py preserves the failed run and creates an independent rerun directory. Its first regression encountered an absent stderr fixture; finish_preparation.py records the corrected synthetic stderr regression and completes preparation. Preserve both as executed history, not a single clean-run recipe.
- fix_terminal_memory_remote.py adds the validated terminal-sample rule.

Validation: sparse PLE old RED/new GREEN and seven adverse cases rejected; memory 10 cases matched expectations; original 52614422 output passed the amended complete gate with 256 valid ranks and 82 terminal-null records excluded. Original failed validation remains preserved. User-approved reruns: 52745397 (OMP comparison), 52746333 (MPI parallel Jacobi). Sequential 52613213 was subsequently cancelled by explicit user instruction.

## WSe2

apply_wse2_policy_remote.py changes the shared gate of queued jobs 52698609, 52698613 and 52698614 without resubmission. Classified warning-only stderr is recorded; fatal and unclassified stderr remain failures. Finite magnetization is recorded without a magnitude cutoff. All-rank exit, completion, finite, SCF convergence, final electron count, memory, topology, restart and complete RT-output checks remain.

Validation: old RED/new GREEN; fatal/unknown stderr, child failure, missing memory, over-capacity, nonfinite magnetization, unconverged SCF and electron-count errors rejected. Full static manifests passed before and after atomic gate replacement, using a temporary old/new-hash verifier bridge. No executable or physical input was changed.

## Access and evidence

Execute remote actions only through the authenticated native Firefox Fugaku OnDemand UI, under applicable human authorization. No SSH CLI, internal communication or credential extraction. No GS data collection. Historical backups, locks and original submission records remain intact. Local calculation records are under OneDrive/calc2026; cloud synchronization is unverified.
