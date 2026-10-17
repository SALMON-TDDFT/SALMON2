# Giant Gamma native SR exchange snapshot

This snapshot records the source used for the native-core, wave-first, summed-return and conjugate-pair implementation, including orbital-distributed learned ACE and its dependencies. The conventional real-space propagation path replaces fixed-LCFO RT in this validated source tree; older unreferenced source files remain in Git history/the checkout but are excluded from its CMake source lists.

Pair reuse applies only to explicitly identical source/target support ACE at Gamma with the supported real even kernel. Unsupported action paths retain the conventional route. Source and target batches are bounded at 8 or 16; conjugate reverse contributions return to their orbital owner. The pair OFF path is retained as the reference.

## Runtime selection

```sh
export SALMON_FACTOR_HISTORY=learned
export SALMON_LCFO_RECONSTRUCTION=cached
export SALMON_HSE_SR_RADIUS_AU=10
export SALMON_HSE_SR_GROUP_FFT=1
export SALMON_HSE_SR_WAVE_FIRST=1
export SALMON_HSE_SR_RETURN_SUM=1
export SALMON_HSE_SR_WAVE_BATCH=16
export SALMON_HSE_SR_PAIR_SYMMETRY=1
```

Use `SALMON_HSE_SR_PAIR_SYMMETRY=0` for the pair reference and batch `8` for the smaller batch reference. The input must retain positive SR tolerance, full-source support, integer occupations and the existing learned ACE exact/history/teacher settings. The production runs use 4 MPI per node and 12 OpenMP threads.

## Existing evidence

Authenticated Fugaku native OOD comparisons 52364942, 52365802 and 52365807 completed with all tested ranks exiting zero and comparison/capacity checks passing. Kfast pair OFF/ON times were 4.33278734/3.62258488 seconds per step at batch8 and 4.14318069/3.34693763 at batch16. These are means of the maximum rank whole-step time over four early steps, not learned steady-state timings. Maximum current difference was 5.4977e-18 and energy difference 2.842171e-14 Ha against same-radius references; this does not establish untruncated SR0 accuracy. RSS+Hugetlb values are sampled rank maxima, not instantaneous/node peaks.

The binary hash and full source SHA manifest are in `source-provenance.json`. Original logs and reproducible analysis are retained in Calc2026 `sr-native-batch-pair-complete-20261017-0940`. No new tests, simulations or post-completion investigations were performed for this commit. Existing production jobs/binaries remain unchanged.
