# Group5: standard TDDFT handle-owned timers

Reuse s_timer_context for each standard handle slot, with a temporary host bank.
Select the bank before initialization or handle work and park on every returned
path, including rejected initialization and close. Keep the running-timer latch
per slot. Initialization/release continue to use the same numerical routines.
Legacy singleton timing and HSE adapter timing are unchanged. No numerical
checkpoint fields, format/revision change or Python module is required.

Expose an additive timer protocol and rank-local completed counters. Paused
intervals and sibling work are excluded; restore never erases time already spent.
Use deterministic native timer-bank tests and alternating real standard handles
against prior independent runs. Do not collect GS data; reuse existing GS read-only,
keep diagnostics in calc2026 and run no CI. Old runtimes retain operations; only
new diagnostics require the capability. Threaded calls remain unsupported.

## Acceptance

PBE rgrid2 and TB-mBJ orbital2: two standard handles versus two old independent
controls,32steps, bit-exact currents.64 sibling counter comparisons per case;
memory rollback and close/reused-generation/failure checks pass. Existing GS
read through references, no GS copies/SCF/disk electronic checkpoints. The
protected step exposed a resolved-path/relative-path identity mismatch; accepting
the exact resolver path fixes admission without weakening unrelated file checks.
Source and installed638 tests:637pass1skip. Deterministic native timer fixture
passes. JP/EN/AI updated; no CI. Fast HSE timer ownership remains separate.
