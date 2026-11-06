# Automatic FFT consensus and dense-k integration plan

Goal: resolve rank-local FFT auto-selection mismatch, then integrate and commit the most effective supported dense-k development implementation without changing remote jobs or unrelated worktrees.

1. Reproduce missing synchronization interface and a one-rank pre-plan error. Add actual MPI kernel/oracle tests before implementation.
2. Add an optional pre-plan consensus callback. Aggregate trial failure and per-layout worst-rank timings, retaining the 1e-5 noise floor and 5% selection margin. Preserve explicit layouts and mismatch rejection. Route local pre-plan errors through the same agreement.
3. Compare committed worker baseline, validated 2026-11-01 workspace/FFT reuse, and 2026-11-04 batched contraction sources with identical dense-k synthetic inputs, independent direct operator oracle and maximum-rank whole-call timings. Preserve source hashes and previous Fugaku RT outcomes.
4. Integrate selected changes only; retain unrelated dirty changes and default/input/science protections. Document scope-dependent performance, extra memory, and experimental opt-ins.
5. Full source/object prohibition audits, MPI/serial builds, numerical/rejection/cache/lifecycle regressions, repeated automatic-layout application starts, focused RT validation, independent review, final commit. No push or remote executable replacement.

## Final selection after human review

The human requested Fugaku evidence and one supported implementation, with uneven block work accepted. Original 53272164 rejects reuse/packed as a performance upgrade; monitoring53550165 rejects batched contraction. Retain the established dense-k implementation already merged with upstream; archive both experimental candidates and do not add their switches. Remove only the invalid global strip spread<=1 requirement; retain quotient/remainder balancing, exact coverage and capacity checks. Commit FFT consensus, uneven acceptance, tests and provenance.
