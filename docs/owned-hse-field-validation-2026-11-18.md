# Owned fast-PR HSE field admission — 2026-11-18

This slice adds a read-only uniform A/c little-group validator to the independent
fast-PR symmetry owners. It does not admit symmetric TDDFT/FDTD or checkpoints.
Fast-PR remains the only supported HSE engine; exchange algorithms are unchanged.

Native status contract: 0 compatible, 1 invalid/unready, 2 incompatible field.
Nonzero A is normalized by its largest absolute component; the Cartesian
rotation residual must be <=1e-10. Zero A still needs a valid initialized map.
No shared symmetry fallback, state mutation, handle copies or MPI in validation.
The TDDFT step hook precedes numerical backup/propagation and participates in the
existing collective guard. Existing symmetric open/checkpoint prohibitions remain.

Local Apple Silicon evidence:

- Fail-first compilation against the preceding modules rejects the two missing
  validator symbols. New audited external source/library passes the native probe.
- Independent C2z, mirror-x and rotated-cell maps; global symmetry poisoning;
  zero/tiny/huge A; direction tolerance boundaries; NaN/Inf; malformed/unready
  maps; stale geometry/operations; wrong/closed owners; original arrays unchanged.
- Existing16-case full/ACE × star/orbits × MPI1/2 processes × OMP1/2 exchange
  campaign passes. Each process uses COMM_SELF and different model/call schedules;
  this is not symmetry-reduced internal-k MPI2. Old full-k reference bytes exact.
- Existing two-rank GS outputs exact. Five-rank PBE/HSE Taylor/PT-CN and seven-rank
  mixed32-step protected numerical failure/retry match prior arrays/ordinary
  output bytes. These are non-reduced profiles.
- Seven-rank field/LD/SBE/PBE/HSE protected direct restore: save16/fail18; all
  original state bytes preserved, continuation and three retried future steps
  exact, zero replayed accepted steps.
- Source and installed-wheel Python suites, source-distribution-to-wheel build,
  native hashes and installed combined/1D fault-retry are retained with evidence.

Evidence directory: `~/.codex/outputs/salmonpy-owned-hse-field-20261118/`.
Tests are small implementation regressions, not material convergence proofs.
Memory budgets are advisory. No GitHub Actions were used.

Next: complete symmetry-map/history/cache persistence, then internal-k MPI and
symmetric mixed evolution/failure/resume. The validator is not a full native
initialization or archive-compatibility certificate.
