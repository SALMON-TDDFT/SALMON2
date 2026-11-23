# Native TDDFT admission/retry — 2026-11-17

## Scope

First subset of native recoverable errors: standard embedded TDDFT validates on the attached material communicator before setting its lifetime latch or entering an RT update. A queryable native protocol distinguishes statuses10 invalid inputs,11 rank disagreement,12 lifecycle,13 capacity and14 pending restore candidate. Numerical update body/order is unchanged. Python restores initialization availability only for declared unanimous pre-setup10/11; known unanimous step rejections raise InputError. Old/unknown binaries remain RuntimeFailure without an unchanged-state claim.

This does not solve deep solver STOP, post-mutation numerical errors, allocation/process/MPI failure, HSE/FDTD equivalents or whole coupled-step rollback. All communicator members must call the same entry. Actual initialization still consumes the one-process lifetime. Failed work directories are retained for diagnosis.

## Evidence

- Old fixed library lacks the new protocol; raw MPI regression fails before implementation. A separate real Python retry regression failed after native admission but before the host latch fix. Neither failing run is counted as a pass.
- New external prepared sources pass manifest/shared-kernel verification. Full explicit Libxc/MPI build with two workers succeeded; final adapter was compiled/relinked in another new directory using unchanged core objects. Existing binaries were not edited. Scalar IEEE checks retained.
- Four ranks split into two independent2-rank groups: one rejects13 raw ABI cases while the other makes valid progress. Cases cover before-init/after-close, nonfinite initialization, invalid capacity/layout, missing file, initial rank disagreement, repeated initialization, pending candidate, nonfinite/stale/rank-inconsistent step input and exhausted capacity. Every live-state rejected call preserves serialized numerical bytes exactly. Discarding the candidate permits progress.
- Two-rank Python initialization rejects unequal dt, releases the unattained lifetime and successfully retries. Eight currents after retry equal uninterrupted old-library currents exactly on both ranks. Installed-wheel raw/Python routes pass the same comparisons.
- Four-rank ordinary mixed FDTD + two LD materials + SBE + Libxc-PBE TDDFT(orbital2):64-step A/E/j and previous-current histories equal the prior fixed library exactly. Existing protected restore late-failure original-byte preservation, continuation, retry and following3 frames pass; zero replayed steps. Budget4GiB is advisory, with four cores and OMP/BLAS1.
- Source and installed-wheel suites: **519 tests passed, one skipped** each. sdist→wheel build/install succeeds. Source tests explicitly select the built LD/SBE library; no import-time compilation.
- Tests caught a stale rt-step-config base hash; updated it to the new verified standard source inventory. Japanese/English and AI guide copies synchronized; whitespace check passed. Read-only review found no blocker. No GitHub Actions used.

## Reproduction and artifacts

`tools/test_native_tddft_admission_mpi.py BUILD GS NEW_OUTPUT [--baseline|--split|--python-init-retry]`, two ranks normally, four with --split. Use an MPI timeout and OMP/BLAS1. The pinned fixture is the local small converged Libxc-PBE Si input; this tests numerical regression, not production optical convergence.

Persistent evidence: `/Users/otobetoshihito/.codex/outputs/salmonpy-native-admission-20261117/`. Includes prepared source, relocated byte-identical native library and receipt, raw/Python records, mixed protected-restore records, suite/build logs and distribution hashes. Paths are local provenance, not installation requirements.
