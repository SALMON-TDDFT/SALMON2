# Fast-PR HSE deep numerical recovery — 2026-11-17

## Contract and limits

Native fast-PR HSE now declares a separate recovery protocol1. Supported scalar fixed-ion periodic HSE06 Taylor4 and learned PT-CN retain accepted electronic arrays, XC payload and the explicit HSE numeric payload before each step. The latter includes ACE endpoints, factor/PT-CN histories, caches and counters, but not FFT plans or MPI owners. Handled operator, convergence, finite-value and history failures are agreed collectively before return. Rollback restores the payload and PP phase data directly, without propagation or exchange-history replay. Accepted A, step, current and energy remain unchanged; ordinary failed-step RT records are not appended. Diagnostic logs may remain.

This is an individual HSE group contract, not rollback of previously advanced field/other materials in a mixed step. The audited internal partition is k-only MPI1/2. Deep initialization, allocation errors, process/MPI failure, uncontrolled third-party STOP and output-device errors remain excluded. Unknown/old errors remain conservative RuntimeFailure. The external checkpoint format is unchanged. Extra backup/checking costs are real; memory budgets remain warning-only. No conventional HSE support was added.

## Verification

Local Apple Silicon, user capacity64GB/18 cores; OMP/BLAS1, explicit build4 workers, at most5 MPI ranks for these runs, advisory48GiB. No remote Actions or fixed-binary modification. Fresh full source preparation/build followed by new-directory incremental relinks; source manifests and shared-kernel checks pass. Production exports no injection hook. A separate module compiled with SALMONPY_RECOVERY_TESTING supplies single-rank fault injection.

The fixture is synthetic8³, two states/four electrons and four k points; it verifies implementation rather than converged physical accuracy. Taylor uses predictor-corrector and exact interval2; PT-CN uses its learned profile. Both histories become active before the warm-step tests.

- Two-rank Taylor/PT-CN failures at first step1, warm exact/prediction steps12/13 and final step32. One group rank receives wave/density/Hartree/XC/endpoint stages20..24, ACE30, history33 and PT-CN32 as applicable. All failure-before/after rank-local numerical checkpoints are byte-identical. Python A/steps/current/energy and ordinary output files remain unchanged; valid retry succeeds.
- Production finite input1e200 exercises deep overflow checks. Production PT-CN input1e100 triggers an actual nonlinear rejection with infinite residual and returns safely through the solver cleanup path. Four ranks form two independent two-rank groups: one rejects/retries while the other continues normally.
- Twelve matched trajectories (eight fault cases plus four independent-group trajectories) match the old fixed fast-PR binary exactly for32 accepted currents and energies, and accepted current/energy numeric output records. No cross-decomposition bit-equality claim.
- Normal128-step FDTD2 + two field-local LD materials + SBE1 + HSE2: A/c, E and j are exactly equal to the old binary (maximum difference0). Existing five-rank protected restore still preserves every original native payload and Python state after late staging failure, continues from the live state and produces three exact retry frames; zero restore replay steps.
- Source and installed wheel suites:526 tests passed, one skipped each. Source-distribution-to-wheel build, native hash check, packaging layout test and Markdown/AI JSON synchronization pass. Installed dependencies were reused offline from the existing local test environment; no import-time build/download was introduced.

## Review and negative evidence

Fail-first Python test observed RuntimeFailure before protocol handling was added. Old-library protocol absence is rejected. Review caught an exhausted-capacity admission return becoming success after a collective reduction. The raw ABI test reproduced that defect; capacity/checkpoint-failed/nonfinite checks now participate in the collective guard, and the final adapter rejects exhausted capacity. Python had masked this with its own capacity check. One intermediate relink accidentally reused the staged old adapter; its failed raw-capacity run is retained as negative evidence, not counted as a pass. Final relink explicitly compiles the current repository adapter. Final read-only review found no remaining blocker in the supported paths.

## Reproduction and artifacts

`tools/test_fast_hse_deep_recovery_mpi.py BUILD GS NEW_OUTPUT --steps 32 --fault-step 13` uses two ranks and an external timeout. Add `--split` for four ranks/two groups, `--overflow-A 1e100` for production PT-CN rejection, or `--faults` for a separately built injection library. `--capacity-probe` bypasses Python admission after the accepted trajectory to verify raw ABI rejection; it deliberately ends that native lifetime. All destinations must be new.

Persistent local evidence: `/Users/otobetoshihito/.codex/outputs/salmonpy-fast-hse-deep-recovery-20261117`. Inventories retain hashes of rejected payloads, accepted trajectories, logs and build receipts; selected raw before/after pairs are retained. Build records preserve original external source/build locations and a relocated runtime library. They are host-specific evidence, not portable binary releases.
