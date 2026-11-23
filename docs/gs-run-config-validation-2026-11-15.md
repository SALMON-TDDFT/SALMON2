# Explicit GS SCF configuration validation

Local Mac verification on 2026-11-15; no GitHub Actions. This continues the existing shared-input extraction, not new XC coverage.

## Change

`s_gs_run_config` contains maximum SCF iterations, convergence threshold, pre-SCF threshold and consecutive-step requirement. `gs_input_config` captures old inputs; optional `gs_run(config=...)` validates before activating and copies the value into its system. Without an explicit value, the same inputs are captured once. SCF loop/convergence/adaptive-readiness checks consume owned values. Unconverted callers retain a legacy fallback.

Finite/nonnegative threshold, positive step count and bounded positive iteration count checks precede allocation. Explicit pre-SCF also repeats the original semantic restrictions: fresh static hybrid DFT, no restart/MD/optimization, positive DC temperature, no checkpoint/shutdown schedule, and density mixing/convergence metric. The external wrapper reports rejected run configuration as nonzero failure; the pre-existing invalid-periodicity band early return is retained.

Independent review identified the semantic-check bypass and wrapper's ignored status. Both were reproduced and fixed. The old semantic probe reached STOP 10 because checkpoint scheduling was wrongly accepted (`/tmp/salmonpy-gs-config-semantic-red-20261115`). The old wrapper accepted nscf=0 with exit 0 (`/tmp/salmonpy-gs-config-wrapper-red2-20261115`). The updated tool checks explicit rejection and the wrapper's nonzero failure.

## Verification

- Initial contract probe failed to compile before the new type/entry existed (`/tmp/salmonpy-gs-config-red-20261115`).
- Four GS cases: capture settings, invalidate all four legacy controls, run A with explicit settings, mutate caller A settings, run B with a different iteration cap. A's retained trace remains unchanged; both scalar/checkpoint results match the original separate-process reference.
- Numeric invalid cases (zero/overflow iteration cap, NaN threshold, zero step count, negative pre-SCF threshold) and semantic invalid cases (checkpoint schedule, potential metric, restart) are rejected before subsequent valid initialization. Live-context reentry and idempotent close remain tested.
- Four compatibility-wrapper GS cases, eight DC/MPI comparisons, eight TDDFT interleaving cases passed.
- HSE serial/MPI and HSE OFF builds passed source/pre-link/object forbidden-eigensolver audits. An HSE OFF PZ wrapper smoke comparison matches its prior build; it does not test independent PZ contexts.
- Source Python suite: 305 passed. sdist→wheel→fresh environment: 305 tests, one optional plotting skip. Installed GS probe repeats all four cases and wrapper rejection checks. Installed source and JP/EN/AI guides match checkout bytes.

## Evidence

- `/tmp/salmonpy-gs-config-v2-{serial,mpi}-build-20261115.json`
- `/tmp/salmonpy-gs-config-v2-retained-20261115/result.json`
- `/tmp/salmonpy-gs-config-v2-wrapper-20261115/result.json`
- `/tmp/salmonpy-gs-config-v2-dc2-20261115/result.json`
- `/tmp/salmonpy-gs-config-v2-dc4-20261115/result.json`
- `/tmp/salmonpy-gs-config-v2-tddft-20261115/result.json`
- `/tmp/salmonpy-gs-config-v2-nohse-smoke-20261115`
- `/tmp/salmonpy-gs-config-v2-wheel-retained-20261115/result.json`
- `/tmp/salmonpy-gs-config-v2-dist-20261115`
- `/tmp/salmonpy-gs-config-v2-python-tests.log`
- `/tmp/salmonpy-gs-config-v2-installed-tests.log`

The reference remains original fast-HSE commit c9e8ce431f43eefd4a88998134adb4796f8ec837. Physical convergence, different material/grid/XC independence, concurrency and full failure atomicity are not established. Convergence metric, mixing, initialization and other material/parallel inputs remain shared. All native builds used new external directories; aggregate build/execution parallelism was at most four.
