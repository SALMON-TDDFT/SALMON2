# Retained GS context storage validation

Local Mac verification, 2026-11-15. This extends the existing state-isolation checklist without adding XC coverage.

## Implementation

The former typed local work objects in main_dft now reside in internal `s_gs_context`: system/grid/decomposition objects, orbitals, density/potentials, pseudopotential data, mixing/optimization/band/DC storage, XC handle and eligible HSE owner. `gs_run` preserves the numerical statement order and retains the state after output. `gs_close` finalizes XC, unbinds/clears its HSE resources and releases the stored allocatables. The external `main_dft` wraps run/close.

The actual context must be declared TARGET and kept at a stable address because system retains a pointer to its HSE member. Live contexts must not be copied. The internal entry rejects a second run on a live object with status 1; close is idempotent. Shared input, timers, parallel/DC setup and native process-terminating errors remain. This is not a stable handle ABI, general concurrent GS API or complete explicit-input isolation.

## Results

- The new probe failed to compile against the preceding build because the context module was absent (`/tmp/salmonpy-gs-context-red-20261115/build.log`).
- Four tiny scalar-HSE/pre-SCF configurations: A remains live while B runs, then B closes. A's stored trace remains byte-identical before/after each operation. Trace includes orbitals, densities/potentials, occupations, forces, allocated density/potential mixing-history arrays and selected HSE caches. It is not an exhaustive comparison of every field, especially PP/DC internals.
- Both A/B scalar outputs and all reference-generated restart payloads remain byte-identical to original fast-HSE separate-process execution. Restart/trace payloads are deleted after verification.
- The main_dft wrapper separately passes all four serial GS cases with explicit HSE ownership required.
- DC/MPI two- and four-rank layouts, four pre-SCF configurations each, match original outputs. These use the legacy HSE resource route; no independent DC contexts are claimed.
- Eight full TDDFT interleaving/close-reopen cases remain identical.
- Serial and MPI HSE builds and HSE-disabled build succeeded. Source and pre-link/all-object forbidden-eigensolver checks passed. One ordinary PZ smoke test through the HSE-disabled wrapper matches the prior HSE-disabled build; this is regression protection, not PZ context independence validation.
- Source Python suite: 305 passed. sdist→wheel→fresh installation: 305 tests, one optional plotting skip. Installed lifecycle probe repeats all four retained-state cases successfully; installed sources and JP/EN/AI guides match the tested checkout.

Independent read-only review found no numerical-order regression and identified two improvements: explicit TARGET/address contract and tracing mixing/occupation/force state. Both are incorporated. No CI was run; maximum aggregate parallelism was four.

## Local evidence

- `/tmp/salmonpy-gs-context-v2-{serial,mpi}-build-20261115.json`
- `/tmp/salmonpy-gs-context-v2-retained-20261115/result.json`
- `/tmp/salmonpy-gs-context-v2-wrapper-20261115/result.json`
- `/tmp/salmonpy-gs-context-v2-dc2-20261115/result.json`
- `/tmp/salmonpy-gs-context-v2-dc4-20261115/result.json`
- `/tmp/salmonpy-gs-context-v2-tddft-20261115/result.json`
- `/tmp/salmonpy-gs-context-v2-nohse-smoke-20261115`
- `/tmp/salmonpy-gs-context-v2-wheel-retained-20261115/result.json`
- `/tmp/salmonpy-gs-context-v2-dist-20261115`
- `/tmp/salmonpy-gs-context-v2-python-tests.log`
- `/tmp/salmonpy-gs-context-v2-installed-tests.log`

The immutable original fast-HSE reference remains commit `c9e8ce431f43eefd4a88998134adb4796f8ec837`. Tests use deliberately coarse synthetic control-flow inputs and do not establish physical convergence, different-input independence, concurrent execution or checkpoint consumption by TDDFT. Earlier audited binaries/sources were left unchanged.
