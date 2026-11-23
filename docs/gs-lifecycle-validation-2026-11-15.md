# Sequential GS lifecycle and checkpoint validation

Local Mac verification on 2026-11-15. This closes the same-input repeated-GS and generated-checkpoint comparison gaps recorded in `gs-hse-owner-validation-2026-11-15.md`, without extending the independent XC or parallel solver support.

`tools/test_gs_lifecycle.py` compiles `tests/native/test_gs_lifecycle.f90` against an unchanged audited serial build. The probe reads one input, invokes the full `main_dft` twice with distinct output directories, and requires explicit GS HSE ownership for both calls. Each run is compared to a separate execution of the immutable original fast-HSE executable. Pre-SCF disabled and 1/2/3 consecutive-step configurations all pass.

For both repeated calls, `gs_info.data`, `gs_eigen.data`, `gs_k.data` and reference restart payloads are byte-identical. Payloads include `wfn.bin`, `rho_inout.bin`, `occupation.bin`, `info.bin`, atomic coordinates and velocities. Their SHA256 values are retained; generated restart directories are removed after comparison. This is output equivalence, not a TDDFT checkpoint-import test.

The fixed fixture is the earlier tiny Si/HSE control-flow case: one atom, 8³ grid, two k points and two states, loose thresholds. Same-input sequential reuse is verified. Different-input reuse, concurrent GS instances, physical convergence, MPI GS reuse and general native failure recovery remain unverified. No solver changes or full native rebuild were required for this follow-up.

Evidence:

- `/tmp/salmonpy-gs-lifecycle-final-20261115/result.json`
- Build `/tmp/salmonpy-gs-owner-serial-build-20261115.json`
- Reference `/tmp/salmonpy-hse-rt-state-reference-build-20261115.json`, original commit `c9e8ce431f43eefd4a88998134adb4796f8ec837`

The initial probe incorrectly attempted to call the main program's internal `set_basic_flag`; inspection established that ordinary DFT only sets `calc_mode='GS'`. The corrected probe uses this same assignment. No solver failure was observed.

Packaging follow-up: source Python suite 305 passed; sdist→wheel→fresh environment 305 tests, one optional plotting skip, no failures. The installed tool independently passed all four sequential/checkpoint cases (`/tmp/salmonpy-gs-lifecycle-wheel-20261115/result.json`). Installed test sources and JP/EN/AI guides match the checkout. Checkpoint payload deletion was verified for both final test runs. Distribution and logs: `/tmp/salmonpy-gs-lifecycle-dist-20261115`, `/tmp/salmonpy-gs-lifecycle-python-tests.log`, `/tmp/salmonpy-gs-lifecycle-installed-tests.log`. No CI runs were triggered.
