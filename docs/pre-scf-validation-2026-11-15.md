# Owned GS/DC pre-SCF controller validation

Local Apple Silicon/macOS verification on 2026-11-15. No GitHub Actions or network runtime builds were used. This extends the owned evaluation phase recorded in `xc-phase-validation-2026-11-15.md`.

## Implementation and scope

Each `s_dft_system` retains pre-SCF threshold, required consecutive steps, current count and first iteration. Construction captures settings; other SCF entry points initialize only if absent. Each SCF invocation resets count/start without reactivating a completed phase. Local readiness (including first-iteration exclusion), collective reduction, count update, diagnostic, transition, mixing reset and potential rebuild remain in the original order.

GS, orbital solves, DC fragment/total energy and hybrid metadata consult the system phase. DC total-system phase is synchronized before its energy evaluation. HSE-disabled compilation retains a false pending-phase guard and does not depend on the HSE-only module.

This does **not** isolate DC communicators/material swapping, adaptive exchange globals or legacy shared exchange resources. The actual DC/MPI collective route was subsequently executed in the two-rank and four-rank comparisons below. The controller probe separately models rejection from another fragment after reduction. Neither test proves simultaneous independent solver isolation. No public Python API is added.

## Tests and results

- `tools/test_pre_scf_state.py`: two controllers with different thresholds, counts and starting iterations, alternated while legacy globals are poisoned, match independent original counting rules. A modeled fragment rejection resets the count. Invalid settings preserve the existing phase/configuration. The probe failed to compile against the preceding build before the new contract existed; it passes against the new build.
- `tools/test_pre_scf_runtime.py`: four serial GS cases (disabled warmup and 1/2/3 consecutive steps) compared with the immutable original fast-HSE PR executable. Phase/energy diagnostic traces and `gs_info.data`, `gs_eigen.data`, `gs_k.data` are byte-identical. The GS information includes total energy and forces. Generated synthetic restart payloads are removed by the tool.
- These GS fixtures use one Si atom, an 8³ grid, two k points, two states and deliberately loose thresholds. They check control flow and numerical regression, **not physical convergence**. Strict equality here is a same-platform migration check, not a portable accuracy criterion.
- `tools/test_tddft_contexts.py`: eight cases with differing timestep, OpenMP width and trajectory capacity remain identical to separate-state runs. Includes legacy phase poisoning. The common-grid/material guards remain in effect.
- `tools/test_hse_rt_state_runtime.py`: 32 cases, k meshes of 2/3/6/8 points, two timesteps, Taylor/ACE/full/PTCN variants. Current and energy maximum absolute differences from the immutable original PR are zero.
- HSE ON and OFF builds succeeded. Source/object forbidden-eigensolver audits passed. HSE OFF also completed a tiny PZ GS smoke calculation.
- Python source suite: 305 tests passed.
- Source distribution → wheel → fresh environment: 305 tests, one optional plotting skip, no failures. Installed guides, AI JSON, new probes and freshly prepared replacement source bytes match the checkout. Manual synchronization and changed-line checks passed.

## Local evidence

These paths identify this machine's observations; they are not installation requirements or portable commands.

- Build record: `/tmp/salmonpy-prescf-build-20261115.json`.
- Original reference record: `/tmp/salmonpy-hse-rt-state-reference-build-20261115.json`, commit `c9e8ce431f43eefd4a88998134adb4796f8ec837`.
- Controller: `/tmp/salmonpy-prescf-probe-20261115/result.json`.
- GS comparison: `/tmp/salmonpy-prescf-gs-20261115/result.json`.
- Independent RT contexts: `/tmp/salmonpy-prescf-contexts-20261115/result.json`.
- Original RT comparison: `/tmp/salmonpy-prescf-default-20261115/result.json`.
- HSE OFF build/smoke: `/tmp/salmonpy-prescf-nohse-build.log`, `/tmp/salmonpy-prescf-nohse-smoke-20261115/run.log`.
- Packaging/tests: `/tmp/salmonpy-prescf-dist-20261115`, `/tmp/salmonpy-prescf-python-tests.log`, `/tmp/salmonpy-prescf-installed-tests.log`.

The build uses gfortran/OpenBLAS/Libxc on the local Mac, at most two build jobs and a maximum aggregate of four execution threads. Cross-platform and physically converged DC validation remain outstanding.

## DC/MPI follow-up

The same immutable source versions were built locally with MPI enabled, GNU C/Fortran, OpenBLAS, Libxc and FFTW. Source checks and pre-link/all-object PZHEEV prohibition audits passed. Two fragments were run with one spatial rank each (2 ranks) and two spatial ranks each (4 ranks), all with OpenMP/BLAS threads set to one. The cases use two Si atoms, a 16×8×8 total grid, Γ only, four states per fragment, four buffer cells in x, finite electronic temperature and LCFO disabled. Pre-SCF disabled and consecutive-step requirements 1/2/3 each passed for both layouts: **eight comparisons**.

For each layout/case, the modified and original fast-HSE executables produced identical phase/energy diagnostic traces, including DC total energy, and byte-identical information/eigenvalue/k-point files for both fragments. This is a same-layout regression check; we do not assert bit identity between different MPI decompositions. Deliberately loose convergence thresholds keep this a control-flow check rather than a converged materials calculation. No Fortran changes were needed in this follow-up.

Reproduce using `tools/test_pre_scf_runtime.py CURRENT_MPI_RECORD ORIGINAL_MPI_RECORD NEW_OUTPUT --dc-ranks 2` (or `4`). The helper uses standard input, Cartesian atom coordinates, and `write_gs_restart_data='no'`, required by the original spatial-hybrid implementation. Optional `--mpiexec PATH` selects the launcher.

Local setup findings: Homebrew `mpicc` selected Apple Clang and rejected the GNU `-fopenmp` flag during a header check; selecting GNU C explicitly resolved configuration. Sandboxed local socket creation prevented MPI startup; rerunning with local process-network permission resolved it. These were setup failures, not numerical failures.

Evidence:

- `/tmp/salmonpy-prescf-mpi-current-20261115.json`
- `/tmp/salmonpy-prescf-mpi-reference-20261115.json`
- `/tmp/salmonpy-prescf-dc2-final-20261115/result.json`
- `/tmp/salmonpy-prescf-dc4-final-20261115/result.json`

Follow-up packaging also passed: sdist → wheel in `/tmp/salmonpy-prescf-mpi-dist-20261115`, fresh installation in `/tmp/salmonpy-prescf-mpi-installed-20261115`, 305 source tests and 305 installed tests (one optional plotting skip). Installed guide/JSON/tool bytes match the checkout. The installed tool ran from `/tmp` with isolated Python and repeated all four two-rank DC cases successfully (`/tmp/salmonpy-prescf-dc-wheel-20261115/result.json`). The serial four-case regression also passed with the extended tool. No CI runs were triggered.
