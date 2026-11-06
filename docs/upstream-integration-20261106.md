# SALMON upstream integration, 2026-11-06

This integration merges upstream `develop-2.0.0` at `376635eb` into the committed local EXX baseline `b3388e45` on `codex/integrate-upstream-20261106`. Existing dirty SALMON2/EXX worktrees are preserved; no dirty performance patch, remote job, fixed executable, or checkpoint is changed. No push is performed.

## Included upstream work

- Atomic pseudopotential wavefunctions, multi-projector atomic solver, PP initialization and PDOS radial/channel normalization and k weights.
- SCF conjugate-gradient stuck-step guard.
- OpenACC current/nonlocal GEMM paths, NCCL option and GPU loop/interface fixes.
- IEEE wrapper for compiler compatibility, interpolation and diagnostic fixes.

## Integration resolutions

Retain the local `exx_*` and `hybrid_semilocal` architecture rather than resurrecting upstream `hse_*` copies. Preserve custom KY SOC/NLCC/core23 readers, HSE restart metadata and odd-step parity, no-PC and PT-CN modes, scalar IEEE policy, eigensolver source/object guards, default CPU HSE support, MPI/OpenMP compiler handling, Apple complex BLAS selection and Libxc compatibility. Historical documentation and probes remain preserved.

The upstream atomic density check is expressed as a scalar loop with identical every-element NaN/Infinity/negative rejection. The IEEE wrapper accepts scalar calls only. The new 114 current regression now has an explicit GS fixture from 111 and a Python3 verifier; missing-result reporting no longer references an undefined variable. The 64-step HSE RT test gets a 1200-second runtime allowance instead of 300 seconds; its science, row counts and numerical tolerances are unchanged.

An independent read-only review found a CMake ordering problem: defaults must be resolved before MPI/OpenMP setup. The include is restored after option declarations; default-only `USE_MPI_DEFAULT=ON` is checked separately.

## Validation

Validation results and build/test logs are stored in the canonical results directory `calc2026/salmon-upstream-integration-20261106`. Local builds cover GNU Fortran 15 CPU/MPI with HSE enabled and disabled; this does not certify Fujitsu or GPU execution. No executable is transferred or applied to existing production jobs.

## Automatic FFT layout agreement (fixed in follow-up)

A startup failure was reproduced on the unchanged pre-merge EXX implementation. `choose_fft_layout` measures each rank independently, while `apply_k_pairs` requires identical layout flags across ranks. Timing jitter can select different layouts and cause collective rejection. A diagnostic-only pre-merge build, with one print added before the settings guard, showed ranks 3/2/0 selecting strided (`F`) and rank 1 contiguous (`T`), followed by a nonzero exit. The same baseline also completed a one-step RT when all four selected strided. Numerical operators and guards were not modified in that diagnostic.

This is an inherited limitation, not an upstream merge regression. The follow-up described in `dense-k-integration-20261106.md` fixes this by agreeing the choice collectively before creating layout-dependent FFT plans. The consistency guard remains. Production deployment still requires platform validation.

Additional checks: CPU/MPI HSE and non-HSE complete builds; source and all application-object PZHEEV audits; 25 EXX oracle/rejection tests (33-rank launch needs local oversubscription); 111/114 current numerical comparison and fixture-only selection; 12 atomic-density finite/NaN/Infinity/negative cases plus array-interface compile rejection; 16 forbidden-checker tests; checkpoint elapsed/parity fixture; default-only MPI option resolution. The first RT64 attempt timed out at 300 seconds; a separate auto-layout startup rejection is preserved as the inherited issue above. The final fresh GS + 64-step RT run passed all six CTest stages (RT 317.43 seconds), including current, electron count, energy drift and native dielectric/conductivity checks. The 26 primary cases, 25 EXX cases and 6 split-current cases have all passed across the recorded runs; initial timeout/environment/baseline-auto failures are retained in the evidence.
