# HSE/SOC integration review guide

This proposal targets `develop-2.0.0` at `376635ebd939687959021a2c4e87747802401c4a`.
It contains the accumulated HSE/SOC development, not merely the latest FFT fix.
The 2026-11-13 preparation excludes developmental giant-Gamma source-support
ACE and finite-range SR group/wave FFT. Ordinary Gamma full-support exchange
is retained. See [scope and September comparison](hse-pr-scope-20261113.md).
The ordinary upstream tests remain unchanged; the proposed hybrid regression
uses only one small Si system through the existing420 GS /421 RT cases.

## Review order

| Area | Main entry points | Review focus |
|---|---|---|
| Build | `CMakeLists.txt`, `cmakefiles/` | CPU defaults and optional-library linkage |
| Native hybrid exchange | `src/xc/exx_native.f90`, `exx_k_exchange.f90`, `exx_ace.f90`, `hybrid_semilocal.f90` | kernel/semilocal consistency and distributed ownership |
| Symmetry and worker layouts | `exx_symmetry.f90`, `exx_spinor_wave_workers.f90`, `exx_worker_layout.f90`, `exx_fft_consensus.f90` | scalar/SOC reconstruction, exact coverage and collective failures |
| RT and restart | `src/rt/`, `src/io/`, `src/xc/exx_factor_history.f90` | propagation, finite-temperature/occupation constraints and saved-state compatibility |
| Pseudopotential / noncollinear support | `src/atom/pp/`, `src/common/`, `src/symmetry/`, `utility/crps_soc_nlcc/` | KY projector/SOC/NLCC/core-density consistency |
| Inputs and the public regression | `src/io/`, `docs/inputs/`, `testsuites/420_bulk_Si_hse_gs`, `testsuites/421_bulk_Si_hse_rt` | input validation and reproducible GS-to-RT results |

Some optional exchange/functional paths share these dependencies. No unrelated
runtime source module was removed merely to reduce the file count. Input changes
and restrictions are documented under [inputs](inputs/hse.md).

## Deferred standalone Si512 work

The standalone giant-Gamma Si512 owner-exchange snapshots and their later
optimization candidates are outside this proposal. They are preserved unchanged
under `modules/si512_owner` in the original base commit `5926312f` and its
original worktree. No SALMON source, build or test depends on them. Their
standalone validation is preserved; native SALMON/RT/SOC integration remains
separate work. Their exclusion is a scope decision, not a rejection of the
standalone numerical or performance results.

## What is deliberately excluded

Development plans, benchmark/raw-output archives, measurement reports, material
status files, standalone probes, extra regression systems, local job-management
scripts and historical development notes are outside the submitted file diff.
They are preserved in the pre-cleanup commit `a0da2902` and a verified local
archive; they are not runtime dependencies or current execution instructions.
No historical evidence is represented as a current testsuite result.

Workspace/FFT reuse with packed transforms and later batched contraction were
not adopted: Fugaku comparisons did not establish the required performance
advantage. Their code and algorithm switches are absent from the submitted tree.
See [selection evidence summary](dense-k-integration-20261106.md).

## Validation limits

The single Si regression is a smoke/numerical integration check. It does not
certify SOC, every optional functional, all MPI layouts, restart combinations,
Fujitsu or GPU execution. Broader previous local development checks are recorded
as historical evidence, not shipped as additional public tests. This proposal
should remain a draft until maintainers review the broad scientific/build scope.
See [reproduction and results](PR_VALIDATION.md).
