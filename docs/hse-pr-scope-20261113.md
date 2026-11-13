# HSE integration scope — 2026-11-13

The working base is `codex/hse-soc-pr-20261106` at
`5926312f66a19ee6afdd0a6cbd68ac99039933de`, not the September
`hse-merge-prep` branch. The original branches remain preserved.

## Approved scope

Retain ordinary HSE06/PBE0/PBEh, dense-k MPI + OpenMP exchange,
SOC/NLCC, occupied-support ACE, restart and ordinary Gamma calculations.
Exclude the developmental giant-Gamma SR group/wave FFT and
source-support ACE routes. Shared full-support spatial/DC exchange is not
giant-Gamma-specific and must remain available.

Dense-k learned factor history is separate from its Gamma source-support
adapter. Excluding the latter does not authorize deleting the former.
Standalone Si512 owner/candidate modules are not linked into SALMON and are not
evidence of integrated runtime support. They are also omitted from this PR;
all eleven files remain preserved unchanged in base commit `5926312f` and the
original worktree. No source/build/test dependency refers to these snapshots.

## Changes missing from the September PR

| Area | Selected November base relative to September |
|---|---|
| Functionals | Shared EXX classification and kernels; PBE0/PBEh and rVV10 integration |
| Dense k points | Hermitian pairs, crystal orbit reuse, MPI + OpenMP FFT work, worker layouts and services, collective FFT selection |
| Parallel memory | Distributed ACE, spatial/orbital ownership and bounded exchange workspaces |
| GS/RT connection | Conventional and DC mesh paths, PBE initialization and occupation validation |
| Spin | Native SOC exchange, noncollinear support and NLCC/core-density handling |
| Propagation | Expanded PT-CN and factor-history implementation, restart/checkpoint handling |
| Portability | Scalar IEEE checks and accumulated compiler fixes; platform execution still requires validation |
| Upstream | Mainline changes through `376635eb`; these are not all hybrid changes |

A path-based comparison of `src/` against `483058cf` found 75 added,
7 removed and 58 modified files. Renames count as additions/removals.
Presence in source does not establish numerical or capacity validation.

## Build safety

The restored September tree fails the workspace prohibition checker at
`src/gs/dc/lcfo_diag_chefsi_complex.f90`; it must not be built or executed.
The selected November source passes that checker. However, `05359b43`
removed the configure/incremental/application-object guard, so that guard
must be restored before building the new preparation tree.

No GS data is collected as part of this preparation. Existing remote jobs,
fixed binaries and historical evidence are outside the edit scope.

## Other development that is not integrated

The separate `work/kli-sic-development` repository has native KLI/HSEX
screening work through `67f7674` (2026-11-12), including optical-limit
screening, Sternheimer and screened-interaction modules. These modules are
absent from the selected HSE/SOC base. They are outside the approved
HSE/PBE0/PBEh integration scope and need a separate scientific review.
`work/hsex-delta-w` has a stale worktree pointer to an unavailable session
path, so its branch state could not be established; no changes were made there.

## Verification on the reduced tree

- GNU Fortran 15 Release: HSE ON/MPI ON, HSE ON/MPI OFF, and HSE OFF/MPI OFF
  application builds passed. Installed OpenMPI, FFTW, Libxc and OpenBLAS were
  used where applicable. These are CPU checks.
- Source and application-object prohibition checks passed. A synthetic fixture
  confirmed safe build, configure rejection, incremental reintroduction rejection,
  and rejection/quarantine when a stale application object references PZHEEV.
  The fixture never calls the prohibited driver.
- Nine application input fixtures confirmed source ACE, positive/NaN/infinite
  SR tolerance and the four removed namelist controls are rejected. Ordinary
  Gamma input passes validation and stops at a deliberately invalid pseudo fixture.
  No GS/restart wavefunctions were loaded or produced.
- The changed spatial wrapper was checked against an independent periodic direct
  convolution for screened HSE and the existing full-range Coulomb convention.
  Occupied, fractional and zero weights, Hermitian pairs and ordered pairs,
  MPI1/2/4/8 orbital layouts and a mixed spatial/orbital MPI4 layout passed.
  Maximum absolute action error was `4.766466e-14` (acceptance `1e-10`).
- Six collective rejection modes passed: NaN/Infinity weights, negative weights,
  NaN/Infinity orbitals and inconsistent spatial weights. Together with the five
  normal layouts, the synthetic runner passed eleven cases.
- Independent read-only review found no concrete production-code regression.
  The guard file was explicitly added despite the case-insensitive ignore rule.

These checks do not certify complete GS/RT trajectories, compact MLWF fallback,
every restart layout, Fujitsu/A64FX, ScaLAPACK or GPU execution. No new performance
claim or physical-spectrum convergence claim is made. The original PR's historical
Si regression does not substitute for a system test of this reduced tree.

Reproduce with an MPI Makefiles build after the normal guarded configure/build:

```sh
python3 tools/tests/check_eigensolver_cmake_gate.py
python3 tools/tests/check_hse_pr_scope.py /path/to/build/salmon
python3 tools/tests/run_spatial_exchange_scope.py /path/to/build
```

Local test-build SHA256 values (not production acceptance):

| Build under `/private/tmp/` | Executable SHA256 |
|---|---|
| `hse-pr-cleanup-mpi/salmon` | `d1c1c76b5991afd6ef17a723906d5d5f4f916c88e2fba90d3f7e8c41b88685fd` |
| `hse-pr-cleanup-serial/salmon` | `79ab13b3fbd7c58e5f063df23a635561175e0aa5813955d0d5b31f64a5b31596` |
| `hse-pr-cleanup-off/salmon` | `64960b434dbf65ffb769dfb0e08f8c895186a5428c713d8e1e377795b9931f6c` |

Verify the recorded hash before reusing an unchanged test executable. These
temporary paths do not imply long-term binary storage or cloud synchronization.
