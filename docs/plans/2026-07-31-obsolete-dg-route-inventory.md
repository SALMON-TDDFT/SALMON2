# Obsolete DG Route Inventory

Status: removal and production acceptance completed through Task 8 of
`2026-07-31-remove-obsolete-experimental-dg-routes.md`.

The only retained experimental DG production path is:

1. buffer-supported, symmetry-preserving overlapping-Wannier GS;
2. its `SALMON_OW_GS_CHECKPOINT_V3` checkpoint; and
3. generalized-eigenvalue exponential coefficient RT, including restart.

Normal SALMON routes remain protected. In particular, normal DC LCFO and
EigenExa, plus their ordinary Wannier90 consumers, are not DG-removal
candidates.

## Inventory method

The inventory was recorded at plan HEAD `24449148` with:

```text
rg -n '<symbol>' src tests samples
rg --files src/common src/gs/dc src/rt/dg tests/dg samples
```

No existing build directory was present in the supplied worktree, so there
was no configured `<existing-build>` on which to run `cmake --build
<existing-build> --target help`. CMake membership was instead read directly
from the source lists below; the clean overlay build in Task 8 remains the
authoritative configured-build check.

## Public input and dispatch inventory

| Route class | Input surface | Dispatch/consumer boundary | Disposition |
|---|---|---|---|
| direct real-space DG GS | `yn_dg_dc_local_periodic`, `yn_dg_fragment_from_dcdft` | `src/gs/main_dft.f90`, `src/gs/scf_iteration_dft.f90` | remove |
| DG-Fragment RT | `yn_dg_frag`, `yn_dg_fragment_rt`, `yn_dg_hse_ace` | `src/rt/main_tddft.f90`, `src/rt/initialization_rt.f90` | remove DG-only controls and route; `yn_dg_hse_ace` is consumed only by the DG HSE exchange |
| Nodal RT | `yn_dg_nodal_rt` | `src/rt/main_tddft.f90` | remove |
| WPW GS/RT | `yn_dg_wpw_checkpoint_rt`, `yn_dg_wpw_production`, all `yn_dg_wpw_*` controls | `src/gs/main_dft.f90`, `src/rt/main_tddft.f90` | remove |
| mixed-z | all `yn_dg_mixed_z*` controls and `dg_mixed_z_*` selectors | `src/rt/dg/rt_dg_fragment*.f90`, mixed-z diagnostics | remove |
| full-H seed | `yn_dg_full_h_eigen_seed`, `yn_dg_full_h_wannier_band_gauge`, `yn_dg_lcfo_seed_exhaustive_check`, `dg_wavefunction_seed_filename` where DG-only | `src/gs/dc/lcfo*.f90`, RT initialization | remove only after normal LCFO/Wannier90 consumer audit |
| adaptive DG basis | `yn_adaptive_basis`, `yn_adaptive_basis_dg`, `yn_adaptive_basis_update` where DG-only | SCF and RT DG update paths | remove DG promotion/update behavior; do not alter unrelated normal behavior without proof |
| experimental promotion/ExpDiag | `yn_dg_subspace_diag`, all `yn_dg_expdiag_*` controls | `src/rt/dg/rt_dg_integrator_expdiag.f90`, DG promotion branches | remove |

The removal checker carries the exact forbidden input set. The four retained
DG inputs are `yn_dg_dc_overlapping_wannier`,
`yn_dg_overlapping_wannier_rt`,
`yn_dg_overlapping_wannier_rt_restart`, and `yn_dg_length_gauge`.

## Implementation and CMake inventory

Direct GS removal candidates:

```text
src/common/dg_dc_direct_sipg.f90
src/common/dg_ground_state_checkpoint.f90
src/gs/dc/dg_dc_ground_state.f90
src/gs/dc/dg_dc_ground_state_adapter.f90
src/gs/dc/dg_dc_handoff.f90
src/gs/dc/dg_dc_local_basis_ground_state.f90
```

WPW candidates are every file under these prefixes, plus
`src/gs/dc/dg_wannier_pw_scf.f90`:

```text
src/common/dg_wpw_
src/gs/dc/dg_wpw_
src/rt/dg/rt_dg_wpw_
```

DG-Fragment/Nodal candidates are:

```text
src/rt/dg/rt_dg_fragment*
src/rt/dg/rt_dg_nodal_*
src/rt/dg/rt_dg_integrator_*
```

Their current CMake entries are in `src/common/CMakeLists.txt`,
`src/gs/dc/CMakeLists.txt`, and `src/rt/CMakeLists.txt`. Residual generic DG
helpers must be retained only when the accepted overlapping-Wannier route has
a direct consumer.

## Focused test and sample inventory

Removal candidates are:

```text
tests/dg/*dg_dc_direct*
tests/dg/*dg_dc_ground_state*
tests/dg/*dg_dc_handoff*
tests/dg/*dg_dc_local_basis*
tests/dg/*dg_ground_state_checkpoint*
tests/dg/*dg_wpw*
tests/dg/*nodal*
tests/dg/*wpw*
samples/exercise_dg_fragment_C2H2/
samples/exercise_dg_fragment_rt/
samples/exercise_dg_rt_hse_test/
```

The obsolete assets above, together with their Full-H, mixed-Z, direct-AMN,
Krylov, and fragment-operator source checks, have been removed. The checker
now rejects their reintroduction, stale CMake registrations, and historical
documents that present a removed route without a `Historical/removed` notice.

## Protected retained consumers

| Protected surface | Evidence/consumer |
|---|---|
| overlapping-Wannier construction | `src/gs/dc/dg_overlapping_wannier_construction.f90` and its focused MPI fixtures |
| V3 route checkpoint | `src/gs/dc/dg_overlapping_wannier_checkpoint.f90`, Si64 gate, restart hash checks |
| generalized-eigenvalue Exp coefficient RT | `src/rt/dg/rt_dg_overlapping_wannier.f90`; calls `zhegv`; no Crank/conventional RT fallback |
| buffer/core support | `src/common/dg_buffer_window_projector.f90`, `src/gs/dc/dg_dc_buffer_core_faces.f90` while direct consumers remain |
| generalized algebra | removed after the Task 4 direct-consumer audit found no retained consumer |
| normal DC LCFO | `src/gs/dc/lcfo.f90` and `src/gs/dc/CMakeLists.txt` |
| EigenExa | `src/gs/eigen_subdiag_eigenexa.f90` and `src/gs/CMakeLists.txt` |
| Wannier90 | `src/gs/dc/lcfo_wannier_*` normal LCFO sources; remove no symbol without an exclusive-consumer proof |

## Review classification

- Critical: deleting or weakening any overlapping-Wannier GS/V3 checkpoint/
  coefficient-RT requirement, or changing normal LCFO/EigenExa behavior.
- Important: an unlisted executable input, dispatch, CMake source, focused
  test, or sample belonging to a forbidden route.
- Minor: historical prose that is already clearly labelled removed and is not
  linked as executable guidance.

The Task 1 specification review expanded the plan's five example forbidden
flags to the complete source-observed family above. Consumer audit confirmed
that the generically named `yn_adaptive_basis` is nevertheless exclusive to
DG-Fragment RT, and that `yn_dg_hse_ace` is exclusive to the DG HSE exchange.
The code-quality review requires identifier-boundary matching, deterministic
sorted output, reports all findings in one run, skips binary files, and
excludes only the checker itself from token scanning.

Task 8 confirmed the final inventory in a clean archive overlay: the only
source below `src/rt/dg/` is `rt_dg_overlapping_wannier.f90`. The retained
production chain, normal DC protection, hashes, and review disposition are
recorded in `2026-07-31-obsolete-dg-route-removal-results.md`.
