# DG PW-Augmented Basis Audit

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

## Existing Complete Paths

- Input and global controls are already wired through `yn_plane_wave_basis`,
  `n_plane_waves_dg`, `k_cutoff_plane_wave`, `yn_dg_subspace_diag`,
  `dg_subspace_extra_states`, and `dg_subspace_pw_vectors`.
- `rt_dg_plane_wave.f90` can initialize selected plane waves, build fragment-PW
  overlap blocks, fragment-PW Hamiltonian blocks, PW diagonal Hamiltonian terms,
  PW momentum/gradient blocks, and mixed dense Hamiltonian diagnostics.
- `prepare_mixed_basis_startup` prepares `S_mat_frag_pw`, `H_mat_frag_pw`,
  `P_mat_frag_pw`, and PW diagonal/operator data without dense diagonalization.
  The intended startup state is fragment coefficients unchanged and PW
  coefficients initialized to zero.
- `update_density_hamiltonian_stage` refreshes PW mixed operator blocks after the
  density/Hamiltonian stage by recomputing fragment-PW overlap/Hamiltonian data
  and calling `build_mixed_hamiltonian`.
- `rt_dg_fragment_ops.f90` already has `apply_mixed_hamiltonian`, PW coefficient
  row ownership, remote PW row fetch, full PW coefficient cache refresh, and a
  full coefficient view helper that returns fragment and PW coefficient blocks.
- `rt_dg_density_reconstruct.f90` has a mixed/PW density route. It can use
  `mixed_transform`/`coef_mix` when a canonical mixed basis is ready, and also
  has direct `coef_pw_full_cache` materialization for orbital density.

## Stop Gates

- `rt_dg_integrator_derivative.f90` stops when `use_plane_wave_basis` is true or
  `coef_pw` is allocated:
  `DG derivative supports the pure fragment block-sparse route only`.
- `rt_dg_integrator_taylor.f90` stops for the same condition:
  `DG Taylor4-PC supports the pure fragment block-sparse route only`.
- `rt_dg_integrator_rk.f90` and `rt_dg_integrator_aetrs.f90` also stop for PW,
  so the initial implementation should focus on Taylor4-PC unless those
  integrators are explicitly needed.
- `rt_dg_observables.f90` stops for PW:
  `DG-Fragment RT: mixed/PW observable route was removed`.
- `rt_dg_fragment_ops.f90` contains mixed helpers, but the active derivative and
  Taylor paths do not call them because of the stop gates above.

## Missing Matrix Blocks

- The propagation path needs a single mixed-basis derivative apply:
  `H_FF`, `H_FP`, `H_PF`, and `H_PP` must act on `(coef, coef_pw)` in one frozen-H
  operation for each Taylor term.
- Vector-potential coupling must be consistent across fragment and PW sectors:
  existing fragment velocity/momentum blocks, `P_mat_frag_pw`, and PW diagonal
  momentum/kinetic terms need to enter the same convention used by
  `build_mixed_hamiltonian`.
- If the fragment basis is not strictly orthonormal in the active route, the PW
  path needs the same overlap/generalized-basis treatment as the fragment-only
  propagation. This is especially important for f-sum and paramagnetic response.
- The observable route must include PW contributions to `J_para` and any
  nonlocal/current correction consistently with the Hamiltonian used in
  propagation. Otherwise `J_tot,z` can remain biased around the diamagnetic
  current.
- The canonical mixed-basis fields `mixed_basis_ready`, `mixed_transform`, and
  `coef_mix` are not currently a reliable runtime propagation representation;
  the practical first path should propagate raw fragment/PW coefficients and
  use canonical mixed fields only when intentionally enabled.

## Validation Inputs

- Add a short PW smoke RT input from the diamond64 DC-flux case with:
  `yn_plane_wave_basis='y'`, `n_plane_waves_dg=8`,
  `k_cutoff_plane_wave=0.1`, `yn_dg_subspace_diag='y'`,
  `dg_subspace_extra_states=8`, and `dg_subspace_pw_vectors=4`.
- The RED expectation is failure at the Taylor/derivative/observable PW stop
  gates. The GREEN expectation is reaching at least `itt=2` with no PW stop.
- After the smoke run, compare impulse `J_tot,z`. The target physics is initial
  kick current, rapid paramagnetic cancellation toward zero, then damped
  oscillation around zero rather than oscillation around `J_dia`.
