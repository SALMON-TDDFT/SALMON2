# Diamond C64 Response Inputs

## Goal

Prepare a shared test set for linear-response RT-TDDFT on diamond using:

- `full_rt_impulse`
- `dg_gs`
- `dg_rt_impulse`
- `dgpw_rt_impulse`

The crystal model is a `2x2x2` supercell of the diamond conventional cell
(64 carbon atoms total), fragmented as `2x2x2` for the DG/DC route.

## Agreed Comparison Setup

- Use `Gamma`-only sampling for this first comparison.
- Use `PZ` LDA and the existing norm-conserving `C_rps.dat`.
- Use one shared `dg_gs` run as the parent for both `dg_rt` and `dgpw_rt`.
- Use `theory='tddft_response'` and `ae_shape1='impulse'`.
- Use `&analysis` so dielectric-function outputs are produced together with RT current.

## Geometry

- Conventional diamond cell in this repository already uses:
  - `al = 6.72 a.u.`
  - 8 carbon atoms at the usual diamond reduced coordinates
- The new supercell doubles the cell in each direction:
  - `al = 13.44 a.u.`
  - 64 carbon atoms

## Initial Input Choices

These files are meant to be practical starting points rather than fully tuned
production inputs.

- `full_rt`
  - `nstate = 256`
  - `num_rgrid = 24,24,24`
  - `num_kgrid = 1,1,1`
- `dg_gs`
  - `yn_dc='y'`
  - `yn_dc_lcfo='y'`
  - `yn_dc_lcfo_diag='y'`
  - `num_fragment = 2,2,2`
  - `nproc_rgrid_tot = 2,2,2`
  - `nstate_frag = 32`
  - `num_rgrid_buffer = 4,4,4`
- `dg_rt`
  - restart from `dg_gs/data_dcdft`
  - same fragment settings as `dg_gs`
- `dgpw_rt`
  - same as `dg_rt`
  - plus `yn_plane_wave_basis='y'`
  - initial mixed-basis settings:
    - `n_plane_waves_dg = 64`
    - `k_cutoff_plane_wave = 0.75`
    - `dg_subspace_pw_vectors = 8`

## Notes

- The `DG` and `DG+PW` RT directories should see `./data_dcdft/` from the
  `dg_gs` run, for example via symlink.
- These inputs prioritize a coherent comparison path over aggressive
  convergence. Grid, state, and PW settings will likely need a second tuning
  pass after smoke runs.
