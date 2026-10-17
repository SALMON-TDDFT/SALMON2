# HSE feature selection

All runtime HSE/ACE/SR feature switches are in `&functional`. Environment variables
no longer override the input. Existing inputs without these keys retain disabled
experimental routes. SR wave batches are fixed at 16 (the separate existing GPU
batch input is unchanged).

| Input | Default | Values / meaning | Replaces |
|---|---|---|---|
| `exx_factor_mode` | `off` | `off`, `strict` (exact endpoint history), `learned` (history prediction) | `SALMON_FACTOR_HISTORY` |
| `yn_exx_fft_orbits` | `n` | symmetry-equivalent dense-k hybrid exchange FFT reuse | `SALMON_EXX_FFT_ORBITS` |
| `yn_ace_shadow_diagnostic` | `n` | shadow comparisons within surrogate trace | `SALMON_ACE_SHADOW_DIAGNOSTIC` |
| `yn_ace_export_factors` | `n` | raw one-orbital H2 factor diagnostic | `SALMON_ACE_EXPORT_FACTORS` |
| `hse_sr_group_fft` | `off` | `off`, `direct` (group FFT without halo exchange), `halo` | absent / `0` / `1` in `SALMON_HSE_SR_GROUP_FFT` |
| `hse_sr_radius` | `0` | positive group-FFT cutoff radius in the selected input length unit | `SALMON_HSE_SR_RADIUS_AU` (previously always bohr) |
| `yn_hse_sr_wave_first` | `n` | wave-first SR group route | `SALMON_HSE_SR_WAVE_FIRST` |
| `yn_hse_sr_return_sum` | `y` | aggregated SR returns | `SALMON_HSE_SR_RETURN_SUM` |
| `yn_hse_sr_pair_symmetry` | `n` | conjugate-pair reuse for identical source/target sets | `SALMON_HSE_SR_PAIR_SYMMETRY` |

`SALMON_HSE_SR_WAVE_BATCH` is removed; batch16 is fixed. All character selections
are case insensitive, broadcast to all MPI ranks, and recorded in `variables.log`.
The radius follows `unit_system` and is converted to bohr internally. Existing
`hse_sr_tolerance` still controls finite-range convolution; it is not replaced.

Group FFT requires an MPI build, positive radius and SR tolerance. Wave-first
requires group FFT; pair symmetry additionally requires aggregated return.
Dense-k FFT reuse supports HSE06, PBE0, PBEh40 and PBEh40+rVV10 with native CPU
exchange, crystal symmetry, and all q points. Global hybrids require fixed ions,
integer occupations, occupied-support ACE, no WF truncation and k-only MPI.
Use `yn_hse_wannier='n'` (default) on both reference and reuse routes.
Cubic real-space/k meshes and an orthogonal cell are still required. The retained
group must preserve the RT field. Learned prediction remains HSE06-only.
Global hybrids use the existing scalar `pbeh_coulomb_radius` convention (zero
selects the existing automatic radius); the q=0 value is finite. No SR
localization is added to the global-hybrid kernel.
Learned history retains the previous static HSE Taylor4 predictor-corrector,
occupation and layout restrictions. Checkpoint mode IDs remain 0/1/2 for
`off`/`strict`/`learned`, so select the same mode when restarting.

Example SR configuration (atomic units):

```fortran
&functional
 xc='hse06'
 exx_ace_support='source'
 exx_mlwf_norm_fraction=1d0
 hse_sr_tolerance=1d-1
 hse_sr_group_fft='halo'
 hse_sr_radius=10d0
 yn_hse_sr_wave_first='y'
 yn_hse_sr_pair_symmetry='y'
/
```

For a dense-k learned calculation, use `yn_exx_fft_orbits='y'` and
`exx_factor_mode='learned'` with the previously validated k-only MPI, symmetry,
occupied-support and predictor-corrector settings. `exx_factor_n_teach` defaults
to 8, and warmup remains derived automatically.

DC-LCFO reconstruction caches fragment payloads temporarily before real-space RT.
Its memory requirement is retained internally. There is no fixed cache budget,
threshold warning or capacity input; `SALMON_LCFO_CACHE_MIB` is removed. Actual
allocation or I/O failures are still reported.

## PBE initialization of hybrid GS

For a supported fresh static hybrid GS, omitting `exx_pre_scf_threshold` selects
`1d-4`. Three consecutive selected-density residuals below it trigger the target
hybrid stage. This is an initialization threshold, not the final accuracy;
`threshold` still controls final hybrid convergence. Explicit
`exx_pre_scf_threshold=0` requests direct hybrid SCF. The automatic default is
zero for RT/restart, nonhybrids, moving ions, unsupported mixing/metrics,
checkpoint/shutdown or diagnostic modes and DC without positive temperature.
The input log records the resolved value. See [PBE pre-SCF](exx-mlwf.md).
