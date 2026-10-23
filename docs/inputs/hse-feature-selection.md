# HSE feature selection

HSE/ACE/SR settings are in `&functional`. Environment variables do not override
the input. Native k-exchange optimization normally needs no switches: the
omitted `exx_k_algorithm` selects the conservative `auto` policy. SR wave batches are fixed at 16 (the separate existing GPU
batch input is unchanged).

| Input | Default | Values / meaning | Replaces |
|---|---|---|---|
| `exx_factor_mode` | `off` | `off`, `strict` (exact endpoint history), `learned` (history prediction) | `SALMON_FACTOR_HISTORY` |
| `exx_k_algorithm` | `auto` | `auto`, `reference`, `pair`, `orbits`; native k-exchange strategy | three legacy exchange switches |
| `yn_ace_shadow_diagnostic` | `n` | shadow comparisons within surrogate trace | `SALMON_ACE_SHADOW_DIAGNOSTIC` |
| `yn_ace_export_factors` | `n` | raw one-orbital H2 factor diagnostic | `SALMON_ACE_EXPORT_FACTORS` |
| `hse_sr_group_fft` | `off` | `off`, `direct` (group FFT without halo exchange), `halo` | absent / `0` / `1` in `SALMON_HSE_SR_GROUP_FFT` |
| `hse_sr_radius` | `0` | positive group-FFT cutoff radius in the selected input length unit | `SALMON_HSE_SR_RADIUS_AU` (previously always bohr) |
| `yn_hse_sr_wave_first` | `n` | wave-first SR group route | `SALMON_HSE_SR_WAVE_FIRST` |
| `yn_hse_sr_return_sum` | `y` | aggregated SR returns | `SALMON_HSE_SR_RETURN_SUM` |
| `yn_hse_sr_pair_symmetry` | `y` | exact Gamma Hermitian pairs for canonical DC; conjugate pairs for eligible SR wave training | `SALMON_HSE_SR_PAIR_SYMMETRY` |

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

For a dense-k learned calculation, use the compatibility switch
`yn_exx_fft_orbits='y'` and
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

## Native k-exchange algorithm selection

Normal inputs omit `exx_k_algorithm` and all three legacy switches.

| Selection | Behavior |
|---|---|
| `auto` | Eligible scalar crystal FFT reuse first, otherwise scalar Hermitian pairs; full SOC pairs on supported native CPU full-q HSE06; reference on unsupported automatic layouts. |
| `reference` | Disable these k-exchange reuse routes for comparisons. Physical settings and other backends are unchanged. |
| `pair` | Scalar Hermitian pairs without crystal symmetry; full pairs for all four SOC spin blocks, including SU(2) crystal expansion. |
| `orbits` | Existing scalar crystal-symmetry FFT reuse; existing mesh and symmetry restrictions apply. |

The selector does not change physical kernels, precision, occupations, k points,
SR cutoffs, ACE approximation settings or checkpoint formats. In particular it
does not infer a finite SR cutoff from a requested optimization. SOC pair performance and full-run memory require workload-specific measurement; automatic selection is not a speedup or capacity guarantee. The same rules apply independently to GS and RT. RT still
requires a field-preserving subgroup (`nyy` for x, `yny` for y in this fixture).

Explicit `pair`/`orbits` require an eligible fixed-ion, full-support native CPU,
k-only, full-q DFT/RT configuration. Scalar pair requires dense k points, occupied
spin pairs and crystal symmetry off; SOC pair retains the existing SOC HSE06
restrictions. Orbits do not support SOC. An unavailable explicit request fails
instead of silently falling back. `reference` only controls native k-exchange
reuse; it does not replace a selected Wannier, spatial or GPU backend.

In `auto`, eligible scalar crystal reuse wins over the scalar pair flag. Explicit
`yn_exx_fft_orbits='n'` remains an opt-out for reference comparisons.

For canonical Gamma DC (`yn_exx_dc_mlwf='n'`), Hermitian pairs are on by default,
including finite temperature and orbital MPI. Explicit Wannier snapshot output
retains its existing reference backend. Source weights are `rocc/2`; the
reverse action uses the other orbital's weight, so zero/fractional occupations
need no division or approximation. One convolution updates both directions.
`yn_hse_sr_pair_symmetry='n'` selects the ordered-pair reference for this backend.
This switch does not turn on SR truncation, wave-first FFT, or change occupations.
The Gamma log records local pair counts and explicit scratch bytes; these are not
whole-process RSS or a speedup guarantee. SR wave-pair reuse still requires
identical source/target sets and aggregated return.

Old inputs with `yn_exx_fft_orbits`, `yn_exx_k_pair_symmetry` or
`yn_exx_spinor_pair_fft` remain readable as compatibility overrides of `auto`.
Omission is distinguished internally from explicit `n`. With a manual new
selector, only consistent explicit legacy values are accepted. Invalid or
contradictory combinations fail. Avoid mixing both styles in new inputs.
`variables.log` records the request, resolved flags and policy reason. Native
initialization prints `EXX_K_ALGORITHM requested=... selected=...` and the reason.
No environment selector is read.

Example comparison setting (omit for normal operation):

```fortran
&functional
 xc='hse06'
 exx_k_algorithm='pair'
/
```

The `exx_factor_mode`, finite-range SR controls and diagnostics remain separate
because they choose approximation/history behavior or explicit diagnostics, not
an interchangeable native-k implementation. Historical learned/orbits inputs
continue through the compatibility switch; the new manual selectors do not
broaden the certified learned-route combinations.

## Additional SOC exchange MPI workers

For static, fixed-ion, full-support SOC HSE06 with ordinary ACE, CPU/full-q
and an eligible orthogonal crystal-symmetry k layout, `nproc_k` is the number
of SCF owners. Launching more MPI ranks assigns the remaining ranks to bounded
r/rprime exchange tasks; they hold no SCF wavefunctions or ACE factors.
The owners retain their k, SCF, and ACE communicator. Native empty-k guards
are not removed. This is an exchange-only service, not a full empty-k SCF layout.

`exx_worker_tile_divisor=16` and `exx_worker_rprime_groups=1` are conservative
namelist defaults; both can be specified in `&functional`. The latter must divide
the total MPI rank count. Existing saved candidate environment variables are not
read by this integrated implementation. Unsupported helper layouts are rejected
collectively before service entry. No-helper inputs retain the normal path.

This integration preserves the existing numerical kernel. Historical local
helper comparisons passed, while the earlier Fujitsu build/layout-dependent
SCF discrepancy remains documented and does not establish production acceptance
for a newly built executable. Existing jobs/binaries are not replaced.
