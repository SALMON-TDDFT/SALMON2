# HSE feature selection namelist implementation plan

Goal: replace all runtime environment feature selectors in the integrated HSE/DC path with functional namelist inputs.

Approved architecture: keep existing defaults; SR batch is fixed at 16. LCFO cache allocation has no fixed budget or abort threshold. Retain internal byte estimates without routine output. Use off/strict/learned history modes (strict collects exact endpoints, learned predicts). SR group route uses off/direct/halo to preserve the previous absent/0/1 distinction. No environment overrides.

1. Add a small input regression probe and observe rejection of the new keys before implementation.
2. Declare new input variables in salmon_global; add functional namelist defaults, broadcasts, lowercase normalization, dimensional radius conversion, variables.log and checks.
3. Replace environment reads in exx_native, exx_spatial and restart readers with these variables. Preserve checkpoint mode numbering and numerical restrictions. Fix SR batch at16.
4. Remove LCFO cache budget input and abort; keep required-byte estimate internal.
5. Document the input mapping. Validate MPI HSE, serial HSE/non-HSE builds, parser/default/invalid-input checks, exchange and pair oracles and existing regression tests. Keep restart probes outside testsuites.

New switches: exx_factor_mode, yn_exx_fft_orbits, yn_ace_shadow_diagnostic, yn_ace_export_factors, hse_sr_group_fft, yn_hse_sr_wave_first, yn_hse_sr_return_sum, yn_hse_sr_pair_symmetry, hse_sr_radius (input length unit).
