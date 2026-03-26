# DG PW Energy Cutoff Plan

1. Add input state
   - Add `pw_energy_cutoff` to `salmon_global`
   - Add it to `&propagation` in `inputoutput.f90`
   - Set default to a disabled sentinel value

2. Add unit conversion and logging
   - Broadcast `pw_energy_cutoff`
   - Convert it with `uenergy_to_au`
   - Write both cutoff inputs to `variables.log`

3. Update plane-wave initialization
   - In `rt_dg_plane_wave.f90`, resolve the effective cutoff:
     - use `pw_energy_cutoff` when positive
     - otherwise fall back to legacy `k_cutoff_plane_wave`
   - Compute internal `dg_frag%k_cutoff_pw`
   - Log both energy and effective k values

4. Verify
   - Build
   - Smoke check one legacy input and one new-input case
   - Confirm logs show the expected precedence and converted values
