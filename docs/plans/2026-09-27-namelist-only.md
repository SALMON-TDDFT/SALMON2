# Namelist-only HSE and LCFO controls

User rule: all added algorithm/diagnostic controls must be namelist parameters.
Remove SALMON_LCFO_* and SALMON_HSE_* environment reads, without fallback.
Use existing yn_hse_wannier for both native and LCFO MLWF; select LCFO via
yn_hse_lcfo_rt. Input validation must distinguish native-only MPI constraints.
Move direct-WF, radius, ACE/U cadence, FFT settings, seed QR and diagnostics to
&functional. Preserve defaults except radius -1 legacy sentinel becomes 0/full.
Low-level reusable modules accept typed optional arguments rather than importing
application configuration. Update active 3D inputsets; historical benchmark records
remain historical. Unknown old environment values must have no effect.

Verification: real namelist-driven MPI GS/RT, fragment x orbital MPI, ACE/U reuse,
full/radius warning, deliberately conflicting legacy environment; distributed seed
unit tests through typed arguments; HSE standalone tests as affected; full build.
Provide archive-source patch with hashes and current Si archive. No large jobs.
