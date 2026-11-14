# ASE warm start, force scope and distribution migration

Approved order: warm-start and Libxc-PZ/PBE/NC-UPF force integration first; isolated,
non-collinear and Ehrenfest remain subsequent native-interface work, not removed
from the first-release target. User budget is 4 GiB / 4 CPU cores.

## Design

- Keep the existing fresh-SCF path unchanged. `DFT.run(initial_guess=completed_gs)`
  explicitly reuses a verified wavefunction/occupation snapshot. Require identical
  cell, grid, kmesh, ordered species, potential bytes/local channels, electron and
  orbital counts, XC and executable hash. Positions may differ. Do not reuse
  mixing/SCF history or old coordinates. Reset iteration count, recalculate density
  from wavefunctions, and run SCF at new input positions. Keep source immutable.
- `SalmonCalculator(warm_start=True)` opts into reuse from its last successful
  compatible result. Default False preserves the previous numerical path. Model
  incompatibility causes a fresh run; corrupted restart data fails explicitly.
  Each result still has a separate directory and records initial-guess provenance.
- Extend Calculator's allowed XC set to PZ, Libxc-PZ and Libxc-PBE, and scalar NC-UPF
  in addition to KY, only after native comparisons. Do not advertise stress or
  cell optimization. Reuse upstream native forces; no numerical finite-difference
  replacement or silent XC fallback.
- Preserve source history and record migration commits in NOTICE/README. Keep
  runtime/source.json pointed at the pinned upstream source. Verify independent
  repo CI, publish its tested runtime ZIPs as a private prerelease, then replace
  user-facing acquisition links and pin exact artifact bytes. Historical records
  retain their original sources and dates. Private Release access needs permission.

## Verification steps

1. Add rejection tests for incompatible initial guesses, native input flags,
   source immutability, Calculator cache/failure/reuse behavior, XC scope.
2. Run real displaced-Si fresh and warm SCF; compare energy/forces, iteration
   counts and final positions. Audit the Libxc executable and use known KY/PBE
   files with their existing provenance. One calculation/thread at a time.
3. Run source Python tests; build sdist→wheel and installed tests after docs/layout
   changes. Synchronize Japanese, English and AI JSON guidance and samples.
4. Observe all runtime CI jobs, download/upload exact tested artifacts to Release,
   verify asset hashes and run the new local runtime's Si example before changing
   acquisition instructions. Do not replace URLs with unverified guesses.
