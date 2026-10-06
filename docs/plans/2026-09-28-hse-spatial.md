# HSE spatial exchange implementation plan

> 2026-10-06：本記録中の`developer_tests`は当時のローカル開発検証です。GitHub配布から除外しました。通常の回帰試験は`testsuites`を使用します。

**Goal:** Extend the existing Gamma spatial MLWF/FFTW/ACE mesh route to HSE fixed-ion RT, keeping legacy implementations until all their supported modes have replacements.

**Architecture:** Share localization, distributed FFT and ACE; select the screened reciprocal HSE kernel with its analytic zero mode, or the existing truncated Coulomb kernel. DC provides initial mesh wavefunctions only. HSE MD remains gated pending force validation.

**Technology:** Fortran, MPI, FFTW, Python unittest.

1. Extend `developer_tests/651_hybrid_exchange/ace/exchange_driver.f90` and `validate_exchange.py` to compare screened exchange against serial Wannier for omega 0, .11 and .3 on 1/2/4 ranks. Observe missing interface, implement optional omega in `src/xc/hse_spatial.f90`, verify parity.
2. Extend `src/xc/hse_native.f90` dispatch and both action calls for HSE DC-initialized RT. Update `src/io/inputoutput.f90` admission and shared RT timing in `src/rt/{initialization_rt,time_evolution_step,em_field}.f90`. Preserve PBEh-only force and MD permissions.
3. Add HSE DC preparation and fixed-ion impulse/pulse serial vs 2/4-rank integration regression. Run existing PBEh suite, HSE-off build and independent review.
4. Document verified scope and remaining migration work. Do not remove legacy multi-k, finite-support, SCF or projected functionality before replacement validation.

## Validation outcome

- Screened/Coulomb exchange: 1/2/4 ranks, omega 0/.11/.3, two refresh stages; maximum absolute action difference 5.064e-15. Negative screening rejected.
- HSE DC-initialized impulse/pulse: serial vs 2/4 ranks, maximum energy difference 1.235e-13 Ha and RT observable difference 2.874e-14 in output units.
- Three HSE integration tests pass, including retained HSE MD rejection and legacy full-propagator input admission (the noncubic fixture reaches its pre-existing cubic-grid guard; it is not a legacy propagation completion test).
- Existing PBEh Ehrenfest suite: 16 tests pass. USE_HSE ON and OFF builds pass.
- Independent review identified overbroad automatic Wannier activation; restricted to spatial or explicitly selected Wannier HSE. Reviewer confirmed resolution.
- Remaining migration: SCF, multi-k, finite support and HSE forces; legacy routes intentionally retained.
