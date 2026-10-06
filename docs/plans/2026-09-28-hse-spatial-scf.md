# HSE spatial SCF implementation plan

**Goal:** Continue the approved common spatial exchange migration with conventional fixed-ion HSE Gamma SCF.

**Architecture:** Reuse spatial MLWF, screened FFTW and reduced ACE operations already used by RT. Admit only ordinary DFT, fixed occupied spin pairs, orthogonal Gamma y/z pencils, full support. Keep fractional DC, multi-k and finite-support legacy implementations until replaced. No arbitrary memory limits.

**Technology:** Fortran, MPI, FFTW, Python unittest.

1. Add serial/2/4-rank SCF convergence and final-energy comparison in `testsuites/653_functional/test_hse_spatial.py`; run to establish current spatial rejection.
2. Extend `src/io/inputoutput.f90` admission and `src/xc/hse_native.f90` spatial dispatch for HSE DFT, restrict Taylor4 requirement to RT, reject unsupported occupations/restarts/snapshots early.
3. Compare converged energy and eigenvalues, exercise invalid inputs, run existing HSE RT tests and ON/OFF builds. Request independent code review.
4. Document scope and results. Preserve old routes required by outstanding modes and the serial reference.

## Results

- Red test confirmed existing spatial SCF rejection at the k-only Wannier guard.
- H4 density residual below 1e-10 on 1/2/4 ranks; final energy max difference 1.5664e-9 eV and occupied eigenvalues 2.0855e-11 Ha.
- SCF iterations 109/49/476: final-state agreement only, not equivalent convergence histories or a speedup claim.
- Final suite: five tests pass, including ten invalid-input subcases and prior HSE impulse/pulse RT comparisons.
- USE_HSE ON and OFF builds pass; diff whitespace check passes.
- Independent review found omitted final GS checkpoint writes; guarded both final-output paths plus full-grid diagnostic exporters, added tests, reviewer confirmed closure.
- Remaining: fractional-occupation DC SCF, multi-k, finite support and HSE forces; legacy routes retained.
