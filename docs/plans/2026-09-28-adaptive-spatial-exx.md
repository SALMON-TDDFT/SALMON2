# Adaptive spatial EXX implementation plan

**Goal:** Reuse localized WF support in distributed SCF, with per-source radii retaining 99.9% norm and a compact convolution path.

**Architecture:** Keep exact Hermitian exchange for the masked sources Q: K=-sum Q v Q*. Do not mask arbitrary targets or silently discard small pair densities. Norm retention is user-selectable; zero disables adaptive truncation. Obtain periodic centers and radii using local grid rows and scalar reductions. Apply masks only after converged localization; full support otherwise. A source-only mask avoids a non-Hermitian ACE metric. Compact FFT evaluation must equal the global discrete convolution for the SAME masked sources, including G=0. DC uses fragment communicators and cell geometry.

**Tech Stack:** Fortran, MPI spatial pencils, FFTW, existing exx_local_fft, Python/Fortran probes.

User approved automatic norm radii and local integration on 2026-09-28. This implementation first preserves the existing source-only approximation: it does not claim linear scaling from heuristic target-WF pair dropping. SCF and fixed-ion native RT use the shared masked-source operator, with separate energy/current checks. Ionic MD remains guarded.

1. Add a distributed norm-support module and tests: periodic boundary centers, 0.999 retained norm, ambiguous-center protection, rank invariance, full-support limit.
2. Add compact spatial convolution adapter preserving the global kernel, with numerical action/Hermiticity tests against full-grid FFT. Never gather full wavefunction grids. Report compact work and fallback.
3. Add exx_mlwf_norm_fraction (0 disabled; 0<value<1 adaptive). Preserve exx_mlwf_radius compatibility. Add input guards and per-refresh radius/loss diagnostics. Localization convergence enables masking; localization failure restores full support. Prevent SCF final acceptance if adaptive support never becomes valid.
4. Exercise HSE/PBEh SCF serial/spatial and DC fragment parity, compare 0.999/0.9999/full support. Document observed errors and restrictions, review, commit.

SCF mixing history reset and convergence invalidation are implemented. Fixed-ion RT is tested over 16 impulse steps. Outstanding follow-on: long-time RT radius continuity and accuracy; sparse target-pair/ACE scaling; full weak/strong scaling remeasurement. No claim of linear scaling or finite-support ionic-force consistency.

## Completion evidence

SCF and fixed-ion RT adapters implemented and verified; see ../adaptive-exx-ja.md for actual numerical differences and scope. Initial spatial identity-gauge saddle discovered by the micro-impulse test was fixed with projected-periodic-position seeding and a dedicated MPI regression. Integration tests passed: 4 adaptive SCF methods, 3 functional RT tests; 14 EXX input tests plus final adaptive guard rerun; distributed mask/compact/seed probes; 27 combined spatial/orbital exchange probes. Independent read-only reviews found no blocking collective, cache, ACE or switching issues. Performance scaling remeasurement and long-time response accuracy remain separate work.
