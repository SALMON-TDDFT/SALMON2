# Native fast-HSE owner validation — 2026-11-15

Local Apple Silicon / GNU Fortran / OpenBLAS / FFTW / Libxc validation. No GitHub Actions were started. The user requested Actions-budget conservation.

## Change and scope

The explicit-routing profile gives a bound native system its own fast-PR FFT kernel, ACE, trial state and histories. Hamiltonian, XC/total energy, Taylor, PT-CN and restart paths resolve the system owner without a process-global active-pointer switch. The optimized pair dispatch is retained; the standalone generic resource primitive is not substituted. Existing calls retain their default owner.

This is serial, periodic, non-spin-polarized, full-support HSE under common immutable process settings. Different cell spacings with common grid/k counts are tested. MPI/Wannier/symmetry/spin paths are rejected for bound owners. Full TDDFT-solver interleaving, arbitrary per-state global configurations, physical convergence, speedups and a new Python/C API are **not** established. Bound systems/owners must not be copied, called concurrently or outlive their owner.

## Results

- New explicit-owner TDDFT driver: 32 cases, 2 OMP threads, k counts 2/3/6/8, dt 0.001/0.002, exact/learned Taylor, direct Taylor and learned PT-CN. Current and energy output files are identical to the immutable fast reference; all maximum differences are zero. Factor/prediction/solver diagnostics match. An activation marker is required, so a fallback to the default owner cannot silently pass.
- Existing default-owner driver: the same 32 cases at 1 OMP thread, identical outputs/diagnostics versus the immutable reference.
- Native HSE owners: 1 and 2 OMP threads; different-cell ACE and direct actions match separate-process references exactly. Interleaved refresh, energy, release/reopen, double binding rejection and unsupported-call output retention pass.
- A deliberately shared FFT kernel in an external mutation build is detected by the separate-process B comparison. The repository and validated runtime are not mutated.
- Existing explicit RT-state routing and standalone resource primitive probes pass at 2 OMP threads.
- Source tests: 305 pass. Installed wheel tests: 305 run, 1 optional plotting skip, no failure.
- sdist-to-wheel succeeds; installed-wheel source preparation matches every updated completed-file hash.
- Source/all-application-object forbidden-eigensolver audits pass before linking.
- Resource-enabled byte accounting reports unavailable (status 2, bytes -1), rather than omitting opaque FFT-plan memory from a purported total.

## Identities

- Reference executable SHA256: `74cd086c228028bca069351fb175d886fc0979138fdf5b0b80c0d98a973d8df7`
- Changed executable SHA256: `7af9fd0a542d122e8ba2aed32956a82a4566d2b4ab2cfa9f04eadf7739af613f`
- Completed-file manifest SHA256: `9cbe653ebc24d55c34d777230267b35d564a1d259a5b79fc5eb10207a27d068f`
- Explicit regression result SHA256: `2abbe960075a0f17f8eeb01af9d3a69aaeafa043b82326a488c4f68381232044`
- Default regression result SHA256: `01473c576e14810f9731d7dbc008f17445c81db31e300a17ab8fb78b73336df3`

Local evidence is under `/tmp/salmonpy-native-owners-final-*20261115` and `/tmp/salmonpy-native-owners-shared-mutation-20261115`. Development failures were kept separately; they are not approved runtimes. No real GS payload was collected. No commit/push or remote CI was performed for this change.
