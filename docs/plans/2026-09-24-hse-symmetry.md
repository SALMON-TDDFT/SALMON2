# HSE06 and fixed-alpha TDCDFT symmetry implementation plan

Goal: support the Si z-field spatial symmetry subgroup and compare against full-k results, preserving ongoing impulse calculations.
Architecture: propagate/store representative k orbitals and ACE factors. At HSE refresh reconstruct each locally owned full-k star, including fractional translation/Bloch phases and little-group averaging of source projectors, call the existing distributed exchange kernel, then retain representative actions. MPI owners retain their own stars; no replicated global orbital gather. Full-k FFT work remains and speed/memory reduction must be measured, not assumed to be eightfold. Existing symmetry density/current machinery is reused. Runtime field/group compatibility guards are required.
Tech stack: Fortran, MPI, BLAS, FFTW; compiled independent symmetry tests and short native parity tests.

User explicitly authorized symmetry support. Existing TDCDFT/HSE long jobs continue using their running executable; build new code in a separate build directory or preserve executable before rebuild.

1. Independent symmetry map + Bloch reconstruction tests including little groups, translation, reciprocal folding, grid validation.
2. HSE adapter integrates full stars and weighted exchange energy; retain old path when symmetry off.
3. Fixed-alpha TDCDFT enables only compatible external fields and group operations, no evolving screening initially.
4. Generate diamond z-preserving sym.dat; count actual stars including boundary multiplicities.
5. Build, run independent MPI exchange regression; native GS/impulse comparison full vs reduced, restart and unsupported-input guards.
6. Document actual k count, limitations, time/memory, then commit/push validated result.

Implementation update: closed conventional-cell affine group has32 operations
(8 rotations times4 FCC centerings). Use density callback to stream each
little-group source projector; never allocate source orbitals times max_little.
Exact half-shifted counts64->12 and4096->576. Native four-step current parity
HSE1.60e-8 / TDCDFT3.28e-7; full MPI4 vs reduced MPI3 tests uneven stars.
Independent review found missing generic Cartesian isometry guard; fixed with
closed shear-reflection rejection regression. Restart HSE/TDCDFT2+2 vs4 is bitwise equal. Native test suite7 pass+1 opt-in skip; real MPI test passes separately;5 symmetry/metadata guard tests pass.
