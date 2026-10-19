# Spinor pair FFT implementation plan

**Goal:** CPU pair tiles for the existing fixed-ion native HSE06 SOC NLCC rectangular-cell route; preserve all four spin blocks and anisotropic uniform k meshes.

**Architecture:** Keep existing SU(2) star/little-group source expansion. Apply pair tiling on the full expanded source, independently for each cross-spin density block, then extract representative actions as before. This first implementation does not compress real-space point orbits. No physics, occupations, GS/restart layout, learned ACE or resource changes. Default route unchanged.

**Tech stack:** Fortran, MPI Alltoall callbacks, FFTW, BLAS; GNU local and Fujitsu Fugaku builds.

## Tasks
1. Add numerical pair/full comparison probes before production changes. Use rectangular real/k grids, shifted meshes, complex offdiagonal spin density, both FFT layouts, uneven/empty MPI ownership and partial pair passes. Assert mode entry and collective invalid rejection; watch missing API fail.
2. Add a bounded full-pair backend to exx_k_exchange, with collective validation, overflow guards, two FFT layouts, exact existing kernel/phase/normalization, and no source approximation. Reuse existing spinor density callbacks; test differences <1e-11.
3. Add explicit optional pair mode to spinor dispatch and environment selection in native; skip unused row FFT workspace only for supported CPU full-q SOC pair mode. Test SU(2), symmetry GS/RT, ACE/direct comparisons, legacy regressions.
4. Review source diff and results. Make independent source/build bundle; preserve all production binaries and current job inputs. Check live queued/unstarted state before any job-specific wrapper changes, no cancellation/requeue. If no eligible jobs remain do not add production jobs. Record binary/input provenance and all material limits.

## Validation commands
`RUN_MPI_TESTS=1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 python3 -m unittest test_native_pair_spinor test_native_spin_nlcc test_native_exchange test_native_symmetry test_spinor_symmetry_field test_spinor_density_scope`
Build `cmake --build build-local -j4`. Use existing validate_symmetry.py short fixture with pair environment, compare identical saved states. No new WSe2 GS jobs for validation.

## Progress
- Design selected: full-pair spinor backend, no spatial orbit compression; protects SU(2) correctness and keeps extension bounded.

- Implementation committed as 037c301. Missing pair API failed before implementation; pair/full numerical tests passed after implementation. Added 8-rank empty-ownership, partial-pass, mixed route, NaN and oversized-tile validation. 17 regressions passed on final source.
- Independent review: narrowing int64 tile count before bounds check must be fixed; corrected before conversion. No other must-fix findings.
- Final local binary f38bd8b3…: all 3 GS residuals <1e-8, x/y 20-step symmetry and ACE/direct comparisons PASS, same saved GS pair/legacy differences 1.68e-16 and 1.90e-15.
- Fugaku separate build source 037c301, exit 0, AArch64 SHA eb0329041a9ddd121c9004fe6f8a1aaf53bb6097febfdf06e69fa8bddb303998. Original WSe2 and Si binaries unchanged.
- Automatic approval review initially rejected queued-job wrapper replacement for insufficient explicit authorization. User explicitly approved all 3 IDs in this chat; replacement then completed after fresh queued/unstarted checks. Inputs/resources/job IDs unchanged; no cancellation/requeue/new compute jobs. Queue and byte/hash evidence downloaded and verified.
- Limits: full-pair backend retains SU(2) source expansion but does not compress spatial point orbits. WSe2 physical convergence, compute-node runtime/RSS and speedup remain unverified while 3 jobs are QUE.
