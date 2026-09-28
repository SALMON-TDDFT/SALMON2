# PBE0 spectrum implementation plan

Approved design: add xc='pbe0' with 25% unscreened Fock exchange, 75% PBE exchange and 100% PBE correlation. No rVV10. Reuse the PBEh spatial/DC/MLWF/source-ACE routes and Coulomb cutoff. Existing functional names retain their definitions.

Goal: compare PBE0 with PBE, PBEh(40), HSE06 for the existing 32 H2 system.
Architecture: recognize PBE0 in hybrid gates, use shared exchange_fraction() for both Fock and semilocal exchange, preserve distinct LCFO metadata and reject incompatible seeds.
Tech stack: Fortran, MPI, Libxc, FFTW, ScaLAPACK, Python.

1. Add compiled tests for fraction/screening, cutoff and LCFO metadata mismatch; verify failure before implementation.
2. Add PBE0 to the existing hybrid gates; pass shared fraction into semilocal evaluation. Test 75% PBE exchange independently through Libxc.
3. Build in a separate directory, run functional tests and small GS/RT checks including source ACE and incompatible GS rejection. Keep ongoing HSE binary untouched.
4. Prepare PBE0 DC GS (PBE warm-up, no fragment localization), short RT probe, 7000-step impulse and zero runs. Same cutoff 4 bohr, .999 source support, dt .05 au and MPI4 as PBEh. At most three concurrent MPI4 runs.
5. Record provenance and timing, extend common-window spectral analysis to four functionals. Timing includes changing contention; do not infer precise speedup from total wall time alone.

## Execution record

- Definition test failed on the old implementation with `PBE0 screening`, then passed after the changes.
- Three compiled definition/Libxc semilocal checks passed; the full MPI GS/RT source-ACE suite (PBE0, HSE06, PBEh40 and fallback) passed 8 tests.
- The legacy LCFO reader uses Fortran STOP, which may return status zero on rejected input. The integration test checks the explicit metadata/preflight rejection and absence of reconstruction/completion instead of relying on the return code.
- Separate Release MPI/HSE/Libxc/ScaLAPACK build completed in `work/pbe0-spectrum-build`; reuses the same Libxc 5.2.3 installation as the ongoing spectra.
- Independent static review found no important issues.
- PBE0 GS/probe/RT scheduling started in `work/h2-dielectric-4x1x1/run-pbe0-queued.py`; calculation results are pending. Existing HSE jobs are unchanged.
