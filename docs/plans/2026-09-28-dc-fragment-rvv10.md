# Fragment-local DC rVV10 implementation plan

**Goal:** Evaluate rVV10 independently in each buffered periodic DC fragment.

**Architecture:** Reuse exchange_correlation with fragment density, mesh and spatial communicator. Its energy-density output already enters the DC core-only sum. Delete the separate total-density rVV10 path. Hartree and scalar energy/occupation reductions remain unchanged. Native mesh RT after DC reconstruction still uses its full propagation cell.

**Tech Stack:** Fortran, existing distributed FFTW/FFTE rVV10, Python integration tests.

Approved by the user on 2026-09-28. This supersedes the total-density DC rVV10 design. The earlier adaptive-WF localization task remains pending; this correction precedes it.

1. Add a finite-buffer DC test requiring reported fragment FFT dimensions, absence of the old total-density backend, and MPI layout energy parity. Run against the old binary and confirm failure.
2. Enable shared rVV10 evaluation for DC in src/xc/salmon_xc.f90. Remove the redundant global rVV10 block and imports in src/gs/dc/dcdft.f90. Log fragment evaluation dimensions.
3. Build existing MPI/HSE/ScaLAPACK build. Run the new finite-buffer test, one-fragment versus conventional SCF, full-buffer partition equivalence, and existing backend tests.
4. Update docs/inputs/pbeh40-rvv10.md: fragment-periodic approximation, core-only energy, buffer convergence and unvalidated truncated-fragment forces. Record actual results; commit after checks.

Validation must not imply that the derivative of the core-partitioned DC energy equals the full-fragment potential, or certify DC-MD forces. 99.9% WF support is a separate approximation and is not introduced here.

## Validation completed

- Existing MPI/HSE/ScaLAPACK build: successful.
- New finite-buffer fragment-grid test failed on the original binary, passed after the change.
- test_dc.py: 4 tests passed (finite-buffer MPI 2/4, one-fragment conventional parity, full-buffer partition parity, localization rejection).
- test_exx_inputs.py filtered to rvv10: 2 tests passed (FFTW/FFTE and fallback).
- PBEhOrbitalSCF.test_dc_orbital_fractional: passed at MPI 2/8/16, including spatial and orbital partitioning and LCFO spectrum parity. Reported energies -76.17659 eV and charge 4 within 6e-14.
- DCForceTest one_fragment_limit and full_buffer_energy_derivative: 2 tests passed. These do not certify finite-buffer forces.
- Read-only independent review found no communicator, energy accounting, initialization or LCFO propagation blockers.

Total: 9 integration tests passed. No weak-scaling improvement is claimed by this correction; adaptive WF support and actual multi-fragment scaling remain subsequent work.
