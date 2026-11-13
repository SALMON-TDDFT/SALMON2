# Independent ASE band bridge validation — 2026-11-13

Environment: macOS arm64, CPython 3.14, ASE 3.29.0, f90nml, GNU Fortran and LP64
OpenBLAS. External SALMON executable is separately built; it is not in the wheel.

- Existing 8-atom cubic PZ Si GS (12³ grid, 2³ SCF k mesh, 32 orbitals/electrons).
- ASE path GX,MR, 17 exact physical points, 24 bands. Native endpoint padding is
  checked and removed. Disconnected path string and special points survive ASE JSON.
- Independent direct SALMON invocation of the emitted input and k list: bitwise
  identical eigenvalues and Cartesian k (maximum energy difference 0 Hartree).
- Caller Atoms, caller BandPath, original GS bundle and internal snapshot unchanged.
- ASE JSON result read back with identical energies; ASE figure rendered headlessly.
- Existing segmented 17/8-point native workflow regression still gives zero difference
  vs direct execution and rejects deliberately unconverged output.
- Native bounds-checking harness covers skew basis, exact batch division, endpoint
  padding, explicit irregular list and absence of a second reciprocal transformation.
- Python tests cover Cartesian Å and reduced bohr input, skew geometry, atomic order/
  position/cell/PBC mismatch, spin/k/band axis conversion, copied mutable ASE objects,
  and consistent cell/k tolerance with ASE's slightly different Bohr conversion.

The optional module imports no ASE or f90nml until an API needing it is called.
The base dependency remains NumPy; the optional `ase` extra declares ASE/f90nml.
Without those installed, ordinary test coverage runs and four ASE-specific tests skip.

No Calculator, energy/forces/stress API, SCF generation or MPI launcher is added.
Native real-GS integration is verified on this Mac, not inferred for other OSes.
Skew-cell tests exercise geometry/path conversion, not a full skew-cell DFT run.

Reproduction (source checkout and active interpreter):

```
python tools/test_band_path_native.py
python tools/test_ase_bands.py SALMON_EXECUTABLE GS_BUNDLE GS_INPUT NEW_DIRECTORY
```

The real-GS script uses the documented Si fixture. Audit the full executable before
execution. Exact binary identity is recorded in the integration validation.json.

57 tests passed with ASE in the source checkout and in the repaired, separately installed Mac wheel.
Without ASE: 53 passed, four ASE-specific tests skipped.
