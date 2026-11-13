# Band workflow validation — 2026-11-13

Validated on macOS arm64 using GNU Fortran, LP64 OpenBLAS and the separately built
Libxc-enabled SALMON executable. This does not establish native execution support
on Windows, Intel macOS or Linux for this new entry point.

- Existing converged PZ Si GS: 8-atom cubic cell, 12³ real-space grid, 2³ k mesh,
  32 orbitals and 32 electrons. Coarse connection fixture, not a converged band gap.
- First path: Γ–X–M–R–Γ, four intervals/segment, 17 physical points, 24 stored points.
- Second independent path: Γ–X, seven intervals, eight physical/stored points.
- First 24 bands requested; eigenvalue-change threshold 1e-5 Hartree.
- Wrapped/direct execution with exactly the same generated input: eigenvalues and
  Cartesian k coordinates bitwise equal (maximum energy difference 0 Hartree).
- One-iteration run with 1e-15 tolerance: native process exits normally, adapter
  rejects explicit unconverged blocks, writes failed manifest, preserves logs.
- Every original GS file and every owned-snapshot file retains its initial SHA256
  after both successful paths and the failed path.

The native path test compiles the actual band_dft.f90 against minimal type modules
with bounds checking enabled. A skew reciprocal matrix checks basis orientation,
exactly divisible point counts, endpoint padding and absence of a second transform.
This is not a skew-cell full-DFT validation.

The native reporting addition emits reciprocal basis columns and convergence for
each block. It does not alter the orbital solver, density/potential updates or
convergence criterion. `max_iterations` follows the native cumulative Miter budget.

Reproduce in a source checkout:

```
python tools/test_band_path_native.py
python tools/test_band_runner.py SALMON_EXECUTABLE GS_BUNDLE GS_INPUT NEW_OUTPUT_DIRECTORY
```

The integration fixture is the documented 8-atom Si bundle; that script is not a
universal test for arbitrary materials. It needs NumPy and Matplotlib. Select the
active environment's interpreter. Audit the full executable before execution.
See the output validation.json for executable SHA256 and numerical results.

Source tests and repaired, separately installed macOS arm64 CPython 3.14 wheel: 52 tests passed each.
