# Native SBE and reusable initial state

**Goal:** connect the existing velocity-gauge SALMON solver to the A→J protocol, with independent trajectories.
**Architecture:** share types and unchanged propagation/local-current kernels between SALMON and a serial C ABI. Keep MPI reduction in SALMON. Python owns immutable initial arrays and separate density matrices; steps commit only after successful finite results.
**Tech stack:** Fortran, BLAS, ctypes, NumPy, CMake.

1. Add failing tests for native initial-state construction, independent runs, midpoint propagation and endpoint current, validation and text loading.
2. Extract `sbe_types.f90` and `sbe_velocity_kernel.f90`; delegate existing solver routines to them. Retain correction orders 0–3 and nonlocal velocity correction.
3. Add `sbe_c_api.f90` and `sbe.py`. Accept atomic-unit SALMON text exports, preserve the existing closed-shell occupation rule, expose only A/c driving.
4. Link BLAS in wheel build, retain audits. Compare native results with a frozen pre-extraction Fortran oracle and rebuild the full executable.
5. Document fixed initial/run parameters and per-call A. Run Python and packaging tests. TDDFT and length gauge remain separate follow-up work.
