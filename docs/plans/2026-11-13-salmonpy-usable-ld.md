# SALMONpy standalone LD implementation plan

Goal: provide a usable first Python entry point for an isolated isotropic LD material.
Architecture: Python validates explicit atomic-unit input and owns independent P/J state; a C-interoperable Fortran adapter calls the existing coefficient and update kernels. Batch and one-step interfaces share that adapter. Frequency response is an analytic diagnostic of the same continuous model.
Tech stack: Python, NumPy, ctypes, Fortran iso_c_binding, gfortran.

1. Add failing behavioral tests: vacuum, Drude constant field, independent state, batch/step identity, invalid inputs, dielectric response.
2. Add an isolated native library build with source/object prohibition audits. Do not change application kernels.
3. Implement explicit-state batch adapter, Python material/state/results API, fixed dt and E history.
4. Verify tests, existing LD regression, documentation example and whitespace checks.
5. Document atomic units, Fourier convention, initial-condition/time alignment, limitations and actual commands. Retain A-only crystalline interfaces as future work.

This increment covers LD. FDTD Python execution, geometry, pulse utilities, analysis readers and crystalline backends remain subsequent increments; do not advertise them as available.
