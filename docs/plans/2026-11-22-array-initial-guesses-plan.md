# Array initial guesses implementation plan

**Goal:** GS/2 and custom providers feed independent SCF owners without GS files.
**Architecture:** one callback-based local-array seed interface. Fortran requests
its actual local bounds after static initialization, validates returned values
and sampling, copies into its independent seed and uses the established SCF path.
Python adapters own immutable source metadata and arrays or verified archive reads.
**Tech stack:** Fortran iso_c_binding, ctypes, NumPy, existing GS/2 reader.

1. Add failing source/API tests for the typed array-seed ABI and immutable custom
   seed/archive adapters. Validate shape/metadata failures and source independence.
2. Add owned local callback descriptor and scalar finite/sampling validation in
   gs_seed; route create_from_arrays through existing candidate creation/cleanup.
   No global selectors or startup builds. Audit before compiling.
3. Add Python GSInitialGuess arrays/archive factories and callback packing;
   extend NativeGSSession.run dispatch while retaining live-GS route unchanged.
4. Extend existing MPI fixture to compare array/provider route against retained
   seed, release/mutate original buffers, inject provider failures and retry.
   Use PZ/HSE k and Gamma spatial; compare unchanged cold controls. Do not save
   real GS datasets. Use synthetic archives for disk adapter tests.
5. Update paired guides/AI JSON/sample/checklist; run source/wheel suites and
   tracked-only source build. Refresh NLCC binary evidence for changed source.
6. Record limits honestly; no material-physics CI or optional heavy dispatch.
