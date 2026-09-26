# Coefficient Hamiltonian handoff implementation plan

**Goal:** Remove redundant grid/coefficient conversions from the approved direct-WF Taylor4 path without adding a truncation approximation.

**Architecture:** Pass optional input/output LCFO coefficient arrays explicitly through hpsi and the HSE adapter. Native kinetic/local/nonlocal actions still use the grid. The LCFO adapter projects that non-exchange result once, adds distributed coefficient ACE (or halo exchange), and returns coefficients without reconstructing an unused projected grid. Standard callers remain unchanged. Taylor accumulates the four powers directly; no implicit solver or new propagation radius.

**Tech stack:** Fortran, existing BLAS and MPI communication wrappers.

Compared approaches: preassemble the entire projected Hamiltonian (larger assembly/storage change); introduce sparse propagation immediately (new truncation and gauge errors); first fuse the existing coefficient handoff (chosen exact optimization, approved continuation of the previous design). Existing development checkout is reused.

1. Extend unequal-partition complex ACE probe to compare coefficient action against dense ACE for several factor ranks; observe missing API failure.
2. Implement additive coefficient ACE kernel with dimension checks and existing spatial communicator.
3. Add optional paired coefficient arguments through hpsi/HSE adapters, reject incompatible calls, preserve ordinary path. Update Taylor to consume coefficient results.
4. Build with GNU15 vectorization disabled. Run coefficient probe and native direct-WF regressions including MPI2/4, finite radius, ACE/U intervals and half dt, sequentially.
5. Compare old/new direct Taylor binaries on C64/MPI8 and C128/MPI16, core16³/dt.02/R6/nt16. One numerical job at a time, no concurrent builds. Record physical errors and timing together, disclose single samples.
6. Review diff, document limitations, update notebook, commit and push authorized development branch. Sparse storage, R_prop and electromagnetic phase handling remain subsequent work.

Verification ledger: coefficient probe first failed because the new API did not exist. After implementation, distributed reference comparisons and finite-value checks pass. Native direct/default regressions pass; HSE and non-HSE MPI builds succeed. Read-only review found no blocker; its finite-value test request was applied. Direct ACE-failure fallback remains reviewed but not force-tested. Sequential Diamond timing completed: C64 42.083→40.527s; C128 64.000→52.698s. Max current difference1.13e-17, density5.00e-15, printed energy difference0. Single samples and unequal background load limit performance attribution. Notebook and JSON updated.
