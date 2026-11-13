# FDTD / LD implementation plan

Goal: preserve legacy arithmetic, phase order and restart format while extracting LD.
Baseline: upstream 3322d056 plus separately verified PZHEEVD replacement and mandatory guards.

1. Add standalone LD behavioral regression harness; observe missing module failure.
2. Replace forbidden driver, validate three workspace types and allocation collectively; preserve guards.
3. Build and retain baseline executable and checksums, run no GS calculations.
4. Extract LD state, coefficient creation, old-field storage, P/J/Q update; retain dense layout and loops.
5. Redirect FDTD allocation, checkpoint and observation references to owned LD state; preserve field-update order.
6. Test standalone LD against frozen legacy kernels, multi-step FDTD bitwise results, restart, serial/MPI.
7. Document AI usage, concrete available interfaces, limitations and evidence.

Results live in canonical calc2026; no production jobs or fixed binaries changed.
