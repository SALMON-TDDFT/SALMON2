# Selected WF column communication plan

Goal: communicate only fragment-support WF columns without changing U or the spherical integration approximation.
Architecture: cache a column request plan on top of the existing row halo. Exchange requested column IDs once, pack only requested row/column rectangles thereafter, and reconstruct directly from compact frames. Preserve dense transport, core diagnostics, and full-range behavior.
Tech stack: Fortran, MPI Alltoallv, BLAS.

The continuing user authorization covers WF/communication optimization and sequential diamond benchmarking. Reuse the existing dedicated branch. Alternatives (sparse approximate U or distributed initial localization) change more of the algorithm and are deferred.

1. Add a failing MPI test comparing selected-column transport with full halo then column selection. Include rank-dependent ordering, duplicate rows, empty columns/rows, zero local rows, and repeated use.
2. Implement a cached selected-column halo plan and compact reconstruction input; integrate only in MLWF source construction.
3. Run standalone tests, native MPI2/MPI4 trajectory regression and build. Run one numerical job at a time.
4. Measure frozen prior binary versus new binary for diamond C64/MPI8 and C128/MPI16 R6 with identical GS, dt=.02, nt=16, ACE1, OMP1/BLAS1. Compare current, density, energy, and wall/phase timings. Record modest or negative speed results honestly.
5. Update both notebooks, verify diff, commit and push.
