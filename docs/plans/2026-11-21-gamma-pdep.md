# Gamma supercell PDEP Implementation Plan

**Goal:** Add an opt-in production Gamma-only static RPA screening mode for amorphous/water cells, using occupied-only real-space Sternheimer, low-rank optical head/wings, orbital MPI and inner Hpsi OMP, and reusable dielectric-mode guesses.

**Architecture:** `hsex_w_solver='pdep'`, Gamma mesh1³, CPU scalar insulator, full spatial grid, nproc_k1 and orbital MPI. A new adapter gathers occupied projectors, partitions source bands, and owns independent single-vector Hpsi workspaces. The response uses local FFTW density/potential transforms rather than real-grid-by-G phase matrices. The low-rank response includes the three analytic optical vertices and all finite G; Schur/Woodbury identities project its inverse into boxes without a full G² dielectric. Refactor the existing optical cell quadrature into a shared integration routine, preserving its head, wings and Gamma-cell convention. Gamma `sternheimer` gives an explicit full-G reference. Existing finite-q/default paths remain selected as before.

**Tech Stack:** Fortran, MPI wrappers, existing Hpsi/grid velocity, serial FFTW transforms inside band-MPI work, BLAS/LAPACK, existing atomic checkpoint conventions. No production Python.

1. Failing tests for adaptive PDEP modes/warm start and low-rank optical-cell equivalence to existing full response; test NaN/Inf and indefinite dielectric rejection.
2. Add residual-driven restarted subspace solver, bounded O(NG*r) storage, with explicit convergence rejection and caller-supplied initial modes. Preserve old finite-q fixed-subspace API.
3. Extract shared optical quadrature, add factored Woodbury/Schur optical-cell projection, verify against the original full matrix including anisotropic wings.
4. Add Gamma FFT response with occupied-band MPI, analytic optical vertices, true CG residual and band work accounting. No grid Hamiltonian, empty-band sum or grid*G phase allocation. Occupied projector is replicated initially; document this limit rather than claiming distributed-grid scalability.
5. Add protected basis cache read/write (new file, no overwrite), strict cell/grid/G/rank compatibility, but allow changed ionic coordinates. Recompute all current-operator actions; cache is only a guess, not old screening reuse.
6. Connect Gamma input, real/complex GS states and independent Hpsi scratch, standard W output; permit Gamma reference and PDEP while keeping finite-q restrictions elsewhere. Preserve default mode and existing jobs/binaries.
7. Audit source before all builds and linked application objects before new executable use. Test native math and 4MPI12OMP band distribution, small periodic water full-reference/full-rank, truncation/warm-cache and RT reading where feasible. HSE-off build; report Fujitsu/GPU untested.
8. Update SALMON-DOCS, reproducibility tests, record only derived results/logs/code under canonical calc2026 (no GS payload collection). No commit unless separately requested.
