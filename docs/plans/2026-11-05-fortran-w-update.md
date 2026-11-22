# Internal Fortran W Update Implementation Plan

**Goal:** Complete static and pump screening, MPI distribution, I/O and W replacement inside SALMON without production Python.

**Architecture:** Reuse the tested native q response kernel from SALMON-linked modules. Port orchestration and dielectric/box/head processing to Fortran, maintaining reference conventions and rejection thresholds. Integrate with td projection and commit W only after collective validation; propagation orbital count remains unchanged.

**Tech Stack:** Fortran, SALMON communication abstraction, OpenMP, BLAS/LAPACK; Python solely for independent validation.

## Approved scope
The user approved internal integration, fractional occupations/intraband response, guarded q->0 treatment, eta comparison, distributed k work and pair-owned W, checkpoint/restart compatibility, and 4 MPI x 12 OMP verification. Preserve existing production jobs and failed evidence.

## Tasks and gates
1. Add tests for dielectric/box assembly, metallic head, q interpolation and complete W versus Python. Observe failure before implementing new modules. Test equilibrium and pump branches, q=k and q<k, and invalid/nonfinite inputs.
2. Implement Fortran screening modules under src/xc with explicit state/arguments, scalar IEEE predicates and BLAS/LAPACK. Reuse q_response without duplicating maintained kernels. Keep residual, eta and head checks.
3. Implement k-distributed driver with collective error handling and bounded communication, checkpoint I/O. Compare serial and 4 MPI outputs against Python reference. No root full Nk*Na^2 gather.
4. Integrate GS and projection_option=td W updates into SALMON, with documented namelist defaults/validation/logging. Avoid dense Hamiltonian filesystem round trips for internal updates. Preserve HSEX projection Hamiltonian and propagation manifold.
5. Audit source, build standard and MPI variants out of source, test complete static and pump fixtures and restart/rejection. Audit all application objects before transfer. Review specification then quality; resolve findings.
6. On Fugaku compile/audit, run 4 MPI x 12 OMP with /vol0004 specified if validation harness requires Spack. Record numerical errors and timings, not just successful submission. Save calculation records under canonical calc2026; do not collect GS data.

## Validation requirements
Compare response matrices, assembled W and update effect on RT current. Existing Python tolerances are authoritative; no tolerance weakening. Check finite data and NaN/Infinity rejection, empty owners and collective failures. All mandatory eigensolver audits remain enabled. Existing dirty changes must be preserved.
