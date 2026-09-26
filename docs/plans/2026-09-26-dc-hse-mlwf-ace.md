# Fragment-periodic DC HSE, Wannier localization and ACE

**Goal:** Build a new branch from hse-merge-prep that solves HSE within each DC fragment plus buffer under periodic boundary conditions, retaining complex multi-k LCFO.

**Architecture:** The existing whole-system Hartree and DC core assembly remain. Each fragment owns its screened exchange and ACE. A rectangular full-mesh Wannier backend supports fractional one-spin source occupations; retained orbital rotations are unitary and do not replace the occupation matrix with identity. Previous localized periodic orbitals seed a polar-transported MV minimization. Full-support exchange is the numerical baseline; localization convergence and performance are reported, not presumed. HSE mixing remains 0.25 and omega remains 0.11 bohr^-1 by default.

**Tech Stack:** Fortran, FFTW, BLAS/LAPACK, SALMON MPI abstraction, Python/NumPy independent numerical fixtures.

## Decisions

- User selected fragment-level HSE-SCF, not merely DC initialization of conventional HSE.
- Use new local branch dc-hse-mlwf-ace based on 483058c; no push or merge is requested.
- Native worktree tool cannot operate in this projectless chat (Not a git repository); use a sibling worktree of the already downloaded repository.
- Initial parallel scope: independent fragments with k-only decomposition within each fragment. No new GPU, spin, DFT+U or ionic-motion support.
- Generalize the new backend to rectangular fragment grids and 1 x Nk_y x Nk_z meshes; preserve the existing native HSE backend for existing inputs.
- All retained fragment states form the localized subspace. For fractional occupations use Q=Psi sqrt(f/2) U as density factors, distinct from the orthonormal localized basis Phi=Psi U. This exactly preserves the fractional density operator at full support.
- Energy accounting: subtract the complete HSE core expectation from the inferred ionic nonlocal term, then add half that expectation once to DC E_xc. Do not add fragment-cell exchange energies to a total-system semilocal E_xc.
- LCFO requires full exchange action on its basis, not occupied-space-only ACE interpolation. Retain the full fragment source during LCFO assembly.
- Reuse localization state in memory through SCF/RT; restart may reinitialize the gauge (same full-support exchange), with this behavior documented.
- Local-support truncation is not silently introduced: source-only truncation is not a variational SCF functional. First validate the untruncated localized operator, then expose only tested controlled approximations.

## Tasks and tests (execute inline)

1. Baseline: configure/build MPI HSE, run existing native unit tests. Save logs outside source.
2. Write failing compiled tests for rectangular full-mesh exchange with fractional occupations, arbitrary target Hermiticity, weighted density invariance and core exchange energy. Include shifted/shuffled mesh and Gamma. Implement hse_wannier.f90 and independent NumPy reference; require agreement at ~1e-11 and nonpositive exchange metric.
3. Write failing compiled polar transport/MV tests. Port the algorithm from TDCDFT, preserving reciprocal boundary phases. Check unitary transformations, occupied-basis covariance, analytic spread gradient, and transport under a pure occupied rotation.
4. Integrate the backend, occupation-sensitive cache, ACE and inputs with full defaults/broadcast/logging/checks. Preserve old HSE tests. Add a DC-HSE numbered test before relaxing guards; verify initial rejection, then successful SCF, finite core energies and complex LCFO eigenvalues on Gamma and multi-k meshes.
5. Verify one-fragment/no-buffer DC against conventional HSE using the same backend. Verify decomposition parity for multi-fragment/multi-k runs, and full vs ACE actions on construction orbitals. Build HSE-disabled configuration.
6. Document implemented scope, approximations, tests, timings and unresolved scaling work. Independent whole-branch review, fix actionable issues, final diff check and local commit. Do not publish or merge.

## Progress

- Implemented rectangular fractional-occupation fragment exchange, polar/MV gauge optimization, ACE, DC core energy accounting and full LCFO action.
- Added phase continuation through the complex ±pi cut, with a compiled regression. Independent review found no remaining blocker for the full-support baseline.
- Four new compiled unit tests pass; eight legacy native tests including MPI pass.
- MPI DC-HSE case 422 and existing complex PZ-LCFO case 130 pass.
- Six fresh GS/DC comparisons pass, including genuinely smaller fragments, Gamma, decomposition parity and one-fragment equivalence.
- Four-step RT comparisons with intervals 1/10 and a step-2 restart pass.
- HSE-enabled MPI and HSE-disabled serial builds pass. Serial C2H2 calculation and Python3 verification pass (the legacy CTest launcher expects a missing `python` alias).
- Updated local input documentation and SALMON-DOCS on its own dc-hse-mlwf-ace branch (290af91).
- Scope remains the exact full-support baseline: no spatial pair pruning/local Poisson/distributed pair scheduler or large-system scaling claim. Some retained-subspace fixtures reach the localization iteration limit; diagnostics expose this, while exchange remains gauge invariant.
