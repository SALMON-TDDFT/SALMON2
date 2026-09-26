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
- The union of fragment states with positive occupation at any k forms the localized subspace; exactly empty bands are omitted from the source while all retained states remain action/ACE targets. For fractional occupations use Q=Psi sqrt(f/2) U as density factors, distinct from the orthonormal localized basis Phi=Psi U. This exactly preserves the fractional density operator at full support.
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


### Si64 follow-up implementation and verification

- Reused the historical Si64 DG fixture; source hashes and a bounded 128-state
  pilot are preserved under samples/dc_hse/si64.
- Exact positive-occupation source selection now reduces source work without
  changing retained targets/ACE. Tests cover source count changes, equal-count
  index changes, all-zero sources and reactivation, as well as a band that is
  empty at only one of several k points.
- Transport-only refreshes skip six unnecessary grid/band overlap products;
  arbitrary-target full exchange remains invariant.
- Added per-fragment final binary snapshots, Gamma NPZ conversion, offline
  localization with explicit occupation-cutoff provenance, and batched screening
  diagnostics. These are diagnostics, not DC restart or production pair pruning.
- Verification after these changes: 11 unit tests; HSE 422 and complex LCFO 130
  CTest prep/run/verify (6/6); 6 GS/DC comparison cases; 4-step RT interval 1/10
  and restart comparison. Max interval energy/current differences were
  1.02e-14 / 2.47e-16; restart energy difference 2.92e-11 in output units.
- HSE-disabled serial build, C2H2 GS, and its Python3 eigenvalue verification pass.
- Review found no correctness blocker; reduced snapshots do not preserve original
  band indices or total retained-state count, which is now documented.
- Si64 physical convergence and useful spatial screening remain unproven. The
  preliminary frozen-state study found only 0.63% pruning at 1e-2 Ha per-fragment
  bound, and direct screened columns lost Hermiticity. Do not enable this
  approximation in ACE/SCF on the strength of an exchange-energy bound alone.
