# Direct WF coefficient Taylor4: untruncated reference

**Goal:** Propagate localized LCFO coefficients with the existing Taylor4/predictor-corrector; keep the integration method unchanged.

**Architecture:** Accumulate Taylor polynomial in local LCFO coefficient rows; reconstruct grids only to call the existing native Hamiltonian. Initially and after accepted endpoint field refreshes, rotate native orbitals into the current WF frame and reset the occupied rotation to identity. Preserve physical WF anchors, ACE operators and midpoint flow; update coefficient cache keys after rebase. This is an untruncated correctness reference with dense gauge correction still present, not the final sparse production method.

**Tech Stack:** Existing Fortran/BLAS/MPI Taylor4 and HSE adapter.

User explicitly directed Taylor4 over proposed PT-CN on2026-09-27. No PT-CN implementation was made; only a rejected-path test was tried. No implicit iteration or PT-CN solver changes belong in this work.

- [x] RED: native direct-WF flag is ignored by old binary and expected coefficient-Taylor marker is absent.
- [x] Configure opt-in SALMON_LCFO_RT_DIRECT_WF=1 for MLWF/Taylor4 only.
- [x] Add coefficient Taylor adapter reusing hpsi; standard polynomial coefficients and predictor/corrector scheduling unchanged.
- [x] Rebase actual native orbitals at initialization and accepted endpoints; retain physical reference frames and repack valid coefficient cache keys without rebuilding invariant exchange operators. Do not truncate propagation support.
- [x] Verify full-range density/current/energy/Gram and half-dt results against standard Taylor; MPI spatial/orbital partition invariance, U/ACE cadence compatibility and rejection of incompatible settings.
- [x] Review, then Diamond reference comparison. Do not claim speedup before timing.

Remaining after dense-limit reference: sparse coefficient ownership and Hamiltonian/ACE action, independent R_prop convergence, and electromagnetic phase handling with finite-basis leakage/boundary tests. exp(i A.r) is not assumed unitary in a truncated fixed LCFO basis.

Ledger: initial old-binary RED was missing coefficient-Taylor marker. First GREEN comparison revealed init config.h missing; added five-frame Gram assertion (RED observed) and included config.h. Rebase initially invalidated all cache keys, causing extra energy-time ACE rebuilds; now repacks actual rebased native orbitals as the new valid key, leaving absent retained-ACE keys absent. Direct/full/finite-radius+U2/ACE4/MPI2xorbital2 and direct half-dt regressions pass; read-only review finds no blocker. Ordinary-mode transport/cadence/native regression pass. Diamond comparisons completed sequentially; results are recorded below.

Orbital-indexed outputs in direct mode refer to the rebased WFs rather than original GS bands. Total observables remain the comparison targets. This reference adds grid-to-coefficient projection/reconstruction and a Gram gate; no sparse performance benefit is claimed.

Diamond ledger: C64 ordinary38.311/direct43.695s; C128 ordinary56.126/direct60.490s. Both35ACE builds,32Taylor calls and17direct Gram checks; numerical comparisons pass. Single samples; reference is slower, default remains ordinary. No sparse speedup or phase compensation is claimed.
