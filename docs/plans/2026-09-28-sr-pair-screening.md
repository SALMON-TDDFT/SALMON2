# SR pair screening implementation plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Diagnose and optionally omit HSE source/target pairs with a conservative discrete exchange-action error budget, preserving the existing spherical support and FFT machinery.

**Architecture:** Evaluate pair densities in the transported localized target gauge. Bound the omitted action using the actual discrete Fourier multiplier, including G=0; accumulate source bounds per target and a global Frobenius bound. Diagnostic mode retains every pair. Opt-in screening is accepted only through the existing strict Hermitian/positive ACE metric checks; otherwise recompute the unscreened action. Preserve an accepted transported gauge when a later minimization fails; never silently change adaptive RT to full support after loss of transport.

**Tech Stack:** Fortran, distributed MPI pencils/orbital columns, FFTW, existing ACE and Python/MPI regression probes.

## Approved design and constraints

The user approved on this turn: retain spherical WF masking; diagnose before dropping; use pair-density and the SR kernel rather than erfc of center separation; verify action, Hermiticity, ACE, and RT switching. Work continues in the existing isolated feature checkout `pbeh40-rvv10-water-md`. No push is included.

For q = occupation-weighted (possibly masked) source and t = target, rho=conjg(q)*t. With the current FFT normalization, C=F^-1 diag(v_G) F. Let dv=grid cell volume, N=total grid size, lambda=max(abs(v_G)), k0=sum(abs(v_G))/N, krms=sqrt(sum(abs(v_G)**2)/N). Then the pair action satisfies

    ||q C rho||_dv <= min(max|q| lambda ||rho||_dv,
                         ||q||_dv k0 sum|rho|,
                         ||q||_dv krms sqrt(sum|rho|^2)).

These are discrete bounds, not continuous-kernel or center-distance heuristics. Per-pair budget is tau/(Nsource*sqrt(Ntarget)); summing source bounds per target and taking the norm gives a total raw-exchange action bound <=tau. Mixing fraction is not included. The bound is for omitted pairs relative to the same WF mask, not the error from that mask or an arbitrary-target ACE bound. Diagnostic mode will establish conservativeness/usefulness first. No speedup is promised.

## Task 1: Preserve transported localization

Files: src/xc/hse_spatial.f90; src/xc/hse_native.f90; testsuites/651_hybrid_exchange/wannier/projected_seed_probe.f90.
1. Add a regression that first accepts a localized state, then forces minimization failure at an unrealistically tight tolerance; require the transported accepted gauge and status to remain usable. Run test_projected_seed.py --build ../work/pbeh40-scalapack-build; expect failure before changes.
2. Save transported U before minimization. On failed minimization with previously accepted state and successful transport, restore U and retain accepted status; expose a retained-gauge flag while leaving the attempt status visible. Apply to both layouts. Keep existing transport-loss reseeding.
3. Adaptive RT must fail explicitly if no accepted localized gauge can be transported/initialized, instead of toggling to full support. SCF retains its existing warmup/fallback.
4. Repeat probe MPI1/2/4 and adaptive RT fixtures, including the previously failing 4x4x1 benchmark. Record energy width and gauge retention.

## Task 2: Discrete pair diagnostics and optional omission

Files: src/xc/hse_spatial.f90; src/xc/exx_spatial_local.f90; new testsuites/651_hybrid_exchange/wannier/pair_screen_probe.f90 and test_pair_screen.py.
1. Write MPI probe with overlapping/tiny-tail/disjoint localized functions and varying SR omega; compare exact and screened action to reported bound, diagnostic equality, compact/global equality, and zero-tolerance identity. Expect missing screening fields before implementation.
2. Add mode (off/diagnose/on), nonnegative raw-action tolerance, counters, accumulated bound, and timing to spatial state. Compute multiplier norms collectively without gathering grids. Compute q/rho norms using small reductions. Count candidates and estimated bound in diagnose but omit only in on.
3. Add optional skip mask to compact apply; pack surviving global-FFT targets in existing batches. Preserve collective calls for empty partitions.
4. Run probes MPI1/2/4; existing exchange validation covers orbital groups and empty layouts. No change to off-mode numerical path.

## Task 3: Native integration and checks

Files: src/io/salmon_global.f90; src/io/inputoutput.f90; src/xc/hse_native.f90; testsuites/653_functional; docs/inputs/exx-mlwf.md.
1. Add input/regression coverage for exx_pair_screening and exx_pair_tolerance, HSE-only and supported spatial/native routes. Default off, zero tolerance means no finite omission. Moving-ion/restart restrictions follow the existing adaptive route. Reject unsupported combinations explicitly.
2. For enabled modes, use transported localized target columns, then rotate the action back to original orbitals without gathering WFs. Build ACE with unchanged strict Hermitian/positive metric validation; on screened failure recompute all pairs and rebuild. Log attempted candidates, bound, skipped pairs and acceptance/fallback distinctly.
3. Compare off/diagnose/on full-support and .999 HSE SCF/native RT using identical seeds; verify actual omission (or honestly report none), action bound, ACE Hermiticity/interpolation, MPI parity and 16-step current/energy behavior. Finite screening is experimental and opt-in.
4. Preserve inputs/logs and numerical report, run relevant existing tests, obtain one final fresh-context review, address findings, commit. No arbitrary memory cap or SR center-distance cutoff.

## Execution ledger

- Plan/design approved by user; implementation begins from c6103d72.
- Ruling: conservative spectral/kernel row bounds are the first supported estimator. Geometric block-distance refinements are deferred until diagnostic evidence justifies them.
- Ruling: strict ACE acceptance and unscreened recomputation are mandatory; a small energy estimate alone cannot authorize a non-Hermitian ACE input.

- Task 1: MPI1/2/4 retention/reseeding probe passed. Failed minimization preserves the transported accepted U. Adaptive RT transport/initial-localization failure now stops explicitly.
- Task 2: raw diagnostic/omission bounds passed MPI1/2/4, global/compact FFT, omega=.11/.22. Native test showed all finite omissions initially fail strict ACE Hermiticity.
- Ruling: add Hermitian completion D=Psi (Psi^dagger Psi)^-1 (A^dagger-A)/2, A=Psi^dagger W. Reserve half the user budget for pair omission; measure the actual distributed correction norm and require raw bound+correction norm <= user budget before strict ACE validation. This preserves the approved total-action bound without relaxing ACE checks. Singular Gram, over-budget completion, or invalid ACE recomputes unscreened W.
- Review: no production correctness blockers; strengthen native tests to require accepted finite omission and accepted bounds, not attempted omission only. Implemented.
- Regression finding: extending accepted-gauge retention to SCF changed a DC-rVV10 convergence trajectory (MPI4 failed within 1000 SCF iterations). Ruling: retention is RT-only, the originally observed continuity problem; SCF keeps its established retry/fallback behavior. Re-run exact DC regression before completion.

- Final verification: native pair screening 3 tests, adaptive RT 3 tests, EXX input 15 tests, adaptive SCF/DC 4 tests all pass. Spatial exchange probe passed 27 MPI/kernel configurations; 18 positive-omega configurations also passed screening with empty orbital groups. Pair-bound/completion probe passed MPI1/2/4 and both FFT routes at omega .11/.22, including subnormal-scale norm regression; accepted-gauge probe passed MPI1/2/4.
- Review follow-up found accumulated tiny bounds/corrections could underflow when squared. Added scaled distributed norm reductions and red/green positive-tiny-bound tests. Final independent review reports no remaining blocker.
- Final RT-only retention DC-rVV10 MPI4 regression converged; full four-test SCF suite passed. Same-seed 4x4x1/MPI8 RT retained spherical support for all 33 updates; energy width about 2.975e-7 Ha, still above prior 1e-7 criterion. This residual is recorded, not marked accuracy-passing.
- Implementation complete; no performance speedup claimed and no push requested.
