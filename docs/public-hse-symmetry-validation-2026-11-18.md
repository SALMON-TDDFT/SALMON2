# Public fast-PR HSE symmetry and protected mixed connection — 2026-11-18

Result: explicit bounded C/Python connection, direct/protected mixed persistence
and failure/retry tests pass. Robust physical accuracy is NOT established: a newly
observed near-zero-gradient Libxc HSE derivative issue remains unresolved below.
No numerical solver, compatibility patch, gradient floor or comparison tolerance
was changed. The fast-PR profile is the only HSE engine.

## Contract

- Existing `sp_fast_hse_init` remains symmetry-denied. New
  `sp_fast_hse_init_symmetry` has identical arguments, with
  `sp_fast_hse_symmetry_protocol()==1`. `NativeFastHSEBackend(...,
  allow_symmetry=True)` explicitly selects it; the default behavior/identity
  construction remains unchanged. The opt-in is identity-bound, group1 only.
- Only declared unanimous finite34 maps to a state-preserving InputError. Correct
  A before retry. Unknown/old protocols do not establish recovery. Internal-k2,
  learned/PT-CN symmetry, malformed initialization and general reentrancy remain
  unsupported. Actual local A, including radiation from other materials, must
  preserve the retained group. No projection or fallback occurs.
- The ordinary mixed sample takes --hse-symmetry and an explicit equivalent
  --hse-reference-bundle for its serial full-k verification. User GS files are
  never converted. Production native continuation does not replay accepted steps;
  the demo separately reruns a serial reference for verification.

## Local Apple Silicon checks

- Fail-first Python constructor/protocol/recovery tests;17 focused adapter and3
  CLI admission tests pass. Old/missing protocols, default-denied symmetry and
  internal-k2 reject before Python consumes the native lifetime.
- Two independent COMM_SELF groups use different screening values.32 steps with
  a forbidden tiny-x input at4: native/Python values and checkpoint bytes remain
  unchanged; corrected z drive succeeds. Fresh7→32 direct resume has exact
  current/energy tails and final checkpoint bytes, zero elapsed accepted replay.
- FDTD1/2 + two local LD + SBE1 + symmetry HSE1 (3/4 ranks) and additionally
  PBE-TDDFT2 (5/6 ranks). Full/reduced serial response checks use unchanged
  rtol2e-10/atol2e-12. Fault/retry and field-decomposition A/E/j histories match
  exactly. The positive fixture explicitly chooses orbital-phase.23.
- Whole-step late-field/SBE/electronic, native deep numerical faults and new
  forbidden-field admission errors preserve original numerical bytes, Python
  metadata and ordinary data streams. Retry results match ordinary execution.
- Four/six-rank protected restore save16/fail18: all original field/LD/SBE/PBE/HSE
  bytes preserved, normal continuation allowed; three retried future frames exact,
  zero replay. Six-rank fresh CLI16→32 direct resume matches continuous histories.
- Existing non-reduced GS/MPI2 and PBE/HSE Taylor/PT-CN mixed regression remains
  exact. Source/installed suites, source-distribution-to-wheel/native hash checks
  and installed mixed/public-protocol tests are retained with evidence. No Actions.

## Newly observed numerical issue (not fixed)

The original phase0 synthetic case completes, but reduced/full currents differ
about5.7e-10 (above the existing comparison tolerance), even without fault
injection. A smaller step does not resolve it. New-native full-k and old serial
full-k comparisons localize it away from the public ABI or rollback path.

A separate diagnostic build finds initial Hamiltonian and exchange actions equal
within about4e-17; after one step density differs about7e-15 but Vxc differs about
0.016. A pure Libxc5.2.3 C-ABI probe at rho=.008 and grad≈1e-16 yields finite but
large gradient derivatives (~0.0055), whereas a slightly smaller gradient gives
near zero. This explains sensitivity but is not a validated mathematical remedy.
`probe_hse_low_gradient.py` records the exact library hash/version/curve and flags
this regime. Diagnostic exit0 is not a physical pass. Failed and negative results
are retained. No smoothing/clamping, silent GS perturbation, library patch or
relaxed tolerance has been applied. Phase.23 deliberately avoids the flat-gradient
nodes for connectivity testing; it is neither a production remedy nor a material
convergence study. A real-GS mixture that generates incompatible local A is also
retained as a correctly rejected case with original_preserved=True.

Highest next numerical priority: derive/verify a consistent stable low-gradient
HSE derivative treatment. Reduced internal-k MPI and broader initialization
remain later work. Memory budgets are advisory, molecular work paused.

Evidence: `~/.codex/outputs/salmonpy-public-hse-symmetry-20261118/`.
