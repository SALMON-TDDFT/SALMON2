# Owned XC evaluation phase (2026-11-15)

Continues commit `65cc34ce` (k/XC/pseudopotential inputs and scalar grid ownership). Local Apple Silicon/GNU Fortran/OpenBLAS/Libxc/FFTW; synthetic small restart fixtures. No GitHub Actions or network builds. At most 4 total workers; native builds 2, numerical runs 1–2 threads. Outputs are external.

## Scope

`exx_set_phase(system, pre_scf)` gives a system explicit phase ownership; `exx_pre_scf(system)` resolves the common phase used for EXX enablement/thermal helpers and semilocal/rVV10 evaluation. RT initialization owns production phase before initial potential evaluation. Legacy GS/DC state machines and no-argument EXX queries retain their shared-phase fallback; full GS/DC phase-state extraction is not claimed.

Changing phase clears the explicitly bound owner's exchange resources/history, preserving its binding flag and system association. This prevents frozen prediction flags from skipping a required rebuild. Callers returning to production must refresh from current orbitals before applying exchange. An unbound system can select construction-time phase but must not reuse legacy shared exchange caches. No arbitrary concurrent-call or general native failure-atomicity claim.

## Red/green and regression evidence

1. Poisoning legacy `exx_pre_scf_active` after independent owner creation failed the previous audited build with `STOP 4` (action changed). This preceded phase implementation.
2. Review identified a frozen-cache hazard in the first setter. A changed-orbital/frozen-owner transition test failed that implementation with `STOP 27`. The final setter clears owner resources/history and retains binding lifecycle state.
3. The final native owner probe checks independent phase selection, zero exchange action during pre-SCF, unchanged other owner, legacy query fallback, and refreshed changed-orbital action against an independent reference, in ACE and full-action modes.
4. Full RT tests poison the legacy phase immediately after opening A; B and reopened states initialize with the poisoned flag still set. Eight k2/k3 propagation/factor cases compare forward/reverse/recreated raw traces against separate unpoisoned processes. Fixed timestep, OMP layout and capacity also differ as previously tested.
5. Thirty-two ordinary periodic full-driver cases compare current and energy against the original immutable fast-HSE PR. These are RT wiring/numerical-preservation tests, not converged material physics or a full legacy GS pre-SCF transition benchmark.

## Final receipts

| Check | Result |
| --- | --- |
| Native ACE/full owner phase and frozen-history regression, OMP 2 | Passed |
| Full RT contexts with shared-phase poisoning, OMP 2 | 8 cases bit-identical |
| Original periodic RT driver, OMP 1 | 32 cases, current/energy difference 0 |
| Python source tests | 305 passed |
| sdist → wheel | Passed |
| Fresh installed wheel | 305 run, 1 optional plotting skip, no failures |
| Installed source and bilingual/AI guides | Byte-identical to tested checkout |
| Source/all-object eigensolver audits before links | Passed |
| New executable Fortran lines and guide synchronization | Passed |

Executable SHA256: `992eb6b75491713fdbd97269498056d5930bb497905583f7001c20a9c22d6247`.

Replacement manifest SHA256: `2ff06bfd89604cf22b75f807165721fad4502b2ac33537c795d9464ea3489cf8`.

Local receipts: `/tmp/salmonpy-xcphase-final-build-20261115.json`, `/tmp/salmonpy-xcphase-{native,contexts,default}-20261115/result.json`. Package: `/tmp/salmonpy-xcphase-final-dist-20261115`. These machine-local paths are validation provenance, not installation requirements.
