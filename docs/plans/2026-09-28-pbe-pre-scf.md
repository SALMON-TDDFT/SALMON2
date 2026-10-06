# PBE Pre-SCF Implementation Plan

> 2026-10-06：本記録中の`developer_tests`は当時のローカル開発検証です。GitHub配布から除外しました。通常の回帰試験は`testsuites`を使用します。

> **For Claude:** REQUIRED SUB-SKILL: Use executing-plans to implement task-by-task.

**Goal:** Opt-in PBE warmup before native/DC hybrid SCF.
**Architecture:** Separate runtime stage gates all EXX and rVV10; reuse the Libxc
PBE semilocal evaluator with exchange weight one, then restore the target operator
at a collective SCF transition. No basis change or new wavefunction storage.
**Tech Stack:** Fortran, MPI, Libxc, Python unittest, existing ScaLAPACK build.

Spec: `docs/plans/2026-09-28-pbe-pre-scf-design.md`.

### Task 1: Behavioral regression (RED)
Create developer_tests/653_functional/test_pre_scf.py. Check staged HSE/PBEh/rVV10
energies against direct convergence, require transition log and no earlier MLWF,
validate controls and rejection of unfinished warmup. Run against old binary;
Expected: new namelist controls are rejected.

### Task 2: Stage implementation (GREEN)
Modify src/io/{salmon_global,inputoutput}.f90 for controls, defaults, broadcast,
validation and logs; src/gs/main_dft.f90 starts runtime stage before potentials;
src/xc/{hse_native,hse_semilocal,salmon_xc}.f90 gates EXX/rVV10 and selects PBE;
src/gs/scf_iteration_dft.f90 implements consecutive collective readiness,
rebuild, mixing reset and completion guard. Build existing ScaLAPACK tree.
Expected: new tests pass; default direct path unchanged.

### Task 3: DC/MPI regression and documentation
Exercise DC PBEh+rVV10 at MPI2/4 and native MPI2; run existing EXX input and
adaptive SCF regressions. Document exact input semantics/limitations and measured
energies in docs/inputs/exx-mlwf.md. Review changes and commit after verification.

Progress: design approved; no shared interface conflict (runtime stage owned by
main_dft/SCF, read-only consumers in XC/EXX). Existing dedicated branch reused.

## Execution record

- Task 1 RED: new namelist controls rejected by original executable.
- Task 2 GREEN: ScaLAPACK build; HSE/PBEh/rVV10 native staged energies agree
  with direct SCF, same PBE-stage histories; incomplete-stage guards verified.
- Independent review found DC spectral consistency and auto-mixing baseline
  defects. Both fixed; auto-mixing regression failed at switch iteration 80
  before reset, then passed. Last-iteration flag and first-real-residual checks
  also reproduced before fixes and passed afterward.
- User clarified DC GS uses 300 K and common chemical potential. Reused 0 K
  fixture was inappropriate; an attempted zero-T filling extension was removed
  completely. DC staged inputs now require positive specified temperature.
- DC300 full-support MPI2/4 and charge conservation passed. Empty states and
  appropriate initial orbitals are needed to pass weighted-capacity validation.
- Exploratory DC300 fraction=.999 did not converge (direct and staged); it is
  recorded as a limitation, not a success regression. The input/control,
  native .999 and DC300 full-support regression suite is the release check.

- Subsequent user steering: DC fragments should not receive additional MLWF
  localization. Added default canonical full-fragment DC source; explicit
  yn_exx_dc_mlwf=y preserves old comparisons. No gauge/previous allocation in
  spatial canonical source, no gauge operations in full-k source. New controls
  failed against the prior binary (RED), then canonical vs localized full-support
  energies passed, including 300 K charge, MPI space/orbital and full-k layouts.
- Second independent review: no blocking issues in canonical weighting,
  source consumers, LCFO action, snapshots or unchanged native RT.
- Final verification: 30 tests in test_pre_scf/test_exx_inputs/test_pair_screen
  passed; four existing adaptive SCF tests passed with explicit legacy DC
  localization; canonical/gauge probe passed under bounds checking on MPI1/2/4.
  HSE/MPI/ScaLAPACK build succeeded. No remaining review blockers.
