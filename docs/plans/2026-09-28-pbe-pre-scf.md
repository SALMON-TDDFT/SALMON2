# PBE Pre-SCF Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use executing-plans to implement task-by-task.

**Goal:** Opt-in PBE warmup before native/DC hybrid SCF.
**Architecture:** Separate runtime stage gates all EXX and rVV10; reuse the Libxc
PBE semilocal evaluator with exchange weight one, then restore the target operator
at a collective SCF transition. No basis change or new wavefunction storage.
**Tech Stack:** Fortran, MPI, Libxc, Python unittest, existing ScaLAPACK build.

Spec: `docs/plans/2026-09-28-pbe-pre-scf-design.md`.

### Task 1: Behavioral regression (RED)
Create testsuites/unit_pbeh_rvv10/test_pre_scf.py. Check staged HSE/PBEh/rVV10
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
