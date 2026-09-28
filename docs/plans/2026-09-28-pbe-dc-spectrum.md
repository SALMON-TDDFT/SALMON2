# PBE DC-to-mesh response and dielectric comparison

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Compare epsilon_xx for the 4x1x1 fragment /32 H2 cell using PBE and PBEh40, with approximately0.5 eV nominal Fourier spacing.
**Architecture:** User approved routing real LCFO files through the existing real reconstruction, including conversion into complex mesh orbitals. Detect file format collectively across every fragment and both files; reject missing/truncated/mixed headers. Keep complex reconstruction and validation. No change to physical model or time propagator.
**Tech Stack:** Fortran MPI, Libxc, ScaLAPACK/LAPACK, Python numerical analysis.

## Task1: Reader regression and fix
- Add testsuites/unit_lcfo_rt/test_real_response.py, taking a binary and existing32H2 PBE seed. Run one RT step on MPI4 and8, compare currents and test mixed/truncated headers.
- Verify unmodified reader fails on real PBE seed (already reproduced: reference metadata invalid).
- In src/gs/dc/lcfo.f90 choose reader from collective on-disk headers, validating real integer headers before allocation. Allow real reader to populate zwf at Gamma, preserving spin/k-point guards. Reject complex files requested as real orbitals.
- Build dedicated USE_LIBXC/USE_HSE/USE_SCALAPACK executable; run reader tests and existing complex/native tests.

## Task2: Matched response
- Existing user-approved cell, x impulse1e-4, PBE versus PBEh40, no rVV10. Hybrid support.999/source ACE, Coulomb cutoff4bohr preserved.
- PBE DC-SCF at300K (no RT temperature), existing hybrid seed. PBE real LCFO uses LAPACK because ScaLAPACK requires complex orbitals; RT MPI4 spatial decomposition1x2x2.
- Short dt.02/.05 and zero-field probes. T=350a.u. gives nominal2pi/T=.4885eV (8.466fs). Choose stable dt based on validation; analysis grid.05eV is not physical resolution.
- Subtract zero-field current from driven current before SALMON cubic-window Fourier transform. Disclose nonstationarity of DC seed: this is DC-initial-state response, not a converged equilibrium dielectric spectrum.
- Monitor long-time norm/energy/ACE acceptance; do not present unstable trajectories as physical spectra.
- Save inputs, provenance, currents, epsilon, plots and Japanese report. No strong-scaling or Full reference reruns.

Execution: continue inline under current user authorization; do not request another execution-mode approval.
