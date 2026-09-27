# Distributed Gamma mesh exchange implementation plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Connect the approved spatial ACE route to native PBEh real-space Ehrenfest RT with distributed MLWF refresh and full-support exchange generation.

**Architecture:** First support Gamma, occupied-only, orthorhombic full-support mesh RT with x-complete y/z pencils. Local wavefunctions and MLWF source rows stay distributed. Reduce localization links and temporal overlaps; start from identity gauge (full-support exchange is gauge invariant). Reuse gauge_minimize and FFTW pencils for pair convolutions with the same spherical Coulomb multiplier. Construct/apply ACE through spatial reductions. No LCFO time propagation. Initial DC reconstruction retains existing single-orbital global scratch and is explicitly outside the scalable RT stage.

**Tech Stack:** Fortran, existing FFTW pencils, native MPI communicators, Python executable regression.

1. Add a failing test in test_ehrenfest.py comparing pulse energy, electron/ion current and final forces/velocities at MPI1/2/4 (1x2x1 and1x2x2). Observe old k-only rejection.
2. Add hse_spatial.f90: local source/gauge/previous state, global small-matrix localization/transport, batched pair convolution, exact Coulomb multiplier, collective layout validation. Add optional reduction to gauge_transport preserving serial callers. Register module in CMake.
3. In hse_native.f90, permit the bounded spatial RT route, refresh spatial exchange/ACE and energy with comm_r, apply ACE with reductions. Preserve old DC/projected/serial paths. In inputoutput.f90, admit only the supported pencil layout for the new route; preserve all other guards.
4. Run parity tests, independent serial-versus-distributed exchange oracle, pulse work/impulse/water regression, HSE ON/OFF builds and independent review. Document limits, measured evidence and local commit. Do not claim orbital decomposition, finite-radius MD or giant-system timings.


## Execution record

- Old input gate observed rejecting MPI2 before implementation.
- Distributed Gamma source/localization/FFT/ACE path and bounded native admission implemented.
- Serial exchange oracle, spatial frames and polar transport agree at MPI1/2/4; 14 native tests, 8 input tests, 10 Wannier tests and gauge-gradient test passed. HSE ON/OFF builds passed; projected LCFO response regression passed.
- Independent review found no native integration blocker and identified early local-return risk in optional gauge_transport callback. MPI2 negative-dv-on-one-rank test timed out before fix; collective validation fixed it. Oracle and full-pulse native parity rerun passed after correction.
- Ruling: initial identity gauge avoids a full-grid pivoted-QR seed; full-support exchange is gauge invariant. Finite-radius locality accuracy is not claimed. Initial DC reconstruction remains unchanged and its global scratch limit is documented.
- Numeric evidence: docs/results/pbeh40-rvv10/spatial-mesh-validation.json. No push or merge.
