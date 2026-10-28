# Exchange Worker Automatic Layout Implementation Plan

**Goal:** Choose worker geometry from symmetry-reduced r points and full r-prime points, with explicit overrides.

**Architecture:** An integer-only selector models the existing blocked-cyclic loops and scratch limits. The worker calls it after symmetry enumeration; input and service validation accept zero as automatic. Existing positive manual settings remain supported.

**Tech stack:** Fortran, MPI/OpenMP, FFTW, Python unittest/oracle, CMake.

1. Add a failing selector probe and an independent Python exhaustive-assignment oracle. Cover production geometries (164,4096,1280/27648 and 305,8000,2560), small and prime communicators, manual overrides, scratch caps and invalid inputs. Run to demonstrate missing automatic selection.
2. Implement `src/xc/exx_worker_layout.f90`, register in application and worker probe build lists. Use 64-bit counts and scalar inquiries only (selector needs no floating point). Test the module independently.
3. Integrate after orbit enumeration in `exx_spinor_wave_workers.f90`. Make input defaults zero; relax only the zero-value guard in input/service. Record selected geometry and work estimates. Update input documentation patch.
4. Extend scalar and SOC worker probes for automatic/manual equivalence, rank-mismatched settings, NaN/Inf and capacity rejection. Source-audit, compile the affected modules, object-audit, then execute 2/4/8-rank comparisons with existing thresholds.
5. Build the complete MPI application with guard checks enabled; run focused registered probes and input regression. Record exactly which compiler/platform were tested. Fujitsu and remote automatic-binary use remain unverified until independently tested; do not change fixed production executables.
