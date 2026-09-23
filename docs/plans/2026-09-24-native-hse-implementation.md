# Native HSE06 MPI/BLAS Implementation Plan

> Execution: superpowers:executing-plans, inline in the existing TDCDFT feature checkout.

**Goal:** Implement native full-kernel HSE06/ACE with MPI and BLAS, validate Si against the Python reference, then integrate native SCF and RT.

**Architecture:** A standalone Fortran exchange/ACE numerical module, a SALMON orbital/source adapter, explicit semilocal HSE06 and exchange energy accounting, and a selectable self-consistent PT-CN path. Keep trial vectors distinct from occupied source vectors. MPI ownership is explicit and unsupported layouts fail early.

**Tech Stack:** Fortran, SALMON communication modules, FFTW, BLAS/LAPACK, Libxc C interface where compiler-compatible Fortran modules are unavailable.

**Spec:** docs/plans/2026-09-24-native-hse-design.md (approved).

## Task 1: Fortran exchange kernel and reference fixture
Files: src/xc/hse_exchange.f90; samples/hse_mlwf_reference/native_exchange_probe.f90; samples/hse_mlwf_reference/test_native_exchange.py; src/xc/CMakeLists.txt.
- Write a compile/run fixture test which requires the native module and compares full exchange on independent complex targets against DistanceExchange. Observe missing native module failure.
- Implement reusable phase/kernel/FFT workspace with blocked density formation and application via ZGEMM. Use rank/size row ownership and leave collective combination to caller.
- Test tiny shifted/shuffled k meshes, uneven/empty row partitions, Hermiticity and energies; compare actual Si state and time kernel build/action separately.

## Task 2: Native ACE
Files: src/xc/hse_ace.f90; native probe and tests.
- Write failing comparison of occupied interpolation and arbitrary target ACE actions with Python ACE.
- Implement grid-weighted metric, LAPACK Hermitian factorization and two ZGEMM applications; retain sign, Hermiticity and conditioning checks.
- Test linearity and negative metric rejection, MPI invariance of assembled actions.

## Task 3: SALMON adapter and HSE06 semilocal integration
Files: src/xc/salmon_xc.f90; src/xc/hse_native.f90; src/common/hamiltonian.f90; src/common/total_energy.f90; src/common/structures.f90; input/initialization sites discovered by call graph; CMake dependency configuration.
- Add guarded HSE06 input/build selection. Read supported geometry/occupation and distribution from system/info; no hardcoded Si dimensions in numerical modules.
- Source updates explicitly gather/pack occupied orbitals, refresh full action/ACE and exchange energy; hpsi uses a fixed source state to act on arbitrary targets.
- Enable matching semilocal HSE06 remainder through Libxc without accepting unsupported hybrids. Add exchange energy once.
- Validate native ground-state operator/energy against exported state, then SCF agreement; assert unsupported inputs fail and semilocal regressions pass.

## Task 4: Native self-consistent PT-CN and restart
Files: src/rt/time_evolution_step.f90; new src/rt/hse_ptcn.f90; source update and restart lifecycle.
- Translate existing Python nonlinear/ACE algorithm with full-exchange endpoint residual checks, preserving dt and field convention.
- Test same initial state/step against Python, norm and energy drift, interrupted/restarted trajectory equivalence. Never accept unconverged endpoints.

## Task 5: Measurements, review and handover
Files: docs/results/si-hse-native/; samples and user-facing documentation.
- Run native1/4/8 MPI measurements with BLAS threads fixed; report whole step, exchange/ACE/local/communication timings and memory.
- Run required regression checks; request one independent fresh code review, fix substantive findings and re-test.
- Commit tested changes and document capability limits honestly. Keep active Python MPI reference job until native parity and stability established. No k convergence scan, push or merge.

## Ledger
- Design approved by user. Existing branch TDCDFT is already isolated from upstream main; continue here to preserve requested development history and inputs.
- Interfaces: kernel uses contiguous grid/orbital/k arrays, no HSE mixing factor; adapter applies0.25 exactly once. ACE uses same unscaled exchange. Semilocal energy is remainder only. PT-CN owns source refresh; hpsi never derives sources from trial vectors.

- Tasks1–2: Fortran full kernel (FFTW/ZGEMM) and ACE (ZGEMM/ZHEEV) implemented; shifted/shuffled meshes, different targets, uneven/empty row ownership and occupied interpolation pass Python parity tests.
- Task3: HSE06 C-Libxc and native MPI-k adapter integrated. Native Hpsi relative error8.9e-11; energy error6.8e-14Ha. Identified and fixed SALMON rdedd sign convention by failing snapshot test. WPBEH derivative noise under one-ulp density perturbation is3.8–7.2e-10, so cross-language Vxc/action tolerance2e-9 is justified separately from unchanged PT residual1e-10.
- Task4: PT-CN callback solver and native adapter implemented. MPI8 first3steps pass, first-step wavefunction relative difference4.6e-12 and current absolute difference5.7e-15. Restart metadata and energy-component consistency tests in progress.
- Ruling: first MPI support is k-only with replicated occupied source during EXX refresh, distributed ACE/local Hamiltonian for repeated inner steps. Unsupported r/orbital process splits explicitly fail; broad process-grid support remains deferred as allowed by approved design.
- Independent review: no energy/ownership blocker; found frozen-functional option silently ignored, now explicitly rejected. FFT plans and ACE target work reused. Native restart metadata stores physics and exact pseudopotential bytes.
- Tasks3–5 completed: native fresh SCF converges after87updates, independently evaluated energy differs from Python by1.43e-12Ha, projected residual9.82e-7Ha. Wavefunction-only SCF imports must reset mixing history; samples document this. No solver change was needed for that setup issue.
- Final MPI1/4/8 median step times14.577/4.628/2.947s versus matched Python8 13.802s (4.68x). Rank RSS snapshot295–298MiB. Background reference job remained active; no isolated-node or k-scaling claim.
- Native100step pilot completed with electron error5.79e-11, overlap error7.40e-11, full residual<1e-10 and sampled post-impulse energy range2.60e-12Ha. Restart2+1 matches continuous3 bitwise. Six negative restart/input cases rejected.78unit tests and6existing Si/TDCDFT regressions pass; enabled/disabled builds and regressions checked. Fresh independent review found no further blocker.
