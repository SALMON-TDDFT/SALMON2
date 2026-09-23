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
