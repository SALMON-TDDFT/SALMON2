# Native SALMON HSE06 MPI/BLAS integration — design draft

## Objective

Integrate the independently validated HSE06 reference into the TDCDFT branch of SALMON as Fortran code using SALMON MPI infrastructure, BLAS/LAPACK and FFT plans. Validate Si linear response against the Python full-kernel reference and measure whole-step performance. This is a branch implementation, not a claim of upstream acceptance.

## Existing evidence and integration points

- `src/common/hamiltonian.f90:hpsi` applies local, nonlocal pseudopotential, and other orbital operators. Exact exchange must act on arbitrary trial orbitals using a separately held occupied source state; trial vectors must never silently redefine the source density.
- `src/xc/salmon_xc.f90:init_libxc` currently rejects hybrid families. Only explicitly supported HSE06 should be admitted together with the actual Fock operator; accepting all hybrids or evaluating only their semilocal remainder is incorrect.
- `src/common/total_energy.f90` has distinct periodic/isolated energy assembly. Add the weighted one-half exchange expectation exactly once, independently of the semilocal XC energy. Avoid double counting through eigenvalue sums or cached ACE energies.
- `src/parallel/init_communicator.f90` already supplies real-space, orbital, k-point and combined communicators. The adapter must define unique ownership and avoid duplicate occupation/k weights in reductions.
- Existing SALMON code already calls ZGEMM. Python reference matrix products also use BLAS, so language conversion alone is not a demonstrated speedup.
- Current Python8 timing: full exchange about1.45s, PT-CN step8.77s; serial controller still owns local Hamiltonian and ACE inner applications. These remaining costs motivate native integration.

## Alternatives and selected proposal

1. Port only the exchange kernel behind a Python interface: smallest validation step, but leaves the serial solver and duplicated ownership; useful only as a temporary test harness.
2. Native HSE operator plus ACE, ground-state and real-time integration: recommended. Reuse existing SALMON data and MPI decomposition, preserve a full-exchange validation path, and optimize measured costs.
3. Implement a general-purpose hybrid framework for every geometry, spin and process distribution immediately: much wider validation scope than needed for this Si comparison; defer unsupported cases behind explicit guards.

## Proposed implementation sequence

### A. Validated Fortran operator

Add a module for the same sampled screened periodic kernel, finite q=0 limit, twist phases, k-grid ordering, FFT normalization and orbital normalization used by the reference. Initially retain the full kernel. Precompute relative offsets/kernel weights where memory permits; reuse FFT plans and work arrays. Use ZGEMM for density-block formation and application. Stream blocks; do not allocate the full real-space density matrix.

Use a small test adapter to compare Fortran actions, energies, Hermiticity and gauge covariance with exported Python fixtures, including arbitrary target vectors, shifted k meshes, unequal/empty row ownership and MPI1/4/8. This adapter is a verification tool, not the completed native feature.

### B. Native ownership and ACE

Connect source/target packing to s_orbital and explicit communicator ownership. First certify the process layout needed for Si; unsupported layouts fail at initialization. Separate source exchange refresh from repeated trial-vector application. Keep owned ACE factors local to the corresponding grid/k partitions where supported. Form distributed inner products with ZGEMM and reductions, factor the occupied metric with LAPACK (eigensolver for direct reference parity; Cholesky only after conditioning checks), and apply ACE with two matrix products. Retain a full-action residual gate and negative-semidefinite exchange checks.

Avoid a permanent root-worker solver design: exchange setup may initially gather sources, but repeated ACE and local Hamiltonian applications must use SALMON's distributed arrays. Measure gather/reduction cost and replicated storage explicitly. General spatial/orbital decomposition is a later supported configuration unless validated in this implementation.

### C. Native HSE ground state and dynamics

Enable explicit HSE06 selection with the matching Libxc semilocal remainder,25% short-range exchange and omega0.11 bohr^-1. Add exchange to hpsi, energy reporting, and source-state lifecycle across SCF and RT. Converge the native ground state against the exported reference before comparing spectra.

Integrate self-consistent PT-CN/ACE as a selectable propagation path so the matched comparison can retain dt0.32au and the existing full-exchange endpoint residual gate. Do not assume the existing Taylor propagator gives the same stable step or freeze the source exchange across nonlinear iterations. Check gauge handling of uniform vector potential and the physical current against the reference.

Persist accepted physical states through SALMON restart; rebuild derived ACE/FFT caches after restart. Save only accepted endpoints; reject incompatible restart physics and unsupported spin, NLCC, geometry, occupation or k-mesh cases explicitly. Initial certification is fixed-ion, unpolarized Si8, uniform shifted4³ k mesh, cubic12³ grid and no NLCC. Parameterize dimensions in the operator, while keeping capability guards until other cases are tested.

### D. Acceptance and timing

- Agreement of full exchange action and exchange energy with reference within floating-point tolerance, including different targets.
- MPI-size invariance and correct electron/occupation/k normalization.
- Native ground-state residual/energy agreement and short impulse response agreement in current, energy and norms with equivalent initial state/integrator.
- Energy-gradient and field-free drift checks; full-exchange residuals and norm gates remain enforced.
- Restart/resume equivalence and regression tests for existing semilocal calculations.
- Measure EXX build, ACE build/application, local Hamiltonian, FFT, communication, allocations and complete steps at1/4/8 MPI ranks with controlled BLAS/OpenMP threads. Report memory as well as wall time. No speedup claim until measured.
- Use the same4³ mesh; no k convergence study requested.

The currently running Python MPI linear-response job remains the reference while the native implementation is validated. Replace a production trajectory only after parity, stability and performance checks.

## Pending design decision

Approve staged native integration (A–D) with the first certified target limited to the present Si setup, retaining the full exchange kernel and ACE. Extending geometric/spin support or introducing a physical cutoff is a separate change.
