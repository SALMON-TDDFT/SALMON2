# Overlapping-Wannier Polarization and HHG Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Emit occupation-weighted polarization and current from the accepted Si64 coefficient-RT route and use polarization to validate linear response and laser-driven HHG.

**Architecture:** Evaluate V3 position and velocity expectation values on the root-owned dense coefficient-RT state while preserving the existing row-owned checkpoint boundary. Publish deterministic, restart-safe time series, then analyze polarization with a standalone NumPy script using documented windows and normalization. Production acceptance uses the existing Si64 V3 checkpoint and never enters a removed or conventional RT route.

**Tech Stack:** Fortran 2008, MPI, LAPACK generalized eigensolver, Python 3, NumPy, CMake, ScaLAPACK, EigenExa.

---

Every task below requires an observed RED, focused verification, specification
review, code-quality review, resolution of all Critical/Important findings,
and `git diff --check` before its commit. Use `systematic-debugging` for every
unexpected failure and `verification-before-completion` before success claims.

### Task 1: Define observable and output contracts

**Files:**
- Modify: `tests/dg/test_rt_dg_overlapping_wannier_mpi.f90`
- Modify: `tests/dg/run_rt_dg_overlapping_wannier_mpi.py`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Add RED observable fixtures**

Extend the existing two-state nonorthogonal fixture with occupations
`[2.0d0, 0.5d0]` and volume `4.0d0`. Require a public evaluator returning:

```fortran
call evaluate_dg_overlapping_wannier_observables(comm,coefficients,occupations,&
  volume,state,polarization,current,ok,message)
```

Compare `polarization(axis)` and `current(axis)` with direct dense
occupation-weighted references. Add collective RED cases for inconsistent
occupations, nonpositive/nonfinite volume, wrong occupation count, nonfinite
coefficients, and stale/uninitialized state.

**Step 2: Add RED publication/restart fixtures**

Require `overlapping_wannier_rt_observables.dat` to contain one header and
columns `step time Ex Ey Ez Px Py Pz Jx Jy Jz`. Test:

- step zero plus every propagated step;
- no duplicate boundary sample after restart;
- strictly increasing step/time;
- one-shot versus split byte identity;
- rejection of a stale/truncated/mismatched observable file;
- rank-zero-only publication and collective failure propagation.

**Step 3: Run RED**

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/run_rt_dg_overlapping_wannier_mpi.py
```

Expected: contract/compile failure because the evaluator and deterministic
publisher do not exist.

**Step 4: Review the test specification**

Confirm the references use
`-sum_n occupations(n)*real(dot_product(c_n,matmul(O,c_n)))/volume`, not a
fixed factor two, and that current is secondary information only.

**Step 5: Commit the RED contract**

```bash
git add tests/dg/test_rt_dg_overlapping_wannier_mpi.f90 \
  tests/dg/run_rt_dg_overlapping_wannier_mpi.py \
  tests/dg/check_dg_overlapping_wannier_route.py
git commit -m "test(dg): specify coefficient RT observables"
```

### Task 2: Implement occupation-weighted polarization and current

**Files:**
- Modify: `src/rt/dg/rt_dg_overlapping_wannier.f90`
- Modify: `tests/dg/test_rt_dg_overlapping_wannier_mpi.f90`

**Step 1: Implement the minimal evaluator**

Add `evaluate_dg_overlapping_wannier_observables`. Validate the collective
contract before evaluation. On rank zero, compute for each Cartesian axis:

```fortran
polarization(axis)=0d0
current(axis)=0d0
do state_index=1,size(coefficients,2)
  polarization(axis)=polarization(axis)-occupations(state_index)*real(&
    dot_product(coefficients(:,state_index),&
      matmul(state%position(axis,:,:),coefficients(:,state_index))),real64)/volume
  current(axis)=current(axis)-occupations(state_index)*real(&
    dot_product(coefficients(:,state_index),&
      matmul(state%velocity(axis,:,:),coefficients(:,state_index))),real64)/volume
enddo
```

Broadcast the six values and collective status. Reject nonfinite results.
Do not allocate another persistent dense operator copy.

**Step 2: Verify GREEN on all MPI sizes**

```bash
python3 tests/dg/run_rt_dg_overlapping_wannier_mpi.py
```

Expected: PASS on 1, 2, 4, and 8 ranks.

**Step 3: Run affected route contracts**

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/check_obsolete_dg_routes_removed.py
git diff --check
```

**Step 4: Perform specification and code-quality reviews**

Review occupation/volume/sign conventions, MPI agreement, finiteness,
allocation ownership, and absence of forbidden fallbacks. Resolve all
Critical/Important findings and rerun affected checks.

**Step 5: Commit**

```bash
git add src/rt/dg/rt_dg_overlapping_wannier.f90 \
  tests/dg/test_rt_dg_overlapping_wannier_mpi.f90
git commit -m "feat(dg): evaluate coefficient RT observables"
```

### Task 3: Publish restart-safe production time series

**Files:**
- Modify: `src/rt/main_tddft.f90`
- Modify: `src/rt/dg/rt_dg_overlapping_wannier.f90`
- Modify: `tests/dg/test_rt_dg_overlapping_wannier_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Wire V3 occupations and volume**

Pass `checkpoint%occupations` to the evaluator. Derive cell volume from the
validated SALMON lattice and require a finite positive value. Evaluate and
publish step zero before propagation on a fresh run and publish after every
successful step. The electric field columns must be the exact field passed to
the propagator.

**Step 2: Implement deterministic publication**

Use a fixed scientific format with enough precision for restart byte identity.
The header must include route magic, electronic sign convention, atomic-unit
units, volume, basis/operator/observable fingerprints, and column names.
Before restart append, parse and validate the header and final row against the
loaded restart state. Never silently overwrite incompatible evidence.

**Step 3: Run RED/GREEN restart tests**

```bash
python3 tests/dg/run_rt_dg_overlapping_wannier_mpi.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: publication, corruption rejection, and one-shot/split byte identity
all PASS on 1, 2, 4, and 8 ranks.

**Step 4: Production smoke**

Run 2-step field-off and impulse cases from the accepted Si64 V3 checkpoint.
Require finite `P/J`, distinct impulse and field-off polarization, no forbidden
marker, and a dedicated RT restart.

**Step 5: Review and commit**

Resolve every Critical/Important finding, rerun focused verification and
`git diff --check`, then:

```bash
git add src/rt/main_tddft.f90 src/rt/dg/rt_dg_overlapping_wannier.f90 \
  tests/dg/test_rt_dg_overlapping_wannier_mpi.f90 \
  tests/dg/check_dg_overlapping_wannier_route.py
git commit -m "feat(dg): publish coefficient RT polarization"
```

### Task 4: Implement polarization spectrum analysis

**Files:**
- Create: `tools/analyze_overlapping_wannier_spectra.py`
- Create: `tests/dg/check_overlapping_wannier_spectra.py`

**Step 1: Write synthetic RED tests**

Generate deterministic field-off, impulse-response, and laser polarization
signals in a temporary directory. Require:

- strict header/step/time validation;
- field-off subtraction;
- selectable Hann or exponential damping window;
- FFT frequency and eV axes with documented Nyquist/resolution;
- linear response normalized by integrated impulse field;
- HHG `omega**4 * abs(P_omega)**2` and harmonic-order columns;
- known sinusoid peaks recovered within one frequency bin;
- nonuniform/truncated/nonfinite data rejection;
- deterministic TSV and JSON summary output.

**Step 2: Run RED**

```bash
python3 tests/dg/check_overlapping_wannier_spectra.py
```

Expected: FAIL because the analyzer is missing.

**Step 3: Implement the analyzer**

Use NumPy only. Never infer pulse metadata from desired peak locations; require
carrier energy, polarization axis, window, and impulse amplitude explicitly on
the command line. Emit raw spectral values without hidden normalization and a
JSON record of every analysis parameter.

**Step 4: Verify GREEN and review**

```bash
python3 tests/dg/check_overlapping_wannier_spectra.py
python3 -m py_compile tools/analyze_overlapping_wannier_spectra.py \
  tests/dg/check_overlapping_wannier_spectra.py
git diff --check
```

Review FFT normalization, angular-frequency conversion, window coherent gain,
zero-frequency behavior, and reproducibility. Resolve all Critical/Important
findings.

**Step 5: Commit**

```bash
git add tools/analyze_overlapping_wannier_spectra.py \
  tests/dg/check_overlapping_wannier_spectra.py
git commit -m "feat(dg): analyze polarization spectra"
```

### Task 5: Add the Si64 response/HHG production matrix

**Files:**
- Create: `tests/dg/data/si64_overlapping_wannier_rt/input_fieldoff.in`
- Create: `tests/dg/data/si64_overlapping_wannier_rt/input_impulse.in`
- Create: `tests/dg/data/si64_overlapping_wannier_rt/input_laser_weak.in`
- Create: `tests/dg/data/si64_overlapping_wannier_rt/input_laser_hhg.in`
- Create: `tests/dg/run_si64_overlapping_wannier_response_hhg.py`
- Create: `tests/dg/check_si64_overlapping_wannier_response_hhg.py`

**Step 1: Write the RED matrix checker**

Require immutable run directories and SHA-256 manifests. Check field-off drift,
finite P/J, `J` versus centered `dP/dt`, impulse half-amplitude linearity,
x/y/z cubic equivalence, weak/strong laser intensity dependence, harmonic
peak alignment, forbidden-component suppression, time-step convergence,
Nyquist/frequency resolution, checkpoint/restart provenance, and forbidden
route absence. Tolerances must be declared constants with physical units and
reported values; do not weaken them after observing results.

**Step 2: Define the initial matrix**

Use the accepted 8-rank Si64 V3 checkpoint, one OpenMP thread, 1.55 eV carrier,
approximately 10 fs duration, a weak reference intensity, and a nonlinear
case near `1e12 W/cm^2`. Choose `dt` from a documented convergence target and
include a half-step nonlinear comparison. Record exact pulse envelope, CEP,
propagation duration, post-pulse duration, and polarization direction.

**Step 3: Observe RED**

```bash
python3 tests/dg/check_si64_overlapping_wannier_response_hhg.py <empty-result-root>
```

Expected: FAIL listing every missing production row.

**Step 4: Run field-off and impulse controls**

Run field-off, full/half-amplitude impulses, and x/y/z impulses. Analyze
polarization and run the checker. Stop and diagnose if drift, linearity, cubic
equivalence, or `J`/`dPdt` fails; do not start HHG from an invalid linear basis.

**Step 5: Run weak and nonlinear laser cases**

Run the weak reference, nonlinear case, and time-step-halved comparison.
Analyze spectra from polarization and run the unchanged checker. If the
near-infrared nonlinear spectrum is unresolved, report that result and write a
separate reviewed plan before starting the 0.4 eV/30 fs extension.

**Step 6: Review and commit assets**

Review raw logs and spectra for physical and numerical consistency. Resolve
all Critical/Important implementation findings. Commit inputs/runners/checker,
not bulky raw production outputs:

```bash
git add tests/dg/data/si64_overlapping_wannier_rt \
  tests/dg/run_si64_overlapping_wannier_response_hhg.py \
  tests/dg/check_si64_overlapping_wannier_response_hhg.py
git commit -m "test(dg): add Si64 polarization and HHG gate"
```

### Task 6: Final clean-overlay acceptance and results

**Files:**
- Create: `docs/plans/2026-08-01-overlapping-wannier-polarization-hhg-results.md`

**Step 1: Run the complete source/focused suite**

Run removal and retained-route contracts, spectrum synthetic tests, every
overlapping-Wannier metadata/projection/metric/construction/operator/
physical-matrix/solver/SCF/checkpoint/coefficient-RT runner on 1, 2, 4, and 8
ranks, normal LCFO/EigenExa contracts, and `git diff --check`.

**Step 2: Run the clean-first parent-prerequisite overlay build**

Create fresh source from `git archive HEAD`, apply the current task diff, and
overlay only enumerated retained parent prerequisites. Explicitly exclude dirty
or removed WPW/Fragment/Nodal sources. Configure MPI, ScaLAPACK, and EigenExa,
then run:

```bash
cmake --build <overlay-build> --clean-first -j1
```

Expected: `[100%] Built target salmon` with no removed source referenced.

**Step 3: Re-run production acceptance with the overlay executable**

Repeat the required field-off, impulse linearity/cubic matrix, weak laser,
nonlinear laser, restart split, and time-step comparison. Require the unchanged
public checker to PASS.

**Step 4: Perform final reviews**

Review the complete branch against the design. Resolve every Critical and
Important finding and repeat all affected verification. Record unresolved
physical limitations without relabeling them as code-quality findings.

**Step 5: Record results and commit**

Record commands, source/binary/checkpoint/output hashes, runtimes, memory,
field-off drift, linearity and cubic-equivalence errors, `J/dPdt` error,
frequency resolution/Nyquist limit, harmonic locations/intensities, time-step
convergence, forbidden-marker scan, and review disposition.

```bash
git add docs/plans/2026-08-01-overlapping-wannier-polarization-hhg-results.md
git commit -m "test(dg): validate polarization and HHG physics"
```
