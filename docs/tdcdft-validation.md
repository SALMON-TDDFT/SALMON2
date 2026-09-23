# TDCDFT implementation validation — 2026-09-23

Repository: SALMON-TDDFT/SALMON2, branch `TDCDFT`, base `ae7d8953` (`origin/develop-2.0.0`).
This is an implementation/smoke-test record, not a converged Si exciton calculation.

## Implemented

Opt-in PZ + head-only LRC and normalized Proca-type field equation; midpoint propagation;
consistent endpoint kinetic/nonlocal current; initial half step; separate xc output;
versioned shared-file RT restart with current history; parameter/compatibility validation.
The classical Ac_tot history and legacy rtdata.bin layout remain unchanged.
GPU and other unsupported physics/propagation paths are rejected explicitly.

The new model is documented in `samples/exercise_si_tdcdft/README.md`.
SALMON-DOCS is a separate repository: `docs/salmon-docs-tdcdft.patch` contains the matching
change to `source/input_keyword_list.rst`, made against a separate checkout.
No remote repository was modified.

## Build environment

macOS ARM64, Homebrew gfortran/GCC 15, OpenBLAS, CMake; serial and Open MPI builds.
An unmodified base executable was built before source changes for baseline comparison.
Legacy source compatibility needed `-fallow-argument-mismatch` (serial communication stubs)
and `-include stdio.h` (posix.c uses snprintf without the corresponding include).
These were build flags, not unrelated source edits.

Serial configuration used:

```sh
cmake -S . -B /private/tmp/salmon-tdcdft-work \
  -DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc-15 \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran \
  -DCMAKE_BUILD_TYPE=Release -DUSE_MPI=OFF \
  -DCMAKE_Fortran_FLAGS=-fallow-argument-mismatch \
  '-DCMAKE_C_FLAGS=-include stdio.h' \
  -DBLA_VENDOR=OpenBLAS -DCMAKE_PREFIX_PATH=/opt/homebrew/opt/openblas
cmake --build /private/tmp/salmon-tdcdft-work -j 8
```

MPI used `USE_MPI=ON`, `mpifort`, and `OMPI_CC=/opt/homebrew/bin/gcc-15` for both
configuration and build. The latter avoids the Apple Clang wrapper rejecting OpenMP.
Local MPI execution needed permission for process-communication sockets outside the sandbox.

## Tests and evidence

- Unmodified baseline: tests 111 (Si GS) and 112 (Si response), six CTest steps passed.
- Modified serial and four-process MPI: tests 111, 112, 142, analytic field test;
  ten CTest steps passed per build.
- Analytic field update: constant drive, zero coupling, smooth sinusoidal drive,
  damped oscillator with second-order timestep convergence.
- Si integration regression: disabled output matches unmodified executable;
  enabled alpha=0 impulse output is identical to disabled; LRC and Proca feedback are nonzero.
- Restart at step 90 of 180 reproduces current and the full response spectrum.
  Missing/incompatible restart state is guarded by the reader; changed-parameter rejection tested.
- Weak impulse 0.001 vs 0.0005: normalized-current relative difference
  approximately 8.11e-5. dt=0.16 vs 0.08: relative current difference approximately 3.54e-4.
  These short, coarse-grid checks demonstrate numerical consistency, not physical convergence.
- Laser alpha=0 vs alpha=1e-10: continuous current limit verified.
- Negative/NaN coefficients, unknown mode, unsupported propagator and changed restart coefficient rejected.
- Full integration regression also passed with two MPI processes.
- Analysis tests: damped-oscillator peak and amplitude scaling, pump subtraction and time-origin shift,
  rejection of inconsistent time grids/zero impulse/empty peak intervals (three Python tests).
- Python Fourier transform matches SALMON's complex epsilon on the Si test within rtol=1e-10, atol=1e-9.
- Six input templates exercised with shortened 180-step runs; pump/probe timing shortened to include the
  probe. Differential current vanishes before the probe, is nonzero afterward, and gives finite spectra.
- Independent read-only code review found a missing initial half step and a laser zero-coupling
  discontinuity; both fixed and covered by tests. Re-review found no further definite correctness issue.

CTest fixture setup originally omitted the execution fixture used by verification; adding the new
GS-dependent test exposed out-of-order verification. The shared helper now registers that fixture,
so the new RT test actually depends on completed GS execution. No reference tolerances were weakened.

Representative commands (paths can be changed):

```sh
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 ctest --test-dir BUILD \
  -R '111_bulk_Si|112_bulk_Si|142_Si|unit_tdcdft' --output-on-failure
python3 testsuites/142_Si_tdcdft/test_update.py
python3 testsuites/142_Si_tdcdft/test_regression.py EXE BASELINE_GS_DIR BASELINE_RT_DATA
# MPI variant appends: /path/to/mpiexec -n 2
python3 samples/exercise_si_tdcdft/test_analyze.py
```

The legacy verification scripts require a `python` command with NumPy; locally this was supplied
by a temporary shim to python3. The new field-unit test only needs Python and gfortran.

## Limits and next scientific work

No fitted Si Proca parameter set, converged absorption curve, E1/E2 shift, bleaching, or exciton
binding energy has been established. The original sample gamma=0.01 a.u. was a numerical demonstration; see the reference update below.
The small 12³ real-space / 2³ k grid and 180-step tests cannot support such claims.
At this initial test stage, full-duration and grid/dt/time checks remained; see subsequent records below.

For lasers, enabled alpha=0 is the comparison baseline: the new path uses endpoint nonlocal phases,
whereas the legacy path uses midpoint phases in that part of the current. Impulse comparisons are
identical because the classical A is constant. The model-field energy is not added to total-energy
output. Zero-field model stability does not guarantee stability of the coupled electronic system.

Before interpreting pump-induced peaks, verify a no-pump delayed probe against equilibrium response,
probe-amplitude convergence, pump-only subtraction, observation-window convergence and coefficient
sensitivity. Fixed alpha does not describe dynamically evolving screening. Test against literature
or BSE with consistent quasiparticle-gap treatment before drawing quantitative conclusions.


## Follow-up: fixed 4³ k grid

The user removed the k-point convergence requirement. The examples now use 4³ k points,
and completed 12000-step ALDA/LRC runs and their initial analysis are recorded in
`docs/results/si-k4/README.md`. The earlier 2³-grid entries above remain the implementation-test record.

## Williams–Ullrich reference update

Primary reference: https://arxiv.org/html/2501.13290v2, Eq. (26), Sections IV.4 and IV.5.
The existing normalized integrator already propagates this field equation. Si samples now use
alpha=0.2, beta=0 and gamma=1e-4 a.u.; gamma is an unvalidated sensitivity trial, not a Si fit.
Its free-oscillator energy scale is 0.272 eV, compared with 2.72 eV for the old gamma=0.01 example.
No full-duration result for the new gamma is reported. The earlier Dewhurst full run was stopped
when the user changed the reference; its partial output is retained separately.

The optional direct a2/a0 conversion remains available. Serial and MPI builds passed; ten
serial CTest checks passed. Extended serial and two-process MPI regressions passed for direct/
normalized equivalence, old normalized restart, massless limit, both coefficient signs, invalid
inputs, and atomic-unit versus A_eV_fs equivalence. No checkpoint format change was needed.
The unit-test module include directory is prioritized to avoid stale module interfaces in parallel builds.

The paper motivates stabilization through time-averaged forces, not cancellation of the xc force
at each instant. No automatic gamma prescription or spatial counter field has been implemented.
The all-electron 2D model's numerical thresholds are not transferred to pseudopotential Si.
Before pump interpretation, compare gamma=0, 1e-4, 4e-4 at fixed alpha and fixed 4³ k grid,
including observation-window dependence, spurious low-frequency response and secular field drift.

## Instantaneous screening implementation

Opt-in `tdcdft_screening='instant'`: instantaneous vector A,J,E,P estimator, scalar positive
screening closure, field-floor hold, trapezoidal P, 12-column diagnostics and version-2 restart.
No temporal averaging. Fixed mode and strength=0 retain the existing propagation.

Analytic tests cover both quadrature zero crossings, calibrated Drude reference, positive bounded
correction, zero-strength limit, field-off hold and pulse tails with remanent polarization.
Serial CTest: 10/10 passed. Extended serial and 2-process MPI regressions passed for active feedback,
zero-strength equality, active restart, state/parameter guards and existing baseline checks.
Read-only independent review found no definite indexing/restart defect; it identified pulse-tail
contamination of the estimator, which was reproduced, tested and documented.

Six short 4³-k laser runs and the field-threshold sensitivity are recorded in
`docs/results/si-instant-screening/README.md`. They demonstrate numerical operation, not improved
physical screening. In particular, post-pulse alpha depends on the field threshold.

The same-pulse long extension is recorded in `docs/results/si-instant-screening/LONG_RUN.md`:
fixed alpha fails normalization near 5.85 fs; instantaneous mode reaches 23.22 fs with electron-count
error below 3.5e-7 but exhibits residual XC electric field and Axc drift. Initial 1200-step output
prefixes match exactly. The post-pulse external fields are zero. The residual field is largely a
constant E_xc-alpha*P offset from the time variation of alpha; no constitutive-law change was made.

## Polarization-consistent closure

Added `tdcdft_screening='polarization'` while retaining the earlier 'instant' acceleration closure.
The new mode integrates a'=-alpha P at beta=gamma=0 using a second-order explicit Taylor step.
Smooth variable-alpha analytic convergence, fixed-alpha equivalence, and removal of the post-switch
field offset pass. Serial CTest 10/10 and full serial/2-process MPI regression suites pass, including
coupled fixed-alpha limit, active restart and incompatible closure/restoring guards. Independent
read-only review found no blocking sign/order/restart problem.

The same 4³-k pulse completes 23.22 fs; E_xc-alpha P after the pulse is below 4e-17 a.u.,
compared with 3.90e-4 in the former closure. Half-dt comparison through 5.805 fs gives 0.055%
relative current difference but a factor-of-two difference in held alpha. Thus field consistency
is fixed, but the estimator remains threshold/time-step sensitive. Detailed results are in
`docs/results/si-polarization-screening/README.md`. No k-grid scan or precision claim was added.

## K0 and estimator-gate audit

`docs/results/si-k0-audit/README.md` records an amplitude-homogeneity counterexample:
scaling the weak A,E,J,P data uniformly by 316.23 leaves K unchanged but moves the absolute-floor
cutoff closer to the pulse end; the old closure then reduces alpha from 0.2 to 2.77e-5 without
any new nonlinear trajectory. A common relative floor removes this scaling artifact offline.
Optical-response and field-weighted K0 candidates instead alter even the weak response.
No new calibrated K0 or production model change was adopted from this audit. Previous numerical
stabilization must not be taken as independent evidence of physical carrier screening.

### alpha0=1, K0=0 trial

See [bare-K trial](results/si-bare-k/README.md). The existing polarization closure completed the strong-pulse trajectory through 23.22 fs with alpha0=1 and reference=0. Residual A and mean current decreased relative to the previous parameters, but the terminal alpha changed by about 103 times when matching the weak/strong relative field gates. This is a numerical comparison, not a validated carrier-screening calibration. No production code changes.

### Frozen-screening pump-probe support

Added `tdcdft_screen_stop` (default -1): freeze K/alpha after an absolute atomic-unit time;
continue P, current and XC field. Checkpoint v3 adds the cutoff, accepts earlier v1/v2
with cutoff disabled, and rejects changed cutoff on restart. Dedicated serial and 2-rank
MPI integration checks pass for frozen K/alpha, evolving P, pre-cutoff equivalence,
restart, parameter mismatch and v2 compatibility. Full serial TDCDFT regression and
three Fourier-analysis tests pass. Read-only code review found no blocking issues.

The full-feedback alpha0=1/K0=0 probe diagnostic returned alpha to 1 after the kick;
normalized differential currents differed by 50.6% on halving the probe. Production
spectra therefore use frozen pump screening and are not the full model's linear response.

Frozen bare-K pump-probe results are in [the result record](results/si-bare-k-pump-probe/README.md).
Probe80 au, observation880 au, same4³ grid. At probe1e-4 and frozen alpha=.000206291,
pumped vs unexcited same-alpha peak3.53 vs3.54 eV; peak height -33.16%, signed2–4 eV
integral -27.47%. Probe5e-4 vs1e-4 spectral L2 difference4.19% in2–4 eV, but full time
current difference48.69%; no claim of globally converged linear response. Ground alpha1
frozen propagation failed norm atstep3064 and is excluded. The earlier alpha.2 comparison
includes both changed screening and excitation. See record for common-window and amplitude checks.


### Time-dependent occupied Wannier diagnostics

[Result record](results/si-time-wannier/README.md) adds offline checkpoint readers,
projected bond-centered Wannier controls and discrete Marzari–Vanderbilt localization.
The MLWF analysis transports the preceding U into the current Bloch basis as an initial guess every10 RT steps; it does
not yet change the production alpha cadence or define a Wannier-to-alpha feedback.
Nine numerical tests pass, including normalization/density, gauge covariance,
periodic boundary phases, analytic spread gradients and invariant/variable spread
separation. Independent review confirmed the Fourier normalization, boundary phases,
gradient and gauge-dependent Armijo objective. No production Fortran changed.

The same 4³-mesh strong pulse projects about1.79 electrons/cell into the finite
reference conduction window by3.096 fs (1.84 holes/cell). This contradicts a
small-excitation assumption but does not validate the original screening estimator.

All3 x161 MLWF snapshots converge at gradient norm <1e-6, preserving orbital labels.
Final mean MV spreads are2.057679,2.057711,4.936765 Å² for none/weak/strong.
No-pump polar transport eliminates additional localization iterations; a final strong
initial-gauge start agrees with the sequential result within2.4e-13 Å².

### Unclipped UEG-distance alpha diagnostic

[UEG-distance record](results/si-ueg-distance/README.md) evaluates alpha=.2D with no upper
clip on all3 x161 existing snapshots. Seven tests and independent review verify FFT
normalization, full momentum coverage, HS distances, ensemble/current conservation,
unitary invariance and the fixed-rank bound. No production feedback changes.
Weak final alpha=.1999489; strong range .1877118–.2185561, final .2024486. D>1 occurs
with larger MLWF and invariant spreads, so D is not a pure localization measure.
Reference purity/occupation changes matter. A rank16-per-k constraint gives a positive
rest-reference alpha lower bound .013298: the abstract zero endpoint is not reachable
within this finite-grid trajectory class. Physical screening validation remains open.


### Time-dependent ELF

[ELF diagnostics](results/si-tdelf/README.md) evaluate the current-corrected one-spin
ELF on483 stored snapshots. Five analytic/unitary/boost/derivative tests pass;
independent review confirmed weights, Bloch gradients and summed local current.
The symmetric even-grid derivative is checked against signed FFT Nyquist derivatives.
The finite-mesh cold electron gas yields ELF=.500105. Strong-pump mean ELF changes
.599693→.574107, bond-region ELF .925046→.893902; weak changes are negligible.
High-ELF bond regions remain. These are localization diagnostics, not proof of a
metallic response or a new alpha law. Existing write_elf was not used: its current
orbital/k sums do not implement the required periodic current-corrected expression.

### Native ELF-dependent alpha feedback

[ELF feedback](results/si-elf-feedback/README.md) adds opt-in `tdcdft_screening='elf'`,
with Q=<n(ELF-.5)^2>/<n> and alpha=alpha0*Q/Q_initial, no upper clipping.
Occupation/k-weighted, current-corrected ELF uses native finite-difference Bloch
gradients and is sampled every10 steps by default. Held alpha drives the
polarization closure; no external-field gate or temporal averaging is used.
Version4 checkpoints retain the initial normalization and stride. A pre-existing
odd-step restart parity error was exposed and fixed in the main RT loop.

Same4³ Si self-consistent runs with alpha0=.2 give strong alpha minimum.167689,
final.178755 at3.09617 fs. Weak-pulse final alpha=.199999716; no-pump change<1.72e-10.
Small-impulse current differs from fixed alpha by7.28e-8 in relative L2 norm.
Every-step ELF changes strong current by0.1303% versus stride10. Independent
checkpoint postprocessing matches native Q within3e-14. The cold finite-mesh
electron gas leaves alpha3.42e-8 versus the exact continuum endpoint zero.
This is an empirical localization model; residual current remains and neither
long-time stability nor improved exciton peaks is established by this short run.

### ELF-Proca connection for pump-probe stabilization

The user identified beta/gamma stabilization as necessary after extending the
beta=gamma=0 ELF experiment. `tdcdft='proca',tdcdft_screening='elf'` now routes
the measured alpha(t) into the existing oscillator equation
`a''+beta*a'+gamma*a=alpha(t)*J`; `lrc+elf` retains `a'=-alpha*P`.
These are distinct variable-alpha closures, explicitly documented. No derivative
of alpha times P is added to Proca. The normalized coefficient input is required;
direct a2/a0 is rejected with ELF. Existing v4 mode/parameter metadata rejects
incompatible restarts without changing the format.

Serial and MPI4 ELF-Proca tests verify the variable-alpha recurrence, odd restart,
normalization, freeze and parameter guards. Legacy ELF tests, ten selected CTest
checks and the existing full regression suite pass. The failed unstabilized
control and completed common-coefficient stabilized comparison are documented in
[ELF pump-probe](results/si-elf-pump-probe/README.md) and its Proca subdirectory.

With beta=0,gamma=.001, all pump and central +/-probe and half-probe trajectories
reach23.2213 fs with maximum logged norm error1.142e-5 electrons. The2–4 eV
central spectra differ by0.2122% on halving the probe (unwindowed current3.235%).
However, large positive/negative structures and observation-window dependence
prevent identifying a converged exciton peak change. The same-coefficient
unpumped maximum is3.34 eV versus3.36 without stabilization. A pumped local
maximum at3.39 eV is only a diagnostic extremum, not an established exciton shift.
