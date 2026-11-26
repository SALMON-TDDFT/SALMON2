# Distinct field/electronic steps: coupling accuracy — 2026-11-26

## What was tested

The user's correction distinguishes the field timestep from the electronic
one: advance light several times during one electronic update. Initial runs
with equal dt were retained only as a baseline. They cannot isolate the errors
of the two solvers.

Added `salmonpy.coupled.FieldSubcycledStepper` in the existing coupling module.
It requires electronic_dt=m*field_dt, a positive integer m. During m native
Maxwell updates it holds the previous accepted physical charge-current density,
then sends final A/c and receives the new electronic endpoint current once.
This is an explicit hold, not a future-current interpolation or implicit solve.
No application-wide communication was added for each light update. Field-local
LD keeps its native fine-step update; the test below has no LD.

Counters are explicit: after electronic update n, field update is n*m.
The electron's local zero maps to field_dt/2, so the accepted electronic state
and current have physical time n*electronic_dt+field_dt/2. Its propagation
midpoint is(n-1/2)*electronic_dt+field_dt/2, which need not coincide with stored
E. Kinetic/nonlocal phases still use the same midpoint A, computed by the
existing electronic propagator from its endpoint A values. No native physics
was changed or a new solver binary built for this schedule.

Composite failures invalidate the response; this new interface does not claim
whole-cycle protected rollback, unequal-step disk restart, arbitrary ratios or
molecular E/H support. Existing equal-step CoupledStepper remains unchanged.

## Fixed calculation and records

Same existing Si PZ GS, no SCF or GS copies:16^3 electron grid,32electrons,
32states,10.26bohr cube, original SALMON2x2x2 full k mesh. Maxwell40^3 grid,
2bohr spacing, one electronic macro point,8bohr^3 macro volume, x-PML and
transverse periodic boundaries. A closed15^3-cell observation box excludes
source/PML and includes all six faces, with the current mapped across the
field-rank partition. Two field ranks plus two electronic k ranks; one thread
per rank. All three cases end at field time16.002au, about0.387fs.

Identical incident pulse: amplitude.001au, omega.5au, duration8au, z polarization,
plane normal x at x=-16bohr. This is a short weak-drive time-discretization
experiment, not a converged optical spectrum or long-time stability campaign.

| Case | Electronic dt | Field dt | Electronic updates | Field updates |
|---|---:|---:|---:|---:|
| A | .003 | .0015 | 5334 | 10668 |
| B: field-only refinement | .003 | .00075 | 5334 | 21336 |
| C: electronic-only refinement from B | .0015 | .00075 | 10668 | 21336 |

Use manual `tools/test_coupled_timestep_mpi.py`, not required physics CI.
--dt/--steps are ELECTRONIC; --field-steps selects m. The tool uses actual shared
TDDFT state handles and references GS in place through symlinks. Existing GS
files match their original2026-11-23 ground-state manifest hashes after the runs.
The unchanged audited library was checksum-verified on reuse:
`1efb6bd8289a56b319eef15ec01ec217c1e23e16a7ba35e24d5aa3960931bc73`.

## Physical-time alignment and independent balance checks

The Fortran analysis module `tests/native/test_coupled_timestep.f90` explicitly
interpolates reference J at n*electronic_dt+field_dt/2 to the comparison times;
E is subsampled at exact common field times. No extrapolation or silent update
interpolation. Electronic energy is interpolated between the two states
bracketing the common final field time. Electron-cell energy is multiplied by
macro_volume/electronic_cell_volume. Initial GS is assigned field_dt/2 before
pulse arrival; this is not a generic nonstationary initial-state alignment test.

Native Maxwell observations retain field energy change, actual injected-current
work and six-face outward flux. Both discrete-update field energy and ordinary
E/H-centered energy are reported. Separately compute electronic energy gain and
midpoint work from the average field in each electronic interval. Never replace
an actual TDDFT energy increase by the field work to force the balance to close.

The Fortran analysis also independently sums point work using the actual held
current at each light update; it agrees with the native scattered-edge work.
The difference between electronic-midpoint work and held-current field work
splits into an endpoint term and an explicit hold term:

    W_mid - W_point = h*V/2 * E_last·J_last
      + h*V/2 * sum_q[(sum_light E - E_last_in_q)·(J_q-J_previous)]

where h=field_dt and initial E=0. The second term vanishes for m=1 but not for
m>1. This identity is checked, not assumed when reporting the energy defect.
Finite, NaN/Inf, shifted linear-clock interpolation and known held-current work
self-tests pass. New analysis-library source/object prohibition audits pass.

## Measured results

| Refinement | Current relative RMS difference | FDTD absorption-work relative difference | Electronic energy-gain relative difference |
|---|---:|---:|---:|
| A -> B: light only | 0.0030064% | 0.265354% | 0.0030741% |
| B -> C: electrons only | 0.0062498% | 0.530555% | 0.0224931% |
| A -> C: both halved, m=2 | 0.0037799% | 0.266609% | 0.0194183% |

All differences use the finer case as reference. Current comparisons use all
three components and their documented common interval; they do not divide by
individual zero-crossing samples. They are well below the5% screening target,
which is advisory rather than a reason to stop valid calculations.

| Case | Field-transferred work / Ha | Actual electron energy gain at common end / Ha | Balance defect / absorbed work |
|---|---:|---:|---:|
| A | 8.16885915e-8 | 8.21572252e-8 | 0.573683% |
| B | 8.14724015e-8 | 8.21597508e-8 | 0.843659% |
| C | 8.19069627e-8 | 8.21412747e-8 | 0.286071% |

The last column is(actual electronic energy gain - actual field current work)
/abs(field current work). The full physical balance additionally includes stored
field energy and six-face flux and is recorded in JSON; the difference is small
but is not omitted. Native field-discrete closure is near floating-point error.

Light-only refinement does NOT improve the material-normalized closure: the held
current becomes the limiting coupling approximation. A->C halves the defect
approximately; two levels do not establish a formal first/second order. The
TDDFT energy-vs-midpoint-work term also halves from~0.0381% to~0.0190%. Do not
claim that every residual is caused by the Yee half-step placement or that field
refinement alone converges the coupled solution. Prediction/correction would be
a separate scheme requiring its own validation.

No NaN/Inf or numerical failure occurred. Measured electron loss stays around
1e-14, far below1e-5; gain is reported too. Late current amplitude remains finite
and consistent across the tests. This confirms short-interval stability only.
No statement about many-cycle, strong-field, many-point or other-XC stability.

## Distribution and remaining work

Japanese/English guides and AI JSON document both clocks, explicit current hold,
manual commands and the different meanings of field-only and total balance.
1:1 scheduler regression retains original CoupledStepper output arrays exactly.
Targeted source tests pass; extracted-sdist wheel/installed tests are recorded in
this campaign's distribution logs. Heavy physics CI is neither added nor
manually dispatched.

Canonical local output: calc2026/salmonpy-coupled-timestep-20261126.
`field2-qdt`, `field4-qdt`, `field2-half-qdt` are A/B/C; field-only-comparison.json,
electronic-only-comparison.json and both-steps-comparison.json contain the native
analysis. Earlier dt/half-dt and report.json are explicitly equal-dt baseline
records, not independent-rate evidence. Native energy/RT output, accepted-step
histories, norm samples, input files, analysis-library identity, source snapshots
and logs are retained. Cloud synchronization is not confirmed.

Unequal-step whole-cycle direct persistence/rollback and normal mixed-sample
checkpoint integration remain unvalidated. Molecular feedback stays deferred.

## 10:1 follow-up and user-typical steps

D: electronic .003au, field .0003au,5334/53340updates to the same16.002au.
A->D gives current RMS4.80917e-5, electronic-energy difference4.96796e-5,
field-work difference0.425202%, and material-normalized defect1.00634%.
Norm loss7.78e-15; native point/scattered work difference1.33e-23Ha.
This supports the observation that a much finer field step can leave current
and electronic energy almost unchanged, while current hold affects work closure.
No equality with the original A-based Maxwell discretization is asserted.

The user identifies electronic .2au / field .02au as typical. A separate pilot
uses20bohr field spacing so .02au passes CFL, with the same16^3 Si electronic
grid(~.34angstrom) and native default Taylor4. It becomes electronically
unstable: measured Ne reaches32.3883 at update8,137.636 at9 and28723 at10,
then the next numerical update fails with recoverable electronic status21.
A zero-drive control exhibits the same runaway, separating this issue from
field-subcycle coupling. These failed trials are retained as typical-pilot and
typical-dark, with no success record. They do not establish that .2au is invalid
for other grids or propagators. Matching the user's typical electronic grid and
propagation order is pending; the small-dt study is not a substitute.


### Measured electronic0.08 au / field0.02 au case

The fixed Si fixture used fourth-order Taylor propagation, two electronic MPI ranks and two FDTD ranks to field time20 au. Electronic/field steps were0.08/0.02,0.08/0.01 and0.04/0.01 au. Comparing the first and last cases gives a relative current RMS difference of0.148%, a TDDFT energy-gain difference of0.484%, and maximum electron-number loss below1.12e-11. This tests short-run stability and current agreement, not material convergence.

The field absorption-work difference is12.1%. The mismatch between field work and electronic energy gain is28.8% of absorption work, decreasing to12.7% with both steps halved. This is an auxiliary comparison against internal TDDFT energy, not an electromagnetic conservation error evaluated with J·E. Judge absorption convergence from the timestep dependence of J·E itself. Maxwell spacing is20 bohr and electronic spacing approximately0.34 angstrom. These steps are not universal recommendations.


### Field-step refinement with electronic dt fixed (2026-11-26)

Electronic dt was **0.04 au**, with fourth-order Taylor propagation. Field dt was0.01,0.005 and0.0025 au, giving4,8 and16 field updates per electronic update. All runs reached field time20 au (500 electronic updates), using two FDTD ranks and two electronic k-point ranks. The existing Si/PZ GS was referenced without copying or recalculating it:16³ electronic grid,2×2×2 k points,40³ Maxwell grid with20 bohr spacing, one weakly driven material point.

| Field dt (au) | Field updates / electronic update | Absorption (Ha, including macro-point volume) |
|---|---:|---:|
| 0.01 | 4 | 7.27652472e-5 |
| 0.005 | 8 | 7.13166003e-5 |
| 0.0025 | 16 | 7.05918506e-5 |

The primary absorption is **Σ dt V J·E**, using the actual charge current supplied to FDTD and E averaged before/after each field update. The previous accepted current is held between electronic updates. Work computed instead from averaged electronic-endpoint currents, and changes in internal TDDFT energy, are auxiliary diagnostics; neither replaces this absorption.

Adjacent relative current RMS differences decrease from **0.0231% to0.0114%**. Relative absorption differences decrease from **2.03% to1.03%**, normalized by the finer result. Currents are compared at common physical times by explicit Fortran interpolation. Electron-number loss stays below3.50e-13. The absolute discrete electromagnetic energy-balance residual using J·E stays below3.00e-14 Ha.

Both adjacent differences satisfy the5% screening target for this short run. Electronic dt, spatial grids and k sampling remain fixed; their convergence, strong excitation and long-time accuracy are not certified. Internal TDDFT energy differences are not automatically electromagnetic conservation errors. Records: `qdt004-light-scan.json` and its two comparison JSON files.
