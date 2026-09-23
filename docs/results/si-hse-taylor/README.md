# HSE06: SALMON Taylor4 versus PT-CN

The two integrators approach the same current response as the time step is
reduced. This is a time-discretization test on the same physical HSE model,
not another comparison between two implementations of PT-CN.

## Conditions and method

Si8, 32 electrons, 16 occupied orbitals, 12³ real-space grid, full 4³ k mesh,
HSE06 (25% screened exchange, omega=0.11 bohr^-1), identical converged initial
wavefunctions, transverse z impulse of 10^-4 au, fixed ions and occupations.
No k convergence scan or orbital renormalization. Taylor calls the original
`src/rt/taylor.f90` unchanged, with `n_hamil=4` and SALMON's predictor/corrector.
The local and Fock operators are both averaged between the initial and
predicted endpoint. Sources remain distinct from the powers of H applied to
trial vectors. The default semilocal SALMON propagation path is unchanged.

`hse_taylor4` averages the two ACE operators. `hse_taylor4_full` applies the
full screened exchange to every trial vector and averages the actual source
density matrices; it is an independent check of the ACE approximation.
Concatenation with a factor of 1/sqrt(2) implements both operator averages.
Fourth order describes the Taylor polynomial of the exponential. For the
self-consistent time-dependent Hamiltonian, this predictor/corrector and
PT-CN both have second-order global time accuracy.

## Time-step convergence over 8 au (0.1935 fs)

Reference: Taylor4 at dt=0.04 au. Currents are compared at the same 25 times
in every row; no interpolation. RMS error is the norm of the vector-current
difference divided by the reference vector-current norm over those samples.
Density error is the endpoint density difference divided by the norm of the
reference *induced* density, rho(t)-rho(0), rather than the much larger ground
state density. `convergence.json` also records absolute errors, total-density
errors, electron count, overlap errors and optimally gauge-aligned orbital
distances. Raw orbital differences across the two gauges are not used.

| Method | dt (au) | Current RMS difference | Induced-density difference |
|---|---:|---:|---:|
| Taylor4 + PC | 0.16 | 0.01057% | 0.04232% |
| Taylor4 + PC | 0.08 | 0.00215% | 0.00911% |
| Taylor4 + PC | 0.04 | reference | reference |
| PT-CN | 0.32 | 0.25421% | 1.72768% |
| PT-CN | 0.16 | 0.06530% | 0.62835% |
| PT-CN | 0.08 | 0.01656% | 0.17050% |
| PT-CN | 0.04 | 0.00432% | 0.05699% |

Successive differences give an observed current convergence order of 1.991
for Taylor4 + PC and 1.973–1.999 for PT-CN. Second-order Richardson
extrapolation from dt=0.08/0.04 gives a current RMS difference between the
two methods of 0.000580% and an induced-density difference of 0.0470%.
This is an extrapolation estimate, not an exact solution or a proof that
all high-frequency density components are already in their asymptotic regime.

The original dt=0.32 PT-CN result is therefore close in current, but is not
identical in accuracy to a finely stepped Taylor calculation. Tight nonlinear
solver residuals alone do not remove this integration error. For a more
accurate point-by-point current/density comparison, dt=0.08 is preferable to
0.32 within the conditions tested here.

## Longer interval: 32 au (0.7740 fs)

Reference: Taylor4 + ACE at dt=0.08. All three completed runs are compared
at 100 common times; endpoint density is evaluated at exactly 32 au.

| Method | dt (au) | Current RMS difference | Induced-density difference |
|---|---:|---:|---:|
| Taylor4 + PC | 0.08 | reference | reference |
| PT-CN | 0.32 | 0.47298% | 4.31002% |
| PT-CN | 0.16 | 0.14933% | 2.92640% |

The current remains close and improves on refinement. The induced-density
error is larger and does not yet show a clean factor-of-four reduction over
this longer interval; equivalent density accuracy is not established at these
coarse PT-CN steps. No long-interval dt=0.08 PT-CN test was performed here.
Electron-number errors are below 7.4e-10 electrons and maximum overlap errors
below 3.8e-10. Post-impulse energy ranges are below 6.8e-12 Ha.
The reused PT-CN dt=0.32 stability run has an exactly identical current prefix
to the fresh short run. See `longer.json` and `auxiliary_checks.json`.

## Exchange, stability and restart checks

Over the first four dt=0.08 steps, ACE Taylor versus full-Fock Taylor has a
current RMS difference of 0.000099%, maximum current difference 4.44e-12 au,
and induced-density difference 0.00303%. This initial audit is separate from
the time-step convergence test; it does not by itself bound long-time error.
An extended full-Fock audit over 8 au at dt=0.16 gives an ACE/full current RMS
difference of 0.007109% and induced-density difference of 0.003575%.
This checks all Taylor trial vectors against the full operator, with identical
dt and initial state (`full_exchange_audit.json`). The time-step extrapolation
above compares the two ACE-based propagation paths; it does not remove ACE
error or establish an exact full-Fock continuum limit.

The auxiliary Taylor4 + ACE dt=0.32 trial is unstable on this grid: by step8
the total energy is hundreds of Ha above the initial state and the run stops
at step9 on the existing electron-count check. It is excluded from the
accuracy table. Taylor dt=0.16/0.08/0.04 all complete. The PT-CN dt=0.32 run
completes with its full-exchange residual and norm gates active.

Taylor four uninterrupted steps versus two plus restart plus two are bitwise
identical. Wrong Taylor order, missing predictor/corrector, and a restart
which changes between ACE and full-Fock propagator names are rejected.
The standalone midpoint-operator tests and the orbital-gauge comparison test
pass, along with the full 79-test numerical suite and six existing Si/TDCDFT
CTest checks. HSE-enabled and HSE-disabled builds both succeed. Independent
source review found no blocker.

## Reproduction

Build and initialize Si following `samples/hse_native/README.md`. Use one
fixed GS checkpoint for every run. Select `hse_ptcn`, `hse_taylor4`, or
`hse_taylor4_full` in the propagation namelist; Taylor modes require
`n_hamil=4` and `yn_predictor_corrector='y'`. Choose nt=8/dt for the short
convergence series, checkpoint_interval=nt, and out_rt_energy_step=1.
Keep dt, field, geometry and occupations unchanged for a genuine restart.

`samples/hse_mlwf_reference/compare_native_propagators.py` compares completed
native runs against a selected reference and checks that the saved checkpoint
is at the current endpoint. `--export` supplies the common initial native GS
export (created with SALMON_HSE_REFERENCE_EXPORT), `--reference` the reference
run directory, `--output` a JSON filename, followed by the run directories.
The working runs are under ignored `calculations/si_hse_native/compare_*`.
Exact inputs are archived in `inputs/` with hashes in `input_manifest.json`;
restart paths refer to the common local GS fixture. `traces.npz` preserves
current and energy outputs used in the comparisons. `comparison.png` and
`comparison.pdf` summarize the current results. Runs shared the machine with
another reference job (and one batch was paused); their wall times must not
be treated as clean performance benchmarks.

These sub-femtosecond tests assess propagation and physical-observable
agreement. They do not resolve the Si exciton peak or replace the separately
running long optical-spectrum comparison.
