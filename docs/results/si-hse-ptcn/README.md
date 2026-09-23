# Si: PT-CN-ACE short time-step comparison

This extends the independent exported-Hamiltonian HSE06 driver with parallel-transport Crank–Nicolson. It does not enable native SALMON HSE. Si8, shifted4³ k mesh, 12³ primitive grid, impulse A_z=1e-4 a.u., all exchange pairs; FFTW and BLAS use one thread on Apple M5 Pro. No k convergence study.

The implemented endpoint equation is Eq.12 of [Jia and Lin, PT-CN-ACE](https://arxiv.org/abs/1809.09609). It removes occupied-space rotations from the orbital derivative. Inner iterations reuse ACE, while fresh full exchange checks the endpoint residual before acceptance. The full accepted endpoint action is cached for the next step's initial endpoint. This cache is valid for the constant-field driver; a time-dependent laser requires explicit endpoint field handling.

The initial PT derivative norm is 0.000339 times the ordinary Schrödinger derivative norm. This reflects gauge removal, not a reduction in physical response, and does not prove a safe time step. Orbitals in different gauges are compared only after unitary alignment; current and density are compared directly.

## Verification scope

Tests cover grid-weighted PT projection, complex occupied-gauge covariance, stationary occupied subspaces, second-order density convergence in a constant-H model, state-dependent exchange with a deliberately stale ACE, cached versus rebuilt initial endpoint actions, and nonconvergence/input failures. Independent review found no equation, cache or energy-factor blocker.

No normalization is performed inside or after PT-CN. Nonlinear trapezoidal integration is not assumed to preserve endpoint orthogonality exactly. The driver records electron-number and Gram errors and fails above its gates. The full relative endpoint residual threshold is1e-10 and the inner threshold1e-12.

The common final time is **0.32 a.u. =0.00774043 fs**. This is a short numerical pilot, not a spectrum or long-time stability study. The reference is full self-consistent RK4 at dt0.08 (four steps), not an independently time-converged spectrum. Comparisons report finite-step differences. Original-gauge midpoint at dt0.32 is included to avoid attributing implicit-integration benefits solely to PT.

Timings distinguish initial full-action/energy setup from propagation. The reported wall interval includes solver initialization, bootstrap and final diagnostics, but excludes input loading/localizer initialization and final compressed checkpoint writing; propagation time sums measured step walls. These are observed short-run costs, not statistically repeated throughput/scaling tests. Full exchange counts include the bootstrap in the PT driver. The earlier midpoint benchmark additionally performs standalone ACE timing and full endpoint energy evaluations, so compare propagation times rather than its total wall time.

Raw JSON and exact timing/accuracy comparisons are in `validation.json`; regenerate them with `collect_results.py`. Large checkpoint arrays are local under `calculations/si_hse_reference/ptcn_*` and ignored by Git. See `samples/hse_mlwf_reference/README_PTCN.md` for execution details.


## Inner-iteration stagnation handling

The initial zero-field dt0.32 run did not meet the nominal inner goal1e-12. Instrumentation showed residual6.34e-8 initially,1.35e-11 by iteration10, then fluctuations around1.1e-11–3.1e-11 through iteration140. This was not accepted as a completed trajectory. The revised inner solve may stop early only after at least9 evaluations, with residual below half the full tolerance and insufficient improvement between the best residuals in two consecutive four-iteration windows. It then rebuilds full exchange at the unchanged candidate. Acceptance still requires the original full residual <1e-10. Fallback events log both inner and full residuals. Tests cover both valid acceptance and rejection when the fresh full residual remains too large. Positive-impulse sweep timings preceding this fallback used strict inner convergence; their values are retained as measured.


The revised zero-field run triggered two inexact-inner exits. Its first fresh full residual was3.14e-10 and was correctly rejected; after updating exchange the residual became1.54e-11 and was accepted. Electron number was32.0 and maximum Gram error1.12e-13. Thus the fallback does not simply declare convergence at the observed inner floor.


## Measured fixed-final-time comparison

All rows end at0.32 au; times below are propagation only, with loading/bootstrap/final checkpoint excluded.

| Method | dt (au) | Steps | Propagation (s) | Current difference vs RK4 | Relative density difference |
|---|---:|---:|---:|---:|---:|
| Full RK4 |0.08|4|181.28|reference|reference|
|PT-CN-ACE|0.08|4|97.18|0.0246%|8.676e-08|
|PT-CN-ACE|0.16|2|72.05|0.0971%|2.914e-07|
|PT-CN-ACE|0.32|1|37.82|0.3313%|6.854e-07|
|Original-gauge midpoint + ACE|0.32|1|62.58|0.3388%|5.148e-06|

At dt0.32, PT-CN is1.65× faster than original-gauge midpoint and its density difference is7.51× smaller. Current differences remain similar (~0.33%). PT-CN energy change is-7.46e-14 Ha versus9.42e-9 Ha for original-gauge midpoint in this pilot; these small energy values are not a long-time conservation claim. The full RK4 reference is a finite-step calculation; listed differences are not experimental errors.

The dt0.16 PT run uses2.52× less propagation time than RK4, with0.0971% current difference in this short test. The dt0.32 cost ratio is4.79, at0.3313% current difference. These are cost/accuracy pairs, not equal-accuracy speedup claims.

All46 Python tests passed after the stagnation fix. No Fortran changes were made. Production time-step selection and long-time HSE spectra remain unvalidated.
