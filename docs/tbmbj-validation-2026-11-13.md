# TB-mBJ native worker validation, 2026-11-13

The worker now permits native `xc='tbmbj'`, without requiring Libxc.
Only the input allowlist changed. Native automatic predictor-corrector dispatch,
current correction to kinetic-energy density, and automatic/fixed cval behavior
remain unchanged. In particular TB-mBJ executes predictor-corrector even when
`yn_predictor_corrector='n'`. cval belongs to the immutable initial input; only A/c
is supplied at each step.

## Tests on Apple Silicon Mac

- The original worker rejected TB-mBJ before initialization (recorded red test).
- The synthetic integration fixture passed with automatic cval and fixed cval=1.1
  (explicit predictor-corrector enabled for the fixed case). Ordinary pre/post
  current and energy outputs were byte-identical. Worker currents, independent
  states, invalid-input rejection and killed-worker recovery passed.
- The automatic-cval fixture also passed using a Libxc-disabled executable.
- Source and application objects passed the forbidden-eigensolver audit.
- All 34 Python tests passed.

Real GS test: 8-atom Si, lattice 10.26 a.u., 12³ grid, 2³ k points, 32 states and
32 electrons; repository Si_rps.dat, SHA256
`8c51fbf4341de256a868ab446bdb3ed19eef1d1dcc3813cfe5214c07213b8c7c`.
The default Gaussian-orbital start produced NaNs at the second SCF iteration in
the pre-change ordinary SALMON binary. That calculation is rejected, not used as
an initial state. The immediate numerical origin was not instrumented; no change
to the native TB-mBJ algorithm or broad claim about all initialization choices
is made.

A second GS run used the previously converged PZ wavefunction as an initial guess:
`yn_restart='y', read_gs_restart_data='wfn', yn_reset_step_restart='y'`, with
`xc='tbmbj'` and default automatic cval. It was fully reconverged under TB-mBJ:
28 iterations, density residual 3.0353314e-10 (threshold 1e-9). Only that resulting
TB-mBJ GS was loaded for RT, with fresh-RT `yn_restart='n'` as required by the worker.

RT: 500 steps at dt=0.05 a.u., A/c z pulse amplitude .001, frequency .2, duration
20 a.u. Max ordinary/worker current difference: 5.421010862427522e-20 a.u.
Two states from the same TB-mBJ GS were advanced for 64 interleaved steps with
opposite A and matched their separate trajectories. These are functional/regression
tests, not spatial/k/time convergence evidence. The KY file header does not
establish its generating functional provenance. No PAW/USPP or other-material
support is inferred. TB-mBJ is used here as a model potential; the reported native
energy is not a validation of a variational TB-mBJ energy functional or forces.

## Artifacts

Full inputs, outputs, initial-guess seed, summaries, scripts, executable copies,
failed-start compressed log and regression logs are saved outside the repository:
`/Users/otobetoshihito/.codex/visualizations/2026/11/13/01a1d275-1ce2-7121-b8f4-cfd9dbc2d897/salmonpy-tbmbj-seeded-20261113`.
SHA256 values are recorded in the summary and archive manifest.

The runtime distribution proposal is in
`docs/plans/2026-11-13-runtime-distribution-proposal.md`; no optional full-runtime
wheel or preparation/capability API has been implemented by this change.

Environment note: absolute paths in this report identify the historical validation host. They are not portable setup defaults or required user paths.
