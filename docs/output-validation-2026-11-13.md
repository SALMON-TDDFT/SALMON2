# Selectable output validation — 2026-11-13

Implemented opt-in postprocessing archives through
`crystal.respond(..., output_directory=...)` and `TDDFTState.finish(directory)`.
Native files are flushed on clean process exit before copying. A/J histories,
endpoint times, units/conventions, input, initial file hashes, executable identity,
file sizes/checksums and native headers are indexed. Unchanged GS restart files
are omitted; all other worker files are retained, including identical preexisting
outputs, with a content-matches-initial flag. Failed exports and interrupted exports
retain worker files in a recovery directory even after close/GC.

Readers: numerical NPZ arrays (including complex arrays), native rectangular
text tables, scalar positive-grid CUBE, and k/band/spin current blocks. The latter
preserves integer step/index arrays and existing occupation/k-point weighting.
No unit conversion or invented field time location is applied.

## Local Mac tests

PZ and TB-mBJ: 24 steps, dt=.05 a.u., 12³ Si grid, same converged GS states as the
prior validation. Output strides 8, CUBE format. Each case compared ordinary
SALMON with selected outputs off/on, then the process-isolated SALMONpy worker.
The native file drive includes a tail beyond the requested RT duration because
the existing file interpolation excludes its last input knot; actual applied A
was explicitly checked before comparing currents.

- Output selection changed current by exactly zero in both cases.
- Worker/reference max current difference: 4.743384504624082e-20 a.u. in both cases.
- 36 selected native files per case were byte-identical: energy/current tables,
  decomposed current, density/density difference, ELF, microscopic current and
  electrostatic-field cubes. This excludes native post-run pulse/FFT output, which
  the interactive worker does not run; A is saved explicitly in arrays.npz.
- CUBE values were finite and retained 12³ shape. Density integrals at step 24:
  31.999999999997726 (PZ), 31.999999999998394 (TB-mBJ), versus 32 electrons.
- Summed band/k/spin contributions matched the total current within
  5.20e-19 (PZ) and 6.98e-19 (TB-mBJ) a.u.
- Partial-run finish saved exactly the accepted two samples and closed the state.
- Normal native A→J regression, independent states and killed-worker recovery passed.
- Unsupported RT wavefunction checkpoint requests are now explicitly rejected by
  the worker instead of silently disabling the request.
- 42 Python tests passed, including save-failure preservation, interruption,
  identical-input/output content, no overwrite, checksums, complex arrays,
  scalar-grid order and weighted-current block reading.

A local Apple Silicon CPython 3.14 wheel was built and passed through the existing
runtime-bundling/hash-refresh repair step. Installed-wheel validation is recorded
with the artifacts. These changes have not been revalidated on remote Windows,
Intel Mac or Linux CI in this turn. No new separate-node transport, CPU affinity,
shared-state thread safety or public E/H coupling API is claimed.

The available output selections and actual usage are documented in
`salmonpy/OUTPUTS.md`. The archive is not a restart checkpoint; large RT wavefunction
outputs and projection/unfolding remain unsupported. TB-mBJ native energy output
is not a validation of a variational energy functional.

Artifacts (both output archives, exact regression logs, source test script and
local repaired wheel):
`/Users/otobetoshihito/.codex/visualizations/2026/11/13/01a1d275-1ce2-7121-b8f4-cfd9dbc2d897/salmonpy-outputs-20261113`.

Environment note: absolute paths in this report identify the historical validation host. They are not portable setup defaults or required user paths.
