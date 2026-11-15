# Python sample audit — 2026-11-15

Scope: current source checkout on Apple Silicon macOS/Python 3.14. Explicitly
loaded the existing LD/SBE library; no import-time compilation, downloads,
GitHub Actions, or physical convergence studies. Native runs used one thread;
MPI runs used three ranks with a 45-second process-group timeout.

## Execution evidence

All 39 Python example modules (excluding `__init__`) were exercised either by
running their calculation or displaying `--help`. This does **not** mean that
all 39 scientific workflows were executed.

Full lightweight runs, all successful:

- `ld_response`, `pulse_spectrum`, `custom_material` (serial)
- `snapshot_restore`, `material_config`, `field_inputs`, `dipole_radiation`
- `sbe_response`, `field_coupling`, `outputs_roundtrip`
- `configured_ld`, `point_ld`, `geometry_preview`

Additional complete runs:

- `prepare_inputs`: staged `{"coefficient": 2.0}` and verified the bundle.
- `prepare_structure`: ASE cubic Si extxyz, spacing 0.4 angstrom, kmesh 2×2×2.
- `fdtd_layers`: existing managed serial runtime, two LD layers, 160 steps.
- Three-rank `mpi_recorded_material`: exact serial/MPI response comparison.
- Three-rank `mpi_point_identity`: persistent material/point IDs.
- Three-rank `custom_material --mpi`: custom material response.

The remaining input-dependent samples were checked with `--help`, not by
repeating DFT, band, TDDFT, or embedded FDTD/TDDFT calculations. No claim is made
about physical convergence or portability to other platforms. The isolated
response and archive-to-dipole-current examples have separate validation records.
The five-/six-rank defaults were not launched under the four-core budget.

## Fixes

Six MPI examples imported mpi4py before parsing arguments. Without MPI installed,
users could not even inspect required inputs. MPI imports now follow parsing.
`mpi_multiscale_tddft` also preserves argparse's normal help/usage exits instead
of treating them as distributed failures. Its real failure handler retains MPI
job cancellation. No solver stepping or numerical formula changed.

`test_example_cli.py` blocks mpi4py explicitly and checks all six help commands.
Both example READMEs now list the previously missing connection/input examples,
distinguish dipole current generation from actual FDTD propagation, and warn that
some MPI defaults exceed four cores. The LD example no longer requires a checkout
in its instructions. The English six-command count was corrected.

## Local evidence paths (not prerequisites for users)

- CLI/lightweight logs: `/var/folders/ng/vj2yg45d3850_stlzd3tfhlc0000gn/T/salmonpy-example-audit-wel427jy`
- Input/FDTD logs: `/var/folders/ng/vj2yg45d3850_stlzd3tfhlc0000gn/T/salmonpy-extra-examples-cexdjadu`
- MPI logs: `/var/folders/ng/vj2yg45d3850_stlzd3tfhlc0000gn/T/salmonpy-example-mpi-klxhrjy4`
- Targeted CLI/manual tests: `/tmp/salmonpy-example-contract-tests.log`
- Full source tests: `/tmp/salmonpy-example-full-tests.log`

## Fresh wheel verification

Built the sdist and then the wheel with `python -m build --no-isolation`,
CMAKE_BUILD_PARALLEL_LEVEL=2; no downloads or Actions. Installed with `pip install
--no-deps --no-index` in a new venv. Reused the existing dependency environment via
an explicit .pth file; this is not an online dependency-resolution test. Native
LD/SBE were loaded from the installed wheel, without the source-suite overrides.

- Wheel: `salmonpy-0.1.0a5-cp314-cp314-macosx_26_0_arm64.whl`.
- Full installed tests: 409, no failures, two skips (optional plot dependency and
  case-sensitive-filesystem-only test).
- All 39 installed sample modules: 13 lightweight executions and 26 help checks,
  all exit zero, with isolated Python (`-I`) outside the checkout.
- Additional `dipole_radiation --archive ... --cell-volume 8000`: successful
  readback of the earlier native isolated archive and conservation check.
- Public isolated loader against the existing validated E-worker: independent
  states, exact replay continuation and killed-worker recovery all pass.
- This run covers Apple Silicon/macOS/Python 3.14 only. Runtime is a separately
  prepared local binary, not bundled in this wheel or newly published.

Local build artifacts: `/tmp/salmonpy-sample-audit-wheel-20261115`.
Installed venv: `/tmp/salmonpy-sample-audit-installed-20261115`.
Build/test logs: `/tmp/salmonpy-sample-audit-build.log`,
`/tmp/salmonpy-sample-audit-installed-tests.log`,
`/tmp/salmonpy-sample-audit-isolated.log`.
Sample report: `/var/folders/ng/vj2yg45d3850_stlzd3tfhlc0000gn/T/salmonpy-installed-samples-6ovcwfqz/report.json`.
Native report: `/tmp/salmonpy-sample-audit-isolated-20261115/report.json`.
