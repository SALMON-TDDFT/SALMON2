# Si: from an installed wheel to charge current

[日本語](SI_INTRO.md) · Topic: `si_intro` · Paired revision: `2026-11-14.1`

This starter calculates a real self-consistent Si ground state, then propagates
it under a prescribed A/c pulse and checks the charge current against a fixed
reference. The wheel supplies the input, Si pseudopotential, reference values
and provenance. The full solver is a separate prebuilt runtime; this route does
not compile or require a source checkout.

## 1. Install the supplied wheel

Use a Python version and wheel matching your operating system and CPU. This new
sample requires a wheel containing `salmonpy.examples.si_intro`; the earlier
wheel workflow at commit `6f8b36c6` predates this sample. A source change alone
does not publish an updated wheel. Obtain the updated wheel from the distributor.

On macOS or Linux:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install --only-binary=:all: ./SALMONPY_WHEEL.whl
python -m salmonpy.doctor --json
python -m salmonpy.examples.si_intro --help
```

Replace `SALMONPY_WHEEL.whl` with the actual downloaded wheel filename, including
its Python/platform tags. Keep `--only-binary=:all:`: a missing compatible wheel
must not trigger a source build. NumPy is the only Python dependency needed for
this example. Full-solver Windows binaries are not provided by this starter;
Windows LD/SBE wheels do not establish Windows TDDFT support.

## 2. Download the matching prebuilt runtime

Use the exact artifact from [successful runtime run 43541732653](https://github.com/SALMON-TDDFT/SALMON2/actions/runs/43541732653),
commit `6f8b36c6c04f7f1f3be4e68e718a35cd22f6453c`:

| Machine | Artifact |
|---|---|
| Apple Silicon Mac | [salmonpy-runtime-macos-15](https://github.com/SALMON-TDDFT/SALMON2/actions/runs/43541732653/artifacts/14092143545) |
| Intel Mac | [salmonpy-runtime-macos-15-intel](https://github.com/SALMON-TDDFT/SALMON2/actions/runs/43541732653/artifacts/14091993925) |
| Linux x86_64 / glibc | [salmonpy-runtime-ubuntu-22.04](https://github.com/SALMON-TDDFT/SALMON2/actions/runs/43541732653/artifacts/14092432957) |

GitHub may require sign-in to download Actions artifacts. These are branch
artifacts, not permanent public Releases. Their recorded expiry is 2027-02-12;
if unavailable, request a replacement starter with new verified provenance.
Do not substitute a different run or disable checksum validation. Runtime OS
compatibility is checked during preparation (the arm64 Mac baseline is 15.8.2).
These profiles provide PZ/TB-mBJ and one MPI rank, without Libxc. This fixed
starter uses PZ; it does not require Libxc.

Keep the downloaded **outer Actions ZIP** intact. You do not need to unpack the
inner runtime ZIP or edit a catalog. In the next command `RUNTIME_ARTIFACT.zip`
means this outer ZIP, using its actual downloaded name.

## 3. Prepare and calculate

Choose a new directory outside the source checkout. Preparation verifies the
pinned outer ZIP, then checks the inner runtime and its libraries through the
normal runtime installer. It copies the fixed input and pseudopotential locally.
Nothing is downloaded implicitly.

```sh
python -m salmonpy.examples.si_intro prepare new-si --runtime-artifact RUNTIME_ARTIFACT.zip
python -m salmonpy.examples.si_intro run new-si --threads 1 --memory-budget-gib 4
```

The example above records a **4 GiB total-job allowance** and uses one solver
thread; it is not a 4 GiB per-process allowance. Set a budget you have permission
to use. `--threads` is limited to 1–4; other numerical-library threads are limited
by the example. DFT and TDDFT run sequentially. The memory budget is recorded,
not enforced as a hard operating-system limit. Unavailable resource measurements
are reported separately and do not prove budget compliance. Do not scale this
example based on an unmeasured memory estimate.

Preparation and calculation preserve existing directories. A retry uses a new
prepared directory; do not erase failed outputs to make the command succeed.

## 4. Check the result

The calculation report is `new-si/result/report.json`. Read calculation and
numerical-validation outcomes separately from optional resource measurements.
A validation failure is a failure even when the native calculation completed.
`result/gs/native/stdout.log` keeps the GS log; `result/rt` keeps the response
archive. `result/current.csv` and `result/current.npz` contain time, A/c and j.
CSV columns are `time_au, Ac_x_au, Ac_y_au, Ac_z_au, j_x_au, j_y_au, j_z_au`. The current is **physical charge-current density j**, in atomic
units. SALMON's electronic matter current is J = −j; do not negate the returned
current again.

The fixed reference in `examples/si_intro_data/expected.json` contains all 21 time
samples. The reference energy is about **−31.7465390575 Hartree**, and the largest
absolute current component is about **1.47275×10⁻⁶ a.u.** The automatic checks use
an energy absolute tolerance of 1e−7 Hartree and per-component current tolerances
`atol=1e-9`, `rtol=1e-5`. An SCF iteration count is diagnostic, not an exact
acceptance value. The reference was generated with the pinned prebuilt macOS
arm64 runtime; availability of another platform artifact alone does not certify
this starter's cross-platform numerical agreement.

This is eight-atom cubic Si, cell length 10.26 bohr, 32 electrons/32 states,
16³ real-space grid, 2³ k mesh, PZ. The GS residual must satisfy 1e−10. The real-time
run is only 20 steps at dt=0.02 a.u.; its short z-polarized pulse is an interface
exercise. This **is not converged optical material data**. Research use requires
separate grid, k-point, state-count, time-step, duration and pulse checks.

## Included provenance

`examples/si_intro_data/input.json` fixes every starter setting.
`provenance.json` records source and runtime identities and file hashes.
`Si_rps.dat` comes from SALMON's `samples/exercise_04_bulkSi_gs` at upstream
commit `3322d0561c0cad5180679136b3eddbe8b909487c`; its SHA256 is
`8c51fbf4341de256a868ab446bdb3ed19eef1d1dcc3813cfe5214c07213b8c7c`.
The upstream Apache-2.0 license is retained as `LICENSE.txt`. This legacy KY file
has no embedded generator-XC metadata; its provenance is the pinned PZ Si
example, not an inferred metadata certification. GS→RT uses the same saved
physical input and pseudopotential. The fixed starter is not an arbitrary-GS
compatibility validator.
