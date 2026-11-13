# SALMONpy — getting started

Full-runtime setup: prefer a compatible prebuilt binary, with an explicit source-build fallback. See [runtime guide](RUNTIME.en.md). Branch CI artifacts and public releases are distinct.

> Standalone layout: full SALMON sources are no longer bundled. Source-checkout
> runtime setup: `runtime/README.md`; use `tools/prepare_salmon_runtime.py` with
> the pinned external upstream checkout. Historical validation commands below
> describe the old source layout and are not current build instructions.

[日本語](README.md) · Topic: `getting-started` · Paired revision: `2026-11-13.1`

SALMONpy provides independent material-response and analysis components. It is a
preview, not a published PyPI release. The project is intended to move to its own
repository; this branch is not intended for merging into SALMON itself.

## Install and diagnose

Obtain a **repaired wheel** matching your OS, CPU architecture and CPython version.
Python 3.10+ is required, but the wheel must also match your specific interpreter.

```sh
python -m pip install "PATH_TO_COMPATIBLE_WHEEL.whl"
python -m salmonpy.doctor --lang en
python -m salmonpy.examples.ld_response
```

Replace the placeholder with the actual file. Use the Python interpreter belonging
to the installation environment. Wheels do not require a Fortran compiler;
importing never compiles anything. Diagnostics load the bundled libraries but do
not calculate, install software or start the external SALMON executable.

## Choose a task

| Task | Python entry | Requirements |
|---|---|---|
| Permittivity and E-driven P/J | `salmonpy.LorentzDrude` | Wheel |
| Crystalline A-to-J with SBE | `salmonpy.sbe` | Wheel; band/transition data for real materials |
| Pulses, geometry and analysis | `fields`, `geometry`, `analysis`, `outputs` | NumPy; Matplotlib for plots |
| Ground-state DFT and forces | `salmonpy.dft` | Compatible full SALMON and pseudopotentials |
| Crystalline A-to-J with TDDFT | `salmonpy.tddft` | Compatible full SALMON and matching GS |
| Electromagnetic propagation | `salmonpy.fdtd` | Compatible full SALMON |
| Band calculations and plots | `salmonpy.bands` | Compatible full SALMON, GS, Matplotlib |
| ASE structures, bands, energy/forces | `salmonpy.ase`, `salmonpy.ase_calculator` | Corresponding solver plus optional ASE dependencies |

Start with the [shared Python samples](examples/README.en.md). Optional packages:

```sh
python -m pip install "matplotlib>=3.6"
python -m pip install "ase>=3.26,<4" "f90nml>=1.4,<2"
```

Install only what your task needs. Full SALMON is **not included** in the wheel.
Use its audited executable and recorded SHA256. Do not assume a standard SALMON
release contains the worker/status interfaces required by this preview.

```sh
python -m salmonpy.doctor --lang en --executable PATH_TO_SALMON --sha256 VALIDATED_SHA256
```

`hash_matches` means file identity only, not compatible capabilities, restart
format or successful execution. Build requirements and historical platform
validation are in [PACKAGING.md](PACKAGING.md).

## Units, state and supported scope

- Core APIs use atomic units. Frequencies are angular frequencies; do not pass
  eV or fs directly. Crystalline input is **A/c**, with E = −d(A/c)/dt.
- SBE/TDDFT return electronic matter current J; electronic charge current is −J.
  LD is E-driven and returns physical P and J, not FDTD's auxiliary Q.
- DFT cells have row vectors; SBE lattice vectors are columns. ASE uses Å and eV;
  the bridge performs explicit conversion.
- Models and initial states fix physical settings. Each `new_state` creates
  independent evolution; dt is fixed per state. Each `respond` starts afresh.
- Python DFT/bands/TDDFT launches currently use one MPI rank. Do not infer arbitrary
  MPI scaling from native SALMON support. TDDFT states own separate processes.
- DFT: periodic, fixed ions, unpolarized, PZ/Libxc-PZ/Libxc-PBE. TDDFT also supports
  native TB-mBJ from a matching converged GS; midpoint propagation only.
- ASE Calculator: PZ with scalar KY potentials, energy and forces only, no stress.
  Establish grid, k-point and SCF convergence; egg-box effects affect force checks.
- FDTD's public entry is batch execution; a Python single-step E/H API, periodic
  boundaries and oblique incidence are not exposed.

## Budgets and failures

Before a first calculation, decide the purpose, **total allowed memory** and
**total CPU concurrency**. Count all processes, threads and numerical libraries.
Detected hardware is not permission to use all resources. Capacity estimates and
hard memory enforcement are not implemented. Start small and measure; resource
fit does not imply physical convergence. Large spatial outputs are opt-in.

Use new writable output directories; existing results are not overwritten. A
nonconverged GS is an error. Inspect `native/stdout.log` and the result manifest.
Output archives require an indexed arrays.npz and matching array names, shapes,
dtypes and finite values. verify=False skips hashes only, not structural checks.
Opening reads arrays one at a time, so allow memory for the largest array.
An output archive is not a restart checkpoint. Preserve complete GS directories
for reuse and do not edit their hashed native files. Close TDDFT states with
`with`; send/response timeouts terminate or recover the worker as documented.

## Ask an AI assistant

```sh
python -m salmonpy.manual --lang en "I want to run TDDFT"
python -m salmonpy.doctor --lang en --json
```

The AI catalog uses stable identifiers, localized intents and shared technical
rules. English mode localizes onboarding and task titles; detailed legacy task
instructions are still mixed-language. It is not yet a fully translated reference.
`sample_files` points to executable examples. Unmatched queries return a task
catalog and getting-started guidance. Reuse known resource answers.

For source development, see [PACKAGING.md](PACKAGING.md); source-only build/test
commands require a checkout. Restart Python after rebuilding native libraries.

Public desktop visualization: [SALMON_VR.en.md](SALMON_VR.en.md).
