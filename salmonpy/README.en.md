# SALMONpy: material response, electronic dynamics and light propagation in Python

[日本語](README.md)

Use LD, SBE, DFT, TDDFT and FDTD independently or combine them. Structure and pulse
preparation and analysis of existing results can also be used on their own.
This is a development release and is not yet published on PyPI.

| How to use it | Calculations | Preparation |
|---|---|---|
| Calculate material response from Python | LD, SBE, pulses, geometry and analysis | Install a matching distribution wheel |
| Run SALMON in a separate process | DFT, TDDFT, FDTD and bands | Also install a compatible SALMON executable |
| Advance a calculation one step at a time from Python | Couple electronic calculations, FDTD and user-written code | Also prepare the shared library and MPI |

Start with the first row. For additional programs, follow the [installation guide](RUNTIME.en.md).
Prefer a matching prebuilt distribution; explicitly build the bundled sources only when needed.
A separate upstream SALMON checkout is not required.

SALMONpy is a set of reusable calculation components for material electron response and light propagation. Use simple material models or first-principles calculations independently, couple them, or connect your own current model. Input preparation and analysis can also be used alone.

## First calculation: two commands

Obtain a distribution wheel matching your OS, CPU and Python version (Python 3.10+).
SALMONpy is not yet on PyPI. Replace the filename below with the actual supplied wheel;
ask its distributor for the file if it is missing. No compiler is needed.

```sh
python -m pip install --only-binary=:all: "PATH_TO_COMPATIBLE_WHEEL.whl"
python -m salmonpy.quickstart
```

`SALMONpy first calculations: PASS` means both calculations finished:
1000 Lorentz–Drude (LD) updates and 500 semiconductor Bloch equation (SBE) updates
for a synthetic two-band model. Each uses one thread. No additional inputs, MPI,
SALMON executable, plotting package or network access during calculation is needed.

The command prints a new output directory containing `summary.json`, both calculation
logs and SBE results. It normally takes seconds; package download/install time varies.
Use `--output NEW_DIRECTORY` to retain results in a chosen location. Existing directories
are never overwritten. These are installation checks, not converged material properties.
If a calculation fails, read its log; do not automatically compile or install other solvers.

Next, choose a task below or inspect the [examples](examples/README.en.md).
[Si current](SI_INTRO.en.md) additionally needs a matching calculation program; it is
not a prerequisite for completing this first test. For diagnosis only:
`python -m salmonpy.doctor --lang en --json`.

## Choose a task

| Task | Python entry | Requirements |
|---|---|---|
| Permittivity and E-driven P/J | `salmonpy.LorentzDrude` | Wheel |
| Crystalline A/c-to-j response with SBE | `salmonpy.sbe` | Wheel; band/transition data for real materials |
| Fields (j→A) coupled to matter: transverse / longitudinal / thin film / 1D propagation | `salmonpy.maxwell` | NumPy and a matter component; [FIELDS.en.md](FIELDS.en.md) |
| Pulses, geometry and analysis | `fields`, `geometry`, `analysis`, `outputs` | NumPy; Matplotlib for plots |
| Ground-state DFT and forces | `salmonpy.dft` | Compatible full SALMON and pseudopotentials |
| Crystalline A/c-to-j response with TDDFT | `salmonpy.tddft` | Compatible full SALMON and matching GS |
| Electromagnetic propagation | `salmonpy.fdtd` | Compatible full SALMON |
| Band calculations and plots | `salmonpy.bands` | Compatible full SALMON, GS, Matplotlib |
| ASE structures, bands, energy/forces | `salmonpy.ase`, `salmonpy.ase_calculator` | Corresponding solver plus optional ASE dependencies |

Start with the [shared Python samples](examples/README.en.md). Optional packages:

```sh
python -m pip install "matplotlib>=3.6"
python -m pip install "ase>=3.26,<4" "f90nml>=1.4,<2"
```

Install only what your task needs. The SALMON **executable is not included** in the wheel; its sources are bundled for explicit builds.
Use its audited executable and recorded SHA256. Do not assume a standard SALMON
release contains the worker/status interfaces required by these interfaces.

```sh
python -m salmonpy.doctor --lang en --executable PATH_TO_SALMON --sha256 VALIDATED_SHA256
```

`hash_matches` means file identity only, not compatible capabilities, restart
format or successful execution. Build requirements and historical platform
validation are in [PACKAGING.md](PACKAGING.md).

## Choosing a parallel workflow

| Goal | Launch | Available scope |
|---|---|---|
| Light propagation through multiple LD materials | Serial Python, `fdtd_mpi_layers` | Existing spatial MPI FDTD; many local LD states per field rank |
| Your own internally parallel A-to-J model | `mpiexec -n 5`, `mpi_response` | Two material groups, each with two ranks and one material/point |
| Two TDDFT responses from a saved GS | `mpiexec -n 3`, `mpi_tddft` | Two one-rank material groups; independent serial worker per point |

See [materials and MPI FDTD](FDTD_MATERIALS.en.md), [response service](DISTRIBUTED.en.md),
and [commands and outputs](examples/README.en.md). The field role in the latter two
supplies prescribed A. For internal-MPI TDDFT and E/H coupling use the separate
[experimental solver coupling](COUPLED.en.md) and `mpi_coupled_tddft` sample
(two field ranks plus two two-rank material groups = six ranks).
Do not nest the batch FDTD launcher inside an MPI job.

## Units, state and supported scope

- Core APIs use atomic units. Frequencies are angular frequencies; do not pass
  eV or fs directly. Crystalline input is **A/c**, with E = −d(A/c)/dt.
- SBE/TDDFT return the electric (charge) current density j. SALMON's electronic matter current J = −j stays internal (`state.matter_current`).
  LD is E-driven and returns physical P and J, not FDTD's auxiliary Q.
- DFT cells have row vectors; SBE lattice vectors are columns. ASE uses Å and eV;
  the bridge performs explicit conversion.
- Models and initial states fix physical settings. Each `new_state` creates
  independent evolution; dt is fixed per state. Each `respond` starts afresh.
- The `salmonpy.dft`, `bands` and process-isolated `tddft` entries use one MPI rank.
  The separate `embedded_tddft` entry supplies internal MPI with the restrictions
  in [COUPLED.en.md](COUPLED.en.md).
- DFT: periodic, fixed ions, unpolarized, PZ/Libxc-PZ/Libxc-PBE. TDDFT also supports
  solver TB-mBJ from a matching converged GS; midpoint propagation only.
- ASE Calculator: PZ/Libxc-PZ/Libxc-PBE with scalar KY/NC-UPF potentials; energy/forces, opt-in wavefunction warm start, no stress.
  Establish grid, k-point and SCF convergence; egg-box effects affect force checks.
- FDTD has a batch entry and a separate restricted E/H update and coupling entry
  in [COUPLED.en.md](COUPLED.en.md). Periodic boundaries are available where documented; do not infer oblique-incidence support.

## Budgets and failures

Before research calculations, decide the purpose, **total allowed memory** and
**total CPU concurrency**. Count all processes, threads and numerical libraries.
Detected hardware is not permission to use all resources. Memory advice is available for documented calculations; it warns rather than enforcing a hard limit. Start small and measure; resource
fit does not imply physical convergence. Large spatial outputs are opt-in.

Use new writable output directories; existing results are not overwritten. A
nonconverged GS pauses the workflow with a warning and residual history for user review.
Inspect `native/stdout.log` and the result manifest before continuing.
Output archives require an indexed arrays.npz and matching array names, shapes,
dtypes and finite values. verify=False skips hashes only, not structural checks.
Opening reads arrays one at a time, so allow memory for the largest array.
An output archive is not a restart checkpoint. Preserve complete GS directories
for reuse and do not edit their hashed solver files. Close TDDFT states with
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
commands require a checkout. Restart Python after rebuilding Fortran/C libraries.

Public desktop visualization: [SALMON_VR.en.md](SALMON_VR.en.md).

For practical electronic-feedback calculations, start TDDFT with 1D propagation/layered structures. Recommend 3D with SBE or a user-written current solver. This is a cost-based recommendation, not a capability restriction; see [dimension selection](COUPLED.en.md).

## calling program API and state snapshots

See [public API and errors](API.en.md) and [in-memory snapshots](SNAPSHOTS.en.md).
LD, SBE and the four Python fields support fixed-dt snapshot/restore since 0.1.0a5.
For saving FDTD, TDDFT or MPI calculations to disk and restarting them, follow the [checkpoint guide](CHECKPOINTING.en.md).

Custom A-driven models can use one serial/MPI contract and an explicit conformance checker: [authoring guide](MATERIALS.en.md), `python -m salmonpy.examples.custom_material`.

See [LD.en.md](LD.en.md) for persistent LD state ownership and memory.
