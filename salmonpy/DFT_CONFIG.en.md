# Material configuration → DFT → TDDFT

[日本語](DFT_CONFIG.md) · Topic: `dft_material_config` · Paired revision: `2026-11-15.1`

`builtin_material_registry()` explicitly includes backend `salmonpy.dft`, model
`dft`, adapter version `1`. This resolves local structure/potential references into
the existing DFT API, with no runtime launch, acquisition or compilation. ASE is
an optional dependency required for this resolver; LD/SBE registration stays usable
without it.

## Recipe

```json
{
  "schema": "salmonpy.materials/1",
  "materials": {
    "film": {"backend": "salmonpy.dft", "model": "dft", "parameters": {
      "structure": {"path": "structures/cell.xyz", "format": "extxyz", "index": 0},
      "pseudopotentials": [{"path": "pseudo/Si_rps.dat", "atomic_number": 14, "local_channel": 2}],
      "xc": "pz", "states": 32, "kmesh": [2,2,2],
      "grid": {"spacing": {"value": 0.4, "unit": "angstrom"}}
    }}
  }
}
```

Numbers illustrate the schema, not converged parameters. Structure, potentials,
XC, states, kmesh and grid are required. `index` defaults to 0. Optional `electrons`
must equal the neutral valence sum read from the chosen potentials. Grid is exactly
`{"counts":[16,16,16]}` or `{"spacing":{"value":0.4,"unit":"angstrom"}}`;
width also accepts three values and bohr. Width resolution uses STRUCTURE.en.md's
positive axis-aligned periodic-cell limits; explicit counts retain the existing
DFT geometry contract. No rank-dependent rounding. Unknown keys are rejected.

Paths must be local relative references staying inside base_directory. Structures
use explicit CIF/extxyz and preserve the supplied supercell, subject to CIF symmetry
expansion. Mixed/partial occupancy, charge/spin and nonperiodic ASE conditions are
rejected under existing adapter limits. Potentials use existing KY/NC-UPF parsing,
unique species mapping and neutral electron validation. This resolver additionally
requires the NC-UPF header element to match atomic_number. KY element metadata and
PP functional/XC provenance are **not certified**; choose verified compatible
potentials. Supported XC names remain pz/libxc_pz/libxc_pbe, not all Libxc or HSE.

Resolved parameters retain source references, normalized settings, electrons,
`structure_record`, `grid_resolution` and `file_sha256`. Files are hashed around
loading and compared again at model creation. A changed file requires explicit
new resolution. Counts-only reports give counts; spacing reports also preserve
requested/actual widths. No large wavefunction data go into JSON.

## Execute and preserve identity

```python
from salmonpy.builtin_materials import builtin_material_registry, material_from_catalog
from salmonpy.config import load_materials
from salmonpy.dft import SCFSettings
from salmonpy.runtime import load_runtime

registry = builtin_material_registry()
catalog = load_materials("inputs/materials.json", registry=registry)
material = material_from_catalog(catalog, "film", base_directory="inputs")
engine = load_runtime("runtime", required=["dft", "tddft", "pz", "single-rank"])
gs = material.run(engine, "new-gs", scf=SCFSettings(max_iterations=300, tolerance=1e-9), threads=1)
with gs.tddft_initial_state(executable=engine, max_steps=100, threads=1) as initial:
    with initial.new_state(dt=0.02) as state:
        j = state.step([0.001, 0.0, 0.0])  # A/c in atomic units; charge current j
        identity = state.material_id, state.definition_digest
```

Runtime, SCF stopping controls, thread count, timeout and warm-start choice are run
arguments, separate from the fixed material definition. Select required runtime
features for your chosen XC (including Libxc when needed). Existing DFT checks,
convergence report and verified GS→TDDFT loader remain in use; the wrapper does
not relax feature/physics limits. The fundamental crystal path is transverse,
A/c-driven. Two new_state calls use the existing independent process workers.

`material.run` returns `ConfiguredGroundState`; `.result` is the original
GroundStateResult, and other properties are delegated. `initial_guess=previous_gs`
accepts either a configured result or the original result, unwrapping before the
existing warm-start compatibility checks. `.tddft_initial_state(...)` returns a
ConfiguredMaterial holding the original TDDFTInitialState; its context closes that
owner. Model contexts with no close method are no-ops. Use normal state contexts.

Material ID and definition digest persist in the live DFT→GS→TDDFT→state path.
`record()` on material/GS returns identity for explicit result metadata. The digest
identifies fixed material input, not SCF execution, GS bytes, pulse or dt. Existing
GS/TDDFT manifests carry their own data hashes. This is **not** a generic JSON
TDDFT backend or a loader selecting catalog.initial_states recipes. Disk identity
persistence, point/run IDs, MPI assignment and arbitrary saved-GS selection remain
separate work. Passing `.model` or `.result` to other APIs transfers responsibility
for preserving that identity to the caller.

## Offline Si example

```sh
python -m salmonpy.examples.configured_si RUNTIME_DIRECTORY NEW_OUTPUT_DIRECTORY
```

Uses the packaged Si starter, explicit local runtime and staged local structure/
potential files. It runs the short existing GS/TDDFT regression, checks interleaved
opposite-drive workers and an exact replay, compares against the starter's tolerated
reference, and records material ID and current outside the input bundle. One thread
per process, at most two TDDFT workers. Memory is not hard-limited. This is a short
connection check, not grid/k/time-converged optics. No network or runtime build.

Explicit ID selection of a pinned saved GS/1 is now available separately via
`tddft_from_catalog`; see [INITIAL_STATES.en.md](INITIAL_STATES.en.md).
Generic material construction itself does not select an initial state.

The separate [PointRun layer](POINT_RUNS.en.md) now carries serial point/run IDs
and saves result mappings. The limitations above still apply to automatic
propagation by the original API alone, MPI and restarts.

## Developer-only standard GS model inputs

The standard runtime adds internal Fortran type salmonpy_gs_model::s_gs_model.
load(path,status) retains a model-owned parser and post-coordinate input values;
run(status) activates them for sequential legacy SCF, and clear() releases input
and scratch. A rejected unsupported replacement preserves the existing model.

Scope: periodic full-k unpolarized fixed-ion PZ/TB-mBJ/Libxc-PZ/PBE. The host first
initializes MPI, comm_attach, math/parallel infrastructure and timers. All subgroup
ranks call the same model in matching order. Runs switch/restore caller input and
unit views, calc_mode, actual symmetry switches, modulo tables and XC scratch.
The model-only call retires persistent communication caches before local buffers
expire and releases scoped child communicators after each call. The default run(status) call retains no GS orbitals. Optional retained storage
and Python handles are now available through [retained GS](NATIVE_GS.en.md).
Simultaneous GS/TDDFT, threads and HSE/SOC/+U/DC GS remain separate work.

Status0 means returned normally,1 means empty model/occupied child scope,2 means
unsupported conditions for this internal entry. SCF nonconvergence remains a
warning and outputs the last iterate. General parser/STOP/allocation failure
recovery is not promised. Relative references use the caller's working directory.
Python DFT.run continues to use its existing separate-process execution path.

Developer probe: tools/test_standard_gs_models.py NEW_BUILD OLD_BUILD NEW_OUTPUT
--reference-sha256 SHA --size 2 --split spatial. Supply explicitly audited matching
runtimes; the tool audits probe sources and all linked application objects before
use. It compares A→B→A with independent old-executable controls for8³/12³ grids,
different cell/atom/electron/XC inputs, checking energies, residuals, iterations,
eigenvalues, k output, restored symmetry and repeated communicator counts.
It sets write_gs_restart_data='no', generating/collecting no GS restart bundles.
Four-iteration diagnostic inputs establish implementation parity, not converged
material properties.
