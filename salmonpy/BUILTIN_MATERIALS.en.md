# Built-in LD and SBE material configurations

[日本語](BUILTIN_MATERIALS.md) · Topic: `builtin_materials` · Paired revision: `2026-11-15.1`

`builtin_material_registry()` explicitly creates a new MaterialRegistry with
`salmonpy.ld`/`ld` and `salmonpy.sbe`/`sbe` resolvers, adapter version `1`.
No registration, resolver or model factory initializes native states, fetches data,
or compiles. `new_state` needs the existing installed LD/SBE library.

```python
from salmonpy.builtin_materials import builtin_material_registry, material_from_catalog
from salmonpy.config import load_materials
registry = builtin_material_registry()
catalog = load_materials("inputs/materials.json", registry=registry)
material = material_from_catalog(catalog, "film", base_directory="inputs")
with material.new_state(dt=0.02) as state:
    print(state.material_id, state.definition_digest)
```

Choose input/state arguments according to the model: LD remains **E-driven** with
`initial_field`; SBE remains **A/c-driven** with `initial_A` and `time0`.
The wrapper forwards numerical methods, current, snapshot/restore and context
cleanup without changing existing contracts or update order. It does not turn
LD into an A-driven model. Use the existing explicit LD coupling adapters for that.

## LD parameters

```json
{
  "schema": "salmonpy.materials/1",
  "materials": {
    "film": {"backend": "salmonpy.ld", "model": "ld", "parameters": {
      "units": "atomic_gaussian", "plasma_frequency": 0.5,
      "strength": [1.0], "resonance": [0.3], "damping": [0.1]
    }}
  }
}
```

All five fields are required. Units are explicit: atomic Gaussian with angular
frequencies, passive isotropic poles and background permittivity 1. Existing
LorentzDrude validation applies; unknown fields/other units are rejected. Numbers
are illustrative, not measured material constants. No automatic eV/fs conversion.

## SBE parameters

```json
{
  "schema": "salmonpy.materials/1",
  "materials": {
    "silicon": {"backend": "salmonpy.sbe", "model": "sbe", "parameters": {
      "units": "atomic", "directory": "sbe", "sysname": "si",
      "bands": 2, "kpoints": 1, "electrons": 2,
      "lattice": [[4,0,0],[0,4,0],[0,0,4]]
    }}
  }
}
```

This is a schema example, not actual silicon band data. Required fields are shown;
optional `active_bands` defaults to all bands, `correction_order` to 0 and
`nonlocal_correction` to false. Lattice **columns** are vectors in bohr; do not pass
ASE row vectors in Å without explicit conversion. `directory` is relative to the
input base (`.` allowed), and sysname is a single filename stem. The three local
files are `si_k.data`, `si_eigen.data`, `si_tm.data` in that directory. Existing
SBE single-spin/atomic-unit export requirements, closed-shell occupations (export occupation columns are ignored) and
matrix shapes apply; see [SBE.en.md](SBE.en.md). No binary checkpoint or generic
NPZ loader is implied.

Resolution reads/validates arrays and hashes all three files before/after loading.
The resolved parameters contain `file_sha256`. Model creation reloads and rechecks
against the catalog. Changed files are rejected, not silently accepted as the
same material. This is initialization work, not a per-step hash. File sizes and
memory are not automatically budgeted; constructing a catalog and later a model
requires loading again. Neither pass initializes native handles.

## Identity and independent states

`material_from_catalog` returns `ConfiguredMaterial`: `.material_id`,
`.definition_digest`, `.model` (the original LD/SBE model) and `.record()` with
identity fields. `.new_state(...)` returns `ConfiguredState` carrying the same
immutable identity fields while retaining the existing state methods. Use separate
calls for separate states and close them. Other existing public model properties
are delegated. The original `.model` can be passed to callers requiring its concrete
type, but then that caller must preserve the catalog identity separately.

Snapshot tokens remain the existing same-model/fixed-dt in-memory contract, not
portable checkpoints. Default initial state is zero-history LD or closed-shell
SBE. The catalog's optional `initial_states` source/prepare recipes are **not
selected or consumed** by this factory. There is no implicit excited-state load.

Custom resolvers can still be added to the registry, but this factory dispatches
the documented built-in backend/model pairs and known adapter version. No automatic
imports from backend names. Generic saved-TDDFT adapters, general capability schemas,
point/run IDs, MPI placement, and automatic output/restart identity propagation
remain separate milestones. Output helpers do not automatically inspect wrappers;
explicitly save `material.record()` alongside your results.

## Offline preparation and example

Use the same registry with `prepare_inputs`; list all three SBE text assets under
the paths used by the recipe. Resolvers run against the staged directory, so missing
references fail before publication. After transfer, verify the bundle, reload
`materials.json` with this registry, and create the model using that bundle's base.
The catalog alone does not guarantee physical convergence or correct real-material data.

```sh
python -m salmonpy.examples.configured_ld new-ld-run
```

This synthetic two-material example prepares an input bundle, checks independent
states, and saves identity/current data in response.json **outside** that input
bundle. No full SALMON runtime or network access is needed.

DFT configuration and GS→TDDFT identity forwarding are available in [DFT_CONFIG.en.md](DFT_CONFIG.en.md); LD/SBE drive contracts are unchanged.

Explicit ID selection of a pinned saved GS/1 is now available separately via
`tddft_from_catalog`; see [INITIAL_STATES.en.md](INITIAL_STATES.en.md).
Generic material construction itself does not select an initial state.

The separate [PointRun layer](POINT_RUNS.en.md) now carries serial point/run IDs
and saves result mappings. The limitations above still apply to automatic
propagation by the original API alone, MPI and restarts.
