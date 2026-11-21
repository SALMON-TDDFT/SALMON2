# Per-material configuration catalog

[日本語](CONFIG.md)

`salmonpy.config` is a configuration-only preview. It accepts Python dictionaries
or UTF-8 JSON with schema `salmonpy.materials/1`, not the planned full simulation
schema. It does not create a solver, read a structure/GS, download data, initialize
MPI, or certify capabilities. No built-in calculation implementation is registered automatically.
Run `python -m salmonpy.examples.material_config` for a solver-free example.

## Entry points

- `MaterialRegistry()`: explicit local registry.
- `registry.register(name, version=..., models=[...], resolve_parameters=...)`:
  register a trusted Python validator. Duplicate names/models are errors.
- `resolve_materials(data, registry=..., base_directory=...)`: Python JSON values;
  base directory must be explicit. Tuples, arrays, Path objects and custom objects
  are not JSON values; convert them explicitly before calling.
- `load_materials(path, registry=...)`: same validation, with paths relative to the
  input file's directory. No dynamic imports from JSON.
- Returned `MaterialCatalog.materials`, `.initial_states`, `.definition_digests`
  are recursively immutable. `.as_dict()` returns a detached JSON-safe record
  with schema `salmonpy.materials-resolved/1`. This record is not a recipe or
  restart and cannot be passed back to `load_materials`.

## Input contract

```json
{
  "schema": "salmonpy.materials/1",
  "materials": {
    "film": {"model": "lab:response", "backend": "local",
             "parameters": {"coefficient": 2.0}, "label": "First layer"}
  },
  "initial_states": {
    "film-ground": {"material_id": "film", "prepare": {"kind": "equilibrium"}}
  }
}
```

`initial_states` and `label` are optional. Each initial-state entry has an existing
`material_id` and exactly one object-valued `source` or `prepare`. These objects
are recorded as recipes only: their keys, external file existence, hashes and
physical compatibility are not validated yet. They must not be treated as
prepared initial states. Placement, MPI assignments and other unknown envelope
keys are rejected. IDs preserve case and Unicode; empty IDs, edge whitespace,
control characters, duplicate JSON keys and non-finite numbers are errors. IDs
are not filesystem paths.

## Custom resolver

The callable has signature `resolve_parameters(model, parameters, *, base_directory)`
and returns a plain JSON object. The calling program explicitly trusts this Python code;
it is not sandboxed. The callback must validate unknown parameters, units,
ranges and combinations, normalize units/defaults, preserve provenance, and return
resolved values. Do not launch solvers, allocate states or modify global settings
here. Inputs/outputs are detached; resolved definitions do not alias caller data.
Schema/ID/reference errors are detected before any callback. Later callback
failure may follow earlier callback execution; callback side effects are not
rolled back. ValueError/TypeError are wrapped in InputError with the material
path and chained cause; other exceptions, including OS errors, propagate.

There is no schema-derived AI parameter generation yet. Registering a resolver
is not a declaration of MPI, checkpoint, threading or multi-state support. Use
[MATERIALS.en.md](MATERIALS.en.md) for the existing executable model contract;
see the built-in adapter guide below for LD/SBE configuration-to-model binding.

## Definition hashes

SHA256 covers calculation implementation name/version, model and resolved parameters, excluding
material ID and label. Canonical encoding v1 uses sorted keys, compact JSON,
UTF-8 and finite numbers (`1` and `1.0` remain distinct). calculation implementation versions must
change when resolution semantics change. File contents and callback code are
not automatically hashed. Resolvers must include data hashes for content-aware
identity and portable references instead of absolute machine paths. These are
configuration hashes, not physical validation or cross-platform result promises.

## Separate structure and grid adapters

ASE Atoms and local CIF/extxyz files are supported through the separate STRUCTURE.en.md adapter. Supercells/complex structures belong in external files; retain the file's
cell and atoms without implicit repetition, rotation or standardization. Record
source, format, frame index where needed, content hash and unit convention.
MateriApps is a discovery portal; data are obtained from the selected provider
(e.g. COD/Materials Project) under its own access rules, then read through an
explicit adapter. Structure data do not determine XC, pseudopotentials, charge,
magnetism or convergence parameters automatically.

An explicitly registered resolver can use grid_from_spacing for per-material maximum electronic spacing with an explicit unit. For supported orthogonal periodic cells, start from N_i = ceil(L_i/h_i),
then increase N_i to meet the selected calculation implementation's minimum/alignment constraints.
Report requested spacing, adopted counts and actual L_i/N_i. Do not silently
alter a finalized grid when MPI resources change. Counts remain an expert option;
counts and spacing together are rejected by the ASE DFT/TDDFT entries to avoid two sources of truth.
Nonorthogonal/isolated cells require specific to each calculation implementation geometry validation. FDTD
spacing is separate. This rounding is not a convergence guarantee. Direct DFT still requires counts; ASE DFT/TDDFT entries also accept the explicit width/unit alternative described below.

## Offline structure references

Do not require database queries/downloads from compute nodes. Stage structure/
supercell files and all referenced data on an allowed connected calling program, recording
source, version and SHA256 before transfer. Use local references in calculation
inputs; verify hashes and relative paths afterward. Separate discovery through
MateriApps or other services from job execution; missing files must not trigger
automatic runtime downloads. See the offline compute-node section of
[RUNTIME.en.md](RUNTIME.en.md) for wheels and runtimes as well. Use the separate structure adapter described in STRUCTURE.en.md; the catalog does not collect or verify files automatically.

Explicit local collection and post-transfer verification are available in [PREPARATION.en.md](PREPARATION.en.md). Structure interpretation and automatic dependency discovery remain separate.

Standalone structure/spacing helpers and an explicitly registered example are now available: [STRUCTURE.en.md](STRUCTURE.en.md). The catalog itself does not automatically load structures.

LD/SBE explicit registration and model construction are available in [BUILTIN_MATERIALS.en.md](BUILTIN_MATERIALS.en.md).

Explicit ID selection of a pinned saved GS/1 is now available separately via
`tddft_from_catalog`; see [INITIAL_STATES.en.md](INITIAL_STATES.en.md).
Generic material construction itself does not select an initial state.
