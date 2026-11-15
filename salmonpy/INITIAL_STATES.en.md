# Select a saved ground state by ID

`tddft_from_catalog` selects a local SALMONpy-generated `salmonpy.gs/1` ground
state for an explicitly configured DFT material. It does not run SCF, download,
compile, or start a native worker. `new_state(dt=...)` starts an independent worker.
Install optional ASE/f90nml and prepare a verified local DFT/TDDFT runtime first.

Add this entry alongside `materials` in a `salmonpy.materials/1` recipe:

```json
"initial_states": {
  "Si-ground": {
    "material_id": "Si-starter",
    "source": {
      "format": "salmonpy.gs/1",
      "path": "saved",
      "manifest_sha256": "REPLACE_WITH_SHA256_OF_GROUND_STATE_JSON"
    }
  }
}
```

The path is relative to the explicit base directory, not a material ID or URL.
Use a complete saved GS directory, including `ground-state.json` and all files
listed in its `files` object. Preserve the material's structure, pseudopotential,
license and provenance files too. `prepare_inputs` can collect these with an
explicit destination-to-source file mapping; verify the bundle after transfer
with `verify_inputs`, using its separately recorded manifest digest. See
[offline preparation](PREPARATION.en.md). Compute-node execution needs no network.
Native restart portability is not promised: this selector requires the exact
executable SHA256 used for GS creation and an otherwise compatible runtime host.

```python
from salmonpy.builtin_materials import builtin_material_registry
from salmonpy.config import load_materials
from salmonpy.initial_states import tddft_from_catalog
from salmonpy.runtime import load_runtime

engine = load_runtime(runtime_directory, required=('dft', 'tddft', 'single-rank'))
catalog = load_materials(bundle / 'materials.json', registry=builtin_material_registry())
with tddft_from_catalog(catalog, 'Si-ground', base_directory=bundle,
                       executable=engine, max_steps=20, threads=1) as model:
    with model.new_state(dt=0.02) as first, model.new_state(dt=0.02) as second:
        j = first.step([0., 0., 1.e-5])  # A/c, atomic units; second is unchanged
    metadata = model.record()
```

The loader checks the pinned manifest, native file hashes, SCF completion,
ordered geometry, cell, electrons, states, grid, k mesh, XC, pseudopotential bytes,
actual scientific GS input and copied snapshot. Reordered or transformed cells
are not silently converted. Malformed metadata and mismatches are rejected.
The digest anchors bytes, not scientific convergence or authenticity on its own.
The short Si starter is a connection regression, not converged optical physics.

`material_id`/`definition_digest` and `initial_state_id`/`initial_state_digest`
are retained on the returned model and its states. The latter digest is the GS
manifest SHA, not a dynamic checkpoint. Save `record()` explicitly with results.
Point/run IDs, automatic disk identity and MPI dispatch remain separate work.
Generic material construction does not select an initial state. `prepare` recipes,
arbitrary SALMON restart directories and GS/2 are not supported by this selector.

The runnable sample accepts an NPZ with scalar `dt` and finite `A` of shape
`(steps+1, 3)` in atomic units, with `A[0]=0`:

```sh
python -m salmonpy.examples.catalog_tddft inputs/materials.json Si-ground \
  /path/to/local/runtime waveform.npz /new/output/directory
```

It writes current arrays and ID metadata; current is charge-current density j = -J_SALMON. Use [FIELDS.en.md](FIELDS.en.md) for signs and timing. Each state has
its own process and snapshot storage; this is not shared-memory TDDFT or a memory
budget guarantee. Keep thread and concurrent-worker counts within your resources.

The separate [PointRun layer](POINT_RUNS.en.md) now carries serial point/run IDs
and saves result mappings. The limitations above still apply to automatic
propagation by the original API alone, MPI and restarts.
