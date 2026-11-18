# Standard TDDFT Python shared models

[日本語](SHARED_TDDFT.md) · Topic: `shared_standard_tddft` · Revision: `2026-11-18.1`

`NativeTDDFTBackend.open_model(initial, dt, context)` returns a
`NativeTDDFTModel`. One model reads an immutable GS/input/PP bundle once and owns
one work copy and one native attachment. `model.new_state(point_id)` creates an
independent electronic state from the original GS. States own wavefunctions,
density/XC, A/current/time, output streams and restore/protection payloads.
A point ID labels a state; it is not a native handle or a material definition.

Use an explicitly built/verified runtime with standard handle/model protocol1.
Imports do not build, download, initialize MPI or fetch inputs. Reusing an
unchanged validated runtime checks its recorded SHA256. The host owns MPI and
calls all subgroup operations collectively on the main thread. Keep the input
bundle/work files unchanged, with work_root outside the initial bundle.

```python
from salmonpy.embedded_tddft import NativeTDDFTBackend

# initial = provider.prepare(context), with verified local GS/PP/input files.
# context.solver_comm is the host's material subgroup communicator.
backend = NativeTDDFTBackend(build_record, work_root,
                            decomposition=(1, 2, 1, 1, 1))
model = backend.open_model(initial, dt=0.003, context=context)
a = model.new_state("film/0")
b = model.new_state("film/1")
try:
    ja = a.step([1e-4, 0.0, 0.0])
    jb = b.step([0.0, 1e-4, 0.0])
    # Rank-local filenames must differ between MPI ranks.
    path = checkpoint_root / f"a-rank-{context.solver_comm.rank}.bin"
    a.checkpoint(path)
    resumed = model.new_state("film/resumed")
    try:
        resumed.checkpoint(path, restore=True)
        # Same next A gives the same current; no accepted steps are replayed.
        ja2 = a.step([2e-4, 0.0, 0.0])
        jr2 = resumed.step([2e-4, 0.0, 0.0])
    finally:
        resumed.close()
finally:
    b.close()
    a.close()
    model.close()
```

`checkpoint_root` must exist; save refuses an existing file. Direct load requires
a fresh state, matching binary/inputs/dt/capacity/MPI partition and authenticated
files. Standalone native files are host-managed, not portable archives; use the
normal mixed sample's aggregate checkpoint for provenance checks. Closing the
model prevents new states but lets live states finish. Close every state and model
explicitly before host MPI finalization; closing one sibling leaves others live.
An open idle model can create another state. This default API retains one model
lifetime per process. Use the [explicit registry](TDDFT_REGISTRY.en.md) for
sequential switching between different-input models; concurrent threads are unsupported.

The supported profile is periodic, full-k (no symmetry reduction), unpolarized,
fixed-ion PZ, TB-mBJ (automatic/fixed c), Libxc-PZ or Libxc-PBE, transverse
middlepoint propagation. **dt, capacity, initial A=0, decomposition and immutable
material input are shared.** Different materials require different MPI groups.
x/y/z and rotating A are accepted. A is Cartesian A/c, dt atomic time, and Python
returns charge-current density j in atomic units. The handle result already has
that sign; do not negate it again. RT text columns retain SALMON conventions.
Only per-owner RT/current and energy streams are supported here; optional global
analysis streams require the existing separate-state route. Native handle and
legacy initialization must not be mixed in one process. Initialization, allocation,
I/O and MPI failures do not gain a general recovery guarantee.

## Ordinary multiple-point placement

Opt in to session ownership with `shared_model=True` and set
`assemble_native_layers(..., points_per_group={"si-tddft": 2})`. The material
binding's backend is this standard backend; `ranks_per_point["si-tddft"]` is the
size of each shared material group and must equal the decomposition product.
The planner keeps point IDs, layer order/thickness and LD regions unchanged;
it groups compatible electronic points instead of removing points. Default
one-point groups and `create_state` remain unchanged. ResponseSession uses
`open_group` to own/close each model and every state it returned.

The installed ordinary mixed sample has a two-point TDDFT layout:

```sh
mpiexec -n 4 python -m salmonpy.examples.mpi_mixed_native \
  STANDARD_BUILD SERIAL_RUNTIME INITIAL_BUNDLE NEW_OUTPUT \
  --layers mixed_shared_tddft_layers.json --tddft-ranks 2 \
  --tddft-points-per-group 2 --steps 128 --dt .003 --protect-steps
```

Copy the JSON from the installed `salmonpy/examples/` directory, or use its full
path. This uses field1 + SBE1 + shared-TDDFT2 ranks, with two local LD materials.
Omit grouping to use six ranks (two independent TDDFT groups). The accompanying
`mixed_shared_tddft_hse_both_layers.json` also allocates two HSE points; choose
`--hse-points-per-group 2` and compatible HSE arguments to share both separately.
No cross-material model sharing is implied. Resource budgets warn instead of
stopping a computable calculation; actual rank-plan mismatch remains an error.

Use `--direct-save-at 63` and a fresh invocation's `--direct-resume CHECKPOINT`
with the same topology to save/restore field, LD, SBE and every TDDFT point without
history replay. Protected restoration stages all owners before publication and
retains original states on a handled staging failure. Each TDDFT state's current
metadata is validated before native staging. Numeric restart does not merge old
text streams. See [coupled](COUPLED.en.md) and [distributed execution](DISTRIBUTED.en.md).

## Reproducible implementation checks

`tools/test_python_shared_tddft_mpi.py BUILD INITIAL_BUNDLE NEW_OUTPUT --size 2`
requires six MPI ranks (three groups). Two shared states follow distinct 3D A
waveforms against two independent legacy controls; fresh states load direct files
and continue, invalid NaN/Inf/complex A is rejected, and sibling/idle lifetimes
are checked. Add `--rgrid` for spatial instead of orbital MPI splitting.

Local checks cover PZ, Libxc-PZ/PBE and automatic/fixed TB-mBJ at MPI1, orbital2
and spatial2. The ordinary LD/SBE/two-TDDFT mixture compares independent/shared
execution and direct resume, including a late rank-local restore failure/retry.
These small existing fixtures verify ownership and numerical parity, not
material convergence. No new GS data is collected and no Fortran kernels change.
See [native C contract](STANDARD_HANDLES.en.md) for lower-level entry points.

## Different material inputs: current boundary and remaining work

Different TDDFT inputs can run in **separate MPI material groups**. Each group can
share multiple states of its own immutable material. Add `--tddft-materials FILE`
to the ordinary `mpi_mixed_native` sample with `--layers` to declare arbitrary
TDDFT material IDs, instead of only the built-in `si-tddft` ID. Example:

```json
{
  "schema": "salmonpy.native-tddft-materials/1",
  "materials": {
    "silicon-pbe": {
      "bundle": "./pbe-initial",
      "input_file": "rt.inp",
      "decomposition": [1, 2, 1, 1, 1],
      "points_per_group": 2
    },
    "silicon-tbmbj": {
      "bundle": "./tbmbj-initial",
      "decomposition": [1, 2, 1, 1, 1],
      "points_per_group": 2
    }
  }
}
```

Paths are relative to the catalog, or absolute local paths. Prepare complete
GS/PP/input bundles on a login/management machine and transfer them before a job;
no internet, downloads or compilation are required by the loader. Default input
is `rt.inp`, decomposition is `[1,1,1,1,1]`, width is1. Width1 uses a separate
legacy state; width>1 explicitly uses the shared profile above. Decomposition is
(k, orbital, x, y, z); its product sets ranks per material group. Different IDs
never share a native model, even if their input bytes happen to match. All groups
use the coupling dt; this does not introduce subcycling.

Layer `material_id` selects one catalog ID. The provided
`mixed_distinct_tddft_layers.json` places SBE1 + PBE2 + TB-mBJ2 points, a gap and
two local LD materials. With width2/inner MPI2, launch6 ranks (field1+SBE1+PBE2+
TB-mBJ2); width1 requires10. Copy `distinct_tddft_materials.json` alongside your
local bundles and edit references before use:

```sh
mpiexec -n 6 python -m salmonpy.examples.mpi_mixed_native \
  STANDARD_BUILD SERIAL_RUNTIME DEFAULT_INITIAL NEW_OUTPUT \
  --layers mixed_distinct_tddft_layers.json --tddft-materials materials.json \
  --steps 128 --dt .003 --protect-steps
```

The positional `DEFAULT_INITIAL` remains the compatibility fallback for the
built-in `si-tddft` or HSE ID. It is not read for custom IDs. A catalog entry for
`si-tddft` overrides that fallback; unused default flags/paths are not restart
conditions. Per-material input hashes, initial-state identity, MPI plan and memory
advice use each actual bundle. `sbe-point`, `hse-point`, `drude`, `lorentz` remain
reserved built-ins. A shared catalog can include unused declarations. Unknown
fields, duplicate keys and invalid decomposition reject before native setup;
resource estimates stay advisory. The existing aggregate direct save/resume
and protected failure/retry work across all material groups.

**Same-process different-input models are available through the explicit registry.**
See [registration/switching](TDDFT_REGISTRY.en.md) and `examples/tddft_models.py`.
Inputs, initialized controls, optimization tables and XC scratch switch per model;
wavefunctions, clocks and checkpoint resources remain state-owned. Existing default
`open_model` and ordinary layout contracts are preserved.

Remaining work includes automatic heterogeneous registration from ordinary material-ID
layouts, full GS contexts and specialized SOC/+U/DC resources. Attachment, timers
and synchronous compatibility views remain process infrastructure; concurrent
threads are unsupported. Existing Si fixtures with different cells, atom counts,
16³/8³ grids, k counts and five XC choices establish implementation parity, not
multi-element physical accuracy or convergence.

## Model-owned initialization input (2026-11-18)

The standard native runtime now owns227 initialization values in `s_model_inputs`:
cell/grid/k points, PP/XC references and coefficients, electron/occupation/mixing,
propagation/dt/capacity/parallel settings, output choices, input units and atom
coordinates/species/names. Allocatable inputs are deep copies; cloning or changing
the legacy parser arrays cannot change the captured model. An uncaptured record
cannot be activated. Unit conversion tables are reconstructed from stored unit
selectors; already converted physical inputs are not converted again.

A native `rt_model` owns this record and initial A. `open_initial_state(owner,...)`
receives its owner explicitly, installs a synchronous compatibility view for the
existing initialization kernels and assigns the context a borrowed model pointer.
Sibling states do not own/deallocate that model; an idle open model can still
create states from the original inputs. Numerical wavefunctions, live MPI/Libxc
owners and RT checkpoint payloads are not copied into this input record.

The owned-input record, explicit initialization boundary, parser workspaces and
defaults now support explicit registration/selection and Python model handles.
The default API profile, numerical order and restart format remain unchanged;
the registry adds a separate ABI. See the registry guide above.

Developer validation uses `tools/test_model_input_ownership.py SOURCE NEW_OUTPUT
--fc gfortran`, compiling the actual record and a two-input deep-copy/resize/release
driver. `tools/test_standard_tddft_handles_mpi.py` compares the new library against
unchanged old controls; `test_standard_checkpoint_finite_mpi.py` checks NaN/Inf in
each restored A component with a rank-local corrupt shard and valid continuation.
The latter retains scalar IEEE inquiries on every checked element. New runtimes
require source/object prohibition audits; unchanged validated binaries use recorded
SHA256. No new GS fixtures are required. Cross-compiler validation is recorded
separately; a local GNU build is not a Fujitsu compilation certificate.

## Model-owned parser workspace and defaults (2026-11-18)

`inputoutput.read_input(path, parser)` accepts an `s_input_parser` owned by each
standard native model. It owns namelist status, coordinate-block flags, scratch
connections, a fresh `defaults` record and final `values`. The original default
assignment order/formulas and physical conversion order are retained. Defaults are
created after that input's unit selection; preceding atom arrays and coordinate
mode are excluded before capture. Effective worker dt/capacity/decomposition are
captured in final values after worker setup and coordinate reading.

The owned route uses unnamed Fortran scratch connections, keeping each split
input alive until consumed, then closing/resetting all handles. It does not
replace preexisting `.namelist.tmp`/coordinate/reentrance files. Parser read status
and coordinate flags agree across MPI ranks; unit numbers remain root-local.
`parser.log_path` chooses a diagnostic log (default `variables.log`); callers
reading several inputs in one directory can set distinct log paths. Input/file
references retain the existing working-directory semantics. No download or
compilation is triggered by reading input.

Legacy `read_input()`/`read_input(path)` callers retain a deliberate legacy
workspace and named temporary files: the MD velocity reader still consumes those
coordinate files later. This preserves batch/MD behavior. The new record covers
selected standard scalar-TDDFT model fields; it is not a complete FDTD/SBE/DC
configuration owner. Parsing and numerical initialization still use a synchronous
legacy value/unit view, so this is sequential input isolation, not concurrent
thread entry or heterogeneous live propagation in one process. Existing error
policy is unchanged; general parser I/O/MPI failure recovery is not added.

Developer test `tools/test_parser_ownership.py BUILD_JSON NEW_OUTPUT --fc mpifort
--mpiexec mpiexec --ranks 2` uses the actual Fortran parser, with matching compiler
module files from the build directory in the receipt. It checks A→B→A, reused
workspace, different units/atom counts/Cartesian-external-reduced-absent sources,
omitted XC/k/pulse/mixing defaults, two interleaved scratch sources and closed
connections. Diagnostic logs match fresh-process legacy parsing exactly; existing
named temporary-file sentinels remain unchanged. Repeat at MPI1/2. No GS data is
needed. Native old/new XC/MPI/direct/mixed comparisons verify numerical parity.
