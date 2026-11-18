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
An open idle model can create another state. There is one native electronic model
lifetime per process, not arbitrary different-input models or concurrent threads.

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
