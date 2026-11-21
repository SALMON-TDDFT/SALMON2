# SBE state ownership and persistent solver handles

[日本語](SBE.md)

Keep using SBEInitialState/new_state/step(A). A is Ac=A/c; public current is charge
current density j. matter_current retains SALMON matter current J=−j. Initial
arrays use Fortran band/k layout and atomic units, with even closed-shell electron
count. Correction order0–3 and nonlocal velocity correction are fixed per model.

```python
# model is an existing SBEInitialState
with model.new_state(dt=0.02, initial_A=(0., 0., 0.)) as state:
    current = state.step([0.001, 0., 0.])
    saved = state.snapshot()
    density = state.density       # explicit dense export only when requested
    state.step([0.002, 0., 0.])
    state.restore(saved)
```

Model matrices, weights, band differences and occupations are prepared once per
initial model and shared by its independent states. Each state owns density,
previous A and current. A normal step sends A3 and receives J3; model matrices no
longer copy each call. density/snapshot explicitly export the dense matrix.
The existing velocity-gauge calculation routines, arithmetic order and outer-k OpenMP directives
are unchanged. This does not newly enable OpenMP in wheel builds.

## Loading and numerical method

```python
from salmonpy.sbe import SBEInitialState
model = SBEInitialState.from_salmon(
    directory, sysname, bands=bands, kpoints=kpoints, electrons=electrons,
    lattice=lattice, correction_order=0, nonlocal_correction=False)
```

Lattice vectors are columns. The text loader accepts single-spin orthogonal
k/eigen/tm exports with explicit esp[a.u.], rejecting eV/nonorthogonal formats.
Array construction accepts a right-handed cell, (band,k) energies and
(band,band,xyz,k) matrix elements. Occupations are2 in each first electrons/2
band, independent of exported occupations. active_bands defaults to all bands.
Four sequential Taylor terms use midpoint A, then current uses endpoint A.
No exact unitarization, adaptive dt or norm correction is added. Check dt,band/k
convergence and physical model validity separately. Length gauge and external
excited-state restart loading remain unsupported; this initial loader builds GS density.

## Failure and lifetime

Compute on a trial density, check every real/imaginary density component and
current with scalar IEEE predicates, then swap accepted buffers. Numerical
rejection preserves density,A,current and clock. Snapshot restore validates
before committing; existing same-model/fixed-dt/time-origin compatibility remains.

close releases the calculation state and is idempotent. Closed states reject step,
snapshot,restore and density export; the last clock/current remain readable.
Use context managers or explicit close. GC provides supplemental cleanup for
existing callers; live states keep their solver model alive. solver model cache
lifetime follows the Python initial model.

Used initial models remain data-only for asdict,pickle and deepcopy. Process-local
handles are excluded and rebuilt at the destination. This is not a contract to
pickle evolving states or opaque snapshot tokens for restart.

Models/states are not thread-shareable; solver registry calls are serialized.
Existing bind_material point-MPI distribution remains available.

## Memory and measurements

For N=nbands²*nkpoints, two complex128 densities consume32N bytes per state:
accepted data plus rollback workspace. One solver read-only model copy is held
in addition to Python initial arrays, approximately104N+16*nbands*nkpoints+
8*nkpoints bytes. Structs,small vectors,calculation routine workspaces and snapshot temporaries
are excluded.

Resident memory need not decrease: many idle states retain extra trial buffers.
Check both copy costs and memory budget. There is no arbitrary machine-wide cap;
integer size overflow and explicit allocate failures are rejected. Automatic
workspace allocation failure inside unchanged calculation routines is not guaranteed recoverable.

Manufactured600step,3repeat isolated-process cases reduced modeled explicit copies
by about9×, with measured speedups around1.02–1.03×. Current/density traces were
byte-identical in every comparison. Copy estimates follow explicit source
assignments, not hardware profiling/network traffic. Peak RSS is the whole worker
high-water mark; small differences are not conclusive memory improvements.
These tests establish neither material physics nor large-system performance.

## Building and scope

Matching wheels contain the new compiled calculation routine; import never builds it. Source
checkouts with old binaries require an explicit rebuild and never silently fall
back to stateless stepping. To preserve an existing fixed binary, build/install
separately and use tools/test_source_sbe_handles.py LIBRARY for validation.
LD handle migration, multiple TDDFT called directly from Python contexts and a public SBE ABI1 plugin
remain separate work. Internal sp_sbe_* symbols implement the Python adapter,
not the public material ABI1 query table.
