# Fast-PR HSE and symmetry

[日本語](FAST_HSE_SYMMETRY.md)

Fast-PR HSE is the primary HSE development and validation path. Symmetry-aware
representative k points, star reconstruction, orbit reuse and ACE are central to
that path. Conventional HSE is not a supported replacement engine. Standard
SALMON symmetry reduction has different arithmetic and is not used to redefine
the fast-PR algorithm.

| Layer | Current scope |
|---|---|
| Fast-PR exchange core | Independent symmetry maps/phase caches/FFT owners; tested scalar periodic CPU, full q, one/two internal k ranks |
| Symmetry/orbit/ACE comparison | Reduced4 k representatives vs independently enumerated full8 k points; C2z and mirror-x groups, different cells/omega/sources |
| Internal Fortran TDDFT | Explicit opt-in, one initialization/process/internal-k1/2, Taylor/PT-CN including learned modes; direct/staged resume tested |
| Python, C and FDTD interfaces | Explicit new-protocol allow_symmetry opt-in; groups1/2 Taylor/PT-CN, mixed tests; default remains denied |
| Checkpoint | Internal exchange V2 and explicit internal TDDFT direct/staged resume tested; bounded public protected mixed resume tested |
| Symmetry internal MPI k2 | Explicit new capability; full/reduced comparison and direct resume verified |

## Independent exchange owners

The internal `exx_bind_system` core can now admit symmetry within the bounded
one/two-k-rank configuration. Each system retains its own k-point configuration and each
HSE owner retains its own symmetry map, phases and FFT resources. The full-grid
check uses that system's ngrid, primitive cell and spacing, without selecting from
shared num_rgrid. Dimension/count mismatch or an unsupported owner layout rejects
before refresh/action mutation when using the declared rt_status entry. Other
invalid symmetry operations/atoms/weights may still terminate solver execution.
The input checking result is not a general numerical-failure recovery guarantee.
The core admits Taylor/ACE with off/strict/learned factors and PT-CN with off/learned factors. New combinations require the history capability described below; full-action remains factor-off test-only.

Create a fresh owner when changing cell/material/group settings; mutating live
owned configuration is not a supported way to change the model. Update source
orbitals through the existing refresh contract. Do not intrinsically clone FFT
plans or communicator owners. Automatic GS binding policy remains unchanged.

Local tests use independent C2z and mirror-x owners with different cells,
screening parameters and orbitals, interleaved refresh/action, grid/layout
rejection and close/reopen. Full-action and ACE, ordinary star expansion and fast
orbit reuse pass OMP1/2 with one/two MPI processes. MPI processes use different model values and call schedules; outputs are checked
for finiteness as well as agreement. These are small exchange-core
regressions, not a converged material result or reduced-k TDDFT validation.
The solver probe is packaged under
`_runtime_source/tests/native/test_hse_symmetry_resources.f90`; compilation is
explicit and requires matching audited fast-PR modules/library.

## Owned field validation

The internal `symmetry_check_field(map, A, status)` and
`exx_validate_vector_potential(system, A, status)` now check uniform A/c against
the retained Cartesian operation matrices, without reading shared symmetry
settings. Status is 0 for compatible, 1 for invalid/unready input, and 2 for a
broken field little group. For nonzero A, the maximum component normalizes the
vector before a 1e-10 directional tolerance is applied; arbitrarily small
amplitude does not excuse a forbidden direction. Zero A still requires a valid
map. These internal Fortran entries do not mutate state or communicate.

The TDDFT adapter invokes the owner check before numerical backup/propagation
and agrees guard failure collectively. Tests cover independent C2z/mirror-x
owners, a rotated cell, tiny/huge fields, invalid values, stale configuration
and closed owners. This checks the field/map boundary, not all symmetry data
or solver initialization. Default opening a calculation process or the C interface and coupled-checkpoint prohibitions
remain unless the explicit public opt-in below is selected; this is not unrestricted support. The explicit
solver probe is `_runtime_source/tests/native/test_hse_symmetry_field.f90`.

## Internal symmetry checkpoint

The existing internal `hse_checkpoint_write/read/capture/install` now support a
value-only symmetry map/cache for the admitted scalar CPU/full-q/internal-k1/2
exchange core. Non-reduced owners still write V1 byte-for-byte. Symmetric owners
write same-binary V2: retained map, calculation routine identity, source/ACE/history values and
optional complete phase cache. The caller authenticates the archive and model;
this raw stream is not a standalone portable interchange format.

Initialize an independent matching owner/operator before reading (operator
construction, no accepted time-step replay). The reader compares all static map
values and trusted shapes before reading the bounded optional cache. It stages
before installation; validate-only and staged reads leave the live state intact.
A truncated stream, differing map/phase partition or malformed cache rejects.
Installation moves value arrays only, leaving live FFT/C/MPI owners untouched.

The solver probe `_runtime_source/tests/native/test_hse_symmetry_checkpoint.f90`
checks full/ACE, star/orbits, cache on/off, OMP1/2 and MPI1/2 processes using
independent COMM_SELF cores. Restored action/energy and subsequent refresh/action
checkpoint bytes match continuous core execution, including reopening an owner.
This is exchange-core persistence; the separate internal TDDFT tests are described below.
Ordinary/default initialization through Python or C stays symmetry-denied; the explicit opt-in is described below. Current learned/PT-CN symmetry, internal-k MPI and ordinary mixed failure/resume coverage are described below.

## Explicit internal TDDFT validation

An internal Fortran `rt_open(..., allow_symmetry=.true.)` opt-in now connects the
exchange core to the complete TDDFT checkpoint path. It is bounded to scalar
periodic CPU/full-q, one/two internal k ranks, Taylor/ACE or PT-CN and a retained
group preserving the actual A/c at each step. Full-action is test-only. One initialization per process is
allowed; initialize a fresh process for each independent TDDFT lifetime.

Density/current symmetry tables are still process-common and immutable. Their
flags, operations, lattice and inverse matrices are checked against owned inputs
before step/current/checkpoint and again before a staged checkpoint commit.
Drift rejects without applying the candidate; after correcting settings the
candidate can be committed. This does not promise heterogeneous same-process
contexts or controlled recovery from malformed initialization assets.

The actual-A owner check returns internal status34 for a broken field little
group, before backup or propagation. An internal caller can correct A and retry.
The new public protocol below now declares this status; do not infer it for old libraries.
The ordinary init/serial worker remain denied; the separate public opt-in enables the bounded configuration.

The solver `test_hse_symmetry_tddft.f90` probe and companion
`prepare_fast_hse_symmetry_fixture.py` / `test_fast_hse_symmetry_tddft.py` tools are
packaged under `_runtime_source/tests/native` and `_runtime_source/tools`.
Compile explicitly against matching MPI fast-PR modules/library, then run the
harness with `PROBE PSEUDO NEW_OUTPUT --ranks 2 --threads 2`. Processes attach
independent COMM_SELF groups, not reduced internal-k MPI2. Synthetic8^3/full8-k
versus reduced4-k C2z/mirror-x cases compare charge current/energy; odd/even direct
and staged resume reproduce internal traces and complete final checkpoint bytes
without replaying accepted steps. These are implementation tests, not a GS or
optical-convergence recipe; no download, SCF or build is performed by the fixture.

Symmetry archive copies mark nonserialized overlap caches dirty. Live flags and non-reduced V1 bytes are unchanged.

## Public native/Python and mixed opt-in

Use `NativeFastHSEBackend(build_record, work_root, allow_symmetry=True)`.
The default False selects existing `sp_fast_hse_init`. True requires
`sp_fast_hse_init_symmetry` and `sp_fast_hse_symmetry_protocol()==1`; the new init
has the same `(const char *path, double dt, int capacity, double J[3])` arguments.
The opt-in is included in calculation implementation identity, admits one/two internal k ranks, and never
changes the supplied operation matrices. Missing/unknown protocols reject before
acquiring a solver initialization. Learned/PT-CN symmetry requires the history capability; full-action remains test-only.

Only unanimous finite status34 plus the specific protocol becomes a Python
`InputError` preserving the accepted state. Correct A before retrying. A reduced
model must match the actual local FDTD field, not only the incident polarization.
Do not project A, drop tiny components or switch engines automatically. For a whole
loop, use `CoupledStepper.step_protected`/`--protect-steps`, and retry only when
`StepRecoveryError.original_preserved` is true and cleanup succeeded.

The ordinary `mpi_mixed_native` sample accepts `--hse-symmetry` and requires
`--hse-reference-bundle FULL_K_BUNDLE` with an explicitly supplied equivalent full-k
state for the existing serial fast-PR comparison. No arbitrary GS conversion is
performed. Select `--hse-ranks 1` or `--hse-ranks 2`. The layout/LD/SBE/semilocal-TDDFT APIs and
direct save/resume flags remain the same. Tested five/six-rank mixtures contain
field1/2, two field-local LD materials, SBE1, PBE-TDDFT2 and symmetry HSE1.
Odd/fault/protected restore tests and fresh direct resume preserve trajectories
without accepted-step replay. Fixtures are synthetic; compare material/field
compatibility and convergence separately.

## Stable near-zero-gradient XC evaluation

Libxc5.2.3 HSE06 mixes WPBEh(0) − 0.25 WPBEh(ω) + PBE correlation.
Separate differentiation of the E1 and log terms in WPBEh(0) cancels large terms
at nearly zero gradients. The original phase0 synthetic full8/reduced4-k case
showed a current difference about5.7e-10; smaller time steps did not remove it.

The current complete fast-PR source evaluates the same formula in u=s² for
reduced gradients s≤1e-4, combining the cancellation before differentiation.
Discarded E1-series terms are below double precision. Screened exchange uses its
analytic quadratic limit only below roundoff in u; correlation retains Libxc.
Points outside this region retain original Libxc output bytes. No density,
gradient or field floor/smoothing, tolerance relaxation or conventional HSE is
used. Spin scaling and fixed-core density/gradient contributions are checked.

The unchanged original phase0 case now passes full/reduced current and energy
comparisons over32 steps, FDTD coupling, failed-step recovery, protected restore
and fresh-process direct resume. Ordinary PBE arrays remain byte equal; ordinary
fast-HSE Taylor/PT-CN and mixed arrays differ only at roundoff. Synthetic numerical
tests do not certify physical convergence for a real material.

Explicit verification (no build/download):

```sh
python tools/test_hse_low_gradient.py BUILD_RECORD COMPILED_NATIVE_PROBE NEW_OUTPUT
mpiexec -n 2 python tools/test_fast_hse_low_gradient_mpi.py \
  BUILD_RECORD FULL_K_BUNDLE REDUCED_K_BUNDLE NEW_OUTPUT
# Replay the original local FDTD field: --drive ORIGINAL_ARRAYS_NPZ --point-index 1
```

Compile `tests/native/test_hse_low_gradient.f90` explicitly against matching new
`hybrid_semilocal.mod` and library. The checks cover80-digit energy differences,
density/screening/branch boundaries, spin/core paths and ordinary-gradient bytes.
See `runtime/NOTICE.md` for formula provenance/licensing. Pinned sources and old
binaries are not modified. Explicitly rebuild **both** the full-k serial reference
worker and MPI library from the new source. Symmetry protocol1 alone does not
identify this numerical fix.

`probe_hse_low_gradient.py BUILD_RECORD NEW_JSON` intentionally still evaluates
raw Libxc. Its `low_gradient_sensitivity_detected` can remain True with a fixed
wrapper. Do not interpret that diagnostic as wrapper success/failure: retain the
new verification results and source/binary hashes. Historical phase.23 integration
fixtures remain; this fix is verified on the original phase0.

Next: broader initialization/partitions and converged real-material tests. Fast-PR-only policy,
warning-only memory budgets and paused molecular work remain unchanged.


## Reduced-k internal MPI (two ranks)

Two internal ranks additionally require `sp_fast_hse_symmetry_kmpi_protocol()==1`.
Old libraries reject before acquiring a solver initialization. The scalar periodic
CPU/full-q configuration is unchanged. Each rank owns contiguous
representatives and their star phases. Empty ownership, real-space/orbital
partitioning remains outside this validation; local representative histories are covered below.

Archives store rank-local phase caches rather than replicating all k points per
rank. Restore compares these to independently initialized ownership and the
pinned MPI rank count/partition. **Use the same MPI decomposition on resume;
checkpoint repartitioning is unsupported.** Consensus uses MAX for recoverable
status34, avoiding multiplication by the number of ranks.

```sh
mpiexec -n 4 python tools/test_fast_hse_low_gradient_mpi.py \
  BUILD_RECORD FULL_K_BUNDLE REDUCED_K_BUNDLE NEW_OUTPUT --group-size 2
mpiexec -n 4 python tools/test_fast_hse_symmetry_mpi.py \
  BUILD_RECORD REDUCED_K_BUNDLE NEW_SAVE --mode save --group-size 2 --checkpoint CHECKPOINT
mpiexec -n 4 python tools/test_fast_hse_symmetry_mpi.py \
  BUILD_RECORD REDUCED_K_BUNDLE NEW_RESUME --mode resume --group-size 2 --checkpoint CHECKPOINT
```

Two independent two-rank groups use different screening parameters. Fresh7→32
resume matches continuous current/energy traces and final checkpoint bytes without
accepted-step replay. Initialization still constructs the required operators.
The unchanged phase0 weak local-field case compares full8/reduced4 k points with
OMP1/2. The ordinary mixed example accepts `--hse-ranks 2`. Old-library rejection
uses `--mode deny-old-k2` against an old build only. Earlier COMM_SELF tests remain
historical evidence; these numerical checks do not certify physical convergence.


## Symmetry-file input checking before initialization

Python checks `sym.dat` before loading solver code, acquiring a lifetime or
creating a work copy. The canonical format has three records per operation and
four finite real values per record. D exponents, commas, `2*0` repeats and inline
`!` comments are supported. No rewriting or missing-value substitution occurs.
The Fortran disabled-axis filter is applied before checking retained signed
permutations, identity, duplicates modulo lattice translations and group closure.

A malformed file on one rank produces collective `InputError` for the entire
solver group. No solver initialization has been acquired: repair the original input,
construct a new material provider and retry. Editing an already prepared
provider's inputs is still rejected by existing hash checks.
The geometry input checking below checks cells, atomic/species mapping and grids.
Generated k-stars, PP references/plain KY and supported legacy GS structural
checks are also available. Physical GS identity, all non-KY contents and general
solver STOP recovery remain outside this scope.

Add `--admission-fault truncated`, `nonfinite` or `nonclosed` to the preceding
save command to corrupt only the last rank of each group and verify rejection,
repair and subsequent start. Normal input/numerical formulas remain unchanged.

Axis flags must contain one to three lowercase `y`/`n` characters; unsupported forms are rejected rather than silently normalized.


## Cell, atom and grid input checking

For periodic fast-PR input with or without symmetry, `NativeFastHSEBackend.create_state`
checks effective geometry before solver loading/lifetime/work-copy creation.
No ASE/f90nml dependency is added. Every solver rank must enter the collective
call. Repair rejected source input, construct a new provider and retry.

- Cells follow solver priorities: al, diagonal al_vec1..3, num_rgrid*dl, and detect off-axis vectors. Finite positive lengths/volume and a consistent inverse are required. The corrected solver core also admits non-axis-aligned, permuted and left-handed nonsingular cells.
- Grid counts/spacing and solver integer bounds are checked. Resolution follows solver int((length+1e-10)/spacing), including its spacing-to-cell tolerance; no ceiling, repair or input rewriting.
- One explicit inline or bundle-local Cartesian/reduced atom source is required. Counts, finite coordinates and species IDs1..nelem are checked. The existing fast-PR adapter remains au-only.
- With symmetry, check Cartesian cell isometries, integer grid maps (corner residual/index bounds, no full-grid allocation) and unique species-preserving atom permutations.

The focused numeric reader supports real/integer literals, D exponents, repeats
and positive index/section assignments. Null or unsupported forms reject rather
than being guessed. Atom records use canonical label/x/y/z/species fields. File
name case is preserved; references escaping the bundle reject. Inputs stay intact.

Use `--admission-fault cell|cell-metric|cell-singular|atoms|grid|grid-translation` in the explicit
probe to corrupt one rank of each group and check collective rejection without
lifetime acquisition, repair and normal start. Current/energy traces and direct
resume remain unchanged. Supported k/PP/legacy-GS structural checks do not
authenticate physical GS identity, every PP format, all solver initialization
paths or physical convergence.

Old binaries hit ACE initialization STOP on axis-permuted fixtures. The corrected
complete source fixes this; rebuild explicitly rather than modifying a fixed binary.
Primitive axes and orbital records are not reordered.

## General cells and x/y/circular polarization

The external field is a Cartesian three-vector A/c. Attached-calling program, zero-impulse
`tddft_response` initializes at A=0 and checks the actual vector at every step.
Ordinary batch initialization keeps its previous checks. No fixed z polarization,
field projection, tiny-component removal or automatic symmetry change is imposed.

| Field | Example retained Cartesian symmetry |
|---|---|
| x polarization | C2x |
| y polarization | C2y |
| z polarization | C2z |
| Circular xy | mirror-z preserving the entire xy plane |

The group must also preserve the actual cell and atomic species/positions. A field
incompatible with that reduced model is refused before updating, preserving the
accepted state. Use a matching full-k initial state when no suitable group exists;
arbitrary GS data are not automatically converted and groups are not changed mid-run.

General cells use solver column-vector lengths for grid resolution and the full
cell matrix for phases, k maps and Cartesian rotations. General-cell symmetry
requires `sp_fast_hse_general_cell_protocol()==1`; retained groups not preserving
z additionally require `sp_fast_hse_symmetry_field_protocol()==1`. Missing protocols
reject before solver initialization acquisition. Old full-k and ordinary z-preserving
diagonal-cell paths do not require unnecessary new protocols.

Rotated orthogonal, left-handed, xy-sheared and z-permuted cells are compared using
full8/reduced4 k points, internal MPI1/2, two independent material groups and own
direct7→32 resume. The same checks cover x/y/circular-xy fields. Synthetic orbital
regressions do not certify real-material convergence.

```sh
python tools/test_fast_hse_general_cell_mpi.py \
  NEW_BUILD FULL_K_BUNDLE REDUCED_K_BUNDLE NEW_OUTPUT --old-build OLD_BUILD
# Optional --roundoff-reference-build ADMISSION_REFERENCE compares nonbinary
# diagonal lengths against old calculation routines with only their admission rejection fixed.
```

Resume still requires the same binary/model/MPI partition; checkpoint encoding
and arrays are unchanged. Arbitrary space groups and
broader real-space partitions remain separate work.

## Symmetry learned histories and PT-CN

With `allow_symmetry=True`, use existing `hse_taylor4` off/strict/learned factors
or `hse_ptcn` off/learned factors. New combinations require an explicitly built
library with `sp_fast_hse_symmetry_history_protocol()==1`. Ordinary factor-off
Taylor does not require this additional capability.

Taylor uses its existing `exx_factor_exact_interval` teacher cadence. Learned
PT-CN requires **`exx_factor_exact_interval=1` and `yn_predictor_corrector='n'`**,
as in the existing solver method. Parameters/methods are not automatically changed.
Choose them before creating a state.

Histories belong to rank-local representatives. Existing numeric payloads retain
coefficients, teachers, ready/prediction flags, accepted counters and PT-CN dt;
direct resume replays no accepted steps. Tests verify actual prediction use and
identical history diagnostics after restore. This adds neither a portable archive
format nor unrestricted multiple TDDFT contexts per process.

```sh
python tools/test_fast_hse_symmetry_history_mpi.py \
  NEW_HISTORY_BUILD FULL_K_BUNDLE REDUCED_K_BUNDLE NEW_OUTPUT
```

This developer test explicitly derives synthetic test inputs from existing
factor-off fixtures. It covers internal MPI1/2, full8/reduced4 k points, split3/7/13
before/after training, two independent screening groups, direct resume and
state-preserving incompatible-field rejection/retry. It does not convert arbitrary
real-material GS data. Learned Taylor/PT-CN also retain all saved state/history
after post-training numerical overflow; corrected retries match an independent
continuous trajectory bit-for-bit.

For the ordinary FDTD mixed sample, supply matching reduced/full-k bundles with
the same propagator and learning conditions. No mid-run method switching occurs.

Optional phase-cache presence is restored as optimization state; a freshly initialized owner need not already have that allocation. Saved cache ranges are checked against trusted calculation routine MPI phase ownership before allocation. Static map, shape, finiteness and truncation validation remains active.
