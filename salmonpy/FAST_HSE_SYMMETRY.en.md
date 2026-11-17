# Fast-PR HSE and symmetry

[日本語](FAST_HSE_SYMMETRY.md) · Topic: `fast_hse_symmetry` · Revision: `2026-11-18.5`

Fast-PR HSE is the primary HSE development and validation path. Symmetry-aware
representative k points, star reconstruction, orbit reuse and ACE are central to
that path. Conventional HSE is not a supported replacement engine. Standard
SALMON symmetry reduction has different arithmetic and is not used to redefine
the fast-PR algorithm.

| Layer | Current scope |
|---|---|
| Fast-PR exchange core | Independent symmetry maps/phase caches/FFT owners; tested scalar periodic CPU, full q, one internal k rank |
| Symmetry/orbit/ACE comparison | Reduced4 k representatives vs independently enumerated full8 k points; C2z and mirror-x groups, different cells/omega/sources |
| Internal Fortran TDDFT | Explicit opt-in, one initialization/process/internal-k1, factor-off Taylor; direct/staged resume tested |
| Python/native-C/FDTD | Explicit new-protocol allow_symmetry opt-in; group1 factor-off Taylor, protected mixed tests; default remains denied |
| Checkpoint | Internal exchange V2 and explicit internal TDDFT direct/staged resume tested; bounded public protected mixed resume tested |
| Symmetry internal MPI k2 | Not enabled by this change; two test MPI processes each own COMM_SELF cores |

## Independent exchange owners

The internal `exx_bind_system` core can now admit symmetry within the bounded
single-k-rank profile. Each system retains its own k-point configuration and each
HSE owner retains its own symmetry map, phases and FFT resources. The full-grid
check uses that system's ngrid, primitive cell and spacing, without selecting from
shared num_rgrid. Dimension/count mismatch or an unsupported owner layout rejects
before refresh/action mutation when using the declared rt_status entry. Other
invalid symmetry operations/atoms/weights may still terminate native execution.
The admission result is not a general numerical-failure recovery guarantee.
The new core scope requires exx_factor_mode='off' and Taylor/ACE or the test-only
full-action variant; learned/PT-CN symmetry still rejects before mutation.

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
The native probe is packaged under
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
or native initialization. Default worker/native-C opening and coupled-checkpoint prohibitions
remain unless the explicit public opt-in below is selected; this is not unrestricted support. The explicit
native probe is `_runtime_source/tests/native/test_hse_symmetry_field.f90`.

## Internal symmetry checkpoint

The existing internal `hse_checkpoint_write/read/capture/install` now support a
value-only symmetry map/cache for the admitted scalar CPU/full-q/internal-k1
exchange core. Non-reduced owners still write V1 byte-for-byte. Symmetric owners
write same-binary V2: retained map, kernel identity, source/ACE/history values and
optional complete phase cache. The caller authenticates the archive and model;
this raw stream is not a standalone portable interchange format.

Initialize an independent matching owner/operator before reading (operator
construction, no accepted time-step replay). The reader compares all static map
values and trusted shapes before reading the bounded optional cache. It stages
before installation; validate-only and staged reads leave the live state intact.
A truncated stream, differing map/phase partition or malformed cache rejects.
Installation moves value arrays only, leaving live FFT/C/MPI owners untouched.

The native probe `_runtime_source/tests/native/test_hse_symmetry_checkpoint.f90`
checks full/ACE, star/orbits, cache on/off, OMP1/2 and MPI1/2 processes using
independent COMM_SELF cores. Restored action/energy and subsequent refresh/action
checkpoint bytes match continuous core execution, including reopening an owner.
This is exchange-core persistence; the separate internal TDDFT tests are described below.
Ordinary/default Python/native-C init stays symmetry-denied; the explicit opt-in is described below. Learned/PT-CN symmetry,
reduced internal-k MPI and ordinary symmetric mixed failure/resume remain open.

## Explicit internal TDDFT validation

An internal Fortran `rt_open(..., allow_symmetry=.true.)` opt-in now connects the
exchange core to the complete TDDFT checkpoint path. It is bounded to scalar
periodic CPU/full-q, one internal k rank, factor-off Taylor/ACE and a retained
group preserving z. Full-action is test-only. One initialization per process is
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
The ordinary init/serial worker remain denied; the separate public opt-in enables the bounded profile.

The native `test_hse_symmetry_tddft.f90` probe and companion
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
The opt-in is included in backend identity, admits one internal k rank, and never
changes the supplied operation matrices. Missing/unknown protocols reject before
acquiring a native lifetime. Learned/PT-CN symmetry and reduced internal-k2 stay
unsupported; full-action remains test-only.

Only unanimous finite status34 plus the specific protocol becomes a Python
`InputError` preserving the accepted state. Correct A before retrying. A reduced
model must match the actual local FDTD field, not only the incident polarization.
Do not project A, drop tiny components or switch engines automatically. For a whole
loop, use `CoupledStepper.step_protected`/`--protect-steps`, and retry only when
`StepRecoveryError.original_preserved` is true and cleanup succeeded.

The ordinary `mpi_mixed_native` sample accepts `--hse-symmetry` and requires
`--hse-reference-bundle FULL_K_BUNDLE` with an explicitly supplied equivalent full-k
state for the existing serial fast-PR comparison. No arbitrary GS conversion is
performed. `--hse-ranks 1` is required. The layout/LD/SBE/semilocal-TDDFT APIs and
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

Next: reduced internal-k MPI and broader initialization. Fast-PR-only policy,
warning-only memory budgets and paused molecular work remain unchanged.
