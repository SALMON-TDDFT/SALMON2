# Fast-PR HSE and symmetry

[日本語](FAST_HSE_SYMMETRY.md) · Topic: `fast_hse_symmetry` · Revision: `2026-11-18.2`

Fast-PR HSE is the primary HSE development and validation path. Symmetry-aware
representative k points, star reconstruction, orbit reuse and ACE are central to
that path. Conventional HSE is not a supported replacement engine. Standard
SALMON symmetry reduction has different arithmetic and is not used to redefine
the fast-PR algorithm.

| Layer | Current scope |
|---|---|
| Fast-PR exchange core | Independent symmetry maps/phase caches/FFT owners; tested scalar periodic CPU, full q, one internal k rank |
| Symmetry/orbit/ACE comparison | Reduced4 k representatives vs independently enumerated full8 k points; C2z and mirror-x groups, different cells/omega/sources |
| Python TDDFT/FDTD connection | Symmetry remains rejected; existing non-reduced fast-PR connection is documented in COUPLED |
| Checkpoint | Internal exchange-core V2 supports map/cache; symmetric TDDFT/direct-protected coupled resume still prohibited |
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
or native initialization. Symmetric TDDFT opening and coupled-checkpoint prohibitions
remain; symmetric time evolution and restart are not certified. The explicit
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
This is exchange-core persistence, not symmetric TDDFT propagation or FDTD resume.
The Python/native TDDFT open prohibition stays in force. Learned/PT-CN symmetry,
reduced internal-k MPI and ordinary symmetric mixed failure/resume remain open.

## Next connection work

Before admitting reduced-k states into TDDFT/FDTD, implement and verify all of:

1. The owned actual-A check is implemented; verify it in symmetric time evolution once persistence is ready.
   Axis flags alone do not certify the final operation matrices. Arbitrary A(t)
   can break a reduction valid for a fixed polarization.
2. Internal core map/history/cache persistence is implemented; connect it through direct and protected
   save/restore without cloning resource handles or replaying accepted steps.
3. Validate representative/star partitioning for internal k MPI, plus failure,
   rollback and resume in the ordinary mixed sample.

The existing adapter/checkpoint prohibitions stay in place until those tests
pass. Do not remove them merely because a core comparison passes, silently
change the field or use conventional HSE. Memory advice remains warning-only;
physical convergence belongs to the user's chosen material/calculation.

Full-action is a test-only comparison inside the same fast-PR sources. It does
not add a supported alternative engine or change production propagation profiles.
