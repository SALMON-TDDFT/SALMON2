# Molecular MPI layout-changing direct restart — 2026-11-26

## Implementation

The caller-owned molecular MPI interface now supports a different orbital/spatial
layout at direct restart. Fortran reads the saved rank payloads into a candidate,
checks ownership coverage/overlap, and copies wavefunction buffers, densities and
potential histories to the new owning ranges. Isolated-grid bounds use the actual
centered indices, including negative indices. Replicated whole-domain Vbox is
copied as a whole, not incorrectly treated as a local density grid.

The candidate receives the new topology's static Poisson/multipole/FFT work
resources, while retaining saved accepted Hartree/potential histories and the
last Poisson iteration count. Layout-specific magnetic vector potentials and
projector phases are prepared from saved E/H coefficients. This is data movement
and derived-resource preparation: no Hamiltonian propagation, SCF calculation or
elapsed time step is replayed. New workers still initialize from existing GS
files before direct loading. No MPI/file/Libxc handle is deserialized from an old
rank into a new rank.

The checkpoint wrapper is `salmonpy.molecular-mpi-checkpoint/2`; the molecular
library protocol is3. The identity binds initial input/library/dt/capacity,
separately from the recorded source decomposition. No old-wrapper/protocol
fallback was added. Same-layout restart retains its direct-load path. The joint
sample permits changed molecular ranks/grid while retaining the field layout.
No per-step communication was added. Repartition currently reads donor files
sequentially on each target rank; large-filesystem scalability remains unverified.

## Build and evidence

Final library SHA256:
`d2da24af566630a23efa09812d2ed40b855f71e2cfda438c9dc7ebc37957ecfa`.
Source inventory SHA256:
`47531303883faa8c36fefa76b70ba55269448730843d6b08554d6e52de4f910e`.
Tracked-file-only working-tree export, preparation and complete CPU MPI/Libxc
build passed. Mandatory source and all-application-object eigensolver prohibition
checks remained enabled before build/link. This export contains uncommitted edits,
not the previous committed HEAD. Scalar IEEE inquiries are retained.

Early failed attempts are preserved separately. The initial read failure exposed
an existing serializer assumption that max_steps was initialized; the donor now
receives the matching current RT configuration before reading. The next failure
exposed a periodic-grid lower-bound assumption; isolated bounds are now used.
A small audited Fortran reader localized these issues without propagation.

The sample also reports/aborts unexpected failures before collective context
cleanup. Previously its error could enter communicator cleanup while the molecular
root waited for another command. This change does not promise general library
rollback or recovery from process/MPI failure.

## Real calculation

Existing PZ C2H2 GS, fixed ions, spinless CPU, one Apple Silicon Mac, one
OpenMP/BLAS thread per rank, dt=0.003 atomic time. FDTD uses two ranks, a 40^3
field grid and spatial E/H on a 24^3 patch at an off-grid molecular position.
One-way coupling: no molecular current return or self-vacuum reference. No new
GS was computed. Save at molecular48/field49, exit the application, then restore
and continue to molecular96/field97.

| Source (orbital,x,y,z) | Target | Target molecular ranks | Max dipole difference (e bohr) | Max electron-number difference |
|---|---|---:|---:|---:|
| 1,1,1,1 | 2,1,1,2 | 4 | 1.24e-13 | 3.66e-13 |
| 2,1,1,2 | 1,1,1,1 | 1 | 1.24e-13 | 3.66e-13 |
| 2,1,1,2 | 1,2,1,2 | 4 | 1.27e-13 | 6.58e-14 |

Each resumed history was compared against an independent continuous execution
in its target layout, using the same final library. All response arrays pass
atol=2e-12, rtol=2e-10. Both ordinary molecular streams also agree within the
corresponding tolerances (energy uses atol=2e-10 eV). FDTD E/H, clocks, step
indices and its two ordinary output files are exact. MPI reduction-order
roundoff is expected; changed-layout molecular output is not byte-identical.

All restorations report `{step:48, updates:0}` and finish with
`{step:96, updates:48}`. Caller MPI survives detach. Same-layout mixed4 restart
matches every array and all four ordinary output files bit-for-bit. Normal
one-rank, mixed4 and spatial4 continuous histories match their pre-change
controls bit-for-bit.

Raw results, `comparison.json` and retained failure/build logs are under canonical
calc2026 in `salmonpy-molecular-repartition-20261126`. Earlier v1/v2/v3 directories
are intermediate failure evidence; `final-*` contains the passing v4 runs.
`v5-*` contains final-build verification after replacing a declaration-initialized
restore selector with an explicit call argument. Final histories match v4
bit-for-bit; no new per-process selector state is introduced. No
cloud synchronization or physical convergence is asserted.

## Scope and next work

Verified CPU PZ fixed-ion spinless C2H2 on one Mac for the listed transitions.
FDTD decomposition must remain unchanged in this joint example. Other XC/spin,
large layouts, multiple nodes, filesystem/I/O scaling, several molecular states
per process and general composed failure recovery remain separate. Projection
states, active tau-XC operator payloads and stored microscopic-current arrays are
not supported by this molecular repartition path; it rejects them explicitly.
Molecular feedback/radiation remains deferred under the current one-way scope.

## Installed-wheel verification

Built an sdist, extracted it and built the wheel from that source. Installed it
into an independent Python environment outside the checkout. All54targeted
runtime-source, near-field, isolated-TDDFT and AI-manual tests passed on source
and installed wheel. The installed six-rank sample loaded the source-generated
mixed4 checkpoint into spatial4 and continued without replay: restore48/0,
finish96/48. Every array and all four ordinary files are bit-identical to the
source execution of the same changed-layout restart, which agrees with target
continuous execution within the reported tolerances. Caller MPI remains alive.
Guide/JSON synchronization and diff whitespace checks passed. No heavy remote
physics CI was dispatched; distribution is locally built, not publicly released.
