# Molecular two-level MPI verification — 2026-11-26

## Execution and ownership

One caller-launched MPI application partitions four ranks into two FDTD ranks and
two molecular orbital ranks. `run_molecular_mpi` attaches only the supplied molecular
communicator. It neither initializes/finalizes MPI nor launches another process.
`MolecularMPIController` sends commands root-to-root; molecular commands and full
spatial E/H frames are broadcast within the molecular group. The shared frame
packet is filesystem-backed and read by the molecular root once per transfer.
No application-wide collective is added to each propagation step. Final caller
MPI communication succeeds after molecular detach.

The molecular channel configuration belongs to the blocking call frame. Its
module pointer is bound only during that dedicated-rank call. This is one lifetime
per process, not a general reentrant/multiple-state molecular library interface.

## Build identity

- macOS Apple Silicon; four MPI processes; one OpenMP/BLAS thread each.
- Source inventory SHA256: `4b0fe7d4317f88b6e1cc8a012922e1948a822d50ab4ca8abc4dcb89f72edea3f`.
- Embedded library SHA256: `e0a37cdf269bba13dc9ad424750f38a910ba97dfbb2cb33a9eb82d7fb5aa7441`.
- Molecular MPI protocol: 2. CPU MPI, Libxc enabled, HSE disabled for this test.
- Tracked-file-only working-tree export, source preparation, complete build and
  application-object audits passed. The export contains these uncommitted edits;
  it is not a claim about the previous committed HEAD. Mandatory eigensolver
  prohibition checks remained enabled throughout the build and link.

## Actual propagation

Existing C2H2 GS data were reused; no new GS was computed. The test uses PZ,
fixed ions, spatial E/H and an off-grid molecular position. A 40^3 field grid
supplies a 24^3 patch, dt=0.003 atomic time. One-way coupling has neither molecular
current feedback nor a self-vacuum reference.

- The separate-process one-rank regression matches the previous executable
  bit-for-bit in all nine observables over four steps.
- The new one-rank host route matches that reference bit-for-bit. Two orbital
  ranks differ by at most 7.11e-15 over the initial four-step comparison.
- Live two-field/two-molecular-rank propagation runs 96 molecular steps. Field E/H
  match the earlier one-rank molecular control exactly; molecular observables
  differ by at most 2.05e-14 (comparison atol=2e-12, rtol=2e-10).
- Source execution saved at molecular step48/field step49, terminated, and resumed
  in a new four-rank application. All complete history arrays and four ordinary
  output files match continuous execution bit-for-bit. Restoration adds zero
  updates; continuation adds exactly48.

## Installed distribution

Built an sdist and then a wheel from the extracted sdist. Installed that wheel in
an independent environment and ran outside the checkout. All53 targeted runtime
source, near-field, isolated-TDDFT and manual tests passed. The wheel uses its own
lightweight kernels; the molecular embedded library is supplied explicitly.

The installed sample created a fresh joint checkpoint with checksum manifest at
step48/49, then a separate MPI launch restored it and reached step96/97. Every
history array and these ordinary streams are byte-identical to the source
continuous control:

- `field/obs0_info.data`
- `field/obs1_at_point_rt.data`
- `molecule-work/C2H2_rt.data`
- `molecule-work/C2H2_rt_energy.data`

Installed metadata reports restore `{step:48, updates:0}`, final
`{step:96, updates:48}`, and `caller_MPI_alive:true`. Raw outputs and comparison
records are stored under the canonical calc2026 directory in
`salmonpy-molecular-mpi-20261126/{v2-continuous,v2-first,v2-resumed,installed-first,installed-resumed}`.
No confirmed cloud synchronization is asserted. Earlier intermediate builds and
failed harness attempts remain separate; they are not final runtime evidence.

## Scope and remaining work

Verified: spinless fixed-ion CPU PZ, one/two orbital ranks for the molecule, two
field ranks on one Mac, matching binary/input and unchanged restart decomposition.
No physical-convergence claim is made. Molecular spatial decomposition, multiple
nodes, large rank counts, multiple states per process and density/magnetization
exports through this channel remain unverified. Unexpected errors abort the
example; arbitrary MPI-failure recovery and composed atomic rollback are not
claimed. Molecular feedback/radiation remain deferred by the current one-way scope.
