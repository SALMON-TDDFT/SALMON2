# Multiple independent material points — 2026-11-17

This campaign extends the existing configured mixed-material path to repeated
material IDs and larger field-local LD volumes. No solver, public API or numerical
update order changed. The reusable layout is
[`mixed_multipoint_layers.json`](../salmonpy/examples/mixed_multipoint_layers.json).

## Configuration

Local Apple Silicon host:64GB/18 cores; sequential jobs, OMP/BLAS1,16 MPI ranks,
48GiB advisory budget. No remote Actions or resource hard caps.

- 64³ Maxwell, field grid2×2×1:four field ranks.
- SBE:four independent points, one rank per point.
- Libxc-PBE TDDFT:two independent points, two orbital ranks per point.
- Fast-PR HSE learned PT-CN:two independent points, two k ranks per point.
- Drude and Lorentz:1,152 cells each, retained by the owning field ranks; no
  dedicated LD ranks. Total constitutive cells2,304; electronic point count8.

Each electronic grid cell owns a different evolving state and logical point ID,
while sharing the declared material definition and initial-state assets. Per-point
solver-group allocation remains explicit. This is the current one-point-per-group
architecture, not arbitrary multiple native contexts in one process. The geometry
has finite transverse widths, not an infinite uniform film. Electronic grid/k
precision is unchanged; three-dimensional TDDFT here is a connection fixture.
Practical TDDFT multiscale use remains1D; molecular work remains paused.

## Results

All runs completed64 accepted steps except the deliberate save-only segment:

| Run | Restored step | New native steps | Maximum loop seconds | Sum of rank peak RSS GiB |
|---|---:|---:|---:|---:|
| continuous | 0 | 64 | 3.533 | 2.350 |
| protected | 0 | 64 | 3.994 | 3.046 |
| save | 0 | 31 | — | 3.045 |
| resume | 31 | 33 | 2.086 | 3.049 |

- Uninterrupted and protected A/E/j/independent serial-reference arrays match
  exactly, and every ordinary `.data` stream has the same bytes. Every material
  point also passes the existing independent serial-response comparison (SBE
  exact, electronic current rtol2e-10/atol2e-12).
- Points at different axial positions have different drives and currents despite
  sharing material definitions: maximum A/c/current differences for selected
  same-material pairs are SBE1.7311e-5/5.4103e-7, PBE1.7496e-5/4.8702e-7 and
  HSE1.3300e-5/1.0367e-7 in the existing atomic-unit conventions. This is consistent
  with distinct states, not broadcast copies of one material response.
- Direct save at31 and a fresh MPI job restoring that state perform only33 new
  native steps, with all64 response-array frames exactly equal to uninterrupted
  execution. The mixed point/material/LD identity record is unchanged and all
  eight logical point IDs are distinct. No accepted history is propagated again.
- Warm step13 failure is injected in only the last HSE point group after native
  success. The remaining electronic groups and field have their own completed
  attempts; the global transaction restores every original numerical payload
  hash, Python state/reference/clock and ordinary stream bytes. A real native
  field-overflow rejection is checked separately. Explicit retries and all final
  arrays/ordinary output streams equal uninterrupted execution exactly.
- Final source and fresh installed-wheel suites:545 tests, one skip each.
  Source-distribution→wheel/native SHA256 checks and guide/JSON synchronization
  pass. The installed wheel contains the new JSON and supports the same16-rank
  direct restart using its packaged layout.

Loop time includes MPI waits and protection; setup/output/serial reference work
is excluded. These are single local runs, not a timing comparison with statistical
confidence. Memory is a sum of rank peaks, not simultaneous total RSS or a heap
bound. Instrumentation and protected checkpoint restore consume extra storage and
memory. Resource-measurement failure is separate from scientific test failure.

Direct restart preserves full response-array history. New RT files record the new
segment; this test does not claim merging or restoring historical ordinary RT
streams. Same inputs/binaries/MPI topology remain required. Unhandled I/O/process/
MPI/STOP failures and unrestricted native reentrancy remain outside the contract.

## Reproduction and evidence

Use the [paired native stepping guides](../salmonpy/NATIVE_STEPPING.en.md), select
`examples/mixed_multipoint_layers.json`,64³,field2×2×1,16 total ranks,TDDFT2/HSE2
ranks per point and64 steps. Keep the separate PBE/HSE initial-state bundles.
Save with `--direct-save-at 31`; in another fresh output directory use
`--direct-resume PREVIOUS_OUTPUT/direct-checkpoint` and retain `--steps 64`.

Persistent local evidence:
`/Users/otobetoshihito/.codex/outputs/salmonpy-multipoint-20261117/`.
Campaign/commands, input layout, per-run outputs, direct checkpoint, fault harness,
file inventories, source/installed logs and built distributions are retained.
Builds use the unchanged production receipts from the native protected-step
campaign. Local paths record provenance, not user installation requirements.
