# Molecular spatial MPI validation — 2026-11-26

## What changed

The one-way `molecular_fdtd_mpi_response` sample now accepts
`--molecule-grid NX NY NZ`. For molecular rank budget N, spatial product S must
divide N; the orbital division is N/S. Both the solver group and controller receive
`decomposition=(1,N/S,NX,NY,NZ)`, which is recorded in output metadata and bound to
the direct checkpoint identity. Default `1 1 1` preserves orbital decomposition.
Invalid products are rejected before MPI initialization. This option changes the
parallel layout, not the physical grid. No numerical Fortran/C routine changed.

## Environment and input

Apple Silicon Mac, one node, one OpenMP/BLAS thread per rank, caller-owned MPI.
FDTD uses two ranks; the molecule uses one, two or four dedicated ranks. Existing
C2H2 GS data are reused, with PZ, fixed ions, no spin, dt=0.003 atomic time and
96 molecular steps. FDTD uses a 40^3 grid and supplies spatial E/H on a 24^3 patch
at an off-grid molecular position. This is the same LD/illumination setup as the
orbital-MPI validation. No new GS calculation was performed. Molecular feedback
and self-vacuum subtraction are absent.

The unchanged library SHA256 was verified before reuse:
`e0a37cdf269bba13dc9ad424750f38a910ba97dfbb2cb33a9eb82d7fb5aa7441`.
Its audited source inventory is
`4b0fe7d4317f88b6e1cc8a012922e1948a822d50ab4ca8abc4dcb89f72edea3f`.
No new solver build or source/object audit was required for this unchanged binary.

## Independent execution comparison

Each decomposition ran in its own MPI application with the same existing GS,
field source and library. The reference has one molecular rank; all cases retain
the same two-rank field group. All response arrays pass atol=2e-12, rtol=2e-10.
Maximum absolute differences over the complete 96-step history:

| Spatial divisions | Molecular ranks | Dipole (e bohr) | Electron number | Total energy (eV) |
|---|---:|---:|---:|---:|
| 2,1,1 | 2 | 1.10e-13 | 3.61e-13 | 8.99e-12 |
| 1,2,1 | 2 | 1.06e-13 | 3.61e-13 | 8.99e-12 |
| 1,1,2 | 2 | 1.27e-13 | 3.66e-13 | 8.02e-12 |
| 2,1,2 | 4 | 9.50e-14 | 3.49e-13 | 3.98e-12 |

Electronic-change differences are at most1.27e-13 e bohr. Field E/H, clocks and
step numbers are exact. Total energy is checked separately from the ordinary
molecular energy stream with atol=2e-10 eV and rtol=2e-10. Small differences reflect
changed floating-point reduction order; no physical-convergence claim is made.

## Direct joint restart

For spatial layouts (1,1,2) and (2,1,2), separate jobs save field/LD and molecular
state at molecular48/field49, terminate, then restore in fresh MPI applications
with the same layout and capacity96. Native rank-specific checkpoint payloads
are loaded directly. No previous step, wavefunction propagation or GS is recomputed.

Both restored histories match their respective continuous controls bit-for-bit in
every response array. All four ordinary streams are byte-identical:

- `field/obs0_info.data`
- `field/obs1_at_point_rt.data`
- `molecule-work/C2H2_rt.data`
- `molecule-work/C2H2_rt_energy.data`

Counters immediately after restoration are `{step:48, updates:0}`; after
continuation `{step:96, updates:48}`. Caller MPI remains alive after molecular
detach. Raw results and `comparison.json` are under canonical calc2026 in
`salmonpy-molecular-spatial-20261126`. No cloud synchronization is asserted.

## Scope

This closes spatial propagation and same-layout direct restart validation for the
listed small CPU PZ cases. Changed-layout restart, combined orbital/spatial
partitioning, multiple nodes, large rank counts, SOC, general distributed failure
recovery and multiple molecular states per process remain separate work. The
existing single-rank process worker and one-way coupling semantics are unchanged.

## Distribution checks

Source targeted tests:54passed. Built an sdist, extracted it and built the wheel
from that extracted source, then installed into a separate Python environment.
All54targeted tests passed outside the repository using the installed package.
The installed six-rank sample restored the source-generated four-spatial-rank
checkpoint in a new process launch. Every array and the same four ordinary
streams match the source continuous run exactly; actual counters are again
48/0 immediately after restoration and96/48 after continuation. Caller MPI
survives. This also exercises the installed `--molecule-grid` option. Human
Japanese/English guides, AI JSON and their synchronized copies were checked.
No remote CI or physical-convergence campaign was launched.
