# Molecular orbital/spatial MPI verification — 2026-11-26

## Tested layouts

The existing `molecular_fdtd_mpi_response` sample already derives the orbital
count from the molecular rank budget divided by the spatial product. No new
solver operation, numerical routine or parallel-layout adapter was necessary.
Both cases were executed using the previously installed wheel outside the checkout,
with one OpenMP/BLAS thread per rank:

| Molecular decomposition (k,orbital,x,y,z) | Molecular ranks | Field ranks | Total |
|---|---:|---:|---:|
| 1,2,1,1,2 | 4 | 2 | 6 |
| 1,2,2,1,2 | 8 | 2 | 10 |

Commands use respectively `--molecule-ranks 4 --molecule-grid 1 1 2` and
`--molecule-ranks 8 --molecule-grid 2 1 2`, with `--field-ranks 2 --steps 96 --ld`.
The outer caller launches MPI once. Each molecular group owns its orbital/space
communication, without a nested launcher or molecular child process. Both groups
return with caller MPI alive.

## Reference and binary identity

Conditions are unchanged from the spatial campaign: one Apple Silicon Mac,
spinless fixed-ion PZ C2H2, existing GS reused, dt=0.003 atomic time, 40^3 FDTD grid,
24^3 spatial E/H patch and off-grid molecular position. No new GS was computed.
One-way coupling has no molecular current return or self-vacuum subtraction.

The independent one-rank control is the previously completed
`salmonpy-molecular-spatial-20261126/independent-one` run. That control uses the
same input, binary, field source and two-rank field layout. It was reused rather
than silently presented as another new reference calculation.

The unchanged calculation-library SHA256 was verified:
`e0a37cdf269bba13dc9ad424750f38a910ba97dfbb2cb33a9eb82d7fb5aa7441`.
Source inventory SHA256:
`4b0fe7d4317f88b6e1cc8a012922e1948a822d50ab4ca8abc4dcb89f72edea3f`.
No new executable was built; existing source/object audit evidence remains tied
to that unchanged library. Documentation updates are not part of the older wheel
used for these numerical runs.

## Independent-run agreement

All response arrays over96steps pass atol=2e-12, rtol=2e-10. Total energy from
the ordinary stream is checked separately with atol=2e-10 eV, rtol=2e-10.
Maximum absolute differences:

| Orbital x spatial | Dipole (e bohr) | Electron number | Total energy (eV) |
|---|---:|---:|---:|
| 2 x 2 | 1.24e-13 | 3.66e-13 | 8.02e-12 |
| 2 x 4 | 1.15e-13 | 3.49e-13 | 5.01e-12 |

Electronic-change differences are <=1.32e-13 e bohr. Field E/H, clocks and step
numbers are exact. This is software agreement, not material convergence or a
scaling benchmark.

## Direct restart

Each layout ran continuously to molecular96/field97. An independent job saved
field/LD and all molecular rank payloads at molecular48/field49, exited, and a
fresh job restored the same decomposition and continued to96/97. Every response
array is bit-identical to its corresponding continuous run. These four ordinary
files are byte-identical for both layouts:

- `field/obs0_info.data`
- `field/obs1_at_point_rt.data`
- `molecule-work/C2H2_rt.data`
- `molecule-work/C2H2_rt_energy.data`

Actual counters are `{step:48, updates:0}` after direct restoration and
`{step:96, updates:48}` at completion. No elapsed propagation or SCF recomputation is performed for restoration. Fresh
workers first initialize from the existing GS, then directly load the saved
numerical state; this setup is not a new GS calculation. Caller MPI survives detach.

Results and `comparison.json` are stored under canonical calc2026 in
`salmonpy-molecular-mixed-20261126/{oz-continuous,oz-first,oz-resumed,oxz-continuous,oxz-first,oxz-resumed}`.
No cloud synchronization is asserted. Expanded existing parameterized layout and
manual-search tests;54targeted source tests pass. Bilingual guides and AI JSON
are synchronized. No heavy remote physics workflow was dispatched.

## Remaining scope

Verified only the listed CPU PZ layouts on one Mac. Larger orbital counts, other
spatial layouts, multiple nodes/site MPI, performance/I/O scaling, SOC, changed
restart decompositions, multiple molecular states per process and general composed
failure recovery require separate work. Molecular feedback/radiation remains
outside the current FDTD -> molecular TDDFT scope.
