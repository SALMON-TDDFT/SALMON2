# Fixed 4³ k-point Si analysis

User instruction: use 4x4x4 k points; do not perform a k-point convergence scan.

Local run directories (generated files are deliberately not versioned):

- `gs/`: PZ ground state, 8 Si atoms, 32 electrons, 12³ real-space grid, 4³ k points.
- `alda/`: impulse response, no macroscopic xc field.
- `lrc/`: identical GS and impulse, macroscopic LRC alpha=0.2.

Each directory contains its exact inputfile and outputfile. RT settings are dt=0.08 a.u.,
12000 steps (about 23.22 fs), impulse=0.001 a.u., z polarization, transverse geometry.
Both RT runs reuse `gs/data_for_restart`. Executable source: commit 84c18671,
MPI gfortran/OpenBLAS build described in docs/tdcdft-validation.md.

GS reached the requested 1e-9 density criterion at iteration 34 (8.52e-10).
Spectra and analysis are saved under `docs/results/si-k4/` after the runs complete.
