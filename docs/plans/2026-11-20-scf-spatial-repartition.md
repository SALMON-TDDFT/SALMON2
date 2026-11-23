# Ordinary-XC spatial SCF restart repartition

Keep the explicit same-runtime/input restart opt-in. Add native spatial protocol1
and interior grid_bounds metadata; allow pure r-grid layouts with scalar periodic
PZ/Libxc-PZ/PBE, not k/orbital mixed layouts or HSE. Fortran reads old rank payloads
and maps physical grid cells, all six wavefunction/CG buffers, potentials, density,
Fourier density and all density/potential mixing histories. Rebuild periodic halo
values by copying accepted interior cells, never by SCF or eigenstate evolution.

New topology already has fresh static grid/FFT/communication resources from the
existing zero-SCF bootstrap. Retain those topology-specific static tables and dead
FFT/XC scratch instead of copying incompatible partial-transform caches. Accepted
XC potentials, density, mixing age and histories remain saved values. Publish only
after global success and complete/nonoverlapping interior coverage.

Verify real/complex1↔2 spatial groups, exact accepted interior cut values, zero SCF
bootstrap and continued source-initialized reference tolerance1e-9/1e-10. Compare
physical interiors, not stale old halo buffers. Inject late donor read failure,
retain original state and source data, then continue. Validate existing k/orbital
and HSE paths, source/installed tests, revised existing sample, paired guides,
AI JSON and C header. No new Python module or routine CI.
