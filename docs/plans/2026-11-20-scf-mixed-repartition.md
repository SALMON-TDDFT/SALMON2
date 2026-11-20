# Ordinary-XC mixed MPI SCF restart migration

Keep run(restart=...,repartition=True) explicit. Add optional mixed-repartition1;
ordinary PZ/Libxc-PZ/PBE only. Source and target may use any compatible k/orbital/
r-grid factorization; physical inputs, assets, material ID and runtime stay fixed.
Existing pure-layout entry points remain available and HSE mixed migration is
outside this step.

Combine saved grid_bounds/orbital_bounds/k_bounds into each rank's five-dimensional
interior ownership box. Validate non-overlap/coverage and actual array bounds.
Fortran maps all six wavefunction/CG arrays by physical grid, orbital and k index,
including periodic halo copies. Replicated density/potential/Fourier-density and
mixing histories are read once per old spatial domain (old orbital/k owner starting
at1), avoiding duplicate accumulation. Static resources and dead FFT/XC scratch
belong to the new topology. Publish only after complete wave/scalar coverage and
a global success vote; no SCF/eigenstate/density recomputation.

Extend existing tests/tool and sample, not module/test counts. Exercise2↔4 and
4↔8 layouts spanning k/orbital/spatial dimensions, exact accepted interiors,
zero bootstrap and source-initialized continuous reference at1e-9/1e-10, late
missing donor preservation/retry, existing pure/HSE controls, installed wheel,
JP/EN/AI/header and checklist. No routine CI.

Acceptance: source 4→2,2→8,8→4,two-atom4→8,real-Gamma4→2,4→4 and16→8
passed. Exact saved physical interiors, zero bootstrap SCF, continuation within
unchanged1e-9/1e-10. Independent/previous-runtime parity and pure-HSE spatial
regression passed. Installed-wheel8→4 and sample4→2 passed. Source/installed637
tests each636pass/1skip; code/docs/native source parity checked. Evidence:
calc2026/salmonpy-scf-mixed-record-20261120. No routine CI. Mixed HSE remains open.
