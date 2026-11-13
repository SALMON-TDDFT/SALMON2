# SALMONpy non-LD entry points

Approved direction: task-level APIs plus explicit low-level state, immutable
model/evolution settings, A-only crystalline material response, legacy equations.

First usable increment:
- Isolated grid, axis-aligned box/ellipsoid geometry and immutable boundary
  settings. Preview at SALMON material-grid coordinates. PML occupies eight cells
  in the current solver; do not expose an unsupported adjustable thickness.
- SALMON Acos2 vector-potential pulse sampling and composition, with A meaning
  SALMON Ac=A/c, E=-dAc/dt. No silent conversion between crystal A and LD E.
- Read numbered SALMON time-series headers, retain labels/units, transform with
  explicit Fourier sign/window/time normalization. No claim that raw FFT is HHG.
- FDTD input generation and execution through an explicit validated SALMON binary,
  new run directory, no shell, fixed immutable calculation, full invocation record.
  Start with isolated impulse-driven dielectric/LD objects and point observations;
  do not introduce a replacement FDTD solver or pretend there is a step API.
- Define A-driven backend/initial-state ownership for TDDFT/SBE in documentation.
  Native reusable TDDFT/SBE contexts require separate state extraction; no fake
  backend or synthetic J result is acceptable.

Tests before code: shape coordinate semantics, pulse derivatives/composition,
Fortran-exponent data reading, Fourier phase/scaling, input ownership, FDTD input
and run failure handling. Integration: compare actual FDTD wrapper outputs with
same explicit namelist executed directly, including emitted shape grid. Package
pure-Python helpers in existing wheels and exercise across the established CI.
