# Transverse periodic plane flux (2026-11-21)

Explicit local E cross H area integration in Fortran, before time integration
or spectral flux ratios. reflection_transmission area_weights is opt-in and
the old uniform path is unchanged. Plane fluence is area-averaged for compatibility;
energy_balance checks and multiplies the declared area once. surface_flux reports
signed total power/cumulative energy and supports additive owner-local patches
without introducing MPI. Ratios require the complete plane; no rank-local ratio
averaging. Acquisition/collocation/directional separation remain explicit.

Validation: uniform limit, unequal area weights, nonnegative/zero masks, scalar
IEEE rejection including zero-weight nonfinite samples; phased propagating order
with zero mean fields but full normal power; partial-plane power/energy additivity.
Real scalar lossless dielectric grating in native E/H FDTD:160x32x4, spacing1,4,2,
dt.0025,2400 steps, periodic y/z, PML x, bar20x32x8, epsilon2.25, sourceomega10/tw2.
Matching vacuum reference,32 y samples spanning full128x8 cell with z invariance.
Raw fields collocated using adjacent x planes and periodic y. Final stored
energy from electromagnetic-only checkpoint, H recentered using the exact TE
Yee update outside PML; endpoint half volumes. Work in dispersive poles is0,
and instantaneous dielectric storage is counted in field energy.

Incident .0043216978325 Hartree; reflected .0002425803222; transmitted
.0040625968990; final storage .0000169812593. Signed residual
-4.6064796894e-7, relative -1.0658958280e-4 (~0.0107%). Crossing averaged
fields gives transmitted .0036258785857, missing about10.1% of incidence.
Weak/finite-time spectral tails can give negative bins: warnings retained,
not clipped. This is not general3D optical convergence or multi-node validation.

No solver update or production binary changed; SHA verification reused the
validated executable. New observation kernel and C reference are source/object
audited. Actual physics remains manual, outside required CI. Canonical records
are calc2026/salmonpy-periodic-surface-energy-20261121; no cloud-sync assertion.

Final checks: source and installed wheel676 tests,1 optional skip; sdist→wheel
passed; stored native grating NPZ example exactly matches source and wheel.
Legacy uniform-port report and all output arrays exactly match preceding wheel.
Logs/results/summary retained in the canonical record directory.
