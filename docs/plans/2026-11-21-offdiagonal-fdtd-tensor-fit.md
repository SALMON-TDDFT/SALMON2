# Off-diagonal Yee constitutive coupling and TDDFT tensor fitting

User requested general reciprocal off-diagonal FDTD, rotated waveplate validation, unchanged isotropic/diagonal results and a TDDFT impulse tensor-to-LD path.

The legacy isotropic/diagonal path remains exact. For each full-tensor material, let D be its diagonal edge coefficient with the existing edge ownership mask, S the adjacent-edge mean, and C the normalized material correlation matrix (unit diagonal). Apply G = sqrt(D)[I + S^T(C-I)S]sqrt(D). Since S is a contraction and C is PSD, G is PSD; for C=I it is exactly the previous diagonal operator. Coupling is within one material's edge subspace, so its pole resonance/damping commute with G. The implicit E operator is background G plus2*pi*dt times pole-C2 operators; use diagonal-preconditioned CG with global reductions and paired mean/transpose-mean halos. P/J/Q use the same projected pole drive. No component is silently dropped. Finite unconverged solves warn; invalid arithmetic preserves native state via the existing stage20 guard. Expose tolerance/budget on FDTD.

Checkpoint: all dynamic histories remain existing E/H/P/J/Q. Immutable tensor context is owner-held and compared/bound in full-tensor checkpoint headers, while scalar/diagonal binary paths stay unchanged. Validate source/object eigensolver audits, Git-only source inventories and installed wheels; no real-material physics in CI.

Analysis: three independent impulse runs supply simultaneous charge-current and total-E histories. Native matrix inversion per frequency yields sigma=J E^-1 and epsilon=I+4*pi*i*sigma/omega, with missing/rank-deficient bins masked, not invented zeros. The reciprocal passive tensor LD fit uses shared poles and PSD factorized strengths, native residual/normalization and optional SciPy optimization. Preserve antisymmetric/noisy input residuals and warn instead of silently symmetrizing.

Acceptance: rotated principal-axis synthetic slab incl.Fresnel multiple reflections; mesh refinement; old scalar and diagonal bitwise controls; MPI1/2 and direct field checkpoint replay-free continuation. Synthetic three-drive tensor recovery and tensor-fit roundtrip; existing TDDFT output adapters provide data but do not claim a new GS/TDDFT material validation.


Implemented this iteration: native three-drive matrix optical conversion and native
PSD-strength tensor-fit objective; public dielectric_tensor and
fit_lorentz_drude_tensor with fixed initial background and optional SciPy;
synthetic mixed-drive inversion and passive negative-off-diagonal fit tests;
ld_fit --tensor sample; paired LD guides and AI catalog entry.

Completed follow-up (2026-11-21): native mean/transpose-mean reciprocal tensor operator,
material-owned EPS/F coefficients and pole drive, diagonal-preconditioned implicit E,
PEC fixed-edge projection, warning-only finite iteration exhaustion and guarded
numerical retry; full-tensor direct checkpoint binds raw matrices, rates, edge
ownership, pole weights, grid spacing/dt and solve options before commit. Scalar
and diagonal checkpoint streams remain unchanged. Added rotated lab-frame native
waveplate diagnostic and --angle-degrees example. The existing native recovery
tool now accepts --tensor rather than adding a separate test family.

Validation: Git-only739-file build, original source/object prohibition audits,
source and installed-wheel658 tests (one optional skip each), wheel/runtime source
byte parity. Rotated quarter/half-wave Stokes errors0.003049/0.005626; half grid
quarter error0.000760 (4.01-fold reduction). Old/new scalar/diagonal bitwise equal;
rotated two-pole dispersion E/H differences below3e-15. 3D cyclic-axis covariance
with all three off-diagonal pairs below2e-21. MPI2 direct restore, one-rank-only
material ownership, late-file restore failure, model mismatch and numerical retry
all preserve the accepted state and match continuous A/E/H and raw observations
bitwise; no replay. Finite maxiter1 run warns and continues. Installed three-drive
synthetic impulse -> epsilon tensor -> passive fit -> native full FDTD completes;
fit weighted RMS9.28e-6. No new real-material TDDFT/GS validation, GPU/Fujitsu or
multi-node claim. No remote CI triggered. Final record and retained audited local
runtime in canonical calc2026/salmonpy-offdiag-final-record-20261121.

Remaining limitations: antisymmetric/nonreciprocal magneto-optics, nonzero native
Bloch phase for this operator, generic material rotation convenience API, field
checkpoint MPI repartition and separate A-only maxwell.py anisotropic backgrounds
are not added by this work. These are explicit scope boundaries, not missing
components silently projected to zero.
