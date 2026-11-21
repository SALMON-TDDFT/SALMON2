# Reference-frequency subpixel smoothing

User authorized a dedicated Fortran/C module connected to FDTD using the input central wavelength. This is an explicit frequency-local dispersive approximation, not Meep's epsilon-infinity-only treatment and not a whole-band second-order guarantee.

Dedicated native module: planar box clipping and analytic primitive level-set normals/fractions; complex local interface correctors at omega_ref; real symmetric energy-congruence lift of each SPD background/PSD LD residue. For two phases use continuity of E_t and D_n to compute T_m, then epsilon_eff_inf=sum f_m Re(T_m^H epsilon_inf_m T_m), F_eff_q=f_owner Re(T_owner^H F_q T_owner). Preserve pole resonance/damping, hence passive causal rational response; it matches the complex planar homogenized epsilon at omega_ref. Different interface cell normals/fractions require cell-owned tensor coefficients, not a palette with a material-ID cap.

Connect initialization-owned per-node/edge coefficients to the existing paired-mean tensor operator. Do not change Maxwell curls, implicit solve, E->J->P/Q->H order or timestep formulas. The coefficient lookup must support spatially varying tensors and globally identified oscillator-rate classes; report allocation overhead. Legacy mode default off is byte-identical. Avoid adding a new per-step solver or frequency-dependent update.

Geometry: boxes/ellipsoids from existing inputs. Integrate a locally planar cut of each Yee dual voxel, rather than sample a fixed subgrid. Corners/multiple interfaces/unsupported material constitutive terms/reference singularities warn and fall back cellwise to the original material, not an invented value. Infer omega_ref only from an explicit carrier; impulse requires an explicit positive reference. Metadata records reference, scope, fallbacks and coefficients.

Tests: independent tangential arithmetic/normal harmonic mean; endpoints/no contrast; absorbing and metal/vacuum complex matching, reciprocity/passivity; geometric plane clips and smooth surfaces; legacy off parity, homogeneous enabled parity, planar/curved grid refinements and reference-band diagnostics; MPI and direct checkpoint roundtrip. Scalar IEEE only. Source/object prohibition audits, inventory updates/Git-only build, installed wheel checks. No GS and no remote physical CI.


## Completed local verification — 2026-11-21
- Dedicated Fortran module and native E/H coefficient preparation connected;
  matching V2 profile required. Shared source copies byte-identical. All740
  current solver inventory assets match the audited Git-only prepared/build tree.
- Opt-in LD interfaces only, normal carrier inference and explicit impulse
  reference. SBE/TDDFT/user response states are never transformed; native response
  neighborhoods excluded. No-LD/unsupported backgrounds warn and retain baseline.
  Keep existing invalid excitation diagnostics rather than asking a disabled
  smoother for a reference.
- Scalar complex interface matching, passivity, fractions/endpoints and shared
  adjacent-layer boundary tested. Synthetic axial-E infinite cylinder reference:
  radius5bohr, wp20, damping4, reference14au, T10au, fixed dt0.0005au,
  h2/1/0.5bohr relative complex center-E errors0.02035998565,
  0.00465516998,0.00109107080. This tangential-field case shows near-quadratic
  refinement, not a broadband/general-geometry theorem. C reference uses outgoing
  cylindrical Bessel functions and native Fourier quadrature.
- Homogeneous two-pole LD enabled/disabled results identical. Prior off-diagonal
  runtime versus new disabled path bitwise equal. Old/new unsmoothed curved case
  also compared. Rotated30deg quarter waveplate Stokes error0.003049, synthetic.
- MPI2 two LD layers/three global poles/one excluded response point: continuous
  and direct restored400-step A/E/H bitwise equal. Missing/NaN/wrong-ID assets
  rejected1113 collectively; late checkpoint/model failures preserve state;
  overflowing finite trial rejected, retry succeeds; caller MPI remains usable.
  Layer boundary fallback count0, local match error2.22e-16.
- Source and installed-wheel suites663 tests, one optional skip each. Wheel MPI
  controls and direct restore also exercised. Updated paired FDTD guides, AI JSON
  and existing waveplate example with --subpixel; preview starts no solver or dir.
- No GS/electronic physical run, remote CI, GPU/Fujitsu or multi-node claim.
  Evidence/runtime: canonical calc2026/salmonpy-subpixel-final-record-20261121.

Remaining boundaries: equal-rate pole coalescing/memory optimization, complex
  shapes beyond boxes/ellipsoids, multi-interface corner accuracy, field MPI
  repartition and general broadband correction are not implemented by this option.
