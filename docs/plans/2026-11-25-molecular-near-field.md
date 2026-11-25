# Molecular near-field exposure and FDTD expansion (2026-11-25)

1. Expose existing SALMON prescribed dipole sources through an explicit isolated-worker mode E3N and public native_near_field=True. Keep E3 unchanged. Preserve SALMON input spatial profile, source locations/vectors/radius and pulse envelopes; caller E remains an additional uniform field. Zero source counts select the unchanged uniform path. CPU/PZ/fixed-ion/single-rank scope remains.
2. Reject invalid native dipole counts (0..2), nonfinite positions/vectors/radius, nonpositive radius and unsupported/impulse source envelopes before near-field initialization/time stepping. Bind mode and source input to snapshots and disk manifests. Runtime evidence must explicitly validate the near-field route; no old-binary capability inference.
3. Synthetic kernel comparison with unchanged set_vonf/calcVbox and tests for all scalar IEEE values. Reuse already saved C2H2 GS if possible without collecting new GS; compare source-worker with native batch, zero-source baseline, independent histories and direct restart. Preserve source/object eigensolver audits before all new builds/links.
4. Add FDTD-to-molecule expansion after first exposing native sources: define expansion origin, real/complex Ylm normalization, radial dependence and units, accepted step numbers and field component locations. Scalar electric potential expansion is appropriate only for the selected quasistatic/longitudinal component; general transverse fields need consistent vector coupling. Recover regular potential coefficients from fields without claiming unresolved spatial data. Fortran/C evaluates profiles and expansion coefficients; Python owns dispatch and provenance. Initial one-way coupling; later adjoint-consistent current/multipole deposition for two-way coupling and work balance.
5. Update Japanese/English guides, AI JSON, examples and the current checklist. Do not add real-material jobs to CI; use local explicit evidence. Do not touch the quantum design document owned by another chat.


## 2026-11-25 steering: E/H is the general near-field interface

The user pointed out that near-field mode should receive E,H directly. Adopt spatial E(r,t),H(r,t) as the general contract.
Ylm is a representation/transfer method, not a conversion of all fields into a scalar potential. The implemented ElectricPotentialSeries is explicitly a quasistatic electric-only projection. Actual propagating-field probes show large residuals and are not physical-accuracy evidence. Retain all failed/high-residual trials.
Next work must align E/H positions and step indices, preserve original field times, map H to the appropriate B for the molecular region, and implement consistent electric and magnetic multipolar interactions. Do not silently drop H or declare general near-field support from the current scalar route.

## E/H multipolar implementation (authorized 2026-11-25)

1. Add NearFieldSeries in near_field.py. Fit the three Cartesian components of E and H independently in a real regular-solid-harmonic basis including l=0. Fix center during fitting (explicit or Fortran sample centroid) and radius; carry them with coefficients, units and sample step indices. No curl-free projection of E. Preserve direct construction from coefficients.
2. Extend the common Fortran solid_harmonics kernel with vector fitting/evaluation and prescribed-field radial-gauge potentials: V_e=r.dot(integral E(s*r) ds), A/c=-r cross integral s*H(s*r)ds/c in molecular vacuum. This includes E1/E2 and orbital magnetic coupling; spin/SOC is outside the current scalar-PZ molecular endpoint. Keep A^2 and the symmetrized kinetic operator. Use ion-centered projector phases with five-point line integration, exact for the supported polynomial order.
3. Add explicit near_field=True worker mode with prepared E/H input. Own coefficient history, geometry, vector-potential grid and projector phases in s_rt, serialize them and restore by copying, without electronic evolution. Preserve the ordinary and quasistatic-only routes.
4. Validate independent algebra, vector-field fitting, magnetic kinetic action, uniform E equivalence, independent molecular states and direct process-exit restart with nonzero H. Reuse the permitted C2H2 GS in the existing campaign. No new SCF/GS collection. Preserve eigensolver audits and scalar IEEE checks; no heavy CI.
5. Update Japanese/English instructions, preparation/execution samples, public API and AI JSON. State magnetic prescription, vacuum units, step indexing and the remaining two-way radiation/physical-convergence work explicitly.


## Implemented and checked

NearFieldSeries/NearFieldFrame, fixed fit geometry, from_yee_samples, prepared and live step_near_field input, orbital magnetic/diamagnetic kinetic coupling and ion-centered nonlocal phases are implemented. Direct state includes coefficients, magnetic arrays and projector phases. Molecular process-exit restart matches continuous execution without additional evolution. The preceding plan statements about missing E/H entry were superseded by this implementation. Spin/SOC, radiation feedback, full Maxwell diagnostics and physical convergence remain outside the tested scope. See the latest checklist and Japanese/English ISOLATED_TDDFT guides.

## Current scope: saved-field input, arbitrary molecular positions

The user separated complete FDTD feedback from this task. Molecular centers may be off-grid;
convert existing spatial E/H to TDDFT input without snapping molecular positions. The live coefficient
API remains a prescribed-field input, not evidence of a complete two-way FDTD coupling.

Pending Tip / Bow Tie study: use rounded tips and rounded triangular bow ties, vary apex radius,
gap, molecule distance/offset, fitting-region radius and polarization. Scan lmax=0,1,2,3,4,6,8 with
adequate independent sample counts (at least (lmax+1)^2 and full rank). Keep separate fit and held-out
points, report absolute RMS/max and relative RMS for E and H, then molecular-response convergence.
First establish FDTD grid convergence. Do not infer required local expansion order from the antenna's
far-field scattering multipole order. The existing 27-point transient test does not establish this order.
The harmonic basis lacks general radial Helmholtz terms: inspect residual plateaus before extending L.
No large FDTD or new GS campaign has been launched for this study.

## Higher-order saved-field scan

Collected 729 vacuum-side points of the existing 40^3 LD-sphere/FDTD transient,
83 frames, dt=.003 au. Fit 486 points and evaluate 243 held-out points; fixed off-grid
center (5.17,4.77,5.31) bohr and radius15 bohr. All fitting/evaluation/RMS arithmetic
uses the checksum-verified Fortran kernels. Tool: tools/scan_near_field_orders.py.
Held-out whole-trace relative RMS E/H for field lmax=0,1,2,3,4,6,8:
97.21/96.96%, 86.24/85.03%, 70.59/68.43%, 68.22/66.08%, 62.08/59.53%,
49.82/47.27%, 42.43/40.23%. At lmax8 absolute RMS is E8.04266e-4 E:au,
H6.89575e-4 H:au. No order reaches5%. This is a broad-region transient diagnostic,
not a Tip/Bow Tie or molecular-response convergence result. Field lmax2 already
includes electric octupole potential terms; field lmax3 is not the definition of
electric octupole in this component-wise field convention. Next isolate fit-volume
and missing radial dependence before inferring a required physical multipole order.
Record: calc2026/salmonpy-c2h2-nearfield-20261125/eh-order-scan.json.

## Molecular-volume / radial-dependence diagnostic

Implemented standalone Fortran complete-polynomial least squares (SVD rank check, scalar IEEE,
application object prohibition audit before linking), used only as a diagnostic; solver unchanged.
A separate r^2 oracle is recovered to maximum absolute1.11e-15, whereas lmax2 harmonic fitting
has46.8% held-out RMS error. This independently exposes the excluded radial term.

Existing729-point LD/FDTD data, off-grid center(5.17,4.77,5.31)bohr, local cube sides
2.12/4.23/6.35/8.47 angstrom; deterministic spatial split75% train/25% held out.
Compare field degrees1,2,3,4,6,8 at steps40,71,82, not an entire-trajectory claim.
At step71, side4.23angstrom, degree2 E/H harmonic errors23.00/21.78%, complete polynomial
18.88/19.54%; degree3:30.62/27.84% versus7.43/6.30%; degree4:20.09/20.56%
versus0.000460/0.000444%. Degree4 absolute RMS E1.7485e-8/H1.5112e-8au for complete
polynomials, versus7.6358e-4/7.0001e-4au for harmonics. Side2.12angstrom degree2:
harmonic19.64/15.73%, complete0.000584/0.000581%. Degree8 cannot be assessed in the
smallest volumes because of sample/rank limitations. Larger-volume degree8 remains poor
for harmonics (side8.47angstrom step71 E50.58%/H48.23%).

Interpretation: shrinking the volume helps but cannot supply missing radial curvature. Increasing
harmonic order alone is not a reliable remedy. These extremely small polynomial residuals are
on an existing Cartesian grid with only3/5/7/9 coordinate levels; polynomial interpolation of
that discrete data is NOT evidence of continuum accuracy, sub-grid resolution, or Tip/Bow Tie
or molecular observable convergence. No new GS or TDDFT run. Diagnostic does not change the
production expansion/Hamiltonian. Next requires independent denser-grid reference or analytic
Maxwell fields before connecting a radial-complete representation to molecular dynamics.
Results: calc2026/salmonpy-c2h2-nearfield-20261125/eh-radial-volume-scan-recorded.json.

## Radial-complete TDDFT input implementation

Approved by the current user request: replace the near-field E/H basis with r^(l+2k)Y_lm,
total degree d=l+2k<=L<=8. Order d ascending, l descending in steps2, m ascending.
P=(L+1)(L+2)(L+3)/6, max165. No compatibility branch: fieldV2, direct diskV5,
near kernel protocol2 and executable-bound near-EH receipt2 identify the new contract.
Fortran fitting, radial-gauge potential integrals (degree+1/degree+2 denominators),
projector line phases and all saved coefficient histories use the same representation.
Five-point Gaussian line integration remains exact up to the maximum degree9 vector potential.
Off-sample analytic r^2/r^4/r^6/r^8 checks and quadratic E/H potentials were added.
The previous harmonic diagnostic records remain historical; pending actual molecular restart validation.

### Radial-complete validation completed locally

Fortran-only400 fit /53 independent off-grid points: degree2/4/6/8 E/H maximum absolute
errors3.4521e-16/4.5797e-16/7.2858e-16/1.2967e-15au. Electric potential maximum1.1701e-15,
A/c maximum6.8441e-19. Existing GS C2H2,20 updates, nonzero r^2 E/H, center(.2,-.3,.1)bohr:
prepared/live exact, independent controls exact, process-exit direct checkpoint at7 and updates8..20
all nine observables exact, restoration0 extra evolutions; electron loss fraction1.1191e-14.
Normal uniform path still matches batch within3.86e-12e*bohr. No new SCF or physical-convergence claim.
Records in calc2026/salmonpy-c2h2-nearfield-20261125/radial-molecule-validation and
radial-complete-basis-crosscheck.json; earlier angular-only diagnostic records are retained.

## Exponential-factor comparison (Fortran diagnostic)

Compared exp(-kappa*z) times the radial-complete basis against the same linear coefficient count
without the factor; known z direction, normalized cube. Fortran uses729 regular-grid fitting points,
129 independent off-grid tuning points,257 separate final points, SVD with rank checking.
Degrees0/1/2/3/4/6/8; decay ratios.5/1.7/2/7.3/8 plus quadratic control. Rate selection is
per degree on tuning points over0..24 in steps.25; no final-set selection. Known-rate comparisons
have the same3P linear coefficients; selected-rate comparisons add1 shared parameter.

336 combinations passed. Uniform preservation, scalar NaN/Inf rejection, rank-deficient rejection,
unchanged predictions on failure and exact known-rate exponential recovery are checked in Fortran.
For exp(-8z)[sin(8x),0,cos(8x)] Laplace field, final relative RMS for degrees2/4/8:
polynomial94.73/62.88/7.183%; supplied-rate factor81.30/25.13/.3143%; selected-rate
factor81.09/25.16/.3101%. Degree8 selected kappa8.5 differs from true8: selected rate is
not uniquely a physical parameter. For two decay lengths at base8, degree8 error1.676%
falls to.01085% with selected factor (kappa12.75). Quadratic control selects0 and remains
roundoff exact. Off-candidate true rates1.7 and7.3 are included, not just rates on the search grid.

Conclusion: promising optional representation for steep directional decay; not a universal
radial exp(-r) rule, not a replacement default. This diagnostic changes no TDDFT Hamiltonian,
field file, checkpoint or runtime. Dynamic H/Maxwell consistency, tip/bow-tie fields and molecular
observables are untested here. Next exponential connection needs weighted radial-gauge integrals,
projector phases, parameter ownership and direct-restart validation.
Records: calc2026/salmonpy-c2h2-nearfield-20261125/exponential-basis-expanded.csv and
exponential-basis-comparison.json. No new GS or CI campaign.

## Representation-independent coupling and optional exponential factor

Current approved step: keep polynomial default and add explicit decay_vector in inverse bohr,
field(r)=exp[-decay_vector.dot(r-center)] times the existing radial-complete basis. Direction/rate
are supplied after preparation; no automatic physical decay-length inference. Separate a common
Fortran coupling module (electric/magnetic path integrals) from basis fitting/evaluation.
Implement stable exponential path moments, retain exact polynomial arithmetic for zero decay,
and use a checked high-order path quadrature for magnetic nonlocal-projector phases.

Fixed center/radius/order/decay belong to the model and checkpoint, not per-step overrides.
FieldV3/directV6/kernel3/near-EH receipt3 identify the contract, no older-format fallback.
Test analytic off-grid exponential fields, zero-decay identity, generic sampled-field coupling,
line-phase accuracy, invalid scalars and output preservation; reuse current permitted C2H2 GS for
prescribed/live equality, independent histories and process-exit restart with zero extra evolution.
No new GS, automatic downloads, heavy CI or speculative full two-way feedback.

### Representation/coupling separation implemented

New shared near_field_coupling.f90 computes representation-independent E and sH path integrals,
and the scalar/vector potential map. Public couple_path_samples accepts field samples and caller
quadrature, with numerical integration in Fortran. The solid-harmonic provider computes exact
exponential path moments without small-argument cancellation; magnetic line phases use16-point
composite quadrature. Zero decay preserves polynomial arithmetic. Fixed inverse-bohr decay_vector
is owned with geometry in frames/checkpoints. FieldV3/directV6/kernel3/near-EH receipt3.

Independent Fortran32-point reference: maximum moment relative error4.21e-15 over arguments
0,1e-12,1e-6,.5,2,8,20,21,100,-1,-8,-30 and degrees0..9; maximum magnetic line-phase
absolute error8.47e-22. Invalid path samples preserve outputs; NaN/Inf and rank checks pass.
Production MGS fitting on729 grid points,257 new off-grid points in an8bohr cube with known
source-free electrostatic exp[-8(z+.5)]*[sin(8x),0,cos(8x)] field: degree8 polynomial E/V
relative RMS6.159%/7.336%; exponential-factor E/V .4003%/.1687%. This is an analytic spatial
accuracy test, not Tip/Bow Tie FDTD or time-dependent Maxwell/molecular physical convergence.
C2H2 prototype20 updates with fixed(.04,-.03,.12)bohr^-1 envelope: prescribed/live and
process-exit direct restart7->20 all9observables bit-exact, extra evolutions0, electron loss0.

Updated Japanese/English guide, AI structured contract, public API and preparation example.
Remaining: independent actual tip/bow-tie fields and molecular convergence, generic spatial
interpolation provider, full E/H Maxwell diagnostics and separate two-way radiation design.

### Final tracked-source and installed-wheel verification

Final source inventory32e8175ce571b895abac80db4baa844b83bdfd9f2ff12cb828b62b56a3c35d64;
tracked-only source preparation/build and application-object prohibition audit pass. Source tests50pass.
Installed wheel77pass/2skip. Its Fortran analytic test repeats the independent strong-decay E/V
and moment/phase checks. Final C2H2 drive varies both E and H coefficients over time; center and
(.04,-.03,.12)bohr^-1 decay remain fixed. Prepared/live histories, independent controls, process-exit
restart7->20 all9observables and restored output streams match exactly;0extra evolutions on restore.
Electron loss3.3751e-15. Mismatched decay is rejected before a step. Zero-decay regression and the
ordinary path are bit-exact against the previous validated results (Fortran maximum difference0).
No new SCF, CI dispatch or cloud-sync claim. Records: calc2026/salmonpy-c2h2-nearfield-20261125/
exponential-coupling-validation.json, exponential-coupling-analytic-wheel.log,
exponential-final-wheel-molecule-validation/, exponential-final-wheel-polynomial-regression/.

## Direct Cartesian-grid E/H input (current approved step)

Implement a third field provider, with no Ylm fit: regular3D grid, trilinear interpolation,
explicit per-component E/H origins and common positive spacing. GridFieldSeries owns field
histories and metadata. molecule_position is the FDTD position of the TDDFT integration center;
center is that same point in molecular coordinates. Fortran performs the translation, packing,
Yee half-step averaging, interpolation and path integration. No extrapolation or implicit wrapping.

Split center-to-electron rays at all component-grid faces;3-point Gauss integration is exact
for the resulting piecewise-trilinear E and sH integrands. Use converged line integration for
magnetic nonlocal phases. Save grid data/geometry with the electronic state; direct binaryV7,
field fileGRID1. Prepared histories use step(E_extra); live grid-streaming and two-way feedback
remain separate. Existing fitted provider and arithmetic stay unchanged.

Verify off-grid staggered linear fields, smooth-field grid refinement, cell crossings, bounded
failure/no output changes; use existing permitted C2H2 GS and actual saved LD/FDTD fields for
continuous versus process-exit direct restart. No new SCF or heavy CI. Update guides/JSON/example.

Direct-grid core implemented and locally checked: component-origin-aware trilinear interpolation,
exact3-point integration after splitting every crossed cell face, Fortran packing/half-step average,
source-to-molecular translation and prepared-history TDDFT provider. Grid metadata/data are owned
and included in the bounded stream codec; isolated directV7 and embedded directV4 identify the
expanded codec. Existing fitted field formats remain unchanged. Native analytic affine E/H and
potentials have maximum absolute error2.08e-17, including staggered origins and cell crossings.
Smooth analytic field refinement h=.5/.25/.125bohr gives E relative RMS6.602/1.455/.3285%,
V3.418/.6464/.09014%; no physical Tip/Bow Tie convergence claim. Invalid sources and queries
preserve outputs. Actual LD/FDTD64^3 sampling completed with4624 sites and82 E/H frames;
physical field time origin.303au is retained, source library reused unchanged by hash.

### Direct-grid molecular check completed

Native grid input at FDTD position(21.2,.7,.9)bohr, TDDFT center0, coarse patch origin
(7,-15,-15)bohr, spacing2bohr, counts16x17x17. Existing C2H2 PZ GS,20 updates dt=.003au.
Uniform direct-grid agrees with batch within3.86e-12e*bohr. Actual saved LD/FDTD history:
raw/public direct process-exit restart7->20 all9observables bit-exact, no extra evolution,
output streams exact and independent state unchanged; electron loss0. Field time origin.303au
is retained in the input comment/model identity; it is provenance, not a molecular-history replay.
Source records: calc2026/salmonpy-c2h2-nearfield-20261125/direct-grid-molecule-validation/.
No new SCF, physical-convergence claim, automatic FDTD launch or complete two-way connection.
