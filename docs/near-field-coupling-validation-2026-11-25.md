# Molecular field representation and coupling validation (2026-11-25)

The optional directional factor exp[-decay_vector dot (r-center)] multiplies the existing
radial-complete real solid-harmonic basis. Decay is in inverse bohr and fixed with center,
radius and degree. Zero decay retains the ordinary polynomial representation. FieldV3,
direct numerical checkpointV6, near-field kernel3 and executable-bound receipt3 identify
this contract; no old-format fallback is added.

Fortran separates field representation from molecular coupling. The representation supplies
integrals of E and sH along center+s*r. A common map builds electron scalar energy and A/c.
The polynomial/exponential provider uses stable path moments. The Python function
`couple_path_samples` dispatches caller-supplied field samples and quadrature to the same
Fortran map; it does not start TDDFT or interpolate an FDTD grid.

## Independent spatial and integral checks

Fortran reference quadrature versus exponential moments (degrees0..9, arguments from-30
to100 including zero and small arguments): maximum relative difference4.21e-15.
Magnetic projector line integration versus independent quadrature and closed-form uniform-field
moments: maximum absolute phase difference8.47e-22. NaN/Inf and rank-deficient inputs
are rejected; failed sampled-path coupling retains output buffers.

Production fitting uses729 Cartesian sample points and257 separate off-grid points in an8bohr
cube, with known electrostatic Laplace field exp[-8(z+.5)]*[sin(8x),0,cos(8x)] in normalized
coordinates. At total degree8, relative RMS errors are:

| Representation | E | Electron scalar potential |
|---|---:|---:|
| Polynomial |6.159%|7.336%|
| Exponential factor, known decay |0.4003%|0.1687%|

The scalar potential is checked against the analytic potential difference from the expansion
center, not only against a second fit. This is a local analytic-field test, not a Tip/Bow Tie
or time-dependent Maxwell validation. No physical decay-length inference is claimed.

## Molecular history and direct restart

Existing C2H2 GS, PZ, fixed ions, unpolarized CPU1rank/1thread, dt=.003au,20updates.
Center(.2,-.3,.1)bohr and decay(.04,-.03,.12)bohr^-1; quadratic radial terms and distinct
E/H time histories. Prepared and live coefficients, independent controls and direct restart
at update7 after process termination match through update20 for all9observables. Restoring
requires0extra evolution calls. Saved output-stream continuation checks pass. A per-step
decay mismatch is rejected before updating state. Electron loss fraction3.3751e-15.

Ordinary uniform-field and zero-decay polynomial runs are bit-exact against the preceding
validated results; Fortran maximum absolute differences are0. No new SCF was performed.

Git-tracked-only source preparation/build and application-object eigensolver audit pass.
Source tests50pass; independently installed wheel77pass/2skip. The installed calculation
library also passes the Fortran spatial/integral test. CI was not dispatched.

Details are in the configured calc2026 results store under
`salmonpy-c2h2-nearfield-20261125/exponential-coupling-validation.json`,
`exponential-coupling-analytic-wheel.log`, `exponential-final-wheel-molecule-validation/`
and `exponential-final-wheel-polynomial-regression/`. Local OneDrive storage is not
confirmed cloud synchronization. Actual Tip/Bow Tie FDTD fields, molecular-observable
convergence, direct grid interpolation and bidirectional feedback remain separate work.

The electric path integral follows the multipolar treatment discussed by
[Takenaka, Taketsugu and Iwasa, JCP154,024104 (2021)](https://doi.org/10.1063/5.0034933).
That paper omits magnetic interactions; the orbital magnetic implementation and exponential
compression tested here are not claimed as results of that paper.
