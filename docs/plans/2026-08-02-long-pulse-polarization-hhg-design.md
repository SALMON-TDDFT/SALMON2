# Long-Pulse Polarization HHG Design

## Scope

Validate the accepted overlapping-Wannier V3 coefficient-RT route with a
fixed-ion Si64 calculation whose cell is too small for quantitatively
converged harmonic intensities.  Polarization is the primary observable and
the HHG spectrum is computed only as `omega^4 |P(omega)|^2`; current remains a
secondary consistency diagnostic.

## Laser and propagation

Use the existing 1.55 eV `Acos2` pulse in length gauge, extended to ten
optical cycles.  Continue propagation for two optical cycles after the pulse
so the total record is about twelve cycles.  The generalized-eigenvalue
exponential propagator permits a production step of `dt=2.0 a.u.` while still
sampling one carrier cycle about 55 times.  Repeat the same physical interval
with `dt=1.0 a.u.` as the focused time-step check.

The coefficient driver continues to pass zero vector potential to the RT
Hamiltonian.  No simultaneous length- and velocity-gauge coupling is allowed.

## Acceptance policy

Retain implementation-level gates: field-off stability, finite observables,
V3 provenance, polarization/current sign and qualitative consistency,
amplitude linearity, exact buffered-fragment symmetry, Cartesian covariance,
and time-step convergence over the plotted harmonic band.

Treat absolute harmonic powers and strong/weak power ratios as diagnostics,
not pass/fail thresholds.  For the small Si64 cell, the physical acceptance is
qualitative spectral morphology: identify whether H2 and H4 are local peaks,
dips, or slopes and require a recognizable H3 local peak in ideal Si.  A fixed
displacement may reduce the exact group and may permit even-order structure;
it is not interpreted as lattice dynamics or dephasing.

## Figure

Generate a tracked plotting script and reproducible PNG/PDF artifacts from
the polarization-derived TSV spectra.  Use a logarithmic vertical axis, plot
the ideal x/y/z spectra together, mark integer harmonic orders, and show a
second panel comparing ideal and fixed-displaced x polarization when both
accepted matrices are available.  Record pulse length, `dt`, checkpoint and
binary hashes in the accompanying JSON evidence.

