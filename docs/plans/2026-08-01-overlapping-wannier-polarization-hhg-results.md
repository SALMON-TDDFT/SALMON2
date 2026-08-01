# Overlapping-Wannier Polarization and HHG Results

## Implementation and numerical verification

Coefficient RT now publishes occupation-weighted polarization and current in
`overlapping_wannier_rt_observables.dat`.  Spectra are calculated only from
field-off-subtracted polarization.  Current is retained as a secondary
consistency observable.  The synthetic analyzer contract verifies linear
response normalization and HHG as `omega^4*abs(P_omega)^2`.

The focused MPI fixture passed on 1, 2, 4, and 8 ranks.  The final clean-first
archive overlay was configured with MPI, ScaLAPACK, and EigenExa and produced
binary SHA-256
`28d88710fcd4e30885459197c02b52bd31856a210b4e24d4b321573c255c1cf9`.
The reused V3 manifest remained
`1a38778c27350ae07173ba891e94b8a9fa0b67f5ce223c83ca821b77e8f0269c`.

For `dt=0.2 au`, 2200-step field-off and x-impulse calculations gave:

- field-off polarization drift: `1.421e-12`;
- field-off maximum current: `2.185e-12`;
- impulse maximum `|Delta Px|`: `2.298e-3`;
- transverse/longitudinal polarization ratio: `4.04e-2`;
- centered `dP/dt` versus current relative RMS difference: `7.46e-3`;
- frequency resolution: `0.3884 eV`.

For 1.55 eV, approximately 9.7 fs x-polarized pulses at `dt=0.5 au`, the
maximum polarization changes were `5.258e-6` at `1e8 W/cm^2` and `5.239e-4`
at `1e12 W/cm^2`.  Their ratio, `99.64`, closely follows the electric-field
amplitude ratio of 100.  The strongest spectral bin was harmonic order 1.10,
within one frequency bin of the fundamental.  No resolved odd-harmonic
plateau was present; the calculation remained essentially linear at the
tested intensity.  Therefore it is not accepted as a demonstrated HHG
response.

## Critical physical-provenance finding

The artifact historically named the accepted "Si64" gate is not silicon.
Its tracked `atom.dat` labels all 64 atoms as `C`; the input uses `izatom=6`,
`C_rps.dat`, and 256 valence electrons.  Consequently the observed dominant
linear-response strength near 69--83 eV cannot be judged against silicon
optical data, and no Si linear-response or Si HHG acceptance claim is made.

A genuine Si64 V3 checkpoint, built with a silicon pseudopotential and
auditable Si input provenance, is a mandatory parent prerequisite for the
remaining Si cubic-axis, amplitude-linearity, time-step-halving, and HHG
acceptance matrix.  Reusing or relabeling the C64 checkpoint would invalidate
the requested physical verification.

## Genuine Si64 physical RED

A genuine 64-Si checkpoint was generated from the tracked silicon
pseudopotential in a 20.52-bohr cubic cell.  Conventional DC SCF converged at
iteration 85 and overlapping-Wannier SCF at iteration 79.  The final density,
coefficient, orthogonality, and minimum-metric-eigenvalue diagnostics were
`1.680e-8`, `2.904e-8`, `1.155e-14`, and `0.2773`, respectively.  Its V3
manifest SHA-256 is
`b001dad6da8edd17358b02ef0c1d21d9530d4ec7926222558b0b9f797e581124`.

The eleven-case matrix used `dt=0.5 au`, 800 steps, and matching half-time-step
runs with `dt=0.25 au`, 1600 steps.  The 1.55-eV laser cases used approximately
10-fs pulses at `1e8` and `1e12 W/cm^2`.  Numerically stable diagnostics pass:

- field-off polarization drift: `4.741e-13`;
- x-impulse amplitude linearity relative RMS: `1.132e-3`;
- current versus centered `dP/dt` relative RMS: `1.706e-2`;
- x-laser transverse/longitudinal ratio: `2.06e-2`; and
- time-step-halving HHG-spectrum relative RMS: `2.002e-4`.

The physical symmetry gates intentionally remain RED.  Cartesian impulse
amplitudes are `1.466e-3`, `1.543e-3`, and `1.619e-3`, a 9.4-percent spread
against the 5-percent cubic gate.  The `(H3+H5)/(H2+H4)` power ratios are
`0.830`, `0.762`, and `0.719`, below the required odd-harmonic dominance of
three.  The strong/weak third-harmonic power ratio is about `1.32e4`, below the
nonlinear gate of `2.0e4`.  These failures must not be hidden by loosening the
acceptance thresholds.

The retained construction currently verifies only fragment translations;
Hermiticity checks on `X/V` do not establish crystallographic vector
covariance.  Exact instantaneous crystallographic subgroups must therefore be
detected per buffered fragment and used to validate/project retained scalar
and polar-vector matrices.  For displaced zero-point or finite-temperature
snapshots the exact subgroup, including a possible C1 result, must be retained
instead of forcing ideal parent symmetry.
