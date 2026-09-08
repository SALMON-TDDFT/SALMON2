# Overlapping-Wannier Polarization and HHG Design

## Goal

Validate the accepted Si64 overlapping-Wannier coefficient-RT route through
linear response and laser-driven high-harmonic generation. Polarization is
the primary physical observable. Current is emitted only as a secondary
consistency diagnostic, and spectra are derived from polarization.

## Production boundary

The calculation must remain inside the accepted chain:

1. buffer-supported, symmetry-preserving overlapping-Wannier GS;
2. accepted `SALMON_OW_GS_CHECKPOINT_V3`; and
3. generalized-eigenvalue exponential coefficient RT.

It must not enter conventional RT, DG-Fragment, Nodal, WPW, adaptive-basis,
or normal checkpoint routes. Normal DC LCFO plus EigenExa remains unchanged.

## Observable definitions

For occupied coefficient columns collected in `C`, the density matrix is
`rho = C C^dagger`. With the V3 row-owned metric, position, and velocity
operators, the primary polarization and secondary current are evaluated as
collective expectation values:

```text
P_alpha(t) = -sum_n f_n Re[c_n^dagger X_alpha c_n] / volume
J_alpha(t) = -sum_n f_n Re[c_n^dagger V_alpha c_n] / volume
```

The occupation weights `f_n` come from the V3 checkpoint and must match the
coefficient columns exactly. The implementation must not silently assume a
fixed spin multiplier for arbitrary inputs.
The sign convention is the electronic convention and must be written in the
file header. Values at step zero and after every propagated step are emitted
by rank zero in one deterministic text file. Restart continuation must append
without duplicating the split boundary sample and must reproduce one-shot
output byte-for-byte.

The observable evaluation uses row-owned operator blocks and collective
reductions. It must not assemble persistent replicated observable matrices on
every rank. Hermiticity, finiteness, state dimension, provenance fingerprints,
and volume must be checked before publication.

## Spectral analysis

Linear response uses three independent small impulses. The field-off
polarization is subtracted before analysis. A documented damping/window is
applied to `Delta P(t)`, and the response is divided by the integrated impulse
field. Linearity is checked by repeating one direction at half amplitude.

HHG uses the polarization acceleration form without numerically differentiating
noisy time data:

```text
I_HHG(omega) proportional to omega^4 |FFT[P(t) - P_background(t)]|^2.
```

The analyzer records the time step, total duration, window, zero padding,
frequency resolution, Nyquist limit, laser carrier frequency, and harmonic
order. It also reports `J(t)` versus a centered finite difference of `P(t)` as
a secondary consistency measure; this diagnostic does not replace the
polarization spectrum.

## Calculation matrix

Start from the accepted Si64 V3 checkpoint and use eight MPI ranks and one
OpenMP thread.

- field-off reference;
- x/y/z impulses with a small integrated kick;
- one half-amplitude impulse for linearity;
- a 1.55 eV, approximately 10 fs pulse at a weak reference intensity;
- the same pulse near `1e12 W/cm^2` for HHG;
- a time-step-halved rerun of the nonlinear case or the shortest sufficient
  converged segment.

If this near-infrared matrix does not produce a resolved nonlinear spectrum,
the next stage is a separately reviewed 0.4 eV, approximately 30 fs mid-IR
matrix. It is not silently substituted into the first acceptance run.

## Physical acceptance

The result is physically acceptable only when all of the following hold:

- field-off polarization drift and current remain below recorded tolerances;
- S norm, charge, energy in field-off propagation, and checkpoint provenance
  remain within the existing coefficient-RT gates;
- impulse response scales linearly with impulse amplitude;
- cubic-equivalent directions agree within numerical tolerance;
- the nonlinear spectrum changes with intensity and is converged against the
  chosen time step and analysis window;
- harmonic peaks align with integer multiples of the carrier within the
  finite-pulse frequency resolution, with symmetry-forbidden components
  suppressed to a quantified ratio;
- `J(t)` is consistent with `dP/dt` after accounting for the documented sign
  and volume conventions;
- no forbidden-route marker appears in any production log.

Failure of any gate is reported as a limitation or a blocking result; analysis
parameters are not tuned after seeing the desired spectrum without recording
the change and rerunning the relevant controls.

## Verification and review

Implementation follows RED-first focused fixtures for collective observables,
restart-safe publication, spectrum normalization, and synthetic known-frequency
signals. Every affected Fortran runner is exercised on 1, 2, 4, and 8 MPI
ranks. Each task includes specification review, code-quality review, resolution
of all Critical/Important findings, `git diff --check`, and a final clean-first
MPI/ScaLAPACK/EigenExa parent-prerequisite overlay build before production
claims are recorded.
