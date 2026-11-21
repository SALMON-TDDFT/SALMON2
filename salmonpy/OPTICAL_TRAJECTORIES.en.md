# Reuse stored optical trajectories

No calculation, new GS, network access or MPI is started by these APIs. All new
transforms, weighted summation and differences use the Fortran observation kernel.
Python selects arrays, schedules cases and records metadata. Keep original data.

## Duration and window sensitivity

```python
from salmonpy.analysis import optical_sensitivity, dielectric_function
result = optical_sensitivity(
    time, charge_current_x, total_E_x,
    frequencies=[0.0, 0.02, 0.04, 0.06],
    sample_counts=[2000, 4000], windows=["none", "cubic"],
    current_convention="charge", field_kind="total", units="atomic_gaussian",
)
reference = result.spectra[result.reference_index]
epsilon = dielectric_function(reference)
report = result.report
```

Inputs are scalar simultaneous components and angular atomic frequencies. Counts
select prefixes of the SAME record, keeping its initial time. Each returned
conductivity uses the same explicit frequencies: no interpolation, zero padding
or hidden simulation. The descriptive reference is the longest prefix and first
requested window. Relative RMS differences use only common valid frequencies;
a zero reference or no common data yields unknown, not a convergence pass.
Results are independent read-only arrays; metadata retains counts, duration, dt,
window, nominal resolution2*pi/(N*dt), Nyquist, masks and input content SHA256.
Selecting many nearby frequencies does not increase the information in a short
trajectory. Out-of-Nyquist and weak/zero-drive bins are masked. DC conductivity
may be valid while DC dielectric stays masked.

Windows use u=n/(N-1): none=1, hann=(1-cos(2*pi*u))/2,
cubic=1-3*u^2+2*u^3. The causal cubic keeps an initial impulse; Hann removes an
impulse at the first sample and can therefore leave no usable drive spectrum.
The transform is dt*sum(window*x*exp(+i*omega*t)), without gain correction or
baseline subtraction, consistent with existing spectrum(). The initial absolute
time phase is retained. The cubic window is a postprocessing choice, not a
change to propagation. Shortened records may truncate the drive or ringdown;
a small difference is not proof of physical convergence, linearity, k sampling,
grid accuracy or sum-rule completeness. No automatic pass/fail threshold.

Requested-frequency cost is O(samples*channels*frequencies), not FFT cost. Use
selected optical bands, not a huge dense grid. Large work counts warn but do not
stop. Inspect duration/window sensitivity before fitting a material or applying
OPTICAL_CORRECTION; do not silently discard bad cases.

## Weighted k contributions

```python
from salmonpy.outputs import read_decomposed_current
from salmonpy.analysis import k_current_contributions
raw = read_decomposed_current("Si_current_decomposed.data")
k = k_current_contributions(
    raw, dt=0.02, time0=0.0, current_convention="matter", units="atomic_gaussian",
    frequencies=[0.0, 0.02, 0.04], window="cubic",
    total_current=simultaneous_charge_current_xyz,
)
# k.current: (sample,k,xyz), k.total: (sample,xyz)
# k.spectrum.values: (frequency,k,xyz), if frequencies were supplied.
```

SALMON's writer already multiplies band occupation AND k weights. Never multiply
either again. Native reduction only sums supplied band/spin entries by their
positive k IDs. Declare the file's sign: matter means J=-j and is converted once;
charge is kept. These are charge-current spectra, not independent dielectric
functions of k points. Relative RMS contributions may exceed one under signed
cancellation and are not probabilities. No missing bands/spins/k points are
inferred; complete-looking labels alone do not prove completeness.

Native step indices times run dt plus time0 give the sampled times. dt is the
INTEGRATOR step, not observation-stride*dt. Optional total_current is already
physical CHARGE current, regardless of the decomposed source sign, at exactly
the same sampled times. Align it explicitly; no resampling or offset inference.
The report gives absolute/relative RMS and max differences. A zero reference
has unknown relative error with the absolute error retained. A mismatch is
reported, never repaired by scaling. This is not a k-point convergence test.
Known non-AU file units must be converted before declaring atomic_gaussian.
Unknown header units require the caller's explicit unit declaration.

## Stored-file sample

Prepare NPZ arrays time/current/E in Gaussian AU with charge current and total E.
The example selects one Cartesian component and keeps optional k analysis separate.

```sh
# Preview: no calculation and no output directory
python -m salmonpy.examples.optical_postprocess trajectory.npz NEW_OUTPUT --frequency-range 0 0.1 64
# Explicit postprocessing only, no solver
python -m salmonpy.examples.optical_postprocess trajectory.npz NEW_OUTPUT --frequency-range 0 0.1 64 --windows none cubic --run
```

Optional --decomposed-current FILE additionally needs --decomposition-dt and
--decomposition-convention. It does not infer time alignment with the NPZ or
silently compare a different set of samples. Numerical correctness is tested on
synthetic signals, including weighted cancellation. No new real-material
accuracy, GPU/Fujitsu or distributed execution claim is made.
