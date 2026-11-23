# CC0 optical measurements and stored LD candidates

This is a small pinned numerical subset of the RefractiveIndex.INFO database,
explicitly distributed under CC0. Listing/loading is offline: no download,
refit, compilation, MPI initialization or SALMON run. Unit and n/k conversions
use the bundled Fortran calculation routine.

```python
from salmonpy.optical_materials import list_materials, load_material
for row in list_materials():
    print(row['id'], row['recommended_fit'], row['fit_band_eV_actual'])
data = load_material('au-johnson-1972')
print(data.metadata['doi'], data.metadata['conditions'])
# Full source samples: data.omega (AU), data.energy_eV, data.nk, data.epsilon
model, report = data.model()  # stored coefficients; never refits
print(report['relative_rms_error'], report['optimizer_converged'])
# Explicit alternative: data.model(poles=8)
```

| ID | Source | Conditions | Candidate poles | Relative epsilon RMS |
|---|---|---|---:|---:|
| au-johnson-1972 | Johnson/Christy 1972 | Room temperature | 6 | 4.11% |
| ag-johnson-1972 | Johnson/Christy 1972 | Room temperature | 20 | 2.62% |
| cu-johnson-1972 | Johnson/Christy 1972 | Room temperature | 6 | 4.60% |
| al-mcpeak-2015 | McPeak et al. 2015 | Evaporated Al, 150 Å/s, template-stripped from Si | 4 | 0.162% |
| si-aspnes-1983 | Aspnes/Studna 1983 | Room temperature, <111>, n-type 2.3e14 cm^-3 | 6 | 4.65% |

The full original wavelength tables are retained. Fits request1.5–6 eV, while
fit_band_eV_actual records the actual included measurement endpoints; wavelength
rounding and discrete points leave a slightly narrower range. These are global,
unweighted complex-epsilon RMS residuals, not maximum pointwise, reflectance or
resonance-position errors. The selected SciPy fits satisfy their local optimizer
termination criteria. This is not a guarantee of a global optimum or physical accuracy.
Keep film, temperature, orientation and doping conditions in view.

The initial metal/Si subset retains66 trials:22 solver,22 cold-start SciPy and22 SciPy trials including
retained coefficients as a seed. Each pole-count record keeps the other trials in
alternative_trials; selection minimizes the weighted fitting objective. Defaults use
the smallest count with global RMS at most5%, recommended_dt*omega_max at least0.01,
and measured-point weighted error at most10%. These are screening heuristics.
Ag now uses the20-pole linewidth-constrained candidate described below. The
old12-pole model (RMS2.74%, local error26.7%) remains an explicit alternative. Al4 now has RMS0.162% with a sane
internal timescale; its previous stiff candidate is preserved in alternative_trials.
Si6 is the smallest passing candidate (4.65%); Si8 gives3.53% and is explicitly selectable.
Inaccurate or unfinished fits warn rather than stopping a runnable calculation.
Field CFL and timestep convergence remain separate.

n+i*k assumes mu=1 and exp(-i omega t). P/J remain dispersive poles; the
instantaneous background is epsilon_infinity. See LD.en.md. Load once and reuse
the immutable model to create independent states.

## Rights, identity and reproduction

Upstream: https://github.com/polyanskiy/refractiveindex.info-database
Pinned commit: e6d15e364f70de1cad8a96a20f75689f2fba8b21
Official distribution terms: https://refractiveindex.info/about (CC0-1.0).
Database citation: M. N. Polyanskiy, Scientific Data11,94 (2024),
https://doi.org/10.1038/s41597-023-02898-2.
Original sources: https://doi.org/10.1103/PhysRevB.6.4370,
https://doi.org/10.1021/ph5004237, https://doi.org/10.1103/PhysRevB.27.985.

optical_data/raw retains byte-identical YAML. CC0-1.0.txt and SOURCE_README.md
retain the upstream dedication/README; catalog.json records paths, URLs, commit,
checksums, references and sample conditions. Original journal articles, figures
and website text are not bundled. Source measurements and derived coefficients,
and the data dedication and code license, remain distinct. Loading checks source
and selected-parameter hashes.

Refitting is an explicit offline operation:

```bash
python tools/fit_optical_materials.py NEW_RESULTS --recorded-settings --optimizer scipy
```

The installed wheel carries this tool under its runtime SDK. --recorded-settings --optimizer scipy reuses each stored trial setting, including the additional Ag10/12/20 trials.
Omit it and supply --poles/--energy-range to explore a new candidate set. --export-catalog is an
explicit replacement of the exported library; no network or solver run occurs.

Refitting requires the optional `salmonpy[fit]` extra; loading stored models does not.

## Oxide thin-film extension

Ten datasets are retained. These are sample-specific scalar n,k films, not generic single-crystal values.

| ID | Sample | Default poles | Global relative RMS |
|---|---|---:|---:|
| sio2-rodriguez-2016 | Reactive e-beam film, substrate heated to300 °C | 1 | 0.0430% |
| al2o3-zhukovsky-2015 |20 nm ALD film on Si, deposition120 °C | 1 | 0.0332% |
| tio2-siefke-2016 |350 nm ALD film, broadband measurements/dispersion model | 8 | 4.47% |
| ta2o5-rodriguez-2016 | Reactive e-beam film, substrate heated to300 °C | 2 | 2.39% |
| hfo2-franta-2015 |112.7 nm e-beam film on c-Si, amorphous/fine-polycrystalline, room temperature | 1 | 1.29% |

The requested band is1.5–6 eV; discrete endpoints are recorded in fit_band_eV_actual.
Al2O3 reaches only about5.88 eV. The fitted band does not cover the full retained table.
TiO2 poles8/10 meet residual screening but exhausted the optimization budget; warnings remain.
Many original optical tables are dispersion-model or Kramers–Kronig reconstructions
from measurements, not unprocessed measured points. n-only datasets are not silently
assigned k=0. There are88 retained trials:66 previous plus22 oxide trials.

## Tensor representation shared by all modes

Relative epsilon[i,j] uses lab xyz and exp(-i omega t). data.epsilon_tensor has shape(N,3,3);
model.permittivity_tensor(omega) returns omega.shape+(3,3), immutable complex arrays.
Scalar n,k sources without directional information store one independent component:
epsilon_ij=epsilon*delta_ij. Off-diagonals are zero; all diagonal entries are equal.
This is an explicit isotropic model assumption, not experimental evidence of isotropy.
The metadata records tensor_representation and axes.

```python
from salmonpy.optical_materials import load_material
from salmonpy.ld import LorentzDrude

data = load_material("sio2-rodriguez-2016")
measured_tensor = data.epsilon_tensor
model, report = data.model()
tensor = model.tensor_config  # LDTensorModelConfig, one independent component
response = LorentzDrude.from_config(tensor)
# FDTD(geometry, materials=(tensor,), dt=..., steps=...)
# LDMaterialProvider(tensor) also accepts the common configuration.
```

FDTD projects the compressed isotropic tensor exactly to the existing coefficients,
retaining arithmetic and update order. The A-endpoint LD adapter accepts the same
configuration but returns dispersive-pole current only; the field owns the instantaneous
background. These scalar source tables do not determine general anisotropy.
Supply independent tensor data for that purpose; see LD and FDTD_MATERIALS guides
for tensor fitting, rotations and off-diagonal field updates.
No supplied component is silently dropped.


## Ag refinement with an explicit linewidth assumption

The twenty-pole Ag default requests1.5–6 eV (actual included endpoints about
1.510–5.981 eV). Relative complex-epsilon RMS is2.6184%, weighted RMS2.4506%,
and maximum measured-point weighted error8.4798%. The unchanged screens are
5% global RMS,10% maximum weighted error and recommended_dt*omega_max>=0.01.
The local SciPy optimizer terminated by ftol; this is not a global-optimum claim.
The full source tables, original six pole-count candidates and their alternatives
remain available. There are now89 retained trials across the library.

All twenty pole damping parameters were constrained to at least0.002 atomic
units. This is an explicitly chosen modeling assumption, not a measured
linewidth. It prevents the almost-zero-width poles found in unconstrained
refinement trials. Those trials fit the supplied points better but can produce
large unmeasured inter-point peaks, so they were not adopted. A10001-point
frequency grid plus every in-band resonance center was evaluated for the chosen
model: finite passive response, maximum Im epsilon about5.8502 in the fitted
band. This checks the model, not experimental behavior between measurements.
Time-domain accuracy and field CFL must still be checked for the actual problem.

Replay the recorded initial model, linewidth floor and optimizer budget explicitly:

```sh
python tools/fit_optical_materials.py NEW_RESULTS \
  --materials ag-johnson-1972 --recorded-settings --optimizer scipy
```

To explore other conditions, use `--minimum-damping VALUE` (atomic units) with
explicit `--poles`, `--starts` and `--max-iterations`. The floor defaults to0;
loading a material never fits or changes a user's model. Existing coefficients
for every other material and all original numerical source files are unchanged.
