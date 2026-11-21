# Long-wavelength optical correction

`salmonpy.analysis.long_wavelength_correction` implements the insulating low-frequency treatment of [Otobe et al., PRB 85, 045134, Appendix A4](https://doi.org/10.1103/PhysRevB.85.045134). Numerical fitting, conductivity conversion and correction run in the Fortran observation kernel. TDDFT propagation, measured currents and Maxwell feedback are unchanged.

```python
from salmonpy.analysis import long_wavelength_correction
# observed: OpticalSpectrum from optical_conductivity or dielectric_tensor.
# Example bounds are angular frequencies in atomic units, NOT electronvolts.
result = long_wavelength_correction(
    observed, fit_range=(0.02, 0.035), cutoff=0.02,
    system="gapped_insulator", units="atomic_gaussian",
)
raw = result.raw_dielectric
corrected = result.dielectric
report = result.report
```

Choose a trusted fitting interval below physical absorption onset; the example bounds are not validated material defaults. Fortran fits Re epsilon=c0+c2*omega² and sets Im sigma=-omega*(Re epsilon-I)/(4*pi) below cutoff, using charge-current conductivity and exp(+i omega t). I contributes only to tensor diagonal components. No extra speed-of-light factor is applied to existing SALMONpy conductivity.

The default `absorption="preserve"` retains Re sigma and Im epsilon, following the paper's real-dielectric replacement. It does not repair absorptive low-frequency artifacts. `absorption="analytic"` additionally fits Re sigma=omega²*(d0+d2*omega²); this is an extra insulating model assumption and permits explicitly modeled DC extrapolation. Default DC validity stays unchanged; the fitted static real dielectric is reported separately.

Opt in only for ordinary zero-DC gapped response. Do not apply to metals, Hall response or photoexcited free carriers. The declaration does not verify the assumption. The method does not infer a gap or a fitting window, impose a static dielectric value, correct a band gap, or guarantee global Kramers–Kronig consistency. The paper's 1 eV region and static value 12.6 are not universal defaults.

Results include independent read-only raw/corrected spectra, validity masks, applied-bin flags, coefficients, residuals and reference information. Unavailable tensor components stay masked. Insufficient fits warn and retain raw data; frequencies at/above cutoff are untouched. Synthetic regression verifies implementation and signs; no new real-material TDDFT validation is claimed.
