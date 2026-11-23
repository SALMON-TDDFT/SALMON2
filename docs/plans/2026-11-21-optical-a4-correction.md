# Observation-only long-wavelength conductivity correction

User requested PRB85,045134 (2012), Appendix A4, implemented entirely in the
Fortran observation/analysis kernel, with no TDDFT propagation/current feedback change.
Read and visually verified supplied PDF page11. Appendix states replacing Re epsilon
below1eV by a second-order frequency polynomial; no recipe is given for its fit band.

Explicit opt-in, declared gapped-insulator use only. User supplies a trusted positive
fit interval below physical absorption onset and cutoff (defaults to fit lower edge),
in Gaussian atomic units. Fortran fits Re epsilon = c0+c2*omega^2 and constructs the
corresponding Im sigma = -omega*(epsilon_real-I)/(4pi), respecting our charge-current
and exp(+i omega t) Fourier conventions. Default preserves Re sigma/Im epsilon exactly
(paper-real-epsilon operation). An explicit analytic absorption extension may fit
Re sigma=omega^2*(d0+d2*omega^2), enforcing the leading A4 insulating limit in both
parts; it must be labeled an additional assumption, not a literal paper prescription.

Return immutable raw/corrected sigma and epsilon, masks, applied-bin flags and native
fit coefficients/counts/residuals/static real estimate. Outside cutoff stays unchanged.
Missing tensor components are not zero-filled. Insufficient data warn and retain raw
component. A DC model value is explicitly extrapolated only in analytic mode;
otherwise keep the original DC validity and expose static Re epsilon separately.

No automatic fitting band, band-gap inference, experimental12.6 preset, Drude removal,
scissor shift, global Kramers-Kronig claim or solver modifications. Metal/free-carrier
and photoexcited response must remain raw unless the user isolates a regular part.
All least squares, masking, corrections and conversion are native Fortran; Python only
dispatches, records and plots. Tensor identity is applied to diagonal components only.

Test synthetic even dielectric plus a confined spurious conductivity pole, constrained
sign/slope, exact untouched high frequencies, preservation vs explicit absorptive
extension, tensor missing-channel masks, insufficient fits and input independence.
Build audited kernels, source and installed-wheel suites. No new GS or real-material
TDDFT run. Preserve pending LD-only subpixel changes; they are a separate workstream.

## Verification (2026-11-21)

Implemented in existing convergence_kernel Fortran observation module, exposed as
analysis.long_wavelength_correction with independent raw and corrected spectra.
Synthetic scalar/tensor regression passes; default absorption preservation and
explicit analytic extension distinguished. Public API and paired guides/AI JSON
updated. Source suite: 663 tests, one optional skip, no failures. No new GS or
real-material propagation performed. Pending subpixel work remains separate.

Installed-wheel suite also completed: 663 tests, one optional skip, no failures.
Validation record: calc2026/salmonpy-optical-a4-20261121/record.json.
