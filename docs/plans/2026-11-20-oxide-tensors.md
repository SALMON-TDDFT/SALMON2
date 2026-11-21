# Oxide optical library and tensor interface

User requested optical oxides, then explicitly requested tensor permittivity.

1. Preserve the approved CC0 numerical acquisition workflow, fixed upstream revision and original bytes. Add SiO2, Al2O3, TiO2, Ta2O5, HfO2 scalar thin-film n/k datasets, preserving deposition/sample/dispersion-model conditions. Do not silently fill missing absorption data or substitute single-crystal axes.
2. Expose LorentzDrude.permittivity_tensor(omega) and OpticalMaterial.epsilon_tensor with shape (...,3,3), relative Gaussian dielectric tensor, xyz lab axes. Present scalar sources as explicitly assumed isotropic diagonal tensors; preserve scalar compatibility APIs. Native Fortran fills tensors and checks every scalar component. General anisotropic fitting/propagation is not implied.
3. Fit retained n/k without a Drude pole through explicit local SciPy trials, calibrate only the actual supplied range, retain poor candidates/warnings and quality metrics. Keep existing five-material coefficients unchanged.
4. Update paired guides, example, catalog/CSV/provenance and AI JSON. Test tensor layout/independence/finite rejection on synthetic data; integrity only for experimental data in required tests. Verify source suite and sdist-to-wheel installation explicitly, no new CI dispatch.

Latest user scope: common tensor representation in ALL modes, superseding the short-lived 3D E/H-only direction. Isotropic data use one independent diagonal scalar and zero off-diagonals; do not silently infer anisotropy. First implementation supplies the exact compressed case and common inputs in standalone LD/FDTD/MPI adapters. Independent three-axis and full off-diagonal native propagation still require solver changes and verification; an isotropic display tensor is not evidence of those capabilities.


## Reciprocal anisotropic LD implementation

Approved continuation: each pole owns a real symmetric positive-semidefinite strength tensor, scalar resonance/damping and shared plasma normalization; the instantaneous background is symmetric positive definite. This supports reciprocal diagonal/off-diagonal dielectric response, not arbitrary gyrotropic/nonreciprocal dispersion. The original one-component path retains its kernel and arithmetic. New native handle coefficients drive vector P/J; all export/set/rollback/checkpoint operations retain owned histories. Native computes tensor epsilon, background D/energy and validity checks using scalar IEEE inquiries. Python dispatches and records only.

Acceptance: rotated diagonal analytical spectra and time-domain scalar controls, negative off-diagonal couplings, invalid tensor rejection, independent state replay, retained scalar results. MPI provider preserves tensor config. Native FDTD needs additional tensor E-update/interpolation work, not a projection into scalar coefficients; unsupported general input must be explicit pending that connection.
