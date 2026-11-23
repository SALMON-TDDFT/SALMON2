# Diagonal tensor native E/H and waveplate verification

User approved native anisotropic FDTD connection and birefringence validation.

1. Preserve scalar/one-component path byte-for-byte. Connect principal-axis diagonal LDTensorModelConfig through explicit media_type ld-diagonal and a hashed input asset, not process-global input arrays. Common resonance/damping per pole already supports independent axes through tensor strengths.
2. Each E component uses its own background/c2 and semi-implicit c4/c5; reuse existing owned P/J/Q and checkpoint arrays, update order and Yee-edge ownership. Full off-diagonal coupling is not approximated by applying a local matrix to staggered fields.
3. Reject old runtimes before execution if they lack the new capability; stage/check assets in both batch and NativeEHField. Validate root read collectively before any use. Preserve scalar source/outputs and reject missing/corrupt diagonal records.
4. Source audits, source manifest/profile updates, tracked-only prospective Git export and new external runtime build/object audits. No modifications to fixed prior runtimes, no GS or new Actions.
5. Compare old/new scalar calculations, diagonal isotropic limit and MPI1/2. Verify quarter/half-wave slabs vs analytic normal-incidence slab transfer including Fresnel/multiple reflections; native diagnostics supply phasors and Stokes. Retain physical verification outside CI.
6. Synchronize paired human docs, AI JSON and samples; run source and installed-wheel checks.
