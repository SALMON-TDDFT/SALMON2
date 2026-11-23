# LD subpixel field-energy observation — 2026-11-22

The observer now uses the reference-frequency subpixel model's actual diagonal update coefficients and spatially varying off-diagonal background interpolation. The initial global subpixel_cross(0) flag selects background interpolation consistently across the field group. No additional communication is needed to discover that flag. Full LD component currents already include the transformed dispersive response. The existing18-component budget and separate constitutive_outward_energy apply unchanged.

## Local verification

40³ grid, spacing2 au, dt0.003 au,64 steps. An off-grid curved ellipsoid of scalar LD material has background epsilon4, plasma frequency20, one zero-resonance pole and damping2. A plane-wave source supplies reference omega20. The observation region cuts the material at x-index22; its upper corner is(27,27,27). One external response-point neighborhood is excluded by the existing smoothing rule.

-544 mixed samples, no geometry/reference fallback samples. Reference match error1.33e-15.
- Maximum absolute modified energy residual7.407680803024543e-13 Hartree; maximum constitutive exchange6.47046913725076e-6 Hartree.
- Grid-aligned interface control: maximum modified residual3.122502256758253e-16 Hartree and zero constitutive exchange.
- Field MPI1/2 and x/y/z partitions: samples and budgets agree at rtol2e-12, atol2e-13.
- Observation/no observation and previous runtime: final per-rank E/H, LD, PML and current checkpoint bytes match exactly.
- Step32 direct checkpoint, process exit and separate MPI job: full history/budgets and final state bytes match continuous execution exactly.
- Existing non-diagonal-background and legacy unit-background controls: histories and final state bytes unchanged.
- Source/installed final wheel:694 tests each, one optional skip each. Wheel direct continuation matches source exactly. Python -I imports the separately installed package; dependency packages are supplied by the validation environment. This is not an offline dependency-installation test.

Current-source preparation and source/application-object prohibition audits pass. The wheel's Python implementation, paired manuals, AI JSON, solver inventory and changed solver file match the checkout byte-for-byte. The field library is separately built and audited; it is not silently compiled on import. No backward compatibility path is added.

The model and propagation are unchanged. Background observation uses the previously implemented owning-group admission/allocation agreements and interpolation halo traffic only when required. No application-wide collective or checkpoint-format change.

This is discrete accounting for the reference-frequency model, not broadband optical accuracy or curved-interface grid convergence. Multiple materials/poles, other reference frequencies and other environments are not established by this single curved-material fixture. Fallback cells, if present in other inputs, are observed using the coefficients actually used in propagation. SBE/TDDFT are not smoothed. Ohmic dissipation remains a separate pending scope.

Results and artifacts: canonical calc2026/salmonpy-subpixel-field-energy-20261122/summary.json and verified-artifacts/. Individual field-case paths are recorded in that summary. Local OneDrive storage is confirmed; cloud synchronization was not checked. No new GS/SCF or optional physical CI was run.
