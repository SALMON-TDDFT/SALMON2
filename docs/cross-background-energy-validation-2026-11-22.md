# Non-diagonal background field-energy verification — 2026-11-22

The current18-component field observer evaluates the old/new cross electric energy and outward constitutive exchange using the solver's actual symmetric background interpolation. The exchange is separate from six-face Poynting flux; both enter corrected residuals. The11-component summary preserves this separation. No old-format fallback is added. Rebuild the field runtime and lightweight kernel together. Derivation: plans/2026-11-22-cross-background-energy.md.

## Local evidence

Field-only40³ grid, spacing2 au, dt0.003 au,64 steps. Background epsilon has diagonal(2,3,4) and cross entries xy=.3, xz=.1, yz=.2. The observation region cuts the material at x-index22. Maximum absolute modified residual:6.088256317683939e-13 Hartree. Maximum absolute constitutive boundary exchange:4.308636920195699e-7 Hartree. This is a discrete consistency check, not physical convergence.

- Field MPI1/2 and x/y/z decompositions: samples and all budget values agree with rtol2e-12, atol2e-13. The y/z decompositions cut the tensor material.
- Observation/no observation and previous solver: final per-rank E/H, LD, PML and current checkpoint bytes match exactly.
- Direct checkpoint at step32, new MPI job, resume to64: all histories/budgets and final state bytes match continuous execution exactly.
- Legacy unit-background and diagonal LD controls: previous observation fields and final state bytes unchanged; new exchange field is zero.
- Source and installed final wheel:694 tests each, one optional skip each. Separate environment and Python -I; dependencies supplied from the validation environment, so this is not offline dependency-installation evidence.
- Final wheel observations, direct restart and state bytes match the source. Wheel implementation, paired manuals, AI JSON, solver inventory and changed solver files match the checkout byte-for-byte.

Current-source preparation and source/application-object prohibition audits pass. Ordinary field updates and checkpoint format are unchanged. Non-diagonal background observation adds two owning-field-group max agreements and two applications of background interpolation, each with six halo exchanges per active material, before the existing budget reduction. No world collective is added. Allocation/nonfinite failures leave output unpublished; unsupported observations warn without stopping propagation. Subpixel and Ohmic observations remain outside this scope.

Results: canonical calc2026/salmonpy-cross-background-energy-20261122/fdtd-summary.json and fdtd-verified/. The report links individual field cases in the existing diagonal-field-energy campaign. Local OneDrive storage is confirmed; cloud synchronization was not checked. No GS/SCF was run and no remote physical CI was dispatched.
