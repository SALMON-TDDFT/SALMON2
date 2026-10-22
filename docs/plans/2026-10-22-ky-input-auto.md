# Automatic controls for KY/SOC inputs

Human scope: the SOC/NLCC/KY reader and CRPS generator only. Preserve explicit physical choices and stopped calculations.

- Resolve omitted `spin`/`yn_spinorbit` from the KY_SOC_NLCC_V1 header, or from one explicitly requested spinor/SOC setting. Legacy scalar inputs retain unpolarized/no-SOC defaults. Resolve before MPI broadcasts, existing compatibility checks, algorithm selection and variables.log. Explicit conflicts fail; no scalar averaging.
- CRPS: infer NA from active projector references, infer Lmax from a complete supplied 3/5/7-channel radius list, default lloc to the resolved Lmax, and derive output name from the configuration basename. Retain all explicit overrides and no-argument legacy Au defaults. Do not choose occupations, radii, local reference or NLCC core radius on behalf of the user.
- RED omitted-spin real reader smoke, missing input resolver unit, and omitted generator controls fail on the old implementation. Verify actual input helper plus full serial/MPI2 reader smoke and byte-identical explicit/automatic W and Se generated tables.
- No production restart, new physical GS/RT calculation, threshold relaxation or existing binary replacement.

## Verified outcome

- GNU checked serial Debug and MPI2 Release: 17/17 focused CTest PASS each (449–454 plus automatic reader smoke). Final coverage includes omitted W/Se lloc_ps and resolved spin/SOC in variables.log.
- Explicit versus omitted controls: all five full-application one-SCF reader configurations have exactly equal energies in serial and MPI2. This is an input-equivalence check, not a converged ground-state or spectrum acceptance.
- W/Se explicit and inferred generator settings produce byte-identical output with the same compiler; standard distributed radial data pass the existing portable numerical comparison.
- Invalid explicit spin/SOC combinations and out-of-range element counts are rejected. Scalar legacy defaults and no-argument generator defaults remain defined.
- Independent review found no important issue; its element-count guard suggestion is implemented and regression-tested. Coding Rules source/object prohibition and diff/Fortran line checks pass.
- Official SALMON-DOCS input_keyword_list.rst was updated in a local checkout; the corresponding patch is included in docs/SALMON-DOCS-ky-input-auto.patch. Apply the zero-context manual patch with `git apply --unidiff-zero docs/SALMON-DOCS-ky-input-auto.patch` inside SALMON-DOCS (using the patch path in this checkout). The SALMON-DOCS change is carried here as a patch; its remote repository is unchanged.
- HSE/LibXC, optional accelerators, ScaLAPACK/EigenExa, and physical 5fs calculations were not tested or restarted. No new environment variable or user-facing control parameter was added to SALMON.
