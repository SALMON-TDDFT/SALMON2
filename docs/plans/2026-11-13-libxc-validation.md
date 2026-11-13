# Libxc PZ/PBE validation plan

Goal: Enable the requested Libxc tests on the local Mac and validate the existing A-driven TDDFT boundary against ordinary SALMON.

Architecture: Preserve every GS/RT numerical routine and its update order. Build a separate OpenBLAS + Libxc 5.2.3 executable (the installed Libxc lacks the legacy f90 interface). If native GS/RT checks pass, extend only the worker's explicit allowlist to libxc_pz and libxc_pbe under USE_LIBXC. Keep other functionals unsupported. Do not infer pseudopotential provenance from a filename.

1. Build the existing source with Libxc; audit source and all objects before execution.
2. Run built-in PZ and Libxc PZ with identical existing Si data, compare GS and RT; run PBE with an independently identified PBE norm-conserving potential when obtainable.
3. Demonstrate the existing worker rejects Libxc PZ/PBE. Extend the allowlist only after that failing integration check.
4. Compare ordinary pre/post-change outputs and worker current for both functionals, predictor-corrector off/on, independent states and recovery using the existing integration runner.
5. Validate real converged GS -> ordinary/native-worker RT with the same potential and functional in each pair. Retain hashes, logs and limitations. Update supported-scope/build documentation.

No new DFT Python API, generic Libxc enablement, GS provenance schema, PAW/USPP support or convergence claim is included in this bounded test. Existing initial bundles remain immutable but do not certify their GS generation settings.
