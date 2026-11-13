# User test starting from a wheel

[日本語](wheel-user-test.md)

Test installation and five official examples without editing SALMONpy: LD,
synthetic SBE, pulse/spectrum, geometry, and output roundtrip. This is an execution
and I/O smoke test. It does not repeat the earlier notebook's additional analytic
solution/conservation checks. Full DFT/TDDFT/FDTD propagation/MPI, real materials,
physical convergence and Jupyter UI are outside scope.

## Give the participant

- A `.whl` from a successful LD wheels CI matching OS, CPU and Python. Extract the
  outer Actions ZIP. Record the source URL, CI run and commit separately.
- This guide and `tools/test_wheel_user_workflow.py`. The runner is a source-side
  testing tool, not a module installed in the wheel.
- Compatible Python, a writable new output location and network access. Fortran,
  CMake, external SALMON and Matplotlib are unnecessary. For offline use provide
  compatible NumPy wheels in a wheelhouse too.

Suggested AI request:

> Install this wheel into a fresh environment and run the five specified examples.
> First ask for missing memory/CPU budgets. Use one numerical-library thread and
> sequential runs. Record installation, calculation, scientific validation and
> resource measurement separately. Never silently switch to source compilation or
> modify SALMONpy. Preserve failures and use new output directories for retries.

## Run

Use Python matching the wheel; replace uppercase placeholders with actual paths.
The output directory must not exist and must be outside the source checkout.
The same Python command works on Windows.

```sh
python tools/test_wheel_user_workflow.py WHEEL_PATH NEW_OUTPUT --memory-budget "2 GiB"
```

Use `--memory-budget unknown` if unknown. This records an allowance, not a memory
limit. For offline installation add `--wheelhouse DEPENDENCY_WHEEL_DIRECTORY`.

The runner creates a fresh venv and installs only binary distributions; there is
no source-build fallback. Normally NumPy is fetched from the package index.
Isolated Python (`-I`), a new working directory and a package-path check prevent
checkout imports. It then runs each official example once in its own process:

1. `ld_response`
2. `sbe_response`
3. `pulse_spectrum`
4. `geometry_preview`
5. `outputs_roundtrip`

Independent examples continue after another example fails. Installation or
package-identity failures prevent calculation. There are no automatic retries,
time limits or hard memory limits. Interrupt an unexpectedly long run and retain
its records.

## Evidence and interpretation

`report.json` records the wheel SHA256, environment, budget and stage outcomes.
Each stage retains `.stdout.txt`, `.stderr.txt` and `.json`; SBE, geometry and
archive examples retain their outputs. `installed_identity.stdout.txt` records
actual package versions and import location.

Calculation exit codes and measurement status are independent. An unavailable
clock does not turn a successful calculation into failure, and a nonzero
calculation exit remains a failure. Peak memory is `not_measured` / null for every
stage, never zero or an invented estimate. Overall success does not establish
compliance with the memory budget or independent scientific validation.

Also retain the initial AI request, AI name/version, manual interventions,
confusing steps and elapsed time including preparation. Keep retries separate
and describe changes. Review personal paths/usernames before sharing; share
reports, logs and results rather than the entire venv.
