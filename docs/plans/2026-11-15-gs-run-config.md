# Explicit GS SCF run configuration

Continue the approved input extraction. A value-only s_gs_run_config holds maximum SCF iterations, convergence threshold, pre-SCF threshold and consecutive-step requirement. gs_input_config captures legacy inputs once; gs_run accepts an optional config, validates before mutation and copies it into the system. SCF loop and threshold checks consume this owned configuration; unconverted entry points retain legacy fallback.

Do not change convergence metrics, update order, collective reduction or physics. Test invalid configs before initialization, poisoning legacy globals, mutation of the caller's copied value, two retained states and original separate-process output. Keep DC/MPI, TDDFT and packaging regressions. Other shared inputs remain explicitly open.

Completed: explicit selected values and invalid input paths verified; original serial/DC/TDDFT comparisons and wheel checks pass. Read-only review's semantic-validation and ignored-wrapper-status findings were fixed and regression-tested. See `docs/gs-run-config-validation-2026-11-15.md`.
