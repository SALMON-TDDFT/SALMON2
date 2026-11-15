# GS context storage extraction

Continue the approved GS/TDDFT isolation checklist. Move main_dft's typed electronic/work objects and HSE owner into an internal Fortran context, with explicit run/close lifecycle. Keep the external main_dft entry as a wrapper. Preserve statement order and output/checkpoint lifetime; move XC/HSE cleanup to close. Do not reinterpret global inputs or allow concurrent calls.

Tests retain A while running B under the same immutable scientific input, compare A's full stored electronic trace before/after B and after B's close, then compare scalar and restart output to the original executable. Also check double-run rejection and idempotent close. Existing serial/MPI and packaging regressions remain mandatory. This is retained state storage, not the final explicit-input GS API.

Completed retained-storage stage: 4 retained-context cases plus wrapper regression, 8 DC/MPI cases, 8 TDDFT cases, HSE OFF regression and fresh-wheel checks. See `docs/gs-context-storage-validation-2026-11-15.md`. Explicit input isolation remains on the existing open checklist.
