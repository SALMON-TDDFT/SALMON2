# Full TDDFT context interleaving

User approved the next stage: preserve the HSE owner work, then test two full electronic states in one process against independent-process references. Keep GitHub Actions disabled for this work.

Use an internal Fortran context containing every argument of initialization_rt/time_evolution_step, its own HSE owner, XC handle and output handles. Read process configuration once; create each state independently from the same synthetic restart. Do not copy an initialized context or switch a global active context. Keep the exact native step and parity/update order. Drive states with different endpoint A vectors; compare current, energy, wavefunctions, densities, potentials and HSE histories with isolated runs. Test reversed order and destruction/recreation while another state remains live.

Initial scope is serial CPU, periodic unpolarized fast HSE, transverse A, fixed ions, common immutable global configuration and zero initial A. This is a development Fortran interface, not a new Python/C API. Reject unsupported modes. Inspect and correct shared initialization/scratch problems exposed by the test rather than treating a lower-level HSE probe as whole-TDDFT evidence. Differences in XC/grid/k count/dt across states remain a later configuration extraction.

Acceptance: a genuinely failing two-state test before fixes; real native step interleaving, separate-process references, comparison with an unchanged native full-driver calculation, close/reopen and both execution orders, source/object audits, local Python/package tests, English/Japanese/AI documentation. Preserve immutable reference binaries and real GS data; outputs remain external. No remote CI launch.

## Verification result

See ../tddft-context-validation-2026-11-15.md. Full-step noninterference, reverse order, close/reopen, immutable-driver pulse comparison and shared-owner mutation detection pass within the common-setting scope. Global configuration/scratch extraction remains the next stage.
