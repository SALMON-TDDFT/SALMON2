# Explicit RT construction configuration

Continue approved input extraction: add a value-only s_rt_run_config for dt,
capacity and the three OpenMP selectors. rt_input_config captures existing input
once; rt_open(config=...) validates and copies it before native allocation.
Legacy delta_t/openmp_mode remain supported but cannot be combined with config.
Native initialization consumes its owned capacity, dt and OpenMP selection. No
global input swapping. Common grid/k/XC and other settings remain guarded.

Tests: A capacity 12, B capacity 6, independent dt/OpenMP, isolated references,
reverse order, close/reopen, rejected over-capacity calls, mutated source config,
poisoned old global nt, unchanged native 32-case reference. Keep audits,
Python/sdist-wheel tests and English/Japanese/AI docs. Local only, no Actions.
