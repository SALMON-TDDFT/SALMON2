# Explicit RT construction configuration — 2026-11-15

s_rt_run_config contains value-only dt, maximum accepted steps and three OpenMP
selectors. rt_input_config snapshots native input; rt_open(config=...) validates
and copies it before allocation. Legacy dt/OMP keywords remain available but
cannot be combined with config. Finite-positive dt, selector names and capacity
1..huge(default integer)-1 are checked. This is not memory estimation/enforcement.

Initialization allocates RT arrays with owned max_steps, including endpoint
slots, and stepping uses owned capacity as both the acceptance bound and final
step argument. It no longer selects capacity from global nt. Common grid/k/XC
and other unextracted input settings remain guarded. No active global switch.

## Evidence

- Red probe against prior build failed on the missing configuration type/entry:
  `/tmp/salmonpy-rt-config-red-20261115/build.log`.
- Eight synthetic full-step cases, OMP=2, k=2/3, four propagator/factor modes:
  A capacity12/dt0.001/orbital, B capacity6/dt0.0015/rgrid. Separate-process B
  uses input nt=6. Current/energy/raw state traces match isolated references
  bit-for-bit for interleaving, reverse order and close/reopen. B rejects attempts
  past step6 without affecting continuing A. The receipt's steps=12 denotes the
  outer driver horizon, not B's accepted steps.
- Invalid zero/overflow capacity and mixed config/legacy keywords are rejected
  before allocation. Tests poison old global nt and mutate source-options
  capacity during calls; accepted state remains independent. dt/selector
  non-aliasing is established by value-copy inspection, not a separate mutation
  test of those fields.
- Original native Acos2 driver comparison passes rtol1e-12/atol1e-14, including
  current sign. Legacy default construction: 32 cases vs immutable fast-PR
  original, formatted current/energy maximum differences zero.
- Source/all-object eigensolver audits pass before linking. Fresh completed-source
  inventory matches; no new executable native lines exceed 132 columns.
- Python source tests:305 passed. Sdist-to-wheel build passes. Fresh wheel:305
  tests, one optional plotting skip, no failures. Installed preparation passes
  audits, and all replacement source bytes match the checkout.
- Read-only independent review found no scoped blocker. Japanese/English/AI
  documents synchronized. No Actions launched; uncommitted at report creation.

## Limits

Local Apple Silicon synthetic serial fixed-ion scalar-HSE only. No physical
convergence, heterogeneous grid/k/XC, MPI, GPU, concurrent host calls or complete
native failure atomicity claim. Force/microAc/AVX/SVE paths retain previous
validation limits. A full native input-object interface remains future work.

## Receipts

- `/tmp/salmonpy-rt-config-context-20261115/result.json`; SHA256 `8b6ec36278b37128159c1bca4c873f098f3bb4d27ba0ce2a81cbe0094e9a7cc1`.
- `/tmp/salmonpy-rt-config-default-20261115/result.json`; SHA256 `63276b12c9dc04205e0cf2e5aa54252914923722fcf0624503fb490690034afb`.
- `/tmp/salmonpy-rt-config-build-20261115.json`; SHA256 `5370f375a5bc66f273cd5e23b7e6661c94a851486d659e82b8751d52f00bdff5`.

Supporting logs: `/tmp/salmonpy-rt-config-python.log`,
`/tmp/salmonpy-rt-config-packaging.log`,
`/tmp/salmonpy-rt-config-installed-tests.log`,
`/tmp/salmonpy-rt-config-wheel-prepare.log`.
