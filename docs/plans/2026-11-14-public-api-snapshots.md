# Public API and in-memory snapshots implementation plan

Goal: implement the user-approved API/version/error contract, then fixed-dt
snapshots for LD, SBE and all four Python field components.

The 2026-11-14 handoff was read. Keep charge current, explicit contexts and
atomic failure. Its full-source import proposal remains superseded by current
AGENTS.md; no upstream source import, C ABI or RT context overhaul here.

Design approved in the preceding review and user's instruction to proceed:
- Version 0.1.0a5, explicit module-level public inventory and compatibility rules.
  Keep root imports lightweight. Package-authored validation/numerical/runtime
  exceptions derive from SalmonpyError and their legacy built-in families.
  OS/dependency/callback errors remain their native exceptions, explicitly documented.
- Opaque immutable StateSnapshot, schema/version checked, in-memory only.
  Restore requires the same model object and matching dt/time origin/capacity
  and sampled field drive; no cross-model, cross-version or disk restart promise.
- Copy all evolution history: LD P/J/previous E/counter; SBE rho/A/cached native
  current/counter/time origin; fields A/counter plus induced/reflected or full
  old/current/new grid and j. Immutable model/drive references stay with the state.
- Validate snapshot metadata and integrity and allocate all restored arrays
  before committing. Snapshot/restore must not advance or invoke a solver.
- Memory scales with copied state; no implicit filesystem/MPI/compiler work.

Steps:
1. Write failing tests for public version/errors and exception compatibility.
2. Add version/errors/API inventory, migrate explicitly raised package errors;
   retain existing RuntimeUnavailable and TDDFT exception names/inheritance.
3. Write replay/independence/rejection tests for LD, SBE and four fields.
4. Implement shared snapshot machinery and per-state history declarations.
5. Test initial and exhausted states, boundary modes, malformed snapshots,
   model/dt/drive mismatches, independent branches and coupled SBE/field replay.
6. Add English/Japanese API and snapshot guides plus runnable example; sync JSON
   using a repeatable tool and validate packaged guide equality/API imports.
7. Run the full Python suite, original numerical regressions, source-distribution
   to wheel build, installed suite and installed snapshot example. No push requested.
