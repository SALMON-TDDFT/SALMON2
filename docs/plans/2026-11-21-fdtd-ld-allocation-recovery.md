# Initial LD state allocation recovery

Bounded group5 R3 continuation after core1118 recovery.

- The shared `ld_allocate` owns18 pole4D and6 current/old-E3D arrays. Optional
  status catches any nonzero STAT locally, returns before zeroing an unallocated
  array and leaves partial state to its owner. Status-absent callers keep their
  former allocation-failure behavior and original arithmetic/initialization order.
- Native initial allocation uses this status and one field-group agreement,
  returns1119 and reuses owned cleanup. Python permits protocol-confirmed,
  unanimous known1119 constructor retry; no budget-based stop or step collective.
- Existing testing-only selectors26/27 fail at LD slots2/20 on one field rank.
  No production selector, fault hook or C ABI change. Both maintained/shared
  kernel copies match; wheel CMake explicitly enables preprocessing for this file.

Measured record: `salmonpy-fdtd-ld-allocation-recovery-20261121/summary.json` in
canonical calc2026. Source/installed GNU macOS MPI1/2 scalar/off-diagonal-tensor
controls, partial-failure/retry and32-step direct continuation pass; all12 new
controls match prior runtime outputs/checkpoint bytes. Core1118 cases retained.
Both Python suites679 run,one optional skip,no failures. Source/object audits
pass and production exports no test fault hook. No new GS or heavy CI dispatch.

Not covered: material/coefficient scratch,tensor auxiliary,subpixel-specific
reallocation,PML/observation allocation,real physical OOM/OS kill or MPI failure.
Further R3 work remains bounded, then R4 and group3 as ordered in the checklist.
