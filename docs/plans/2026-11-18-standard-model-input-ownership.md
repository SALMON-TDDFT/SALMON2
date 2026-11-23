# Standard TDDFT model input ownership

User approved first heterogeneous native ownership stage. Goal: make immutable
model initialization input explicitly owned, retaining legacy numerical kernels
and the single-live-model admission. No new heterogeneous handle constructor yet.

Create a typed model input record for the supported initialization sections,
including independently allocated coordinates/species and input unit strings.
Capture once after read_input/worker configuration. Own that record and initial A
in the native model object; create RT contexts from that explicit owner, restore
the legacy initialization view only at the synchronous compatibility boundary,
and recompute unit conversion tables from owned unit strings. Do not serialize
communicators, Libxc handles, output units or numerical states into the record.

1. Compile/run a focused ownership test (capture two configurations with different
   atom counts/XC/k/grid/PP/dt, switch them, mutate globals, deep-copy and release).
   Missing record must fail before implementation; check source before compilation.
2. Implement typed record and standard model initialization path. Keep default
   propagation and C/Python ABI restrictions unchanged. Common parser scratch,
   optimizer/cache/I/O/timer/framework ownership remains separate work.
3. Update pinned source inventories. Prepare a fresh external tree, audit before
   build, retain configure/incremental guards, audit every new application object.
4. Compare new/old standard native traces at five XCs, MPI1/orbital2/spatial2,
   interleaved owners and direct restart. Include normal mixed independent and
   protected save/retry, update JP/EN/AI/checklist, local-only packaging tests.

No new GS collection; reuse existing fixtures in temporary scratch. Store selected
results under canonical calc2026. No Actions. Budget excess remains advisory.

## Completed verification

- Owned227 configuration fields plus explicit native model/initial-A lifetime;
  actual standalone Fortran deep-copy/resize/absence/clone switching passes.
- Five XC profiles × MPI1/orbital2/spatial2:15 exact old/new native comparisons,
  interleaved contexts, idle recreation and direct restoration.
- Normal PBE/TB-mBJ2-point-each mixture:10 independent vs6 shared ranks,
  128-step continuous and63→128 direct restore/failure-retry exactly match.
- Raw MPI admission/retry and restored-A scalar NaN/Inf tests (12 cases) pass.
- Source/installed616 tests each, one skip; sdist→wheel/source inventory/native
  hashes and installed ownership + protected mixed restart pass. No Actions.
- Source/object audits passed; selected response records and audited runtime are
  stored under the user-designated calc2026 (no GS/native snapshots collected).
  Cloud sync is not asserted. Fujitsu compilation remains unmeasured.

Next: isolate parser defaults/scratch and introduce heterogeneous model registry
and activation without mutating live siblings; complete common derived resources
and heterogeneous Python handles. Current model-owned initialization is not
heterogeneous process reentrancy.
