# Group5: native E/H FDTD directory scopes

Extend native path coverage to Maxwell input/geometry/output/restart operations
and advertise directory protocol3. GS accepts1/2/3, TDDFT selects2/3, FDTD selects3;
older runtimes retain scoped-chdir. Reuse _native_group._directory_call with the
field communicator; no new module or public call arguments. Preserve E/H/LD
update order, current convention, numerical state fields and checkpoint layout.

Failure-first wrapper test and old-protocol fallback. Explicit source/object
audits and builds. Compare actual MPI E/H and local multi-material LD arrays
against audited previous binaries, then run SBE+two-LD layered coupling with
no-chdir guard, failed protected restore/retry and direct resume vs continuous.
Use existing no-GS sample configuration; do not collect new GS data or run SCF.
Check source/wheel, JP/EN/AI/sample guidance and actionable queue; no routine CI.
Main-thread sequential native calls and one field lifetime/process remain; path
scoping alone is not multiple FDTD instances or fatal MPI/OOM recovery.

Accepted: both audited profiles match prior MPI2 field arrays exactly; wheel
MPI4 two-LD/two-SBE old/new and direct resume exact, protected failure retains
all bytes and retry3 steps exact/replay0. Source/wheel637:636pass1skip each.
Evidence: calc2026/salmonpy-fdtd-directory-record-20261120. No GS or CI.
