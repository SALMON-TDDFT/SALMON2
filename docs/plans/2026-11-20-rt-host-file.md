# Group5: preserve host-owned files at TDDFT shutdown

Native standard TDDFT closes unit30 despite never opening it in initialization
or the embedded worker. Remove that unowned close only in the native entry;
leave historical standalone main_tddft behavior unchanged. Owned ofile streams
and MPI children still close normally. No numerical/state/ABI change.

Extend existing I/O native fixture with a test-only host-unit C bridge. Audit
source and objects before linking; unchanged helper use checks SHA. Existing
standard handle probe opens unit30 on each material rank and verifies it stays
writable after last state close and native attachment retirement. Reproduce
old failure, compare new currents to old independent runs, test rollback/
sibling timers; reuse GS by references only, no SCF/GS copies or CI.

## Acceptance

Audited old runtime reproduces last-state host-file disconnection. Corrected
runtime MPI2/group PBE/rgrid matches old independent currents exactly, preserves
host file after all state/attachment closure and retains rollback/timer isolation.
Installed TB-mBJ/orbital2 all3groups preserve all6 host files through shared and
legacy singleton closure. Synthetic IO output parity/directory probes pass.
Source/installed638:637pass1skip. JP/EN/AI updated, no GS copies/SCF/CI.
Evidence: calc2026/salmonpy-rt-host-file-record-20261120; uncommitted.
