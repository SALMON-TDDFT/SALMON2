# Source inventory and archive-module consolidation: first batch

Effective-source inventory is a read-only operation added to the existing
runtime_sources tool, not a new application module. It verifies full tree and
ordered complete-file replacement hashes, rejects invalid additions/order, and
records the supplying file. Actual prepared standard670 and fast-HSE725 trees
match the inventories exactly, and the installed tool produces identical JSON.
Their comparison is573 identical/77 differing/20 standard-only/75 HSE-only paths.
This does not certify semantic parity or consolidate the native trees yet.

GS archive internals move from three private modules into _gs_archive.py, without
forwarding stubs. Public APIs and separate GS-to-TDDFT ownership are unchanged.
Function ASTs match after only a private export rename and removal of a self
import. The existing source-inventory/export tests first failed for the missing
new entry/module before implementation. Existing tests were extended/repointed,
with no increase in unit-method count. Python modules70→68, private19→17.

Verification:

-628 source and628 installed-wheel test methods pass, one skip each.
-14 real-MPI archive cases pass in source and installed environments: k/orbital/
 spatial/mixed tiles, bounds/replica/overlap/gap/metadata/options/allocation/writer
 failures, exact retry and existing-output preservation.
- Real native MPI4 two-group PZ orbital2 import/independent states/current/direct
 no-replay resume passes from source and installed tools. Current SHA equals the
 pre-refactor runtime result. Synthetic GS input only; no SCF call.
- Offline sdist→wheel passes required source/object prohibition hooks. Installed
 native SHA matches its receipt. Old three module files are absent in the fresh
 wheel environment; the new internal module is present. No electronic rebuild.
- Source and installed LD/SBE examples execute. Si entry help is checked; this
 batch does not repeat a real SCF or certify the short Si smoke test's precision.
- Beginner first-three guide and both languages/AI snapshots are synchronized.

Canonical evidence contains logs, manifests/hashes and distribution artifacts,
not GS datasets:
/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/calc2026/salmonpy-consolidation-first-batch-20261119

Local writes verified; cloud sync unconfirmed. No CI dispatch. Existing fixed
runtimes and historical sources remain unchanged. Full fast-PR solver-source
unification and the grid/dt convergence workflow are subsequent work.
