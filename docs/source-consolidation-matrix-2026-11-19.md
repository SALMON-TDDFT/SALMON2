# Effective solver migration matrix

The raw fast-HSE reference is not the current fast-HSE engine. Resolve RT-state
and explicit-routing complete-file replacements before comparing sources.
`python tools/runtime_sources.py` prints verified path/hash/origin inventories
without creating work directories, executing code or changing a source tree.
Its results exactly match separately prepared standard670/fast-HSE725 files.
There are573 identical,77 differing,20 standard-only and75 fast-HSE-only paths.
These are byte comparisons, not a certificate of behavior/physics equivalence.

| Area | Current migration requirement |
|---|---|
| Retained GS / in-memory GS→TDDFT | Standard-only gs_model/storage/seed/embedded_gs and their initialization/input/communication contracts must move; raw fast-HSE cannot replace these APIs |
| Standard embedded TDDFT / failure handling | Retain embedded_tddft, worker and model-input interfaces; reconcile fast-HSE state/solver structures rather than overwrite either version |
| FDTD / LD / SBE | Retain embedded E/H, kernel, shared LD/SBE files; enforce native/solver shared-kernel equality |
| Fast PR HSE / symmetry / PT-CN | Keep current explicit-routing sources and host native_adapters; these are not replaced by legacy standard HSE |
| NLCC / XC / forces | Compare effective force/salmon_xc/semilocal paths and existing validation before porting; file-name or symbol presence alone does not establish that a fix is preserved |
| Checkpoint / MPI / model inputs | Reconcile differing implementations and model/partition identities; keep independent candidate restore and numerical rollback evidence |
| Build / source provenance | One maintained current tree with no-HSE and fast-HSE profiles; preserve origin commits/licenses/inventories and separately pinned reference evidence |

The first Python cleanup merges only GS archive internals. It does not consolidate
these native implementations. Three private modules become one, while public
gs_archive/gs_native/native_gs entry points and _gs_tddft lifecycle remain.
Original function ASTs are compared after only private export renaming and removal
of a now-unneeded self import. No compatibility forwarding modules are retained.

Current source trees are preserved pending native migration and parity checks.
Do not build or execute historical snapshots just because they have a recorded
name/hash. Keep the mandatory source/object prohibition checks on every new build.

The user-approved convergence workflow is separate: low-k PP-radius/3 to about
0.1 angstrom grid scan, then impulse dt scan from about0.2au; numerical failure
stops that trial only and nonbreaking dt differs from value-converged dt.
No real GS/SCF or convergence campaign is started in this cleanup batch.
