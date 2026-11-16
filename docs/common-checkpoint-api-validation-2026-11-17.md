# Common checkpoint capability / compatibility API — 2026-11-17

## Public contract

`salmonpy.checkpointing.CheckpointCapability` is an immutable declaration of versioned format,save/restore modes,actual decomposition and optional caller-selected model/input/runtime identity. Detached JSON conversion validates explicit keys and preserves the absence of an empty identity.

`check_checkpoint_compatibility` accepts a host-owned communicator or None for serial use, ordinary group mappings and explicit field/point declarations. No SALMON context is required. It agrees validation errors,operation,saved records,exact ownership and within-group capabilities before admission. Optional required-file names are checked against caller inventories. It performs no payload I/O,hash authentication,persistence callback,propagation,MPI import/initialization/finalization or native loading. Lifecycle checks and actual SAVE/LOAD remain caller responsibilities.

The native mixed coordinator now bridges its solver-specific lifecycle/ABI checks to this API. Its direct-compatibility/1 representation and native-direct/2 archives remain unchanged. Custom sample adapters may expose a public capability object/mapping; prior format/modes/decomposition attributes remain supported. Generic archive I/O,changed-topology redistribution and rollback are separate work.

## Verification

- Source and installed-wheel suites each pass485 tests,one skip. New coverage includes immutable/detached declarations,model identity mismatch,malformed capability maps,required-file/path checks,numeric-type distinctions in saved metadata,quiet import with MPI blocked,custom sample integration and preservation of custom declarations through the native coordinator. Empty local exceptions retain an error type so they cannot be mistaken for successful validation.
- Public API real MPI3 regression uses ordinary mappings and a caller-owned communicator without ParallelContext or a native runtime. Positive metadata roundtrip passes; ten failures are rejected collectively: unsupported mode,changed partition,operation disagreement,saved-record disagreement,missing inventory,local caller error,plan disagreement,changed model identity,missing declaration and empty local exception. MPI remains owned/live at return.
- Existing three-rank coordinator regression still rejects eight cases with zero LOAD calls, including missing native payloads and within-group format disagreement.
- Packaged custom_checkpoint example implements its own synthetic A-driven recurrence and JSON persistence with SHA256. Save17/load/new15 through32 gives zero current error, independent states and rejection of changed model identity before any load call. No native library or MPI is needed for this user-host example. It is not a material physics model.
- Source-distribution→wheel build succeeds. Installed paired guides,API inventory,AI catalog and packaged example are checked; a dedicated task is found by checkpoint compatibility search.

## Existing native archives and numerical invariance

The final installed wheel reads the archives created before the extraction without rewriting them. Their manifest SHA256 values equal the preceding validation record. The native solver libraries,GS/input bundles and numerical kernels were reused unchanged. Each case contains two field-local LD materials and SBE.

| Native route | Restored step | New calls | Full point-ID,A/c,E,j,serial-reference arrays versus continuous |
|---|---:|---:|---|
| Field1 / fast-PR HSE k2,learned Taylor |127|129|byte-identical through256|
| Field1 / fast-PR HSE k2,learned PT-CN |127|129|byte-identical through256|
| Field1 / Libxc-PBE orbital2 |63|65|byte-identical through128|
| Field2 / serial fast-PR HSE worker1 |127|129|byte-identical through256|

Every rank's communication counters show only the remaining calls. Static initialization and independent serial reference calculations are separate; no elapsed propagation is replayed. This validates metadata refactoring and continuation, not changed-topology portability,native reentrancy,atomic rollback,resource ceilings or physical convergence.

Local environment: Apple Silicon/macOS,Python3.14,OpenMPI and the existing audited solver profiles. Tests use at most four compute ranks with one OpenMP/BLAS thread; wheel compilation uses two workers. No remote Actions. Molecular work remains paused. Read-only independent code review found no blocker; its execution evidence is the separate tests above.

## Reproduce

```sh
python -m salmonpy.examples.custom_checkpoint NEW_OUTPUT
mpiexec --timeout 90 -n 3 python tools/test_checkpointing_mpi.py
mpiexec --timeout 90 -n 3 python tools/test_restart_preflight_mpi.py ANOTHER_NEW_OUTPUT
```

The MPI cases also ran with python -I against the installed wheel. Native examples use the existing COUPLED recipes and saved native-direct/2 fixtures. New user-facing guidance is in CHECKPOINTING.md / CHECKPOINTING.en.md.

Persistent evidence root: `/Users/otobetoshihito/.codex/outputs/salmonpy-common-checkpoint-api-20261117`.

- `verified-comparison-final.json`: exact byte comparisons,old manifest hashes,new call counts and final source hashes.
- `source-tests-final.log`, `installed-tests-final.log`, `package-final.log`, `dist-final/`.
- `public-mpi-installed-final.log`, `coordinator-mpi-installed-final.log`, `custom-installed-final/result.json`.
- `hse-{taylor,ptcn}-resume-final`, `pbe-resume-final`, `hse-serial-field2-resume-final` and their logs.
