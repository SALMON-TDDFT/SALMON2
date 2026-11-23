# Common backend-owned save/restore execution — 2026-11-17

## Interface and ownership

CheckpointTarget binds an explicit capability and zero-argument save/restore callbacks. execute_checkpoint validates declarations and callback availability collectively, then calls only the bound backend functions. Paths,JSON/binary/other formats and in-memory tokens belong to the backend. The return dictionary contains compatibility JSON and local field/states results. Opaque callback values are neither copied nor MPI-serialized; callbacks signal failure by exceptions, not return-value truthiness.

Groups execute their field/point order in synchronized component slots. Groups can own different point counts; an error on one rank stops every rank before later slots. RuntimeFailure reports callback failure, retaining an originating exception as the local cause. Completed callbacks are not rolled back. A collective backend must complete its own MPI protocol before raising; this entry cannot recover a stranded native collective,Fortran STOP or process crash. Hosts retain MPI ownership,hash verification,completion publication and candidate lifetime policies.

The custom JSON example now performs actual save/load through this API. The native mixed coordinator also dispatches field/PML/LD,electronic and SBE snapshot I/O through it, retaining codecs,native-direct/2 metadata and existing lifecycle/inventory checks. It performs planning admission before directory setup and rechecks admission immediately before execution. Completion marker publication errors are now caught/agreed on all ranks rather than leaving peers in a final barrier.

## Acceptance evidence

- Source and installed-wheel suites each pass489 tests,one skip. New tests cover opaque result identity,callback availability/compatibility failures with zero calls,stopping a later point after an earlier failure,and actual backend-owned LD snapshot restore with exact future results. Existing quiet import/admission tests remain passing.
- Real MPI4: field1 plus a one-point rank and a two-point/two-rank group. Save/restore callbacks use host-owned subgroup collectives. Opaque lambda results remain local and restore exactly. A missing callback prevents every callback; a rank3 exception after its group collective produces RuntimeFailure on every rank and stops the next point. MPI is not finalized by the API.
- Real MPI3 admission regression:ten collective negative cases. Existing coordinator preflight:eight rejection cases with zero LOAD calls.
- Real MPI3 publication failure injection: root's completion rename raises OSError. Every rank receives RuntimeFailure and no complete.json appears.
- Packaged custom JSON model:save17,load once,new15→32; current equals continuous exactly,original/restored states are independent,changed model identity rejected before load. This is synthetic interface data,not a physical material prediction.
- Sdist→wheel,installed example/API inventory,paired English/Japanese guides and AI catalog checks pass.

## Actual native save and direct continuation

The installed wheel executed NEW saves through the common entry,then restored them in separate MPI jobs. Existing audited solver binaries and input bundles were unchanged. Each case has SBE and two field-local LD materials. At most four outer ranks and one OMP/BLAS thread per active rank; wheel compilation uses two workers.

| Route | Saved / new calls | Full output arrays against preceding continuous baseline |
|---|---|---|
| Field1 / fast-PR HSE k2,learned Taylor |127 /129 through256|byte-identical|
| Field1 / fast-PR HSE k2,learned PT-CN |127 /129 through256|byte-identical|
| Field1 / Libxc-PBE orbital2 |63 /65 through128|byte-identical|
| Field2 / serial fast-PR HSE worker1 |127 /129 through256|byte-identical|

Equality includes point IDs,A/c,E,j and independent serial-reference arrays. Every rank's communication counters show only remaining propagation calls. Compatibility records equal those from previous saves. A separate HSE-k2 job also restored the pre-existing native-direct/2 archive with byte-identical results. No elapsed propagation replay; static initialization and independent serial validation remain separate. No native numerical source,kernel or update formula changed.

This verifies dispatch and continuation on Apple Silicon/macOS/Python3.14/OpenMPI. It does not certify changed-topology restore,native reentrancy,atomic rollback,resource ceilings,scaling,other platforms or material convergence. Only fast-PR HSE is used; molecular work remains paused. Read-only independent review found no executor blocker; its adjacent completion-publication finding was fixed and fault-injection tested.

## Reproduction and persistent records

```sh
python -m salmonpy.examples.custom_checkpoint NEW_OUTPUT
mpiexec --timeout 90 -n 4 python tools/test_checkpoint_execution_mpi.py
mpiexec --timeout 90 -n 3 python tools/test_checkpoint_publication_mpi.py ANOTHER_NEW_OUTPUT
```

Use the existing COUPLED save/resume recipes for native tests; execution now goes through the common dispatcher without changing CLI arguments or archive codecs. CHECKPOINTING.md / CHECKPOINTING.en.md document callback binding and local return values.

Evidence root: `/Users/otobetoshihito/.codex/outputs/salmonpy-checkpoint-execution-20261117`.

- `verified-comparison.json`: exact arrays,remaining call counts,unchanged compatibility records and source hashes.
- `source-tests.log`, `installed-tests.log`, `package.log`, `dist/`.
- `execution-mpi-installed.log`, `publication-mpi-installed.log`, `admission-mpi-installed.log`, `preflight-mpi-installed.log`.
- `custom-installed/result.json` and user-owned JSON files.
- `{hse-taylor,hse-ptcn,pbe,serial-field2}-{save,resume}` and `old-archive-resume`,with logs.

No remote Actions were run. Remaining work concerns archive interoperability/migration,changed topology,native recovery/ownership and resources; backend payload encoding remains intentionally backend-owned.
