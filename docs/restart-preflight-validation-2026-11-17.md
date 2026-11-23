# Direct restart capability / MPI-layout preflight — 2026-11-17

## Implemented scope

The native mixed example coordinator now writes `salmonpy.native-direct/2` manifests with a `salmonpy.direct-compatibility/1` record: exact group names/ranks/point IDs, per-rank routes/formats,save/restore modes and field/electronic decompositions. Standard electronic layout is (k,orbital,x,y,z); field layout is (x,y,z). HSE k-only MPI and serial-worker routes remain distinct. Output-directory identities are excluded. Original input/library/GS hashes still bind the actual native payload. No native binary format, numerical formula, update order or solver binary was changed.

Before LOAD, every rank validates its capability and fresh candidate clocks, agrees current topology and within-group declarations, and compares the entire record with the saved manifest. Required field/electronic/worker binary filenames must appear in the inventory; existing distributed hashes verify every recorded file. Duplicate electronic/worker filenames per rank are rejected. All admission failures occur before native/state LOAD callbacks, leave candidate counters unchanged and do not invalidate the response as a failed numerical step. Unsupported saves fail before creating the output directory.

Static candidate/resource initialization still precedes sample preflight; there is no propagation, snapshot or restore call in compatibility checking. This is the supplied example's same-binary/environment/topology contract, not a general public archive API or automatic backend certification. Custom direct-state/1 callbacks need explicit checkpoint_format,checkpoint_modes and checkpoint_decomposition. Built-in LD/SBE snapshot routes are admitted; an arbitrary snapshot method does not prove direct persistence and may replay time evolution. Legacy archives without compatibility metadata require their original matching reader. Native internal validation still runs during LOAD; late I/O/allocation/backend failure may invalidate the job after partial state changes. No cross-group atomic rollback or changed-topology restart is claimed.

## Acceptance

Real three-rank MPI with field/material LOAD spies: exact layout accepted; eight collective rejection cases, zero LOAD calls:

1. Same total ranks, different k/orbital partition.
2. Changed group rank ownership/order.
3. Changed saved persistence format.
4. Missing compatibility metadata / legacy manifest.
5. One material rank declares no restore support.
6. Required binary payload absent together with its manifest entry.
7. Missing explicit format on an otherwise callable direct backend; no save directory created.
8. Valid but differing format declarations within one material group; no save directory created.

The missing-payload case reproduced actual premature LOAD-spy calls before the fix; the failed diagnostic log is retained. Unit tests also cover stale/advanced candidates,unknown snapshot-only backends and per-rank native filename collisions. No native runtime is needed for the manufactured MPI rejection test.

Installed-wheel real solver continuation, each within four outer ranks and one OMP/BLAS thread per active rank:

| Route | Saved step | New propagation calls | Full output arrays versus preceding continuous run |
|---|---:|---:|---|
| Field1/SBE1/fast-PR HSE k2, Taylor learned |127|129|byte-identical through256|
| Field1/SBE1/fast-PR HSE k2, learned PT-CN |127|129|byte-identical through256|
| Field1/SBE1/Libxc-PBE orbital2 |63|65|byte-identical through128|
| Field2/SBE1/serial fast-PR HSE worker1 |127|129|byte-identical through256|

Each field owns the existing two-material LD states. Equality includes point IDs,A/c,E,j and independent serial-reference arrays. Communication counters on every rank contain only the remaining calls. No elapsed propagation replay; fresh static initialization and independent reference validation remain separate. Same previously audited native libraries/inputs were reused unchanged. These are continuation/interface regressions, not material convergence, aggregate memory ceilings or scaling certification. Environment: Apple Silicon/macOS,Python3.14,OpenMPI,existing OpenBLAS/FFTW/Libxc profiles.

Source and installed-wheel suites each pass476 tests,one skip; source-distribution→wheel build,guide synchronization and diff checks pass. No remote Actions were run. Molecular work remains paused.

## Reproduction and evidence

Use the direct checkpoint recipes in `salmonpy/COUPLED.md` and `COUPLED.en.md`; saving now records compatibility automatically. The existing `--direct-resume` entry checks it before LOAD.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec --timeout 90 -n 3 python tools/test_restart_preflight_mpi.py NEW_OUTPUT
```

The manufactured test also runs against an installed wheel with `python -I`. Persistent local evidence: `/Users/otobetoshihito/.codex/outputs/salmonpy-restart-preflight-20261117`:

- `verified-comparison.json`: exact array checks,accepted steps,new call counts,compatibility records,manifest and source hashes.
- `*-save/direct-checkpoint/complete.json`, `*-resume-verified/result.json` and corresponding logs.
- `mpi-installed-verified.log`: eight preflight rejections; `mpi-missing-payload-red.log`: reproduced old omission.
- `source-tests-verified.log`, `installed-tests-verified.log`, `package-verified.log`, `dist-verified/`.

Remaining work: generic public persistence/capability contracts,changed-topology portability,native recoverable errors/handle ownership,larger-scale resources and additional platforms.
