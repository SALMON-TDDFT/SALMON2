# Fast-PR HSE outer-MPI/FDTD direct restart — 2026-11-16

Base018b0182. New fast-only command worker, Python material provider/backend and native mixed-sample integration. No standard numerical field/LD/HSE kernels or preexisting binaries were modified. Worker rt_current queries physical initial charge current; rt_step already returns charge j, so no extra sign conversion or background subtraction occurs.

## Topology and interfaces

FDTD uses its native host-MPI field communicator, including unchanged local two-material LD updates. SBE occupies one outer rank. Each HSE point owns one serial fast-PR child on a single-rank material group; its supervising MPI rank waits for replies. This is outer MPI coupling, NOT HSE internal MPI. Tests use field1+SBE1+HSE1 (three outer ranks) and field2(y)+SBE1+HSE1 (four outer ranks), OMP/BLAS1, at most four active computational cores. Extra process IDs do not represent extra allocated MPI ranks; aggregate memory was not measured/enforced.

`tools/build_fast_hse_worker.py` verifies the pinned fast-PR explicit-routing source marker/configuration and unchanged origin binary/objects, audits sources/objects, explicitly compiles adapters and links a new executable. New build record includes worker/origin/source/object hashes. No conventional HSE or comparison baseline engine is added. Worker SHA256:e7de01d0243276780b153f719aebf777dfd17945a63286a3b73dcb674d80c430; fast origin:d0c47cd74c4911cb18fd5f8eb668b65bfac83e3293ce9471f8e4e2aef1422728. Field library SHA256:98f77a50cbfe66f18df7a96c20078a0216b5bb975823e9ba6a2a97ce7fb25627.

Fortran POSIX worker preserves the command stdin via /dev/stdin while the legacy Namelist reader reads inputfile. Flushed prefixed READY/STEP/SAVE/LOAD/CLOSE replies distinguish protocol from native output. Python host cwd is unchanged; import performs no compile/download/MPI initialization. Provider hashes/rechecks source bundle and stages independent work copies. Backend enforces single-rank execution, endpoint A, process locking, timeout/protocol checks and cleanup. Work directories inside input bundles are rejected before any directory mutation.

Direct coordinator admits explicit direct-state callbacks, writes each worker's native internal archive with field/PML/LD/SBE and Python state at one epoch, and verifies file hashes/canonical manifest before loading. Restore checks actual worker step/current, authenticated A versus response history, and resets material/response/adapter clocks consistently. One worker per rank is enforced to avoid filename collisions. Changed worker/configuration/bundle/topology is refused. Cross-group failure atomicity is not claimed; a failed native job is aborted.

## Fixture and rejected attempts

New synthetic real constant/cosine orthonormal orbitals on8³ grid,3k (3,1,1),2 occupied bands/4 electrons,cell8³,one Si atom(.3,.4,.45),existing KY Si pseudo. This time-reversal input has initial j approximately[-1.7311468703452483e-20,7.123912695700652e-24,1.4041606191058251e-24]au; no baseline is subtracted. The zero-field example rejects initial HSE current above1e-12au. It is a numerical connection fixture, NOT converged GS/optical material data.

The initial Gamma learned-mode attempt was correctly rejected by the fast engine; a first PT-CN input without explicit serial parallel settings also failed its existing admission guard. Failed runs are retained and are not evidence. Final fixtures explicitly select nproc_k/nproc_ob/nproc_rgrid=1 and use the supported dense-k learned profiles. No fallback engine or relaxed physics tolerance was introduced.

## Actual installed-wheel acceptance

Apple Silicon/macOS,Python3.14,OpenMPI. dt.001,256 steps,40³ field, two nine-cell LD materials,SBE and fast HSE. Learned ACE Taylor and learned PT-CN. Each case has a continuous job, a save job stopping at127, and a NEW MPI job/workers loading numerical state then stepping128..256.

| HSE | Field ranks | Outer MPI ranks | Saved step | Resumed step calls | Max full-history A/c,E,j errors | Payload bytes |
|---|---|---|---|---|---|---|
| Learned Taylor | 1 | 3 | 127 | 129 | 0,0,0 | 54,321,256 |
| Learned Taylor | 2(y) | 4 | 127 | 129 | 0,0,0 | 56,900,108 |
| Learned PT-CN | 1 | 3 | 127 | 129 | 0,0,0 | 54,321,283 |
| Learned PT-CN | 2(y) | 4 | 127 | 129 | 0,0,0 | 56,900,133 |

Strict np.testing.assert_array_equal plus raw array-byte comparison verifies entire A/E/j arrays against SAME-topology continuous execution. Field1/field2 histories are also bit-identical for each propagator. Every archive hash is checked. Actual resumed worker logs contain exactly129 STEP commands,first128,last256,following LOADED127; there are no elapsed1..127 STEP commands. Static initialization and an independent NEW-worker256-step reference replay are separate validation work. All independent SBE/HSE currents match exactly; clocks/A integration/previous-current feedback checks pass. Only external validation replays are run from0, not the resumed trajectory.

Seven protocol/lifecycle tests cover physical initial current/no sign flip, stale-A/capacity,independent save/load,timeout/malformed replies/reaping,hash/path/profile/rank rejection and nonmutating nested-root rejection, atomic-unit/external-A-only input validation and restored clocks. Two new coordinator tests cover inconsistent A before load and filename collision prevention. Final source and installed-wheel suites each459 tests,1 skip,no failures. Sdist-to-wheel succeeds, public API and JP/EN/AI guides synchronized. Read-only review found no important issue in ownership/counter/A/current restoration. The installed builder also linked a fresh worker successfully; a direct public-backend smoke test changed dt/capacity to.0015/4, saved at1 and resumed one further step, with exact current/energy, restored time.003, spaces in the archive path, unchanged host cwd and no host MPI import. A prior smoke test exposed a missing time property, now fixed and tested. Existing Libxc-PBE innerMPI2 archive continuation also passes the updated coordinator. No GitHub Actions used.

## Evidence and limits

Persistent developer root ~/.codex/outputs/salmonpy-state-inputs-20261116: fast-hse-worker/build.json, hse-coupled-fixture-k3-v2, hse-coupled-fixture-ptcn-v2; hse-mpi-final-{taylor,ptcn}-f{1,2}-{continuous,save,resume} and logs; hse-coupled-final-comparison.json; hse-mpi-complete-dist and hse-mpi-complete-{package,source-tests,installed-tests}.log; hse-mpi-complete-ptcn-resume; fast-hse-worker-installed/build.json and hse-worker-installed-smoke-v2/result.json. Earlier candidate builds are retained separately. Paths are evidence locations, not installation prerequisites.

Public Python backend preview and example checkpoint cover the tested POSIX serial-worker profile with same binary/environment/topology. Inner-MPI HSE/native host-library embedding,Windows/GPU/symmetry/localized/distributed profiles,changed-topology restore,generic public archives,recoverable native errors,resource caps,production scaling and physical convergence remain open. Native historical text-file positions are not resumed; Python response arrays contain the complete saved+continued history. Molecular work stays paused.
