# Native E/H FDTD numerical recovery — 2026-11-17

## Implemented contract

New standard host-MPI libraries declare E/H recovery protocol1 and expose a checked experimental step accepting the current host A alongside j. Each supported step retains a detached allocatable field/PML/LD/current/clock payload in memory; communicator and halo-cache owners remain live. Optional eh_calc gates check E (20), LD (21), H (22), sampled E/H, integrated A and scheduled ordinary observation arithmetic (23), before scientific output. On a handled collective failure the adapter restores the payload directly and returns a known status. Python raises NumericalError without changing A/E/step/failed. No propagation replay, dt change or automatic amplitude correction occurs.

The legacy batch eh_calc path has no guard callback and retains the original arithmetic/order. Native guarded E -> LD -> H ordering is also unchanged. The external checkpoint format is unchanged. TTM, dedicated linear-response, ASE/ART and frequency-observation extensions reject before mutation pending their own side-effect audit. Initialization/allocations, process/MPI errors, uncontrolled STOP and output-device faults remain excluded. This protects a field group, not a whole already-progressed mixed CoupledStepper/ResponseSession step. Diagnostic progress stdout may survive a failed attempt. Backups/checks add runtime and memory cost; budgets remain advisory.

## Verification evidence

Local Apple Silicon, host64GB/18 cores, advisory48GiB, OMP/BLAS1, build4 workers and at most5 MPI ranks; no GitHub Actions. Prepared current verified bundled sources into a new external directory, performed a full standard MPI/Libxc build, then relinked the final adapter in a new directory. Production has no fault-injection symbol. Only a separately compiled adapter has SALMONPY_RECOVERY_TESTING and the one-rank hook. Source inventories and the dependent step-config base identity were updated; shared LD/SBE kernels remain byte-identical.

- Real40³ E/H field with PML and two distinct local LD materials,32 steps. MPI1/2 plus two independent MPI2 groups on four ranks. One rank rejects stages20/21/22/23 at first1, warm12 and final32. All rank-local rejected payloads are byte-identical before/after; Python references/values/clocks/failed are preserved; failed attempts append no observation records.
- Production finite extreme A/current exercises arithmetic overflow at the pre-output gate and permits valid retry. These stress values are not physical inputs. No claim of production timestep/mesh convergence follows from this test.
- All accepted32-step A/E/H arrays and ordinary observation numeric records match the old fixed binary exactly. Three-step observation stride also passes, including unscheduled steps whose old-H output scratch must not be read. Serial/MPI equality holds for this fixture, not a universal reduction-layout claim.
- Four-rank independent groups demonstrate that one group's rejection/retry leaves the other group's normal evolution unchanged. Raw mixed legacy/checked entry calls collectively reject11 with unchanged bytes; raw exhausted capacity rejects13 with unchanged bytes. The wrapper still handles old/unknown failures conservatively.
- Two separate normal128-step five-rank regressions: FDTD2 + local LD + SBE1 + Libxc-PBE TDDFT2, and FDTD2 + local LD + SBE1 + fast-PR HSE2. Old/new A/c/E/j maximum difference0 in both. Both existing protected restore tests pass after late material staging failure: all original bytes preserved, live continuation works, three retry frames exactly match, zero restore replay.
- Installed wheel independent-group MPI fault/overflow path passes too. Source and installed suites528 tests passed with one skipped; sdist-to-wheel build/native hashes, layout test and Markdown/AI JSON sync pass. Local dependencies were reused offline; imports do not compile or access the network.

## Review and negative evidence

The new Python regression initially failed before the checked call/protocol mapping existed. Review found mixed legacy/checked calls could enter different collectives; mode is now agreed before optional-A broadcasts and a real raw ABI test covers it. A subsequent review found no remaining blocker in the supported path. Initial test harness runs lacked the field-leader flag / accessed SampledField instead of its values; those failed harness runs are not counted as solver evidence. The first source suite also caught a stale step-config pinned source identity; its current-base hash was synchronized without weakening verification, and all tests then pass.

## Reproduction

`tools/test_eh_deep_recovery_mpi.py BUILD NEW_OUTPUT` with one/two ranks, OMP/BLAS1 and an external timeout. `--split` selects four ranks/two independent groups; `--faults` requires the separately compiled test adapter; `--baseline` skips recovery probes for old binaries. `--observation-stride 3` covers scheduled/unscheduled output. All destinations must be new.

Persistent evidence: `/Users/otobetoshihito/.codex/outputs/salmonpy-fdtd-deep-recovery-20261117`. Retained production/test libraries, source snapshot/build receipts, accepted arrays/logs, complete hash inventories and selected rejected before/after payloads. Libraries and checkpoint payloads are host/compiler/topology-specific evidence, not portable releases.
