# Native TDDFT deep numerical recovery — 2026-11-17

## Implemented contract

New standard native TDDFT exposes recovery protocol1. Each supported step retains the existing numerical checkpoint payload plus XC cval/operator payload in memory. Context-owned RT enable/error fields default disabled for ordinary batch callers. Collective scalar IEEE checks cover propagated orbitals, density, Hartree, XC potential/energy and endpoint current/energy. Protected propagator configuration errors return status instead of terminating. Handled errors restore accepted arrays, fields, histories, counters and original overlap flags, refresh PP phases, and return NumericalError through Python without propagation replay. Unknown/old errors remain conservative RuntimeFailure.

Supported guarded profile: periodic/unpolarized/fixed-ion middlepoint PZ, automatic/fixed TB-mBJ, Libxc-PZ/PBE. Unsupported SOC,+U,jellium,ionic/force,Gram–Schmidt,projection/unfold and microscopic Maxwell profiles reject before mutation. External checkpoint encoding is unchanged. This is a TDDFT-group guarantee; it does not roll back previously advanced FDTD/other materials in a ResponseSession. Process/MPI failure, allocation errors, uncontrolled third-party STOP, deep initialization and output-device errors remain excluded. Budgets remain advisory; per-step copies/checks add time and memory cost.

## Fail-first observations and correction

Old fixed library fails the new protocol regression. A fake-backend unit test initially returned RuntimeFailure instead of the expected protocol-declared NumericalError. First real overflow run returned safely but exposed a single serialized byte difference: external-LOAD normalization changed a preexisting overlap flag. Transactional rollback now uses the original in-memory flags and refreshes only PP phases. The subsequent strict byte checks pass; the failing run is not counted as evidence of success.

## Local verification

Local Apple Silicon, user-supplied64GB/18-core capacity. OMP/BLAS1, native build4 workers, jobs at most11 concurrent MPI ranks, no remote Actions. Full new source preparation and source/shared-kernel verification; complete owned-step replacements prepare correctly and the native step-config field/Taylor ownership probe also passes; explicit standard MPI/Libxc build and new-directory relink for final adapter. Fault module alone is separately compiled with SALMONPY_RECOVERY_TESTING. Production exports no test injection symbol. Existing fixed binaries are not modified.

- Actual finite input A=1e100 causes overflow during propagation. MPI2 returns NumericalError with A/step unchanged; serialized numerical state is byte-identical and valid retry succeeds.
- Compile-gated fault injection targets only group rank1 at five deep stages20–24, after real kernels have mutated state. Every failure preserves serialized state and the ordinary RT output files; following16 accepted-step currents match the old fixed binary exactly.
- PZ uses an8³ synthetic two-k fixture. Automatic/fixed TB-mBJ,Libxc-PZ/PBE use the existing small16³,2³-k Si regression inputs. All five profiles pass internal orbital MPI2; current and accepted current/energy numeric records match the corresponding old binary exactly. This is numerical/interface verification, not optical convergence.
- PBE also passes k-only MPI2 and real-space MPI2, with corresponding old-layout baselines exact; no claim of bit equality across different reduction layouts.
- Four ranks split into two independent2-rank groups: one receives overflow and one-rank faults while the other makes normal progress. First-step PBE/TB-mBJ initialization-cache cases, and odd/even wavefunction-buffer paths, all pass. Final-step failures immediately before scheduled energy output also preserve the state and yield exactly the old accepted output records after retry. Installed-wheel production and split/fault routes pass, with exact old-baseline currents.
- Normal installed5-rank mixed FDTD2 + two local LD materials + SBE1 + TDDFT2: all128-step response.npz arrays (A/E/j and serial replay) match the old fixed library exactly, maximum difference0. Existing late-staging-failure protected restore test also passes with original-byte preservation, independent continuation and exact3-frame retry; zero elapsed replay.
- Source and installed suites **525 tests passed, one skipped** each. sdist→wheel installation and guide JSON synchronization pass. Read-only review found no blocker.

## Reproduction

`tools/test_tddft_deep_recovery_mpi.py BUILD GS NEW_OUTPUT` with2 MPI ranks, timeout, OMP/BLAS1; --baseline records valid old/new trajectories. --faults requires a separately built test library. --split uses4 ranks with two independent groups. --fault-step 3 exercises odd buffers. --partition selects orbital(default), k or rgrid. All output paths must be new. No automatic compile/download or scientific-parameter changes.

Persistent artifacts: `/Users/otobetoshihito/.codex/outputs/salmonpy-tddft-deep-recovery-20261117/`. Prepared complete sources, production/test library receipts and hashes, build/test logs, per-case byte-comparison inventories, currents/accepted traces, representative raw before/after payloads, mixed results and wheel distribution retained. Paths are local provenance, not installation requirements.
