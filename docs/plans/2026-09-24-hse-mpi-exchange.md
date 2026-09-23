# MPI row-block HSE exchange implementation plan

**Goal:** Speed up the active HSE/TDCDFT response comparison using verified MPI exchange, then switch at an accepted checkpoint without changing physics.

**Architecture:** Rank0 owns the existing PT-CN/ACE driver, local Hamiltonian, history and files. Persistent MPI workers receive exchange requests, share occupied/target arrays, evaluate disjoint primitive row blocks, and reduce the unique-row actions to rank0. Other solver operations remain serial. The mathematical full HSE kernel, dt, field and checkpoint schema remain unchanged.

**Tech stack:** Existing Python/NumPy blocked exchange, OpenMPI5.0.9, mpi4py in an isolated venv using the existing NumPy installation. No native SALMON changes.

## Steps

1. Test row partitions including uneven/empty ownership and invalid rank/size; implement optional partition arguments preserving the serial default. Sum of partitions must reproduce full action, both full and finite radius.
2. Implement persistent MPI exchange service: root-only dispatch, array broadcasts, collective error reporting before reduction, root-only file output. Verify real multi-rank results and failure recovery with tiny synthetic orbitals before using Si.
3. Add optional supplied backend to resumable driver/job. Test incompatible method rejection. Root creates MPI service, workers serve requests; root closes workers on completion/error. No process other than rank0 calls file/analysis code.
4. Benchmark Si exchange with4 and8 processes, compare against serial on the same fixed snapshot; compare a complete resumed PT-CN step against serial from that checkpoint. Record communication-inclusive wall time and errors, select measured better configuration.
5. Independently review collective and restart behavior. Stop the currently authorized serial job with SIGINT only after the MPI path is validated, allowing its existing handler to save the last accepted state. Confirm locks are released, then launch the MPI job at the same target2750. If stopping cannot be clean, retain the last atomic checkpoint and document any replayed steps; do not run concurrent writers.

## Constraints

- No changes to exchange cutoff, ground state, dt=.32, impulse1e-4, 4³ k grid or analysis windows.
- Existing serial job continues while development/tests run. Tests/benchmarks use separate output paths and fixed checkpoint copies.
- MPI has replicated input arrays and per-process working memory. This is exchange parallelism, not whole-solver scaling; report actual wall times rather than promising process-count speedup.
- A failed rank must not leave the root waiting silently in a reduction. Python computation errors are exchanged collectively; unrecoverable MPI failures abort the job.
- Only rank0 may save/plot/analyze. Cooperative shutdown returns workers after stopping commands; filesystem locks still enforce a single active trajectory writer.

## Ledger

- User approved MPI switch. OpenMPI available, mpi4py installation in progress. Initial active serial state was step1135.
