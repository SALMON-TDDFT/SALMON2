# Fast-PR HSE numerical failure recovery

NativeFastHSEBackend serial/k-only MPI1/2, scalar fixed-ion HSE06 Taylor4 and
learned PT-CN only. Per-step electronic numerical backup plus explicit HSE
numeric_payload (ACE/endpoints/history/cache/counters), never copy FFT plans,
opaque owners or communicators. Return known deep failures collectively, restore
without propagation/exchange-history replay, preserve native/Python latches and
accepted clocks/currents/energy. Known protocol enables NumericalError; old or
unknown failures remain conservative RuntimeFailure. No mixed-step transaction.

Reuse stage finite checks and gate test injection separately from production.
Handle PT-CN convergence/Gram/retraction errors through optional status and FFT
cleanup; explicit ACE errors through owner status. Commit history before ordinary
RT output on protected paths so failed attempts append no accepted records.
Preserve batch paths and valid arithmetic; prove old/new exact on matched layouts.
Test actual overflow, one-rank staged faults including ACE/PT-CN/history, warm
histories, odd/even/first/last steps, independent groups, mixed/direct-restart
regressions, source/wheel and JP/EN/AI sync. Host64GB/18 cores; advisory48GiB,
OMP/BLAS1, no Actions. Allocation/process/MPI failures and uncontrolled external
STOP remain excluded. No conventional HSE engine or immutable binary edits.
