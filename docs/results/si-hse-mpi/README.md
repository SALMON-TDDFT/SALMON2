# Si HSE MPI exchange validation

Persistent MPI workers split primitive-grid row blocks of the full periodic HSE06 distance-kernel operator. Rank0 retains PT-CN, ACE, local Hamiltonian, checkpoint and final comparison ownership. No physical cutoff or parameter changed; Si8, shifted 4³ k grid, dt=0.32 au, impulse=1e-4.

## Measurements

| Processes | Exchange action (s) | Whole PT-CN step (s) | Step speedup |
|---|---:|---:|---:|
| 1 | 8.34–9.04 | 28.457 | 1 |
| 4 | 2.328 | 11.557 | 2.46 |
| 8 | 1.451 | 8.774 | 3.24 |

Exchange values are medians of three communication-inclusive MPI repetitions, against one serial action per benchmark: 3.58× for four ranks and 6.23× for eight. Whole-step measurements each advance the same accepted step1150 checkpoint to1151 once. Timings are indicative: the original serial production calculation remained active, and startup, checkpoint compression, and final plotting are excluded from step times. The remaining serial solver limits whole-step scaling.

All exchange actions and final wavefunctions agree bitwise with serial; current, energy, electron number, Gram error and full-exchange residual also agree exactly. Residual=8.84e-12, Gram error=2.42e-9, electron-number error=5.82e-11. See exchange_n4.json, exchange_n8.json and step_comparison.json. Selected exchange workspace is about255 MB/rank (not total RSS); orbitals are replicated.

75 unit tests pass. Real four/eight-rank smoke tests cover distinct source/target shapes, exact numerical parity, worker exceptions and root result-allocation failures, followed by successful service reuse. Unknown failures after dispatch abort the communicator to avoid a hung collective. Independent exception-path review found no remaining blocker.

## Run

Install mpi4py in an environment sharing the tested NumPy and a compatible MPI runtime; set OPENBLAS_NUM_THREADS=1. From the repository root:

```sh
mpiexec -n 8 python samples/hse_mlwf_reference/mpi_run.py \
  calculations/si_hse_reference/export \
  calculations/si_hse_reference/scf/state.npz \
  calculations/si_hse_reference/linear_response \
  calculations/si_hse_reference/comparison_tdcdft \
  docs/results/si-hse-tdcdft-linear
```

The default target is2750 total steps. This resumes a validated checkpoint and then runs the existing matched-window spectral comparison. Use --propagation-only and a separate copied checkpoint for pilots. Only rank0 writes files; existing advisory locks reject concurrent writers. Graceful service shutdown follows coordinated computation errors; a hard MPI abort relies on the latest atomic periodic checkpoint.
