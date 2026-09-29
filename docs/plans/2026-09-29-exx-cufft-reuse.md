# Reuse compact EXX GPU data and plans

User approved moving CPU/GPU transfer outside target batches.

- Own a stateful optional backend in the compact exchange plan. CPU remains default.
- Gather the compact source once per source orbital, before the target loop.
- Prepare GPU support indices, filter and source once per source. Compare retained data and upload only changes.
- Reuse cuFFT plans and buffers while padded dimensions, support size and capacity match. Zero inactive tail columns.
- Each target batch uploads targets and downloads actions. MPI collectives remain on CPU.
- Explicitly release resources on normal/error exits; retain fatal MPI abort for OpenACC runtime failures.
- Lifetime is one exchange-action call; cross-time-step persistence is not part of this change.
- Verify CPU stateful dispatch, collective failures and disabled-backend behavior; provide opt-in GPU parity and reuse-counter tests. NVIDIA execution remains unverified on the development Mac.

## Verification on 2026-09-29

- MPI and non-MPI hybrid builds passed.
- Stateful compact-dispatch oracle passed on MPI 1/2/4: one prepare per source, target tails, source replacement, skipped pairs and prepare/apply failure propagation.
- MPI 1/2/4 negative-release fault injection passed both at reinitialization and the exchange-action boundary. Signed error codes are normalized before maximum reductions.
- Conventional GS-to-RT, input and CPU cuFFT stub/fatal-callback tests: 13 passed; one real-GPU test explicitly skipped.
- Numbered CTest 422–437: 48 preparation/run/verification stages passed.
- GPU fixture now checks individual transfer/plan counters and FFTW parity after unchanged preparation, source/filter/index updates, capacity/geometry changes, tail batches and release/finalization. GPU compilation/execution remains unverified because this machine has no NVHPC/NVIDIA device.
- Independent static review completed; two cleanup-status findings were fixed and regression-tested.
