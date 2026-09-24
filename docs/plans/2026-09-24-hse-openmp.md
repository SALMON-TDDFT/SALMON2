# Native HSE MPI + OpenMP

User approved hybrid implementation and fixed-total-core comparison. Python
production reference explicitly terminated at user request; never restart it.
Use existing TDCDFT checkout. Preserve exchange, ACE and propagation physics.

Parallelize independent local-k BLAS operations and tile packing/kernel loops;
keep MPI collectives and shared FFT plans on the calling thread outside OpenMP
regions (SALMON requests MPI_THREAD_FUNNELED). ACE application uses one small
private overlap workspace per thread. ACE construction stays serial initially:
its eigensolver/error paths are small and measured cost is minor. BLAS remains
single threaded to avoid nested oversubscription. No thread-private orbital
or FFT work copies. FFT itself remains serial per rank; document that limit.

Test independent MPI kernel with OpenMP1/2/4, shifted/shuffled mesh, arbitrary
trials and invalid rank-local inputs; expose actual kernel team size as useful
runtime metadata. Verify missing implementation via failing team-size test.
Then run full unit tests, native Taylor and PT-CN parity/restart, fixed-total16
MPI16x1/8x2/4x4 two repeats with memory monitoring and explicit no binding.
Review races and record actual speed/memory (no assumed hybrid speedup).
Publish implementation and experimental limitations with results.

Final design supersedes the proposal above: production-sized integration found
incorrect results with concurrent external ACE BLAS callers on this local
runtime combination. ACE therefore uses one caller and one overlap buffer.
The initial fixed-core measurements forced BLAS to one thread. Follow-up in
2026-09-24-hse-blas-openmp.md moves every kernel BLAS call outside OpenMP and
removes that force, allowing internally threaded vendor BLAS.
