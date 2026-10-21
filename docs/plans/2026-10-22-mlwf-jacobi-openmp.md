# Jacobi OpenMP comparison (experimental)

Independent worktree based on c8ec8181a. Existing remote jobs and fixed binaries preserved. No push. PZHEEV forbidden in source and every linked application object; audits pass.

## Design

Independent local block gathering, angle selection, packing and unpack/update loops use OpenMP static worksharing. Transfer offsets are precomputed in original destination order; output spans are unique. MPI remains outside parallel regions (communication initializes MPI_THREAD_FUNNELED). Row updates begin after column update completion. Explicit recursive small eigen routines ensure automatic per-call work arrays. Matrix layout, round ordering, tolerances and capacity guards are unchanged.

## Local measured results

Same read-only Si64 elongated 8x1x1 GS, 128 occupied orbitals, 8 MPI. Jacobi time below sums maximum-rank sweep wall times, excluding initialization/functional checks. One measurement per setting; not a scaling claim. 8 MPI x12 OMP can oversubscribe this Mac. Full HSE06 and R10 masked actions evaluated on identical two fixed complex probes; this is cross-gauge comparison, not just same-gauge worker-route agreement.

|Version|OMP|Jacobi sweep sum (s)|Final gradient|Unitarity max error|
|---|---:|---:|---:|---:|
|baseline1|1|0.300085|1.300e-08|1.599e-14|
|baseline12|12|0.355506|1.300e-08|1.521e-14|
|finalomp1|1|0.233396|1.300e-08|1.554e-14|
|finalomp2|2|1.000334|1.300e-08|1.676e-14|
|finalomp4|4|1.627511|1.300e-08|1.577e-14|
|finalomp12|12|3.706707|1.300e-08|1.488e-14|

All 10 sweeps converge at unchanged 1e-7 gradient tolerance. All 8 ranks exit0. Actual team sizes 1/2/4/12 observed. Maximum full-action old/new relative L2 difference: 1.1666828288402034e-15. R10 masked-action maximum: 1.1676359351512204e-15.

The unconditional team launches slow this small workload; do not promote as a speed improvement. At 256 MPI ranks only about two angle solves belong to each rank in a 1024-orbital round. Native comparison is necessary before choosing a production configuration.

Round regression: MPI4 and OMP1/2/4/12, N1/3/8/17/65, block2/32, ordered/reversed ranks, all-skip and NaN rejection PASS. Reviewer found no important race/correctness defect; noted small-loop overhead and compiler-storage risks, addressed by explicit recursive work arrays.

## Native pending

Separate directory parallel-jacobi-omp, 64 nodes/256 MPI/6h. New unit-round gate then same Si512/1024 orbital input: mpi12 baseline, OMP candidate1, OMP candidate12. Each stage strict rank/finite/localization/5GiB RSS+Hugetlb gate. Cross-case spread relative threshold1e-8; full cross-gauge EXX validation currently local only. PJM final checked separately. No automatic resubmit.

GS is not downloaded; OneDrive cloud synchronization not verified.
