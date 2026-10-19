# Elapsed-time checkpoint regression

`verification` compiles the shutdown expression from `main_tddft.f90` and checks
measured elapsed time against enabled, disabled, future, and exact thresholds.
The Fugaku observation was 7.32561350998003 seconds elapsed and a
123.6121805900475 second absolute inter-node offset. Subtracting the offset from
an elapsed duration incorrectly suppresses the checkpoint. The existing SALMON
clock implementation is unchanged; root broadcasts the shutdown decision so all
ranks leave the RT loop together.

Learned source-ACE checkpoints now write format v3 with orbital layout and SR
route metadata. Readers validate both before allocating restored factors. Older
v1/v2 files remain readable for `isize_o=1`; they are rejected for distributed
orbital layouts because those files lack the required ownership metadata.
Transient FFT plans and source caches are rebuilt rather than serialized.

The feature was also validated with 16 local MPI ranks and 32 Fugaku MPI ranks:
normal laser/impulse save at non-exact step 11 and resume to 13, plus time-only
shutdown during warmup at step 1 and resume through step 13. The 32-rank runs
completed with all rank exits zero, finite data, capacity checks, and distributed
SVD checks passing. Maximum resumed current difference was 2.30e-18 a.u. and
energy difference 5.69e-14 Ha. These short checks establish restart consistency,
not production-duration physical accuracy.
