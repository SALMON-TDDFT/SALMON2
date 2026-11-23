# Early native grid/boundary/time retry

### Retrying handled early grid/boundary/time errors

Native E/H already returns1101 for inconsistent cell/spacing/grid-count input,
1102 for incompatible periodic/boundary faces,1103 for an invalid duration/
step-count/time-step specification, and1104 when dt is not below its existing
CFL bound. Python now releases the failed initialization reservation for these
recognized unanimous returns only when initialization recovery protocol1 exists.
The native cold candidate/owned-child cleanup remains the same. Correct the
input and retry with a new work directory on the same caller communicator.

No new input restriction, automatic repair, CFL threshold change, native binary
change or collective is introduced. This aligns Python's error classification
with existing native behavior. Unknown/mixed statuses, missing protocol and
other unhandled failures remain conservative. Rejected directories are retained;
this does not promise rollback of staged or diagnostic files. No automatic job
resubmission or solver execution is added to preparation/checking workflows.

The existing --grid-init-retry fixture now includes these four public-constructor
rejections twice each. Combined with --geometry-admission, GNU MPI1/2 complete35
rejection/retry checks per run, verify parent MPI/child cleanup, then32 steps and
direct checkpoint restore. Field/point/shape/checkpoint bytes match earlier valid
controls. This is a local lifecycle regression, not all initialization/IO-error
recovery or cross-platform/physical-convergence certification.
