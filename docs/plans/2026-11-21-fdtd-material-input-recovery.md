# Material-file cold-rejection retry

### Retrying rejected LD tensor/subpixel input

The existing native readers already return1110 for missing/invalid diagonal or
reciprocal-tensor material sidecars and1113 for a missing/invalid reference
subpixel profile. For these cold errors the native initializer destroys the
candidate and releases its child resources. Python now recognizes these two
unanimous statuses with sp_eh_init_recovery_protocol1, detaches, and releases its
failed lifetime reservation. Correct/stage the material file and retry in a new
work directory on the same caller communicator. Error messages identify the
material input instead of reporting a generic native rejection.

No Fortran routine, input validation, coefficient arithmetic, file format,
communication or native binary changed in this correction. The existing audited
runtime is reused by recorded SHA256. Unknown/mixed statuses or an absent/wrong
recovery protocol do not enable retry. This does not extend recovery to numerical
coefficient errors1111/1114, unsupported Bloch/tensor mode1112, arbitrary allocation,
IO-device, MPI or process failure.

The existing E/H recovery tool accepts --material-input-mode diagonal|tensor|subpixel
and --material-input-admission. GNU MPI1/2 repeat missing, invalid-header, truncated
and nonfinite-file rejection twice each through the public Python constructor;
candidate/attachment resources are released and the parent communicator remains
usable. Valid material initialization then completes32 steps and direct restore.
Compare against the pre-change Python package with the same native binary; this
is an input/lifecycle regression, not a new material-convergence certificate.
