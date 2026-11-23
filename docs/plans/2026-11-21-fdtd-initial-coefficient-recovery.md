# Handled initialization coefficient/mode returns

### Retrying handled coefficient/special-mode initialization errors

Python also recognizes existing native cold returns1111 (nonfinite diagonal/
tensor LD update coefficients),1112 (nonzero Bloch phase with the unsupported
off-diagonal tensor mode), and1114 (invalid subpixel update coefficients).
Only unanimous codes with initialization recovery protocol1 release the failed
lifetime. The native initializer already deallocates its candidate, releases its
send/receive cache before its owned MPI children, and returns; the correction is
Python's classification and diagnostic text, not a native numerical change.

Correct the input/mode and use a new work directory on the same communicator.
This does not implement Bloch/off-diagonal tensor evolution, silently set Bloch
phase to zero, relax finite checks or change coefficient formulas. Files already
written during rejected initialization remain in the rejected work directory;
the rejection is not a promise of filesystem rollback. Unknown/mixed statuses,
missing protocol, trapped floating-point exceptions and allocation/device/MPI/
process failures remain outside this bounded retry guarantee.

--coefficient-admission in the existing E/H tool uses material-input-mode
(diagonal,tensor,subpixel) to exercise1111/1112/1114. Finite but oversized test
rates provoke the existing coefficient guards; a native periodic test input
provokes the existing unsupported Bloch/tensor return. GNU MPI1/2 repeat each
failure and then run valid32-step/direct-resume controls. Reader rejections are
combined in the same fixture. The old SHA256-verified native binary is reused;
no solver rebuild, extra numerical calculation in Python or new collective.
