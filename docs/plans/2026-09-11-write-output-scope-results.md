# write.f90 output responsibility audit: DG cleanup

The DG additions to `write.f90` included stored polarization history, a Fourier
transform with smoothing, and susceptibility/dielectric/conductivity assembly.
Those are not formatting or output responsibilities.

Reference audit found no production call to `write_dg_polarization_data`; its
only external reference was an unused import in `main_tddft.f90`. Thus the only
code that could allocate/populate `dg_polarization_history` was unreachable.
`write_dg_polarization_response_3d` returned immediately on unallocated history.
These were dead remnants, not a used response pipeline to relocate.

Removed both routines, the saved history, unused imports, and the no-op response
call. This deletes 122 lines from `write.f90`. The active DG observable
evaluation and sample-writing routines remain unchanged. Existing files and
saved calculation evidence were not deleted. Code is recoverable from Git
history before this commit (parent `a543148a`).

## Verification

The new absence/retained-route assertion failed before the deletion and passed
afterward. All nine Task 6 numerical-integration contract groups passed, along
with the distributed-v5 architecture and overlapping-Wannier route checks.
The production build exited 0; its log is
`/tmp/si8-route-cleanup-20260911/build-write-cleanup.log`. `git diff --check`
passed. No new physical calculation was run for this dead-code deletion.

## Remaining scope

This does NOT make all of `write.f90` output-only. The conventional SALMON
response/pulse routines still contain Fourier transforms and response-data
construction. Those live calculations need a separate, behavior-preserving
extraction with numerical and output-format regression coverage. No claim is
made that that wider migration is complete. New DG data preparation should not
be added to `write.f90`.
