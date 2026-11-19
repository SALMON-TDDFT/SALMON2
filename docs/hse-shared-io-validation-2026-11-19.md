# Shared initialization-input and RT-output infrastructure

Fast explicit routing now references the existing standard IO CMake list,
model_inputs and step_output files. The byte-identical adapter step_output copy
is deleted. Host builders use the compiled solver object; for compatible older
prepared builds lacking it, the same standard source is staged, compiled, audited
and linked. No new Python module or public method is introduced.

The standard s_model_inputs value record is now available to the fast solver as
preparation for GS/ordinary TDDFT migration. Its capture/activate roundtrip,
independent deep assignment and allocated/unallocated array lifetime tests pass
using actual prepared fast Fortran sources. It does not replace the larger HSE
snapshot or make the ordinary GS/parser/TDDFT registration lifecycle available.

## Serial profile findings and changes

A fresh USE_MPI=OFF build found an unconditional reference from symmetry to the
MPI-only comm_is_attached. A standalone worker then found unconditional child-
communicator management imports/calls. Both are now guarded by USE_MPI. Serial
child-owner metadata starts at zero; no MPI resource exists there. MPI owner
allocation/release and collective rejection are preserved. Symmetry validation
for the actual drive is retained; nothing is projected or silently disabled.
No new IEEE inquiry or eigensolver call is introduced.

Final source builds host-MPI fast HSE, serial fast HSE plus worker, and MPI with
HSE OFF. Required source/object audits, including external adapter objects,
remain enabled. Old fixed sources/binaries are not modified. An earlier attempt
to relink a much older serial object set failed on its absent symmetry_check_field
interface; that unrelated ABI mismatch is retained as evidence, not bypassed.

## Numerical evidence

- Internal MPI2, six ranks total: two different full/reduced HSE models match old
  independent current, energy, learned history, fresh state and direct restart.
  Closing registry does not retire live states; inconsistent selection rejects.
- Four ranks: ordinary LD/SBE/two different full-k HSE model sample128 steps
  matches the previous validated fast-source result. Direct step8 checkpoint
  resumes with120 newly propagated steps, zero replay and exact result.
- Source and clean installed-wheel suites:628 test methods each, one skip;
  no new test method or module added.
- Two-rank native field recovery preserves state on real overflow/admission
  failure and retries successfully.
- Serial Taylor and PT-CN learned:24 steps match old worker current/energy;
  save8/resume16 matches uninterrupted continuation exactly, zero replay.

These are connection/regression tests on existing synthetic fixtures. They do
not establish physical convergence. No new SCF, GS dataset collection or CI run.
Source/installed-wheel test results and source hashes are retained with this
record in calc2026/salmonpy-hse-shared-io-20261119. Local files are verified;
OneDrive cloud synchronization is not independently confirmed.

Remaining work: parser workspace, GS numerical owners, ordinary TDDFT lifecycle
and remaining effective source differences. Current inventory fast732/standard670:
587 identical,70 differing,13 standard-only,75 fast-only. Five alternative HSE
files are not migration targets. Source consolidation and the approved SALMON-
specific convergence workflow remain unfinished.
