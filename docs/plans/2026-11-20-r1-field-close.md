# R1 native field restore-candidate lifetime

sp_eh_close must release pending and step_saved even on a repeated/empty close,
retain existing return codes and run eh_finalize only at capacity. No forced
final analysis on early close, numerical update/checkpoint format or MPI change.
Add optional rank-local native allocation flags for verification, no required
Python API/capability gate.

Existing synthetic I/O fixture executes the actual extracted close body with
allocator/finalization stubs across8 lifecycle combinations, before/after fix.
Existing real MPI E/H probe verifies early/full with/without candidate and saved
transaction payload; field/LD trajectories and ordinary output byte controls.
No GS data, no CI. R2 communicator ownership remains separate.

## Acceptance

Old extracted body fails with pending allocated; corrected8 lifecycle cases
pass in source/wheel. Real MPI2 early/full with/without pending plus saved
transaction cases all free allocation flags. Trajectory/prefix and full ordinary
output bytes match prior accepted output; early-close outputs unchanged.
Source/installed638:637pass1skip; audits and guide/tool/code parity pass.
No GS copies/SCF/CI. calc2026/salmonpy-r1-field-close-record-20261120.
R2 MPI resource retirement remains separate; current changes uncommitted.
