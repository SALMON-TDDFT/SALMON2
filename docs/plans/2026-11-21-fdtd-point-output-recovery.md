# Optional point-observation I/O

Bounded R3 work following vector allocation1120. Optional observation I/O warns
instead of stopping a calculation that can still progress.

- Private Fortran point writer handles unchanged headers and append formatting.
  Check open/write/close IOSTAT,close only owned units,retain the first error,
  warn once per field state/output-owning rank. No new MPI collective.
- Future observations retry append after destination/header repair. Missing or
  partial rows are not repaired/replayed. Diagnostic warning flag is state-owned
  and excluded from the physical checkpoint layout. No production fault hook.
- Consolidated existing fixture adds real header/append directory and read-only
  faults,then repairs the path. No new Python test module or physical algorithm.

Record:`salmonpy-fdtd-point-output-recovery-20261121/summary.json` in canonical
calc2026. Prior runtime fails directory collision. Source and wheel GNU Mac
MPI1/2 cover18 fault/caller-unit cases;normal scalar/tensor outputs and physical
checkpoints are byte-identical to previous controls. During output failure,
fields/checkpoints remain equal;repaired rows match their exact expected suffix.
Zero-replay direct continuation and caller-file/MPI lifetime pass. Both Python
suites679 run,one optional skip,no failures. Source/object audits and tracked-only
preparation/build passed. No GS collection or optional CI dispatch.

Scope excludes other outputs,filesystem-wide failures and close-only fault
injection. Next:remaining supported plane-integral/info/required-output errors
and selected auxiliary buffers,then native contracts and group3 as ordered.
