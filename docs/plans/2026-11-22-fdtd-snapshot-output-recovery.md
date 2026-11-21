# Optional plane-snapshot I/O

Bounded R3 continuation after point/integral/final-info output handling.

- Check snapshot open/write/close IOSTAT. Warn on failure,close owned units,
  retain the first error,and preserve all subsequent plane/component reductions
  and staging-array releases. No added MPI collective or process-global state.
- Keep normal formatting,unit conversion and physical checkpoint layout.
  Existing Fortran source calls accept an optional trailing warning flag.
- Later times write their own files. Missing frames are not fabricated;after
  destination repair,user-driven direct-restored continuation can write that
  same time again. Restore itself never replays the numerical evolution.

Record:`salmonpy-fdtd-snapshot-output-recovery-20261122/summary.json` in canonical
calc2026. GNU macOS MPI1/2 source/wheel:10 real directory/read-only/caller-file
fault cases,repair and direct continuation. Normal scalar/tensor outputs and
physical checkpoints match prior bytes;all90 sampled frames match after explicit
continuation17/18. Caller files and state/MPI cleanup pass. Both suites679 run,
one optional skip,no failures. Source/object audits,tracked-only build and
sdist→wheel pass. Stride6 reduces fixture I/O without changing API defaults.
No GS collection or optional CI dispatch.

Scope excludes staging allocation,close-only fault injection,frequency maps,
required restart I/O and filesystem-wide failure. Permission faults were tested
with a non-root Mac. Remaining R3 work stays bounded,then R4/group3.
