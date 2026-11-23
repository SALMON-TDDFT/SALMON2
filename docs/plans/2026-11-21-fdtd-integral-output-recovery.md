# Optional plane-integral I/O

Bounded R3 continuation after point-output warnings. A failed optional diagnostic
must not abort valid physical evolution or skip peer MPI collectives.

- A private Fortran writer checks integral-header/append open,write and close.
  Keep first error,close only owned units,and bound warnings using field-owned
  state. Missing/partial rows remain incomplete;later appends retry after repair.
- Existing integral loops,unit conversions and MPI sums are unchanged. Root I/O
  failure never skips a subsequent plane reduction. Zero additional collectives.
- Public Fortran source calls retain compatibility via an optional final warning
  flag;no production injection hook/global state or physical checkpoint change.

Record:`salmonpy-fdtd-integral-output-recovery-20261121/summary.json` in canonical
calc2026. Old runtime fails directory collision. Source/wheel GNU Mac MPI1/2:
18 real directory/read-only/caller-unit fault cases,repair and zero-replay direct
continuation. Normal scalar/tensor bytes match prior controls;faulted xy output
preserves yz/xz/point files,fields and checkpoints. Caller file/MPI lifetimes pass.
Both suites679 run,one optional skip,no failures. Source/object audits,Git-only
preparation/build and sdist→wheel pass. No GS collection or optional CI dispatch.

No close-only failure injection or general filesystem-failure certificate.
Remaining R3 outputs:snapshots/frequency maps,info and required restart outputs;
selected auxiliary arrays remain a separate bounded concern. Then R4/group3.
