# Optional batch LR/ASE/ART spectral I/O

Bounded R3 continuation. Keep native unsupported return15 and required-save
policy unchanged;do not bypass protected-step contracts to enable this test.

- Guard LR trajectory/response and ASE/ART table open/write/close IOSTAT. Private
  per-call helpers close owned units and retain first error;local warning flag
  bounds messages. No new module/global selector/MPI collective or physical math.
- Window units/status/names are independent. Retain all energy-loop integration,
  MPI sums and subsequent Fourier/output work even after a write failure.
- Missing/partial optional files are not fabricated or retried after finalization.

Record:`salmonpy-fdtd-optional-spectra-recovery-20261122/summary.json` in canonical
calc2026. Source/wheel GNU Mac MPI1/2:isolated/periodic/background LR,ASE,ART;
20 normal-configuration comparisons match prior bytes,40 real directory/read-only
faults retain all other outputs and caller files. Spectral units and MPI finish
normally. Legacy background-normalization behavior preserved. Both suites679
run,one optional skip,no failures. Runtime/host source/object audits,tracked-only
build and sdist→wheel pass. No GS or optional CI dispatch.

Read-only faults measured on non-root Mac;close-only faults not injected. Batch
validation does not establish native spectral stepping/restart. Next bounded R3
work:required restarts,selected initialization/helper I/O and auxiliary buffers;
then R4/group3 as ordered in the active checklist.
