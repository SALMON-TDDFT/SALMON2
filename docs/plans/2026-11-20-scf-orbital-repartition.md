# Ordinary-XC orbital-only SCF restart repartition

Continue priority4 using the existing rank-local disk format and explicit
run(restart=...,repartition=True) opt-in. Add optional orbital-repartition ABI1;
existing k ABI and default exact-layout semantics remain unchanged. Save global
orbital bounds beside k bounds. Admission allows only pure orbital1↔2 layouts,
ordinary PZ/Libxc-PZ/PBE, identical physical inputs, binary and assets.

Fortran reuses donor decode/staging and maps the orbital axis of wavefunction and
CG buffers, preserving replicated density/mixing/Poisson/XC/residual state and
binding continuation parallel inputs to live resources. No SCF replay or Python
numerical kernel. HSE orbital repartition and mixed/spatial topology remain outside
this step. Verify exact cut slices, zero replay, source-initialized continuous
trajectory agreement at existing1e-9/1e-10, checksum/failure handling, prior binary
normal-path parity, existing k/HSE migration and installed-wheel parity. Extend
existing tests/tool rather than adding another suite/module. Update paired guides,
AI JSON and current checklist after numerical verification.

## Verified

Source and installed637 tests (636pass,1skip); final-build ordinary PBE2→1,
PZ1→2 and Gamma Libxc-PZ1→2/PZ2→1; exact restored cuts and zero replay,
continuation at preset tolerances. Installed Gamma migration and sample passed.
Ordinary/HSE prior-binary SCF and16-step RT exact; previous k/HSE reduced-k
migrations still pass. Final builds source/object audited. Failed attempts that
started before the build receipt existed are retained as setup-failure evidence,
not solver failures. No CI. Record: calc2026/salmonpy-scf-orbital-record-20261120.
