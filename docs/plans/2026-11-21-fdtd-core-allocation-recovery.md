# Native FDTD core-allocation recovery

Implemented bounded group5 R3 work; do not expand this to all legacy modes.

- Material ID/scratch and core split-E/H/coefficient/saved-E/H 3D arrays use
  STAT on native initialization. Negative or positive STAT failures normalize
  to one flag; remaining allocations in that phase are skipped.
- Two owning-field-group startup agreements precede shape/halo work and
  downstream coefficient setup. No per-step or application-world communication.
- 1118 uses existing cache-before-communicator cleanup and protocol-gated,
  unanimous known-error Python constructor retry. Successful allocation and
  zero-initialization order is unchanged. Budget warnings remain advisory.
- Test-only stages24/25 fail after partial material/core-field allocation on
  the last field rank. Selectors reuse the existing testing-only hook; no new
  production global or injection export is introduced.

Measured record: `salmonpy-fdtd-core-allocation-recovery-20261121/summary.json`
in canonical calc2026. Source/test/wheel GNU Mac MPI1/2 controls and direct
resume match the preceding runtime byte-for-byte (5/7 files per control).
External-CUBE MPI2 test/production controls match 7 files. Resource counts return
zero after failure/close; caller MPI remains alive. Both suites: 679 tests,one
optional skip,no failures. Production/test source and object prohibition audits
passed. No GS computation/collection or CI dispatch.

LD pole, PML, observation/auxiliary allocations and OS/process/MPI failures are
not covered. Real physical-memory exhaustion was not induced. Further R3 work
must be bounded to actual supported paths, followed by R4 and then group3.
