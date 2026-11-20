# R0 scalar IEEE compliance

Explicit all-element rank1–4 FDTD accumulators replace array IEEE inquiries.
History g/b/coeff use scalar loops even when metadata already rejects. Keep
precision/expressions/thresholds/collective rejection/stream format unchanged.
Existing IO native fixture extracts real routines and tests finite values,
NaN/+Inf/-Inf at every position and invalid metadata. Existing MPI E/H probe
adds rank-local J/A nonfinite rejection with byte-exact state/output retention.
Compare normal field arrays and HSE current/energy/history with prior controls.
No GS copies/SCF/new real electronic checkpoint collection or routine CI.
Fujitsu compiler is not available; do not claim its compilation.

## Acceptance

Source/wheel native finite/history probes pass. Field MPI2 nonfinite J/A
collective rejection preserves state/output bytes; overflow recovery remains.
Normal source/wheel field A/E/H traces exactly match prior accepted traces.
HSE heterogeneous MPI current/energy/history and memory rollback exact.
Source/installed638:637pass1skip, audits pass. Fujitsu is unavailable/untested.
JP/EN/AI/list synchronized; no GS copies/SCF/electronic disk checkpoints or CI.
Evidence: calc2026/salmonpy-r0-scalar-record-20261120; uncommitted.
