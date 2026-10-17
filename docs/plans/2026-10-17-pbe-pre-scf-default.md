# Automatic PBE pre-SCF default

Approved design: for eligible fresh static hybrid GS, omit exx_pre_scf_threshold to obtain 1e-4 and three consecutive density residual checks. Explicit zero disables. RT, restart, MD/optimization, unsupported density metrics/mixing, checkpoints and diagnostic snapshots do not auto-enable. Existing explicit-positive incompatibility errors remain.

1. Red test: bounded H4 default GS must log a PBE stage, transition, and converged target hybrid.
2. Use an internal omitted-value sentinel, resolving after all input broadcasts/normalization and before dump/check. Do not change the PBE/hybrid stage algorithm.
3. Validate eligible GS, explicit zero/positive, nonhybrid, RT, restart, DC temperature, unsupported modes and diagnostics; rebuild serial/MPI and run existing regressions. Keep restart probes outside testsuites.
4. Update current input documentation and record the exact tests. Preserve unrelated uncommitted global-hybrid symmetry work.
