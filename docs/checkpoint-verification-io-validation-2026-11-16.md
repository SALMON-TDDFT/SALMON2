# Partitioned checkpoint verification — 2026-11-16

Direct restore previously streamed every payload for SHA256 on every outer rank. The internal mixed-example coordinator now partitions sorted payloads round-robin, hashes each file once in total, and gathers all errors BEFORE loading any numerical state. All ranks compare a canonical manifest digest before restoration. A manifest mismatch, bad file/path/hash/inventory or local metadata rejects the entire restore. Payloads must be the same files in one shared directory, not different per-node copies. Metadata/stat operations remain collective/local; only redundant payload checksum streaming is removed.

The sample reports per-rank files_checked and bytes_checked in checkpoint_validation. This is verification I/O instrumentation, not aggregate process memory measurement or a runtime speedup claim.

## Fresh evidence

Installed wheel,4 ranks(field1/SBE1/Libxc-PBE TDDFT orbital2),same authenticated XC library and converged Si fixture from xc-inner-mpi-direct-checkpoint-validation-2026-11-16.md. Restored checkpoint63 and propagated65 steps, with complete A/c,E,j bit-identical to the existing continuous128 trajectory. No native or physics source changed.

- Eight payloads totaling82,191,025 bytes were checksum-streamed once across the job, rather than328,764,100 bytes under fourfold verification. Reduction75% applies to checksum payload reads; native restore and ordinary metadata/frame reads are additional.
- Per-rank streamed bytes:14,353,101;14,353,195;53,436,736;47,993. Each rank verified two files. This is file-count partitioning, not optimal byte balancing.
- Eight direct-checkpoint tests cover prior admission, disjoint exhaustive hash partition, assigned corruption, path escape and inconsistent manifests before mutation. Additional review identified the manifest-agreement precondition, now enforced and tested.
- Final installed/source suites each450 tests,1 skipped,no failures. Sdist-to-wheel and guide/AI synchronization pass. Existing archives remain compatible. No Actions used.

Persistent local evidence under ~/.codex/outputs/salmonpy-state-inputs-20261116: checkpoint-io-final-dist, checkpoint-io-final-{package,installed-tests,source-tests}.log, checkpoint-io-final-pbe-resume/result.json,response.npz and corresponding.log. The previous candidate was retained separately before adding canonical-manifest agreement.

## HSE constraint

The project lead reaffirmed fast-PR HSE only. No conventional HSE or comparison baseline is added as an execution engine. Fast HSE direct persistence remains open: its dedicated propagators, ACE/factor histories and native resources require a separate implementation; the semilocal middlepoint whitelist is not expanded to pretend support.
