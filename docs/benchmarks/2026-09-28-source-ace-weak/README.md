# Source-support ACE weak scaling

All eight weak-scaling shapes are complete:16 historical Full runs and16 source-ACE runs.
`results.json` is a weak-only projection; `source_matrix_complete=false` records that
strong scaling is still pending. Full references were explicitly reused at the user's
request after exact input/seed/output verification. Each historical Full row retains
its original executable hash, commit and source-results hash. Current source-ACE rows
use the top-level binary hash unless a per-row hash is provided.

`summary.json` also compares the previous occupied-vector .999 results from
`../2026-09-28-h2-blocks-before-pairs/results.json`. All three timing series are
historical comparisons across measurement times; Full/previous use the older binary.
The numerical inputs and seed payload hashes match; the old and new source-ACE
executable hashes differ deliberately. Do not describe speedups as same-binary results.

`raw.tar.gz` contains all32 completed inputs/outputs/time series/rank RSS records;
interrupted strong work is excluded. Large GS binary seeds are retained locally and
identified by hashes in results.json.32 H2 and smaller use minimum time from3 repeats,
64/128 H2 use one run; RSS is median of max-rank lifetime peaks for repeated cases.
`weak-scaling.png`/`.svg` plot only the 1D shape sequence. Full8-shape tables and
interpretation are in `../../h2-block-source-ace-weak-ja.md`.
