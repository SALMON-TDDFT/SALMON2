# EXX memory reduction evidence

Measured production commit: `429f9b4f`.

- `results.json`: five final runs with executable/input/seed provenance and historical source-ACE reference rows.
- `summary.json`: three sizes, maximum-rank lifetime RSS and 16-step RT times. Eight H2 has three repeats; 64 and 128 H2 have one each.
- `raw.tar.gz`: final inputs, outputs, timing/RSS records and measurement script. Restart wavefunction files are not bundled; their hashes are recorded in the results.
- `validation.tar.gz`: build, RED/GREEN regression and review-fix logs, final 17-test Ehrenfest result, and completed progress ledger. RED logs intentionally contain failures from before the respective fixes.
- `memory.png` / `memory.svg`: summary figure.
- `manifest.json`: SHA-256 of every other file in this directory.

Before-memory source-ACE measurements are reused from the earlier weak-scaling campaign. Full calculations were not rerun. Pre-review pilot measurements (explicit metric inverse) are excluded. Timing comparisons are historical on a shared machine, especially MPI16; they do not establish a definitive speedup.
