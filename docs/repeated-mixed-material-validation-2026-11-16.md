# Repeated materials in native mixed layers — 2026-11-16

Code:01d43cc1. Installed macOS arm64 cp314 wheel SHA256:e15e20706dc353d3be91e5fa990c778ea9b8519848974236817dbbe6e662e2f7. No kernel/code changes for this check. Apple Silicon, local host MPI, OMP/BLAS1, at most four MPI ranks. No remote Actions used.

Starting with the installed mixed_layers.json fixture, set exactly one layer thickness to4 Bohr (grid spacing2). Pass `--tddft-ranks 1 --rank-budget 4` and launch four ranks. Other parameters remain128 steps,dt.003,40³ field grid and the same Si starter PZ GS. Following layers move one cell to preserve the stated order/gap.

| Repeated layer | Topology | Max replay errors for all three electronic points | Difference between repeated points' A/c | Difference between repeated points' j |
|---|---|---|---|---|
| Si | field1 + SBE1 + TDDFT1 + TDDFT1 | 0,0,0 | 1.7304869559946396e-05 | 4.813352538532136e-07 |
| SBE | field1 + SBE1 + SBE1 + TDDFT1 | 0,0,0 | 1.7407476425437524e-05 | 5.440425116013658e-07 |

Both native LD materials remain on the field rank. Checked that the repeated pair shares material ID, has different point IDs and disjoint MPI groups, receives different histories and produces different currents. The installed sample replays every recorded A/c history in a freshly created serial state; all currents match exactly in these one-rank electronic cases. It also checks previous-current feedback and native field/electronic clocks.

Persistent evidence beneath `~/.codex/outputs/salmonpy-state-inputs-20261116`: `layers-repeated-{si,sbe}.json`, `layers-repeated-{si,sbe}.log`, `native-layers-repeated-{si,sbe}/result.json` and `response.npz`. Those paths are developer evidence locations, not user prerequisites. The Japanese/English COUPLED guides include a portable JSON-edit recipe.

Scope: independent material states in separate MPI groups/processes. Not general in-process native TDDFT reentrancy, physical convergence, whole-session restart or memory-cap enforcement. No new claim about two field ranks plus two internally parallel TDDFT groups, which exceeds the four-core budget. Molecular work stays paused.
