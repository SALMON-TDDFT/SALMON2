# H2 supercell PBEh scaling measurement

> 2026-10-06：以下の`benchmarks`は当時のローカル性能測定です。マージ対象から除外し、測定コードはブランチ外へ保存しています。

User-requested weak shapes: 1x1x1, 2x1x1, 4x1x1, 8x1x1, 16x1x1,
2x2x1, 4x4x1, 2x2x2. One MPI rank per molecule. Strong: fixed 4x4x1,
MPI1/2/4/8/16. Three fresh-process repetitions, common MPI16 case shared.

Run real converged native PBEh40 SCF with 8-bohr cubic primitive cells,
1.4-bohr x-directed H2, 0.5-bohr grid spacing, Gamma, occupied-only,
4-bohr fixed Coulomb radius, full EXX support, MLWF 5/100/1e-7 and rho_dne
SCF threshold1e-10. y/z-only spatial decomposition keeps 4096 owned grid
points per rank in weak cases; full occupied columns still grow. Document
serial/distributed algorithm change, finite-cell/k-sampling effects and
SCF iteration counts, and avoid claims of linear complexity or giant-water
performance. Production implementation stays unchanged.

Pilot the smallest and largest jobs, then execute sequentially to avoid mutual
interference. Capture native max-rank SCF timing, root total, launcher wall time,
child per-rank lifetime RSS, energy, iteration count and localization diagnostics.
Archive inputs/logs and metadata; check strong-case energies before reporting.

Pilot finding: default Broyden alpha_mb=0.75 failed to converge in 500 steps
for 16x1x1. Changing only initial density to pp also failed. Changing only
alpha_mb to0.1 converged in264 steps. Use alpha_mb=0.1 uniformly and rerun
all measured cases; preserve failed/default runs as diagnostics, not scaling data.

Repeat showed alpha0.1 alone is not robust (500-step failure). gauss10 with
alpha0.1 converged in330 steps. Final uniform settings: gauss10, alpha0.1,
SCF maximum2000. Preserve all initialization diagnostics; do not select only
successful timings. Repeat the entire matrix under the final settings.

User steering during execution: other applications contend with MPI16; trust
shorter timings more. Primary scaling metric is minimum observed per-run
SCF time/iteration, with minimum total SCF time secondary. Preserve medians
and full ranges; avoid combining unrelated per-metric minima as one trial.

Completed: all36 final uniform-condition runs converged below1e-10; strong
final energies agree within1.596e-8eV. Inputs, outputs, source snapshots,
metadata, initialization diagnostics, and plots are archived under
`docs/benchmarks/2026-09-28-h2`; report: `docs/h2-pbeh-scaling-ja.md`.
