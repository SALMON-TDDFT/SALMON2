# H2 PBEh impulse RT scaling

> 2026-10-06：以下の`benchmarks`は当時のローカル性能測定です。マージ対象から除外し、測定コードはブランチ外へ保存しています。

Continue the requested H2 weak and strong matrix, now impulse16steps rather
than GS. Keep geometry, grid, Coulomb cutoff, full EXX support, MLWF controls,
MPI layouts and three repeats. Fixed ions; x impulse1e-4au, dt0.02au.
Prioritize minimum times due user-reported contention; memory uses median
max-rank peak RSS, retaining every raw value.

GS wavefunctions were not exported by the previous benchmark and native PBEh
RT reads DC-LCFO seeds. Prepare each shape once with a full-cell single fragment,
no buffer, occupied-only, zero-temperature GS. Native RT reconstructs the mesh
and does not propagate in a retained LCFO basis. Share identical seeds across
strong layouts. Hash and archive payloads. Exclude GS preparation from RT timers
and process RSS. Distinguish RT-loop time, propagation sub-timer, whole-job time,
and whole-RT-process highwater memory (including loading/reconstruction).

Pilot2x1x1 atMPI2 passed full-cell export and16-step mesh propagation. Execute
all8 preparations and36 RT jobs sequentially. Validate16steps, finite data,
post-kick energy width, cross-rank currents/energies, and seed immutability.
Archive complete reproducibility data and report weak/strong performance.

Completed:8 full-cell preparations converged; all36 RT jobs completed16steps.
Strong current/energy max differences1.867e-20au/1.101e-12Ha; maximum
post-kick energy width2.700e-12Ha. Best observed16-step loop time atMPI8
3.1539s versusMPI1 14.412s. Full logs and reusable seed payloads archived
in docs/benchmarks/2026-09-28-h2-rt. No production modifications.
