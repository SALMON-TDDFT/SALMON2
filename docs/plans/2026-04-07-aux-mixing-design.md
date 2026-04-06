# Auxiliary-Field Mixing Design

**Problem**

`r2scan` のような meta-GGA SCF では、`tau -> vtau -> Hamiltonian -> orbitals -> tau` の feedback loop が強く、primitive Si の厳しい条件で 2-cycle plateau が出る。`vtau` mixing は output-side heuristic としては効きが弱く、Wisteria の primitive exact case でも plateau と `SIGXCPU` を解消できなかった。

`TB-mBJ` では `tau` だけでなく `j` も evaluator に入るため、`tau` だけを特別扱いする設計では不足する。

**Observation**

- `tau` と `j` は [salmon_xc.f90](/Users/mizukitani/Documents/DFT/SALMON2/.worktrees/stress-tensor/src/xc/salmon_xc.f90) の [calc_tau_from_orbitals](/Users/mizukitani/Documents/DFT/SALMON2/.worktrees/stress-tensor/src/xc/salmon_xc.f90#L616) で軌道から毎 iteration 再計算される。
- `r2scan` は `tau` を使うが、`j` は使わない。[salmon_xc.f90](/Users/mizukitani/Documents/DFT/SALMON2/.worktrees/stress-tensor/src/xc/salmon_xc.f90#L1498)
- `TB-mBJ` は `tau - |j|^2/(2rho)` を使うため、`j` が必要である。[builtin_tbmbj.f90](/Users/mizukitani/Documents/DFT/SALMON2/.worktrees/stress-tensor/src/xc/builtin_tbmbj.f90#L57)
- 既存の Broyden/Pulay solver は generic な実ベクトルに対して plain inner product を使うため、wrapper 側の block packing で auxiliary field を組み込める。

**Decision**

`vtau` や `TB-mBJ` の有効ポテンシャルは混ぜない。`tau` と必要なら `j` を input-side auxiliary field として mixing state に取り込み、mixed auxiliary fields から毎 iteration fresh に XC operator を再計算する。

**User-Facing Controls**

- `yn_aux_mixing = 'y' | 'n'`
- `tau_mixrate`
- `tau_metric_weight`
- `j_mixrate`
- `j_metric_weight`

`yn_aux_mixing='n'` を default にする。`tau_mixrate=0` は `tau` block を無効化し、`j_mixrate=0` は `j` block を無効化する。`LDA/GGA` や `yn_aux_mixing='n'` の経路では、配列確保・演算・アルゴリズムを従来どおりに保つ。

**Chosen Strategy**

採用案は「SCF auxiliary state を導入し、Broyden/Pulay では `rho` と auxiliary block を一つの history に載せる」方式である。

- `simple` / `simple_potential`
  - `rho` は従来の mixing
  - `tau` と `j` は block ごとに独立の linear mixing
- `broyden` / `pulay`
  - 有効 block だけを pack した composite vector を使う
  - 例:
    `x = [rho, sqrt(w_tau) * tau, sqrt(w_j) * jx, sqrt(w_j) * jy, sqrt(w_j) * jz]`
  - `tau_mixrate` と `j_mixrate` は、各 block の residual/update damping として扱う

これは物理量を加算する意味ではなく、solver に渡す state vector を block-wise に連結する意味である。

**SCF Flow**

`solve_orbitals -> calc_density(rho_raw) -> calc_aux_fields(tau_raw, j_raw) -> mix(rho_raw, aux_raw) -> exchange_correlation(rho_mixed, aux_mixed) -> fresh operator payload -> Hamiltonian`

- `r2scan/libxc mGGA`
  - `tau` block のみ使用
- `TB-mBJ`
  - `tau` と `j` block を自動で使用

**Memory / Lifetime**

- 大配列の `allocate/deallocate` を SCF loop 内で繰り返さない
- [s_mixing](/Users/mizukitani/Documents/DFT/SALMON2/.worktrees/stress-tensor/src/common/structures.f90) に auxiliary history と persistent scratch を持たせる
- `j` は新規の大きい 4D work 配列を毎回作るのではなく、`jx/jy/jz` の 3 block として扱う
- `yn_aux_mixing='n'` なら auxiliary history/scratch は確保しない
- `tau_mixrate=0` / `j_mixrate=0` の block は確保しない

**Restart**

restart 時に auxiliary history は復元しない。restart 後は history 未成熟状態から再構築する。raw `rho/tau/j` と mixed `rho/tau/j` の扱いは通常の初期 SCF と同じにする。

**Code-Level Changes**

1. `tau` / `j` を計算する reusable helper を XC evaluator から切り出し、SCF driver が呼べるようにする
2. [structures.f90](/Users/mizukitani/Documents/DFT/SALMON2/.worktrees/stress-tensor/src/common/structures.f90) の `s_mixing` を auxiliary-state aware に拡張する
3. [mixing.f90](/Users/mizukitani/Documents/DFT/SALMON2/.worktrees/stress-tensor/src/gs/mixing.f90) の wrapper を block-aware にし、`rho/tau/j` の pack/unpack と block metric / damping を担わせる
4. [scf_iteration.f90](/Users/mizukitani/Documents/DFT/SALMON2/.worktrees/stress-tensor/src/gs/scf_iteration.f90) で raw auxiliary field 計算と mixing orchestration を行う
5. [salmon_xc.f90](/Users/mizukitani/Documents/DFT/SALMON2/.worktrees/stress-tensor/src/xc/salmon_xc.f90) の `exchange_correlation` に mixed auxiliary field override を渡し、fresh payload を生成する
6. [inputoutput.f90](/Users/mizukitani/Documents/DFT/SALMON2/.worktrees/stress-tensor/src/io/inputoutput.f90) と `variables.log` 出力を `aux mixing` に更新する
7. 既存の `vtau` mixing helper は削除し、`TB-mBJ` でも `tau+j` override を扱えるようにする

**Why Not Other Options**

- `tau` だけを混ぜる:
  `TB-mBJ` の `j` が残り、設計が再び二系統になる
- `TB-mBJ` だけ stop:
  当面は安全だが、tau-using XC を統一的に扱う設計にならない
- solver 本体を block-aware に全面改造:
  将来はあり得るが、現段階では wrapper 側で十分対応できる

**Risks**

- auxiliary history による memory 使用量増加
- `simple_potential` では potential mixing と auxiliary-field mixing が混在する
- `TB-mBJ` の `j` block は `tau` よりノイズが大きい可能性があり、`j_metric_weight` と `j_mixrate` の default を保守的にする必要がある
- `yn_aux_mixing='n'` の既存経路を厳密に不変に保つ必要がある

**Expected Outcome**

- `LDA/GGA` や非 aux 系には影響しない
- `r2scan` primitive exact case の plateau onset がさらに遅れる、あるいは解消する
- `TB-mBJ` も同じ auxiliary-field mixing framework に載る
- output-side heuristic ではなく、QE/VASP に近い input-side stabilization へ移行できる
