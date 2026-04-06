# Tau Mixing Design

**Problem**

`r2scan` を含む meta-GGA SCF で、primitive Si の厳しい条件に 2-cycle plateau が出る。現状の `vtau` mixing は SCF trajectory を少し変えるが、Wisteria の primitive exact case では plateau と `SIGXCPU` を解消しなかった。

**Observation**

- `rho` は [density_matrix.f90](/Users/mizukitani/Documents/DFT/SALMON2/.worktrees/stress-tensor/src/common/density_matrix.f90) の `calc_density` で計算される。
- `tau` は persistent field ではなく、[salmon_xc.f90](/Users/mizukitani/Documents/DFT/SALMON2/.worktrees/stress-tensor/src/xc/salmon_xc.f90) の内部 `calc_tau` で局所的に作られて捨てられる。
- `vtau` は `tau` と `rho` から XC evaluator が毎 iteration 作る output-side quantity である。
- Broyden 実装 [broyden.f90](/Users/mizukitani/Documents/DFT/SALMON2/.worktrees/stress-tensor/src/common/broyden.f90) は generic な実ベクトルに対して plain inner product を使う。

**Decision**

`vtau` は混ぜない。`tau` を input-side auxiliary field として mixing に組み込み、mixed `tau` から `vtau` を毎 iteration fresh に再計算する。

**User-Facing Controls**

- `yn_tau_mixing = 'y' | 'n'`
- `tau_mixrate`
- `tau_metric_weight`

`yn_tau_mixing='n'` を default にする。`tau` を使わない XC では追加配列と追加演算を発生させない。

**Chosen Strategy**

採用案は「simple 系では `tau` の独立 linear mixing、broyden/pulay では `rho` block と `tau` block の composite vector mixing」である。

- `simple` / `simple_potential`
  `tau_mixed = (1 - tau_mixrate) * tau_old + tau_mixrate * tau_raw`
- `broyden` / `pulay`
  mixing vector を
  `x = [rho(r_1..N), sqrt(tau_metric_weight) * tau(r_1..N)]`
  として pack し、unpack 時に `tau` block を元の scale に戻す

これは `rho + tau` を物理的に加算する意味ではなく、mixing solver に渡す state vector を block-wise に連結する意味である。

**SCF Flow**

`solve_orbitals -> calc_density(rho) -> calc_tau_raw(tau) -> mix(rho, tau) -> exchange_correlation(rho_mixed, tau_mixed) -> vtau fresh -> Hamiltonian`

**Code-Level Changes**

1. `tau` calculation helper を XC 内部 subroutine から外へ出し、SCF driver が呼べるようにする
2. `s_mixing` に `tau` history/state を追加する
3. `simple`, `simple_potential`, `broyden`, `pulay` の各 mixing path に `tau` を組み込む
4. `exchange_correlation` に `tau_override` を追加し、present のときは raw `tau` を再計算しない
5. 既存の `vtau` mixing helper は削除する
6. `checkpoint_restart` と `variables.log` に `tau` mixing state / input を反映する

**Why Not Other Options**

- `tau` simple mixing only:
  最低限の改善は見込めるが、QE/VASP の設計意図から外れる。`tau_metric_weight` も意味を持たない。
- generic auxiliary-field framework:
  将来拡張性は高いが、今回の問題に対しては変更面積が広すぎる。

**Risks**

- `tau` history の memory cost が増える
- `simple_potential` では `Vh/Vxc` potential mixing と `tau` density-like mixing が混在する
- `TB-mBJ` など `tau` を使う他の XC との整合を壊さないように、適用条件を `tau`-using XC 全体で設計する必要がある

**Expected Outcome**

- non-meta-GGA には影響しない
- `r2scan` primitive exact case の plateau onset が遅れる、あるいは解消する
- `vtau` mixing のような output-side heuristic ではなく、input-side stabilization に移行できる
