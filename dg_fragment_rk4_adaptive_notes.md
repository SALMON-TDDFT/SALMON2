# SALMON DG-Fragment RT: RK4・基底更新・検証ノート

## 1. 目的
本ノートは、DG-Fragment RT 実装の以下を整理する。

- 理論的な位置づけ（特に RK4 段内の自己無撞着更新）
- 実装上の主要変更点
- 実行手順（GS/RT、2frag/4frag 比較）
- 得られた結果の解釈
- 実運用時の推奨フロー

対象コード:
- `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/rt_dg_fragment.f90`

---

## 2. 理論的ポイント

### 2.1 RK4 の段内自己無撞着更新
論文準拠の古典 RK4 では、各段で状態を更新し、その段の密度からハミルトニアンを再構築する。

- `k1`: `t`
- `k2`: `t+dt/2`
- `k3`: `t+dt/2`
- `k4`: `t+dt`

重要点:
- 段ごとに `coef -> rho -> (Vh, Vxc) -> H` を再計算しないと、自己無撞着 RT と一致しない。
- 係数更新だけ行って `H` を固定すると、論文の時間発展と異なる。

### 2.2 観測量の正規化
本実装の電流は最終的に `system%ngrid` で割る「体積平均電流密度」。

- 出力 `J` が小さいのは、真空を含む大きいセル平均のため。
- 本質比較には、割る前の `current_sum` も併記する。

### 2.3 ゼロ場ドリフト切り分け
`I_wcm2_1=0` でも電流が残るなら、外場応答ではなく初期不整合（GS精度/GS→RT受け渡し整合）起源の可能性が高い。

推奨評価:
- `ΔJ(t) = J_strong(t) - J_zero(t)`
- `ΔE(t) = E_strong(t) - E_zero(t)`

---

## 3. 実装変更（要点）

### 3.1 RK4 を古典形で明示実装
`time_integrator=3` 時:
- `k1..k4` を古典 RK4 の段で計算
- 最終合成 `y_{n+1} = y0 + dt*(k1+2k2+2k3+k4)/6`

### 3.2 RK4 段内での `rho->H` 再構築を追加
`yn_fix_func='n'` のとき、RK4 各段で段時刻の `A(t_stage)` を使って再構築。

### 3.3 RK4 最終状態の時刻整合
RK4 最終係数作成後、`A(t+dt)` で最終 `rho->H` を再構築し、次ステップ初期条件の整合を向上。

### 3.4 RK4 と基底更新トリガの整合
RK4 では外側更新を抑制していたため、基底更新判定が走らない経路があった。
修正:
- `yn_adaptive_basis='y'` の場合は RK4 でも外側 `update_density_and_hamiltonian` を実行し、トリガ判定を有効化。

### 3.5 観測量の補正
- `nocc` を closed-shell 前提で `nelec/2` に修正。
- 観測量集計時の MPI 平均化と `ngrid` 正規化を明確化。
- IFLOW（界面電流保存）診断を追加（ペア反対称性/総電荷収支チェック）。

---

## 4. 主要検証結果（今回）

### 4.1 厳格 GS 収束
GS 入力を
- `nscf = 1000`
- `threshold = 1.0d-12`
に設定して再計算。

結果:
- 2frag/4frag とも閾値水準まで収束。

### 4.2 ゼロ場 vs 強場（厳格GS後）
`dt=0.0125, nt=40, t=0.5`。

- 強場とゼロ場の差分 `ΔJ` は非常に小さい。
- 現状の主成分はゼロ場ドリフト側。

### 4.3 基底更新強制発火テスト
設定:
- `yn_adaptive_basis='y'`
- `yn_adaptive_basis_dg='y'`
- `basis_update_threshold=1d-14`
- `yn_fix_func='n'`

結果:
- 2frag/4frag とも step 1..40 で基底更新トリガ発火。
- `NaN/STOP/SIGSEGV` なし。
- step40 の `J, E` は更新なしと同程度（破綻なし）。

---

## 5. 実行手順（推奨）

### 5.1 GS（高精度）
1. `inputfile_gs` を作成
2. `&scf` を以下に設定
   - `nscf = 1000`
   - `threshold = 1.0d-12`
3. 実行
   - 2frag: `mpirun -np 2 salmon < inputfile_gs_strict > gs_strict.log`
   - 4frag: `mpirun -np 4 salmon < inputfile_gs_strict > gs_strict.log`

### 5.2 RT 比較
1. 強場入力とゼロ場入力を用意
   - 強場: `I_wcm2_1 = 1.0d12`（例）
   - ゼロ場: `I_wcm2_1 = 0.0d0`
2. 同一 `dt, nt` で両方実行
3. `step40` など同時刻で比較
   - `current_sum`（割る前）
   - 出力 `J`（割った後）
   - `E`
4. 最終評価は `ΔJ, ΔE`。

### 5.3 基底更新有効化テスト
`&dc` / `&dg_fragment` / `&propagation` で以下を確認:
- `yn_adaptive_basis = 'y'`
- `yn_adaptive_basis_dg = 'y'`
- `basis_update_threshold = 1d-14`（発火確認用）
- `yn_fix_func = 'n'`

ログ確認文字列:
- `Adaptive basis update triggered`
- `Global ||ΔH||`
- `NaN|STOP|SIGSEGV`

---

## 6. 運用上の注意

- 低閾値（`1d-14`）は「発火確認用」。本番では更新過多になる。
- 本番閾値は `||ΔH||` 時系列を見て調整（例: `1d-4`〜`1d-2` a.u. を探索）。
- 比較は必ず 2frag/4frag で
  - 同一 GS 精度
  - 同一 `dt, nt`
  - 同一基底条件（cap/PW/状態数方針）
 で揃える。

---

## 7. まとめ
- RK4 段内自己無撞着更新は実装済み。
- 基底更新トリガ経路（RK4時）も有効化済み。
- 低閾値発火テストで破綻なし、値は更新なしと同程度。
- 現状差の主因は外場応答よりゼロ場ドリフト寄り。評価は `ΔJ, ΔE` が妥当。
