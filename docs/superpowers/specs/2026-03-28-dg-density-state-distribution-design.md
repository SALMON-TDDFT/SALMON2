# DG 密度行列再構成：状態分散並列化 設計書

Date: 2026-03-28
Target file: `src/rt/dg/rt_dg_density_reconstruct.f90`

---

## 背景と問題

`calculate_density_from_fragments` の `project` ステージが全体時間の 74% を占め、そのうち 91% が timing breakdown に計上されていない。

### 未計時コードの特定

通常の RT ステップ（`density_phi_cache_valid = .true.`）において、
`t_project0`〜`t_project1` 区間内で以下が計時されていない：

1. **`zgemm(D)` 呼び出し** (line 559): `D = C^H C`、フラグメントあたり 1 回
2. **`comm_bcast(D)`** (line 563): `density_matrix_frag` を `icomm_frag` 全ランクに配信、フラグメントあたり 1 回
3. **n_pw>0 パスの `comm_bcast(coef_blk)`** (lines 615–616, 642): グリッドブロック × 状態バッチごとに実行（重大な通信冗長）

### 実測値（典型的実行）

| パラメータ | 値 |
|---|---|
| `nocc_spin` | ~2200 |
| `nbf_max` | ≤ 71 |
| `isize_frag` | 4 |
| `ifrag_count` | ~68 フラグメント/ランク |
| `n_pw` | 0（通常ステップ）/ >0（PW 混在時） |

`zgemm` が小サイズ GEMM（71×71×2200 complex）で BLAS 効率が低く（ピークの 10–20%）、
フラグメントあたり ~10–30 ms × 68 フラグメント ≈ **~2 s/ステップ** が未計時で消費されている。

n_pw>0 での通信量試算（現状）:
```
(nblocks/4) × (nocc/64) × (nbf×64×16 + n_pw×64×16) bytes ≈ 数 GB/ステップ/フラグメント
```

---

## 設計方針

### 核心：グリッド分散 → 状態分散 への切り替え

| | 現状（グリッド分散） | 本設計（状態分散） |
|--|---|---|
| 各ランクが担当 | グリッドブロックの 1/N | 状態の 1/N |
| `coef` 配布 | ブロック×バッチごとに bcast | フラグメントあたり **1 回** bcast |
| rho 収集 | `send_pack` → `comm_summation` → slot mapping | `Reduce(partial_rho)` to root |
| `isize_frag=1` 時 | 非分散 | bcast/Reduce がノーオプ → 完全等価 |

状態範囲の計算:
```fortran
nocc_per_rank = (nocc_spin + isize_frag - 1) / isize_frag
io_s_frag = dg_frag%id_frag * nocc_per_rank + 1
io_e_frag = min((dg_frag%id_frag + 1) * nocc_per_rank, nocc_spin)
```

---

## n_pw=0 パス（密度行列経路）

### 現状フロー

```
[ブロックループ内、初回のみ]
  root:  coef_occ_weighted ← gather from coef  (未計時)
  root:  zgemm('N','C', nbf,nbf,nocc) → D_complex  (未計時、8× 冗長)
  全ランク: comm_bcast(D_complex, icomm_frag)  (未計時)
[ブロックループ内、毎回]
  全ランク: density_mat_re ← real(D_complex)
  全ランク: dgemm(phi_blk, density_mat_re) → density_tmp
  全ランク: rho_blk = diag(phi_blk × density_tmp)
  send_pack / comm_summation → rho_send  (既存グリッド分散機構)
```

### 変更後フロー

**D 構築フェーズのみ変更。ブロックループは既存のグリッド分散を維持する。**

```
[ブロックループ前、フラグメントあたり 1 回（NEW）]
  root:    coef_re(1:nbf, 1:nocc_spin) ← gather from coef  (実数のみ)
  bcast    coef_re → icomm_frag 全ランク  [通信量: nbf×nocc×8 ≈ 1.2 MB/frag]
  各ランク: dsyrk('U','N', nbf, nocc_per_rank, coef_re(:, io_s:io_e))
           → D_partial_re (上三角のみ, real(8))
  AllReduce D_partial_re → D_re over icomm_frag  [80 KB/frag]

[ブロックループ内、既存のグリッド分散を維持（変更なし）]
  各ランク: density_mat_re ← D_re (新規 real バッファから、変換不要)
  各ランク: dgemm(phi_blk, D_re) → density_tmp
            ※ dsymm への切り替えは対称性の正しい利用が確認できた後に行う
  各ランク: rho_blk = diag(phi_blk × density_tmp)
  send_pack / comm_summation → rho_send  (変更なし)
```

### 速度改善

| ステップ | 現状 | 変更後 | 改善 |
|---|---|---|---|
| D 構築 | `zgemm` complex (8×) | `dsyrk` real + 対称 | **8× 高速化** |
| 状態分散 | root のみ | isize_frag 分割 | **isize_frag× 高速化** |
| bcast(D) サイズ | nbf²×16 bytes (complex) | nbf²×8 bytes (real) | **2× 削減** |
| ブロックループ | 変更なし | 変更なし | — |

D 構築部分の isize_frag=4 での理論値: **32× 高速化**

---

## n_pw>0 パス（状態ループ経路）

### 現状フロー（重大な問題）

```
[ブロックループ内]
  do block_offset = ...
    do io0 = 1, nocc_spin, state_block_size   ← 状態ループがブロック内
      root: fill coef_blk_re/im
      comm_bcast(coef_blk_re, icomm_frag)    ← ブロック×バッチごとに通信！
      comm_bcast(coef_blk_im, icomm_frag)
      dgemm + phase_cache
      do ipw0 = 1, n_pw, pw_block_size
        comm_bcast(coef_pw_blk, icomm_frag)  ← さらに PW 分も通信！
      end do
      rho_blk_accum += |psi|²
    end do
  end do
```

### 変更後フロー

```
[ブロックループ前、フラグメントあたり 1 回]
  root:   coef_re/im(1:nbf, 1:nocc_spin) ← gather
  root:   coef_pw ← from coef_pw_full_cache (1:n_pw, 1:nocc_spin)
  bcast   coef_re/im + coef_pw → icomm_frag 全ランク
          [通信量: (nbf×2 + n_pw×2) × nocc × 8 bytes / frag]

[ブロックループ内]
  do block_offset = ...
    do io0 = io_s_frag, io_e_frag, state_block_size  ← 自分のスライスのみ
      psi_blk = phi_blk × coef_blk + phase_cache × coef_pw_blk
                                     ← bcast ゼロ
      rho_blk_partial += occ_factor × |psi_blk|²
    end do
  end do

[フラグメント終了後]
  Reduce(rho_blk_partial → rho_local) to root over icomm_frag
  [通信量: ngrid_frag × 8 bytes / frag ≈ 40 KB/frag]
```

### 通信量比較（フラグメントあたり、nblocks=100, n_pw=50 の場合）

| | 現状 | 変更後 |
|--|---|---|
| n_pw=0 | bcast(D) 80 KB × 68 + zgemm 未分散 | bcast(coef_re) 85 MB + AllReduce(D) 5 MB |
| n_pw>0 | ~8 GB | bcast(coef) ~300 MB + Reduce(rho) ~3 MB |

---

## `send_pack` 機構の削除

現在の `distribute_project` パスは以下を使用：
- `density_subgroup_send_slot_map`
- `density_subgroup_send_count`
- `density_subgroup_self_ixg/iyg/izg`
- `send_pack / send_sum`
- `subgroup_send_offset`

これらはすべてグリッド分散のためのスロットルーティング。
状態分散では `partial_rho` の Reduce で置き換えられるため、**これらの配列と関連ロジックを削除する。**

削除対象フィールド（`rt_dg_fragment_types.f90` 内の `s_dg_fragment_rt`）:
- `density_subgroup_send_count`
- `density_subgroup_send_slot_map`
- `density_subgroup_self_ixg/iyg/izg`
- `density_block_first_offset`
- `density_block_step`
- `density_block_nblocks`

---

## `isize_frag=1` との後方互換性

```
distribute_project = (icomm_frag /= COMM_GROUP_NULL .and. isize_frag > 1)
```

`isize_frag=1` のとき：
- `io_s_frag = 1`, `io_e_frag = nocc_spin`（全状態を担当）
- bcast はノーオプ（自分自身へ）
- AllReduce/Reduce もノーオプ（1 ランク）
- コードパスは `distribute_project = .false.` と同じ挙動

---

## 実装スコープ

### 変更が必要なファイル

| ファイル | 変更内容 |
|---|---|
| `src/rt/dg/rt_dg_density_reconstruct.f90` | メインロジックの変更（大） |
| `src/rt/dg/rt_dg_fragment_types.f90` | `s_dg_fragment_rt` フィールドの削除 |
| `src/rt/dg/rt_dg_fragment.f90` | 削除フィールドの allocate/deallocate 除去 |
| `src/rt/dg/rt_dg_fragment_soi.f90` | 同上 |

### 変更しないファイル

- `rt_dg_fragment_parallel.f90`（`icomm_frag` の設定は変更なし）
- fragment 間の `rho_send/rho_recv` 通信（`comm` ステージ）は変更なし

---

## タイミング計測の追加

変更と同時に以下の timing 変数を追加：

| 変数 | 計測区間 |
|---|---|
| `time_project_coef_bcast` | coef の upfront bcast |
| `time_project_dmat_build` | dsyrk + AllReduce(D) |
| `time_project_rho_reduce` | Reduce(partial_rho) |

breakdown 出力行に追加して `project` 合計との差が常にゼロになるよう確認する。

---

## 成功基準

1. `project` timing の breakdown 合計が `project` 合計と一致する（未計時ゼロ）
2. `project` 時間が現状の 2.15 s から 0.2 s 以下に削減される（isize_frag=4 時）
3. `isize_frag=1` の場合に数値的に同一の rho が得られる
4. n_pw>0 ケースで同様の高速化が確認される
