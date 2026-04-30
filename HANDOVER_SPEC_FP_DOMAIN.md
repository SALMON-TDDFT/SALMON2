# DG-Fragment RT-TDDFT 引き継ぎ仕様書
**作成日**: 2026-04-28  
**対象ブランチ**: `develop-2.0.0`（ローカル変更あり、未コミット）  
**担当**: GitHub Copilot セッション引き継ぎ用

---

## 1. リポジトリ / ワークスペース

| 項目 | パス |
|------|------|
| ルート | `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2` |
| 主要ソース (変更済) | `src/rt/dg/rt_dg_plane_wave.f90` |
| 参照ソース | `src/rt/dg/rt_dg_integrator_derivative.f90` |
| フラグメント型定義 | `src/rt/dg/rt_dg_fragment_types.f90` |

---

## 2. ビルド

### ビルドディレクトリ
```
<root>/build_manual_mpi_probe/
```

### ビルドコマンド
```bash
cd /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2
cmake --build build_manual_mpi_probe -j4
```

### ビルド設定
- コンパイラ: `mpif90` (gfortran-15 バックエンド, OpenMPI)
- MPI: 有効
- 実行ファイル: `build_manual_mpi_probe/salmon`

### cmake 初期設定（再構成が必要な場合）
```bash
mkdir -p build_manual_mpi_probe && cd build_manual_mpi_probe
FC=mpif90 cmake .. -DENABLE_MPI=ON
```

---

## 3. テストケース

### テストディレクトリ
```
samples/exercise_dg_rt_hse_test/H2/
```

### 主要インプットファイル

| ファイル名 | 用途 |
|-----------|------|
| `inputfile_h2_periodic_20_dg_rt_2frag` | **メイン**: H2 20分子 2フラグメント DG-RT（nt=200） |
| `inputfile_h2_periodic_20_dg_rt_2frag_nt40_smoke` | スモークテスト（nt=40, 高速確認用） |
| `inputfile_h2_periodic_20_gs_2frag` | GS計算（RT前提の基底状態） |
| `inputfile_h2_periodic_20_dg_rt_4frag` | 4フラグメント版 |

### GS データディレクトリ
```
samples/exercise_dg_rt_hse_test/H2/data_dcdft/
```
RTインプット中に `yn_dg_fragment_from_dcdft = 'y'` が指定されており、
ここのGSデータを読み込む（事前にGS計算済みを想定）。

### 擬ポテンシャル
```
samples/exercise_dg_rt_hse_test/H2/H_rps.dat
```
（または同ディレクトリの `PS_H_KY_n.dat`）

---

## 4. 実行方法

### 基本実行（MPI 2プロセス、stdin リダイレクト必須）
```bash
cd samples/exercise_dg_rt_hse_test/H2
mpirun -np 2 \
  /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/build_manual_mpi_probe/salmon \
  < inputfile_h2_periodic_20_dg_rt_2frag \
  2>&1 | tee run_2frag.log
```

> **注意**: inputfile を `argv` で渡すとバナーでハングすることがある。必ず `< inputfile` （stdin リダイレクト）を使うこと。

### スモークテスト（nt=40, 短時間）
```bash
mpirun -np 2 \
  ../../build_manual_mpi_probe/salmon \
  < inputfile_h2_periodic_20_dg_rt_2frag_nt40_smoke \
  2>&1 | tee smoke_nt40.log
```

---

## 5. 実装済み機能: S_fp / H_fp 積分域トグル

### 変更ファイル
`src/rt/dg/rt_dg_plane_wave.f90`

### 動作
環境変数 `SALMON_DG_FP_INTEGRAL_DOMAIN` で積分域を切り替える。

| 環境変数値 | 動作 |
|-----------|------|
| 未設定 / その他 | **core** ドメインで積分（デフォルト）|
| `b`, `B`, `p`, `P`, `f`, `F`, `1` | **buffered** ドメインで積分 |

### ログタグ
- `[FP-DOMAIN]` — 初回起動時のモード表示、`||S_fp||_F`、`||H_fp||_F` の Frobenius ノルム

### 使用例

**Core モード（デフォルト）:**
```bash
mpirun -np 2 ./build_manual_mpi_probe/salmon \
  < inputfile_h2_periodic_20_dg_rt_2frag_nt40_smoke \
  2>&1 | tee run_core.log
grep '\[FP-DOMAIN\]' run_core.log
```

**Buffered モード:**
```bash
SALMON_DG_FP_INTEGRAL_DOMAIN=buffered \
mpirun -np 2 ./build_manual_mpi_probe/salmon \
  < inputfile_h2_periodic_20_dg_rt_2frag_nt40_smoke \
  2>&1 | tee run_buffered.log
grep '\[FP-DOMAIN\]' run_buffered.log
```

### A/B 比較スクリプト例
```bash
#!/bin/bash
TESTDIR=samples/exercise_dg_rt_hse_test/H2
EXE=build_manual_mpi_probe/salmon
INPUT=$TESTDIR/inputfile_h2_periodic_20_dg_rt_2frag_nt40_smoke

# Core
mpirun -np 2 $EXE < $INPUT 2>&1 | tee /tmp/fp_core.log

# Buffered
SALMON_DG_FP_INTEGRAL_DOMAIN=buffered \
mpirun -np 2 $EXE < $INPUT 2>&1 | tee /tmp/fp_buffered.log

# 比較
echo "=== Core [FP-DOMAIN] ==="
grep '\[FP-DOMAIN\]' /tmp/fp_core.log
echo "=== Buffered [FP-DOMAIN] ==="
grep '\[FP-DOMAIN\]' /tmp/fp_buffered.log
```

---

## 6. 関連する既存環境変数

| 環境変数 | ファイル | 機能 |
|---------|---------|------|
| `SALMON_DG_FP_INTEGRAL_DOMAIN` | `rt_dg_plane_wave.f90` | S_fp/H_fp 積分域 core/buffered 切替 (**今回実装**) |
| `SALMON_DG_OP_MIX_TRACE` | `rt_dg_integrator_derivative.f90` | op-mix-norm(M), op-mix-norm(H0) ログ出力 |
| `SALMON_DG_M_BLOCK_AUDIT` | `rt_dg_integrator_derivative.f90` | `[M-BLOCK-AUDIT] raw_norms ff fp pf pp all` ログ出力 |

---

## 7. 未完了タスク / 次セッションへの申し送り

### A/B 実験（最優先）
`core` と `buffered` の両モードで同一ケース（nt=40 スモーク）を実行し、
`[FP-DOMAIN]` ログの `||S_fp||_F` / `||H_fp||_F` 数値を比較する。

期待されること:
- buffered で値が増加（積分域が広いため）
- GS エネルギー・電流値の違いを `[M-BLOCK-AUDIT]` と合わせて評価

### M_fp 積分域整合（次ステップ候補）
`rt_dg_integrator_derivative.f90` の L382 付近で `S_mat_frag_pw` を使う M_fp ブロック構築がある。
S_fp を buffered にした場合、M_fp 側も同域で計算するか検討・実装する必要がある可能性がある。

---

## 8. 関連ドキュメント

- `DG_BUFFERED_RT_SEED_CONTRACT_SPEC_JA.md` — バッファ領域設計仕様
- `ADAPTIVE_BASIS_IMPLEMENTATION_ja.md` — 適応基底実装仕様
- `doc/NOTE_RT_DG_TDDFT.md` — RT-DG-TDDFT ノート
