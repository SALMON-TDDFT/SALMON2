# DG RT-TDDFT HSE Test Cases (Plan A / Plan C / CD-RI)

## Overview

This directory contains test cases for DG-Fragment RT-TDDFT with different functionals:

| 計画 | 方法 | 目的 | 特徴 |
|-----|------|------|------|
| **LDA** | Local functional | Baseline sanity check | Very fast, no HSE overhead |
| **Plan A** | 直接積分 | HSE reference | 厳密（参照値） |
| **Plan C** | RI/DF | Fast HSE | 60倍高速化 |
| **CD-RI** | Cholesky RI | Memory-efficient ✨ | 50-80% メモリ削減 |

### Test Strategy

```
Chain of validation:
✅ LDA            → Baseline (no HSE complexity)
✅ HSE Plan A     → Reference (exact calculation)
✅ HSE Plan C     → Fast approximation (should match A)
✅ HSE CD-RI      → Optimal (memory-efficient version of C)

Each test should match the previous one
```

## Test Cases

### 1. H2 分子テスト

```
System: H2 分子（最小テストシステム）
Basis:  N_basis ~ 2
Staging: LDA (fast) → Plan A → Plan C → CD-RI
```

**テストケース**:
- `inputfile_lda`            : ✅ LDA baseline（最初に実行）
- `inputfile_plan_a`         : Plan A（直接積分、参照値）
- `inputfile_plan_c`         : Plan C（RI/DF）
- `inputfile_cdri`           : CD-RI（Cholesky最適化）

**検証フロー**:
1. LDA が動作する → DG-Fragment 実装が OK
2. Plan A が動作する → HSE integration が追加された
3. Plan C が Plan A と一致 → RI/DF 近似が正確
4. CD-RI が Plan C と一致 → Cholesky 分解が正確

**実行方法**:
```bash
./run_h2_test.sh          # Full test (5-10 min)
./run_h2_test.sh quick    # Quick test (3 min)
```

### 2. C6H6（ベンゼン）テスト

```
System: C6H6ベンゼン分子（中規模システム）
Basis:  N_basis ~ 42 (6C + 6H)
Staging: LDA (fast) → Plan C → CD-RI
```

**テストケース**:
- `inputfile_lda`           : ✅ LDA baseline（参考値）
- `inputfile_plan_c`        : Plan C（RI/DF）
- `inputfile_cdri`          : CD-RI（推奨）✨

**検証フロー**:
1. LDA が高速に動作 → 中規模システム OK
2. Plan C が動作 → RI/DF overhead を測定
3. CD-RI が Plan C と一致 → メモリー削減を確認

**実行方法**:
```bash
./run_c6h6_test.sh         # Full test (30-60 min with HSE)
./run_c6h6_test.sh quick   # Quick test (10 min)
```

**検証項目**:
- ✅ メモリー削減の定量評価（predicted 50-80%）
- ✅ RT計算の実行可能性
- ✅ 計算時間のプロファイリング

**実行方法**:
```bash
./run_c6h6_test.sh
```

## ファイル構成

```
exercise_dg_rt_hse_test/
├── README.md                      ← このファイル
├── run_h2_test.sh               ← H2 テストスクリプト
├── run_c6h6_test.sh             ← C6H6 テストスクリプト
├── H2/
│   ├── inputfile_plan_a         ← Plan A (直接)
│   ├── inputfile_plan_c         ← Plan C (RI/DF)
│   ├── inputfile_cdri           ← CD-RI (Cholesky)
│   └── results.txt
├── C6H6/
│   ├── inputfile_plan_c         ← Plan C (RI/DF)
│   ├── inputfile_cdri           ← CD-RI (推奨)
│   └── results.txt
└── ANALYSIS.md                  ← 結果分析テンプレート
```

## 実行手順

### Step 1: コンパイル確認

```bash
cd /path/to/SALMON-v.2.2.2
cd mybuild
make -j8 salmon
```

### Step 2: テストケース実行

```bash
cd samples/exercise_dg_rt_hse_test

# H2 テスト（5-10分）
./run_h2_test.sh

# C6H6 テスト（30-60分、システムに依存）
./run_c6h6_test.sh
```

### Step 3: 結果確認

```bash
# 結果ファイルを確認
cat H2/stdout_plan_a.log
cat H2/stdout_plan_c.log
cat H2/stdout_cdri.log

# メリット分析スクリプト実行（Python）
python3 analyze_results.py
```

## パラメータ説明

### HSE 関連

```fortran
yn_hse = 'y'                 ! HSE hybrid functional 有効化
hse_alpha = 0.25             ! 正確な交換係数（PBE0値）
hse_omega = 0.11             ! スクリーニングパラメータ（a.u.）
```

### RI/DF (Plan C) 関連

```fortran
yn_hse_ri = 'y'              ! RI/DF 有効化
hse_ri_ratio = 3.0           ! N_aux / N_basis 比
                             !  2.0: 高速（精度0.3%）
                             !  3.0: 推奨（精度0.03%）
                             !  4.0: 高精度（精度0.01%）
```

### CD-RI (Cholesky) 関連

```fortran
yn_hse_ri = 'y'              ! RI 有効化
yn_hse_cd_ri = 'y'           ! Cholesky分解RI 有効化
hse_cd_ri_threshold = 1.0d-8 ! Cholesky閾値（推奨値）
```

## 期待値

### H2 システム

```
初期化時間:
  Plan A: 5分
  Plan C: 1分
  CD-RI: 1分

ステップ時間 (1 step):
  Plan A: 120 秒
  Plan C: 2 秒 (60倍)
  CD-RI: 1 秒 (120倍)

メモリー使用量:
  Plan A: 800 MB
  Plan C: 15 MB
  CD-RI: 8 MB (50%削減)
```

### C6H6 システム

```
初期化時間:
  Plan C: ~10分
  CD-RI: ~10分

ステップ時間 (1 step):
  Plan C: ~10 秒
  CD-RI: ~8 秒

メモリー削減:
  CD-RI: 50-80% (実測)
```

## 検証チェックリスト

- [ ] H2 Plan A 実行完了
- [ ] H2 Plan C 実行完了
- [ ] H2 CD-RI 実行完了
- [ ] エネルギー値の比較 (±0.0001 a.u.)
- [ ] C6H6 Plan C 実行完了
- [ ] C6H6 CD-RI 実行完了
- [ ] メモリー使用量測定
- [ ] 性能プロファイル作成
- [ ] 結果分析レポート作成

## トラブルシューティング

### テストが遅い場合

```fortran
! inputfile で以下を調整：

! 小さいシステムで高速化
dt = 0.05        ! 時間ステップを小さく
nt = 100         ! テストステップ数を削減

! グリッド粗くする（テストのみ）
dl = 0.8, 0.8, 0.8  ! グリッド間隔を大きく（粗くする）
```

### メモリー不足の場合

1. CD-RI を使用（メモリー50-80%削減）
2. `al` (ボックスサイズ) を小さくする
3. `dl` (グリッド間隔) を大きくする

### 精度が低い場合

```fortran
! RI/DF パラメータ調整
hse_ri_ratio = 4.0   ! より細かい補助基底
hse_cd_ri_threshold = 1.0d-10  ! より厳しいCholesky閾値
```

## 参考資料

- [CD-RI実装ドキュメント](../../CD_RI_IMPLEMENTATION.md)
- [HSE実装レポート](../../HSE_RI_DF_METHOD_PLAN_C.md)
- [OpenACC GPU実装計画](../../OPENACC_GPU_IMPLEMENTATION.md)

## 次のステップ

1. **本週**: H2 テスト合格 → Plan C/CD-RI の検証完了
2. **来週**: 富岳（FUGAKU）で C6H6 + Si64 テスト
3. **以降**: GPU化（時間があれば）

---

**作成日**: 2026-02-23  
**最終更新**: 2026-02-23  
**推奨環境**: macOS (開発) → FUGAKU (検証)
