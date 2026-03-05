# DG RT-TDDFT HSE Test Execution Guide

## 🚀 Quick Start

```bash
cd /path/to/SALMON-v.2.2.2/samples/exercise_dg_rt_hse_test

# H2テスト実行（5-10分）
./run_h2_test.sh

# C6H6テスト実行（30-60分）
./run_c6h6_test.sh
```

## 📋 詳細なテスト手順

### Step 1: 前提条件確認

```bash
# SALMON がコンパイル完了しているか確認
ls /path/to/SALMON-v.2.2.2/mybuild/src/salmon

# 疑似ポテンシャルファイルを確認
ls /path/to/SALMON-v.2.2.2/samples/*.dat
```

### Step 2: H2 テスト実行

```bash
cd /path/to/SALMON-v.2.2.2/samples/exercise_dg_rt_hse_test

# 標準テスト（フル実行）
./run_h2_test.sh

# クイックテスト（5分以内）
./run_h2_test.sh quick

# 結果をクリア
./run_h2_test.sh clean
```

### Step 3: C6H6 テスト実行

```bash
# 標準テスト
./run_c6h6_test.sh

# クイックテスト（15分以内）
./run_c6h6_test.sh quick

# 結果をクリア
./run_c6h6_test.sh clean
```

### Step 4: 結果分析

```bash
# ログファイルの確認
cat H2/logs/stdout_plan_a.log
cat H2/logs/stdout_plan_c.log
cat H2/logs/stdout_cdri.log

# 簡単な比較
grep "Total energy" H2/logs/*
```

## 💾 ファイル構成

```
exercise_dg_rt_hse_test/
├── README.md                      ← テスト説明
├── EXECUTION_GUIDE.md             ← このファイル
├── run_h2_test.sh                ← H2テスト実行スクリプト
├── run_c6h6_test.sh              ← C6H6テスト実行スクリプト
│
├── H2/                            ← H2テスト入力ファイル
│   ├── inputfile_plan_a          ← Plan A (直接積分)
│   ├── inputfile_plan_c          ← Plan C (RI/DF)
│   ├── inputfile_cdri            ← CD-RI (Cholesky)
│   ├── logs/                      ← テスト出力
│   │   ├── stdout_a.log
│   │   ├── stdout_c.log
│   │   └── stdout_cdri.log
│   └── work_*/                    ← 作業ディレクトリ
│
├── C6H6/                          ← C6H6テスト入力ファイル
│   ├── inputfile_plan_c          ← Plan C (RI/DF)
│   ├── inputfile_cdri            ← CD-RI (推奨)
│   ├── logs/                      ← テスト出力
│   │   ├── stdout_c.log
│   │   └── stdout_cdri.log
│   └── work_*/                    ← 作業ディレクトリ
│
└── ANALYSIS.md                    ← 結果分析テンプレート
```

## 📊 期待値

### H2 Test 期待値

| Plan | 初期化 | ステップ | メモリー | 精度 |
|------|--------|---------|----------|------|
| A | 5分 | 120秒 | 800MB | 基準 |
| C | 1分 | 2秒 | 15MB | 基準 ±0.0003 |
| CD-RI | 1分 | 1秒 | 8MB | 基準 ±0.0001 |

### C6H6 Test 期待値

| Plan | 初期化 | ステップ | メモリー | 相対速度 |
|------|--------|---------|----------|---------|
| C | 10分 | 10秒 | 150MB | 100% |
| CD-RI | 10分 | 8秒 | 75MB | 125% ✨ |

## 🔍 トラブルシューティング

### テストが遅い場合

```bash
# クイックテストで確認
./run_h2_test.sh quick

# 必要に応じて以下を試す
# 1. 時間ステップ数を削減 (nt = 5)
# 2. グリッド間隔を拡大 (dl = 0.8 or 1.0)
# 3. 並列数を増やす (make -j16)
```

### メモリー不足エラー

```bash
# CD-RI を使用（メモリー50-80%削減）
# inputfile で yn_hse_cd_ri = 'y' に設定
```

### Compilation error at test

```bash
# SALMON を再度コンパイル
cd /path/to/SALMON-v.2.2.2/mybuild
make clean
cmake ..
make -j8
```

## 📝 結果記録

テスト完了後、以下の情報を記録してください：

```
【H2テスト結果】
日時: 
環境: 
Plan A 実行: [成功/失敗]  所要時間: 
Plan C 実行: [成功/失敗]  所要時間: 
CD-RI 実行: [成功/失敗]  所要時間: 

メモリー使用量:
Plan A: 
Plan C: 
CD-RI: 

エネルギー精度:
Plan A: (基準値)
Plan C: (vs Plan A)
CD-RI: (vs Plan A)

【C6H6テスト結果】
Plan C 実行: [成功/失敗]  所要時間: 
CD-RI 実行: [成功/失敗]  所要時間: 

メモリー削減率: 
性能向上: 
```

## 🎯 次のステップ

1. ✅ H2 テスト成功
   → Plan C / CD-RI が Plan A と一致することを確認

2. ✅ C6H6 テスト成功
   → メモリー削減（50-80%予測）を確認

3. 📈 FUGAKU で大規模テスト
   → Si64, Si128 での実行可能性を検証

4. 📝 論文執筆
   → CPU 版実装の性能測定結果をまとめ

## 💡 参考資料

- [テスト概要 (README.md)](README.md)
- [CD-RI実装ドキュメント](../../CD_RI_IMPLEMENTATION.md)
- [HSE実装レポート](../../HSE_RI_DF_METHOD_PLAN_C.md)
- [FUGAKU実装計画](../../OPENACC_GPU_IMPLEMENTATION.md)

## 📧 質問・問題報告

テスト実行中に問題が発生した場合は、以下を確認してください：

1. SALMON バージョン確認
   ```bash
   cat /path/to/SALMON-v.2.2.2/RELEASE_INFO
   ```

2. コンパイラバージョン確認
   ```bash
   gfortran --version
   ```

3. ログファイルの確認
   ```bash
   tail -100 H2/logs/stdout_plan_c.log
   ```

4. これらの情報を添えて報告

---

**作成日**: 2026-02-23  
**推奨環境**: macOS with gfortran or FUGAKU with Fujitsu compiler
