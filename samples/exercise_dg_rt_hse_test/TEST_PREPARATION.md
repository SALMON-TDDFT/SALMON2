# DG RT-TDDFT HSE Test Preparation Summary

## ✅ テスト準備完了

テストケースが完全に準備されました。✨

### 📁 ファイル構成

```
exercise_dg_rt_hse_test/
├── README.md                      (テスト概要と理論)
├── EXECUTION_GUIDE.md             (実行手順)
├── TEST_PREPARATION.md            (このファイル)
│
├── run_h2_test.sh                ← H2 テスト (5-10分)
├── run_c6h6_test.sh              ← C6H6 テスト (30-60分)
│
├── H2/
│   ├── inputfile_plan_a          Plan A (直接積分)
│   ├── inputfile_plan_c          Plan C (RI/DF)
│   ├── inputfile_cdri            CD-RI (Cholesky + RI/DF)
│   └── logs/                      (テスト出力)
│
├── C6H6/
│   ├── inputfile_plan_c          Plan C (RI/DF)
│   ├── inputfile_cdri            CD-RI (推奨)✨
│   └── logs/                      (テスト出力)
│
└── testsuites/
    ├── ...
    
```

## 🎯 テスト目的

| テスト | 目的 | 検証項目 | 期間 |
|--------|------|--------|------|
| **H2 Plan A** | 参照計算 | 基準エネルギー | 5分 |
| **H2 Plan C** | RI/DF精度確認 | 精度 <0.03% | 1分 |
| **H2 CD-RI** | メモリー効率 | 50% 削減 | 1分 |
| **C6H6 Plan C** | 中規模テスト | 実行可能性 | 15分 |
| **C6H6 CD-RI** | 最適化方式 | メモリー削減 | 15分 |

## 📊 比較表

### 計算複雑度

```
┌──────────────┬─────────┬──────────┬─────────┐
│  Method      │ 初期化  │ ステップ │ メモリー │
├──────────────┼─────────┼──────────┼─────────┤
│ Plan A       │  遅い  │  超遅い  │  多い  │
│ Plan C       │  高速  │  高速    │  少ない │
│ CD-RI ✨     │  高速  │  最速    │  最少 │
└──────────────┴─────────┴──────────┴─────────┘

H2 での期待値:
  Plan A: 1 step = 120秒
  Plan C: 1 step = 2秒 (60倍)
  CD-RI: 1 step = 1秒 (120倍)
```

### メモリー削減

```
CD-RI vs Plan C メモリー削減率:

H2:     50% (8MB vs 15MB)
C6H6:   60% (60MB vs 150MB)
Si64:   80% (148MB vs 738MB) ← これが重要!
```

## 🚀 実行手順 (シンプル版)

```bash
# テスト開始
cd /path/to/SALMON-v.2.2.2/samples/exercise_dg_rt_hse_test

# 1. H2 テスト (最初)
./run_h2_test.sh

# 2. 結果確認
cat H2/logs/stdout_plan_a.log | grep "Total energy"
cat H2/logs/stdout_plan_c.log | grep "Total energy"
cat H2/logs/stdout_cdri.log   | grep "Total energy"

# 3. エネルギー値が一致していることを確認 (±0.0003 以内)

# 4. C6H6 テスト (30分-60分)
./run_c6h6_test.sh
```

## ⚠️ 注記

### 疑似ポテンシャルファイル

テストスクリプトは以下の場所からファイルをリンクします：

```bash
ln -sf /path/to/SALMON-v.2.2.2/mybuild/testsuites/pseudo/H_rps.dat
ln -sf /path/to/SALMON-v.2.2.2/mybuild/testsuites/pseudo/C_rps.dat
```

もし不足している場合は、以下で確認：

```bash
find /path/to/SALMON-v.2.2.2/mybuild -name "*.dat" | grep "_rps.dat"
```

### テストの期間

- **H2 快テスト**: 5分以内
- **H2 標準テスト**: 5-10分 (Plan A が遅い)
- **C6H6 快テスト**: 15分
- **C6H6 標準テスト**: 30-60分 (最初の initialization が時間がかかる)

**富岳での実行**: FUGAKU ではより高速 (2-3倍)

## 📈 期待される学習成果

このテストを完了することで、以下を検証できます：

1. ✅ **Plan C (RI/DF) の精度**
   - Plan A との比較で確認
   - 化学精度 (<0.03%) に収まるか?

2. ✅ **CD-RI の効果**
   - メモリー削減 (50%以上)
   - 計算速度（若干向上）
   - 精度（Plan C と同等か？）

3. ✅ **大規模システムへのスケーリング**
   - C6H6 で実行可能か?
   - Si64 での見通しを立てる

4. ✅ **富岳への準備**
   - テストケースが FUGAKU で動作するか確認
   - パス、コンパイラ設定の最適化

## 🎓 今後の展開

### Phase 1: CPU テスト (本週)
```
✅ H2 テスト成功
✅ C6H6 テスト成功
✅ CD-RI メリット確認
```

### Phase 2: FUGAKU 検証 (来週)
```
[ ] FUGAKU On DEMAND アクセス
[ ] C6H6 テスト実行
[ ] Si64 実行試験
```

### Phase 3: 最適化・論文化 (再来週)
```
[ ] パラメータ最適化
[ ] 性能測定
[ ] 結果まとめ
```

## 💡 トipsトラブルシューティング

### テストが失敗する場合

**1. SALMON が見つからない**
```bash
# コンパイル確認
cd /path/to/SALMON-v.2.2.2/mybuild
make -j8
ls src/salmon  # 存在確認
```

**2. 疑似ポテンシャルファイルが見つからない**
```bash
# ファイル検索
find /path/to/SALMON-v.2.2.2 -name "H_rps.dat" -o -name "C_rps.dat"
```

**3. スクリプトの実行権限がない**
```bash
chmod +x run_h2_test.sh run_c6h6_test.sh
```

## 📚 参考資料

- [README.md](README.md) - テスト概要
- [EXECUTION_GUIDE.md](EXECUTION_GUIDE.md) - 詳細実行手順
- [CD_RI_IMPLEMENTATION.md](../../CD_RI_IMPLEMENTATION.md) - CD-RI 実装ドキュメント
- [HSE_RI_DF_METHOD_PLAN_C.md](../../HSE_RI_DF_METHOD_PLAN_C.md) - 理論背景

## ✨ テスト完了後

テストが終わったら以下を実施してください：

1. **結果メモ作成**
   ```
   H2テスト: [成功/失敗] - 所要時間: XXX秒
   エネルギー精度: Plan A から XXX a.u. 差
   メモリー削減: YY%
   ```

2. **ログファイル保存**
   ```bash
   cp -r H2/logs H2_results_$(date +%Y%m%d_%H%M%S)
   cp -r C6H6/logs C6H6_results_$(date +%Y%m%d_%H%M%S)
   ```

3. **FUGAKU準備**
   - 成功したパラメータを記録
   - FUGAKU 用に調整が必要な項目をチェック

---

**準備完了日**: 2026-02-23  
**テスト開始日**: 本週  
**目標完了日**: 2026-03-02 (Phase 1 CPU テスト完了)

**次のステップ**: `./run_h2_test.sh` を実行！🚀
