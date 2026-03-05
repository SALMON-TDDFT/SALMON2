# DG RT-TDDFT テスト実装ガイド（最終版）

## 🎯 テスト全体戦略

```
Step 1: LDA baseline (sanity check)
  ├─ H2 (1分)       → DG-Fragment init OK?
  └─ C6H6 (2分)     → Medium system OK?

Step 2: HSE methods comparison (accuracy validation)
  ├─ H2 Plan A (5分)   → Reference calculation
  ├─ H2 Plan C (1分)   → RI/DF approximation
  └─ H2 CD-RI (1分)    → Cholesky decomposition

Step 3: Medium-scale system (C6H6)
  ├─ C6H6 Plan C (15分)  → RI/DF scaling
  └─ C6H6 CD-RI (15分)   → Memory reduction measurement

Total time: ~45分 (quick test: ~15分)
```

## 📋 準備完了ファイル一覧

### H2 テストケース
```
H2/inputfile_lda            ✅ LDA baseline
H2/inputfile_plan_a         ✅ HSE Plan A (direct)
H2/inputfile_plan_c         ✅ HSE Plan C (RI/DF)
H2/inputfile_cdri           ✅ HSE CD-RI (Cholesky)
```

### C6H6 テストケース
```
C6H6/inputfile_lda          ✅ LDA baseline
C6H6/inputfile_plan_c       ✅ HSE Plan C (RI/DF)
C6H6/inputfile_cdri         ✅ HSE CD-RI (Cholesky)
```

### テスト実行スクリプト
```
run_h2_test.sh              ✅ H2 automated test
run_c6h6_test.sh            ✅ C6H6 automated test
```

### ドキュメント
```
README.md                   ✅ テスト概要・理論
EXECUTION_GUIDE.md          ✅ 詳細実行手順
TEST_PREPARATION.md         ✅ 準備サマリー
LDA_BASELINE_TEST.md        ✅ LDA 専用ガイド
```

## 🚀 実行ステップ

### Phase 0: 前提条件確認（1分）

```bash
# SALMON がコンパイルされているか確認
ls /path/to/SALMON-v.2.2.2/mybuild/src/salmon

# テストディレクトリに移動
cd /path/to/SALMON-v.2.2.2/samples/exercise_dg_rt_hse_test

# ファイル確認
ls -la H2/ C6H6/
```

### Phase 1: LDA Baseline テスト（5分）

```bash
# クイックテスト（LDA only）
./run_h2_test.sh quick

# 結果確認
cat H2/logs/stdout_lda.log | grep "Total energy"
```

**期待値**:
- ✅ 実行完了（exit code 0）
- ✅ Total energy 出力
- ✅ 計算時間 < 5分

**成功判定**: LDA H2 が完了 → 次へ

### Phase 2: H2 HSE 比較テスト（10分）

```bash
# 標準テスト（LDA + Plan A/C/CD-RI）
./run_h2_test.sh

# エネルギー比較
grep "Total energy" H2/logs/stdout_*.log
```

**期待値**:
```
stdout_lda.log:     E = -1.123456 a.u.
stdout_plan_a.log:  E = -1.654321 a.u.  (LDA と異なる = OK)
stdout_plan_c.log:  E = -1.654318 a.u.  (Plan A と一致)
stdout_cdri.log:    E = -1.654320 a.u.  (Plan A と一致)
```

**成功判定**: 
- ✅ Plan C が Plan A と ±0.0003 以内で一致
- ✅ CD-RI が Plan C と ±0.0001 以内で一致

### Phase 3: C6H6 テスト（30-60分）

```bash
# クイックテスト（10分）
./run_c6h6_test.sh quick

# または標準テスト
./run_c6h6_test.sh
```

**期待値**:
- ✅ LDA: 最速（2分以下）
- ✅ Plan C: RI/DF overhead（15-20分）
- ✅ CD-RI: 約 Plan C と同じ速度（slightly faster）
- ✅ CD-RI memory: Plan C の 50-70%

**成功判定**:
- ✅ すべてのテスト完了
- ✅ CD-RI メモリー削減確認

## 📊 期待される結果

### H2 での期待値

```
【LDA】
  初期化: 30秒
  10ステップ: 1秒
  合計: <1分

【HSE Plan A】
  初期化: 5分
  ステップ: 120秒/step
  合計: ~20分 (slow!)

【HSE Plan C】
  初期化: 1分
  ステップ: 2秒/step
  合計: ~1.5分 (60倍高速)

【HSE CD-RI】
  初期化: 1分
  ステップ: 1秒/step
  合計: ~1分 (120倍高速)
```

### C6H6 での期待値

```
【LDA】
  初期化: 1分
  100ステップ: 10秒
  合計: <2分

【HSE Plan C】
  初期化: 10分
  100ステップ: 1000秒 (~17分)
  合計: ~30分

【HSE CD-RI】
  初期化: 10分
  100ステップ: 800秒 (~13分)
  合計: ~25分
  メモリー: Plan C の 50-70%
```

## 🎓 各テストの検証目標

| テスト | 検証内容 | 成功条件 |
|--------|---------|--------|
| H2 LDA | DG-Fragment 基本動作 | 実行完了、reasonable エネルギー |
| C6H6 LDA | 中規模システム対応 | 実行完了 |
| H2 Plan A | HSE 積分正確性 | 実行完了 |
| H2 Plan C | RI/DF 近似精度 | Plan A と ±0.0003 一致 |
| H2 CD-RI | Cholesky 精度 | Plan C と ±0.0001 一致 |
| C6H6 Plan C | RI/DF スケーリング | 実行完了、メモリー測定 |
| C6H6 CD-RI | メモリー削減効果 | 50% 以上削減確認 |

## 💡 トラブルシューティング

### テストが遅い場合

```bash
# グリッド粗くする（テスト用）
# inputfile で以下を修正：
dl = 0.8, 0.8, 0.8  (default 0.5)
nt = 5  (default 10)
```

### メモリー不足

```bash
# CD-RI を使用（メモリー 50-80% 削減）
# または system サイズを小さく
al = 15d0, 15d0, 15d0  (小さい box)
```

### コンパイル済みか不確か

```bash
# 再コンパイル
cd /path/to/SALMON-v.2.2.2/mybuild
cmake ..
make -j8
```

## 📈 次のはfase

テスト完了後:

1. **結果まとめ** (1時間)
   - エネルギー比較表
   - 性能測定表
   - メモリー削減率確認

2. **FUGAKU 準備** (1日)
   - テストケース転送
   - FUGAKU 環境セットアップ
   - 大規模テスト計画

3. **論文執筆** (1週間)
   - CPU 版（Plan C / CD-RI）の実装＆性能
   - FUGAKU での計算結果
   - GPU化の見通し

## 🎯 最後のチェックリスト

```
テスト準備:
  ☑ SALMON compiled
  ☑ Test files created
  ☑ Scripts executable
  ☑ Pseudo-potentials available

実行:
  ☑ LDA H2 successful
  ☑ LDA C6H6 successful
  ☑ HSE H2 all methods successful
  ☑ HSE C6H6 successful

検証:
  ☑ Energy values reasonable
  ☑ Accuracy: Plan C ≈ Plan A
  ☑ Accuracy: CD-RI ≈ Plan C
  ☑ Memory reduction: CD-RI < Plan C

Next:
  ☑ Results documented
  ☑ FUGAKU test planned
  ☑ GPU strategy confirmed
```

---

## 🚀 開始コマンド

```bash
# テストディレクトリへ移動
cd /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/samples/exercise_dg_rt_hse_test

# 推奨実行順序：
# 1. クイックテスト（15-20分）
./run_h2_test.sh quick && ./run_c6h6_test.sh quick

# 2. 標準テスト（1-2時間）
./run_h2_test.sh && ./run_c6h6_test.sh
```

---

**準備完了**: 2026-02-23
**推奨開始**: 本日または明日
**目標完了**: 2026-02-28（Phase 1 CPU テスト完了）

**次のステップ**: `./run_h2_test.sh quick` を実行！ 🎉
