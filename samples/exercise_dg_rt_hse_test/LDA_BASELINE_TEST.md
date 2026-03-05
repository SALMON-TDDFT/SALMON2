# LDA テスト（ベースライン検証）

## 🎯 目的

**DG-Fragment RT-TDDFT 実装の sanity check**

LDA は最もシンプルな関数で、以下を確認:
- DG-Fragment initialization が correct
- Time propagation が numerically stable
- Output (energy, dipole, etc.) が reasonable
- **その後 HSE（より複雑）を安全に追加**

## ✅ テスト実行順序

### 推奨フロー

```
1️⃣ LDA H2 テスト      (1分) ← 最初に実行！
   ✓ 実行成功?
   ✓ Output ファイル生成?
   ✓ エネルギー値 reasonable?

2️⃣ LDA C6H6 テスト    (2分)
   ✓ 中規模システム OK?
   ✓ メモリー reasonable?

3️⃣ HSE Plan A (参考値を得る)
4️⃣ HSE Plan C (RI/DF の精度確認)
5️⃣ HSE CD-RI (メモリー削減確認)
```

## 📊 LDA テストケース

### H2 (H2/inputfile_lda)

```fortran
&functional
  xc = 'lda'              ! Local Density Approximation
/
```

**特徴**:
- ❌ HSE functional なし
- ❌ RI/DF approximation なし
- ✅ Just basic DG-Fragment initialized
- ✅ Super fast (< 1 min for 10 steps)

**期待値**:
- 初期化: ~30秒
- 10ステップ: ~1秒
- **合計: < 1分**

### C6H6 (C6H6/inputfile_lda)

```fortran
&functional
  xc = 'lda'              ! LDA only, no GGA
/
```

**特徴**:
- ❌ No hybrid functional
- ✅ Medium-scale system (12 atoms)
- ✅ Tests fragment handling

**期待値**:
- 初期化: ~1分 (DC-LCFO construction)
- 100ステップ: ~10秒
- **合計: < 2分** (very fast!)

## 🚀 実行コマンド

```bash
# テストディレクトリに移動
cd samples/exercise_dg_rt_hse_test

# LDA baseline テスト（3分で完了）
./run_h2_test.sh quick  # This includes LDA
./run_c6h6_test.sh quick

# または標準テスト
./run_h2_test.sh        # LDA + Plan A/C/CD-RI
```

## ✅ 成功の判定基準

### LDA H2 テスト

```
✅ Exit code = 0 (no errors)
✅ stdout.log に "converged" または "Success"
✅ Total energy printed (e.g., "Total energy = -1.123456 a.u.")
✅ Output ファイル存在:
   - Jxyz*.out (current)
   - *.txt (log files)
```

### LDA C6H6 テスト

```
✅ Exit code = 0
✅ 100 timesteps completed
✅ Total energy reasonable (negative, matches LDA theory)
✅ Computation time < 2 minutes
```

## 📈 LDA から HSE へのパス

LDA テスト成功後:

### HSE Plan A (最も正確 = 参照値)
```fortran
&functional
  xc = 'pbe'
  yn_hse = 'y'           ! Enable HSE
  hse_alpha = 0.25
  hse_omega = 0.11
  yn_hse_ri = 'n'        ! Full integral (no approximation)
/
```

### HSE Plan C (高速)
```fortran
&functional
  xc = 'pbe'
  yn_hse = 'y'
  yn_hse_ri = 'y'        ! RI/DF enabled
  hse_ri_ratio = 3.0
/
```

### HSE CD-RI (最適)
```fortran
&functional
  xc = 'pbe'
  yn_hse = 'y'
  yn_hse_ri = 'y'
  yn_hse_cd_ri = 'y'     ! Cholesky decomposition
  hse_cd_ri_threshold = 1.0d-8
/
```

## 💡 Troubleshooting

### LDA テストが失敗する場合

| 症状 | 原因 | 対策 |
|-----|------|------|
| Segmentation fault | メモリー不足 | グリッド粗くする (dl=0.8) |
| Slow initialization | DG-fragment overhead | 正常（LDA initialization） |
| NaN in energy | Numerical issue | グリッド間隔を適切に |
| Compile error | Macro undefined | make clean && make -j8 |

### 詳細ログ確認

```bash
# ログファイル確認
cat H2/logs/stdout_lda.log | head -50
cat H2/logs/stdout_lda.log | tail -50

# エネルギー抽出
grep -i "energy" H2/logs/stdout_lda.log
```

## 📝 記録

テスト完了後、以下を記録:

```
【LDA H2 テスト】
実行日時: YYYY-MM-DD HH:MM
結果: [成功/失敗]
所要時間: XXX 秒
エネルギー: XXXX.XX a.u.
メモリー: XX MB

【LDA C6H6 テスト】
実行日時: 
結果: [成功/失敗]
所要時間: XXX 秒
エネルギー: XXXX.XX a.u.
メモリー: XX MB
```

## 🎓 学習ポイント

LDA テストで以下を理解:

1. **DG-Fragment 実装の正確性**
   - LDA = sanity check
   - 基本機能が correct か?

2. **Performance baseline**
   - LDA: 最速 baseline
   - HSE: LDA + exchange integral overhead

3. **Scale understanding**
   - H2 1sec/step → C6H6 0.1sec/step (LDA)
   - Which means HSE takes 60-120x longer

## 🔗 参考資料

- [README.md](README.md) - テスト全体概要
- [EXECUTION_GUIDE.md](EXECUTION_GUIDE.md) - 詳細実行方法
- [TEST_PREPARATION.md](TEST_PREPARATION.md) - 準備サマリー

---

**重要**: LDA テスト成功 = HSE テストも成功する可能性が高い！

**次のステップ**: `./run_h2_test.sh quick` で LDA baseline を確認 🚀
