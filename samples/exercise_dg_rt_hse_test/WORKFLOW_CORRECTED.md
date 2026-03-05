# DG-Fragment RT-TDDFT: 正しい実行フロー完成

## ✅ 修正内容

ユーザーの指摘通り **GS → RT** の2段階ワークフローに修正しました。

## 📋 新しいテス ファイル構成

### H2 テストケース（修正版）

```
Ground State (GS) 計算:
  H2/inputfile_lda_gs         ← LDA GS（基底生成、一度だけ実行）

Real-Time (RT) 計算群:
  H2/inputfile_lda_rt         ← LDA RT（LDA GS 出力を使用）
  H2/inputfile_plan_a_rt      ← HSE Plan A RT
  H2/inputfile_plan_c_rt      ← HSE Plan C RT (既存)
  H2/inputfile_cdri_rt        ← HSE CD-RI RT (既存)
```

### 実行スクリプト

```
run_h2_test_proper.sh         ← 新规（GS→RT ワークフロー自動化）
  Usage:
    ./run_h2_test_proper.sh gs    # GS のみ実行
    ./run_h2_test_proper.sh rt    # すべての RT 実行
    ./run_h2_test_proper.sh full  # GS + すべての RT
```

## 🚀 実行方法

### 推奨フロー

```bash
cd /path/to/SALMON-v.2.2.2/samples/exercise_dg_rt_hse_test

# STEP 1: GS 計算（基底生成、一回だけ）
./run_h2_test_proper.sh gs

# STEP 2: RT 計算群（複数実行）
./run_h2_test_proper.sh rt

# または全実行
./run_h2_test_proper.sh full
```

### 出力ファイル

GS 実行後（STEP 1）:
```
H2/
├── nutilus_frag_001.dat    ← これを RT が読込
├── nutilus_frag_002.dat
├── Info_frag_001.txt
├── ene.out
└── stdout_lda_gs.log
```

RT 実行後（STEP 2）:
```
H2/
├── stdout_lda_rt.log
├── stdout_plan_a_rt.log
├── stdout_plan_c_rt.log
├── stdout_cdri_rt.log
├── Jxyz_0001.out           ← 各 RT が生成
└── dipole.txt
```

## 💡 ワークフローの論理

```
【LDA GS】
  Input:  H2 molecule coordinates
  Output: Fragment basis (nutilus_frag*.dat)
          ↓ ここから下で再利用

【LDA RT】
  Input:  Fragment basis (nutilus_frag*.dat from GS)
          LDA functional
  Output: Energy, current, dipole (LDA)

【Plan A RT】
  Input:  Fragment basis (same nutilus_frag*.dat)
          HSE direct functional
  Output: Energy, current, dipole (HSE direct - slow, reference)

【Plan C RT】
  Input:  Fragment basis (same nutilus_frag*.dat)
          HSE RI/DF functional
  Output: Energy, current, dipole (HSE RI/DF - 60倍高速)

【CD-RI RT】
  Input:  Fragment basis (same nutilus_frag*.dat)
          HSE Cholesky RI functional
  Output: Energy, current, dipole (HSE CD-RI - memory efficient)
```

## 📊 計算時間見積

```
GS 計算:        ~5-10 分（一度だけ）
LDA RT:          ~30 秒
Plan A RT:       ~5 分
Plan C RT:       ~1 分
CD-RI RT:        ~1 分

合計（full）:    ~15 分
```

## 🎯 検証フロー

```
✅ GS 実行
   └─ nutilus_frag*.dat 生成確認

✅ LDA RT 実行
   └─ エネルギー値確認

✅ Plan A RT 実行
   └─ エネルギー値記録（参照値）

✅ Plan C RT 実行
   └─ Plan A との比較（精度確認）

✅ CD-RI RT 実行
   └─ Plan C との比較（メモリー効率確認）
```

## 📝 次のステップ

1. **即座に実行**
   ```bash
   ./run_h2_test_proper.sh full
   ```

2. **C6H6 用も同様に作成**
   - inputfile_c6h6_lda_gs
   - inputfile_c6h6_lda_rt
   - inputfile_c6h6_plan_c_rt
   - inputfile_c6h6_cdri_rt

3. **大規模システム（Si64）へ向けた準備**
   - 富岳（FUGAKU）でのテスト計画

---

**重要**: これ以降のすべてのテストは **正しい GS→RT ワークフロー** に従います ✅
