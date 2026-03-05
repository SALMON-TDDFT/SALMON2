# DG-Fragment RT-TDDFT 計算フロー

## 🎯 正しい実行順序

DG-Fragment RT-TDDFT には **2段階** の計算が必要です：

```
┌─────────────────────────────────────────┐
│ STEP 1: Ground State (GS) Calculation   │
│                                         │
│ calc_mode = 'GS'                       │
│ → DC-LCFO fragment basis 生成           │
│ → nutilus_frag*.dat, Info_frag*.txt    │
│                                         │
│ Output: Fragment orbital basis files    │
│ Time: ~5-10 分 (LDA)                   │
└─────────────────────────────────────────┘
                    ↓
          (save output files)
                    ↓
┌─────────────────────────────────────────┐
│ STEP 2: Real-Time (RT) Calculation      │
│                                         │
│ calc_mode = 'RT'                       │
│ use_dg_fragment = 'y'                  │
│ → GS output (nutilus_frag*.dat) を読込 │
│ → Time evolution 実行                   │
│                                         │
│ Output: Jxyz*.out, dens.txt, etc.      │
│ Time: 数秒～数分 (system size など)    │
└─────────────────────────────────────────┘
```

## 📋 テスト実行手順

### Phase 0: GS 計算（基底生成）

```bash
cd /path/to/SALMON-v.2.2.2/samples/exercise_dg_rt_hse_test/H2

# 1. LDA GS 計算（必須）
SALMON < inputfile_lda_gs > stdout_lda_gs.log

# 2. 出力ファイル確認
ls -la nutilus_frag*.dat
ls -la Info_frag*.txt
ls -la ene.out
```

**期待値**:
- ✅ 完了（exit code 0）
- ✅ Fragment basis files 生成（nutilus_frag*.dat）
- ✅ エネルギー収束メッセージ

### Phase 1: RT 計算（時間発展）

```bash
# 3. LDA RT 計算
SALMON < inputfile_lda_rt > stdout_lda_rt.log

# 4. 出力確認
ls -la Jxyz*.out
ls -la dipole.txt
```

**期待値**:
- ✅ 完了
- ✅ Jxyz*.out (電流密度)
- ✅ dipole.txt (双極子モーメント)

### Phase 2: HSE テスト

```bash
# 5. LDA GS の出力を保存（再利用）
mkdir -p backup_lda_basis
cp nutilus_frag*.dat backup_lda_basis/
cp Info_frag*.txt backup_lda_basis/

# 6. HSE Plan A RT 計算
#    (LDA GS base + HSE RT functional)
SALMON < inputfile_plan_a_rt > stdout_plan_a_rt.log

# 7. 出力確認
grep "Total energy" stdout_plan_a_rt.log
```

## 🔄 高度なワークフロー

### Option 1: LDA GS をすべてのテストで再利用

```
LDA GS (実行一回)
    ├─→ LDA RT      (LDA basis + LDA functional)
    ├─→ Plan A RT   (LDA basis + HSE direct)
    ├─→ Plan C RT   (LDA basis + HSE RI/DF)
    └─→ CD-RI RT    (LDA basis + HSE CD-RI)
```

**利点**:
- GS 計算が1回のみ
- すべての RT 計算の条件が統一
- 計算時間大幅短縮

**注意**:
- GS basis は LDA
- RT functional は各々異なる（問題なし）

### Option 2: 各方式で GS を生成（厳密）

```
LDA GS       → LDA RT
Plan A GS    → Plan A RT
Plan C GS    → Plan C RT
CD-RI GS     → CD-RI RT
```

**利点**:
- 完全に自己無撞着
- GS と RT で同じ functional

**欠点**:
- 計算時間が 4倍
- HSE GS は非常に遅い

## 💡 推奨アプローチ

**LDA GS をすべてで再利用** (Option 1) が最実用的:

```bash
# GS: 一回だけ実行
SALMON < inputfile_lda_gs

# RT: 複数回実行（GS 出力を再利用）
SALMON < inputfile_lda_rt
SALMON < inputfile_plan_a_rt
SALMON < inputfile_plan_c_rt
SALMON < inputfile_cdri_rt
```

**理由**:
1. ✅ GS 計算が1回（時間大幅短縮）
2. ✅ Fragment basis は物理的に reasonable（LDA でも低い励起エネルギーは OK）
3. ✅ RT 計算での functional 違いは明確に測定可能
4. ✅ SALMON の通常的な利用法

## 📂 ディレクトリ構成（推奨）

```
H2/
├── inputfile_lda_gs       (GS: 一回だけ実行)
├── inputfile_lda_rt
├── inputfile_plan_a_rt
├── inputfile_plan_c_rt
├── inputfile_cdri_rt
│
├── stdout_lda_gs.log      (GS 出力)
├── nutilus_frag*.dat      (GS basis - 全 RT で共有)
├── Info_frag*.txt
├── ene.out
│
└── results/               (各 RT 計算の出力)
    ├── Jxyz_lda_*.out
    ├── Jxyz_plan_a_*.out
    ├── Jxyz_plan_c_*.out
    └── dipole.txt
```

## 🏃 クイック実行フロー

```bash
cd H2

# Step 1: GS 計算（一回）
salmon < inputfile_lda_gs > stdout_lda_gs.log

# Step 2: RT 計算群（複数実行）
for model in lda plan_a plan_c cdri; do
    echo "Running RT for $model..."
    salmon < inputfile_${model}_rt > stdout_${model}_rt.log
done

# Step 3: 結果比較
grep "Total energy" stdout_*_rt.log
```

## 🎯 各テストケース間の関係

```
inputfile_lda_gs
    ↓（generates nutilus_frag*.dat）
    ├→ inputfile_lda_rt
    ├→ inputfile_plan_a_rt
    ├→ inputfile_plan_c_rt
    └→ inputfile_cdri_rt
```

所有している fragment basis (nutilus_frag*.dat) は、
すべての RT 計算で自動的に読み込まれます。

---

**結論**: 
- GS は LDA で一回実行
- RT を複数実行してさまざまな functional をテスト
- これが最も実用的で時間効率が良い

**次のステップ**: RT 計算用ファイルのみを作成
