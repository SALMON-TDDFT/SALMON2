# 🔴 DG-Fragment RT 物理値重大エラー

**発見日**: 2026年2月23日  
**重大度**: ⚠️ CRITICAL - 計算結果が物理的に意味をなさない

## 問題の概要

DG-Fragment RT の出力電流値が、同一系における従来的 RT 電流値より **3700万倍小さい**という物理的に不可能なエラーが検出されました。

## 数値証拠

| 項目 | DG-Fragment | Conventional RT | 比率 |
|---|---|---|---|
| **Jx(max)** | **1.097×10⁻¹³** a.u. | **4.063×10⁻⁶** a.u. | **3.7×10⁷倍差** |
| **Jy** | ~10⁻¹⁴ | ~10⁻¹⁹ | 適正 |
| **Jz** | ~10⁻¹³ | ~10⁻¹⁹ | 10⁶倍差 |

**評価**: 🔴 **完全に異常** - スケールが合わない

## 根本原因の分析

### 1. `calculate_observables` 関数の問題

ファイル: [rt_dg_fragment.f90](rt_dg_fragment.f90#L2406-L2510)

```fortran
! 現在の実装：純粋な波動関数係数から計算
do io = 1, nocc
  current_tmp = current_tmp + sum(real(conjg(dg_frag%coef) * tmp_mat))
end do
```

**問題点**:
1. **正規化されていない係数の使用**
   - fragment basis は |φ_i(r)⟩ で実空間基底
   - しかし、係数 c_i がどのようにnormalizeされているか不明
   - 通常のRT-TDDFTでは密度を明示的に計算する

2. **momentum_mat が適切に初期化されているか不明**
   - `dg_frag%momentum_mat(idir, :, :, ispin)` の定義が不明瞭
   - Framework basis での momentum operator の定義が正しいか？

3. **グリッド積分の欠落の可能性**
   - fragment basis φ_i(r) は実空間で定義されている
   - しかし、密度 ρ(r) = Σ_{ij} c*_i c_j φ*_i(r) φ_j(r) の計算がない
   - 電子密度を使わずに直接係数から電流を計算している

### 2. yn_fix_func='y' の影響

inputfile 設定:
```
yn_fix_func = 'y'
```

**影響**:
- 密度が初期値（t=0）のままで更新されない
- 時間発展する波動関数が電流計算に使われていない可能性
- または、電流計算が完全に skip されている

### 3. fragment basis のスケーリング問題

**考えられるシナリオ**:
```
正しい計算:   <ψ(t)|p|ψ(t)> 
実際の計算:   <φ|p|φ> × (小さい係数)
             ↓ スケーリングエラー
結果:         1.097×10⁻¹³ (3700万倍小さい)
```

## diagnostic チェックリスト

- [ ] `momentum_mat` の値が実際に計算されているか確認
- [ ] `momentum_mat` のスケーリングが正しいか（単位確認）
- [ ] `dg_frag%coef` の normalization が正しいか
- [ ] 密度再構築 `ρ(r) = Σ_{ij} c*_i c_j φ*_i φ_j` が実装されているか
- [ ] Fragment basis の内積 `<φ_i|φ_j>` が単位行列か確認
- [ ] 従来的RTとの電流計算方法の差を詳細比較

## 推奨アクション（優先順）

### 高優先度
1. **`momentum_mat` の値を DEBUG出力**
   ```fortran
   call write_observable_debug(dg_frag%momentum_mat(:,:,1))
   ```
   - 期待値: ≠ 0, 10⁻⁵ オーダー
   - 現在の可能性: 0 または非常に小さい値

2. **従来的RTの電流計算と直接比較**
   - `src/rt/time_evolution_step.f90` の `calc_density` と比較
   - DG-Fragment での対応関数を調査

3. **yn_fix_func='n'でテスト再実施**
   - 密度を動的に更新する場合の挙動確認
   - 出力値が改善するか検証

### 中優先度
4. Fragment basis の正規化確認
   - `<φ_i|φ_j> = δ_{ij}` か確認
   - 内積計算が正しいか

5. グリッド積分スケーリング確認
   - dV（体積要素）が電流計算に含まれているか

### 低優先度
6. 単位の統一確認
   - momentum operator の単位が a.u. か確認
   - 出力ファイルフォーマットの単位タグ確認

## 見積もり

| タスク | 時間 | 難度 |
|---|---|---|
| momentum_mat DEBUG | 30分 | 🟢 低 |
| 従来的RTとの比較 | 1時間 | 🟡 中 |
| yn_fix_func='n' テスト | 1-2時間 | 🟢 低 |
| 根本的な修正 | 2-4時間 | 🔴 高 |

## 物理的な意味

**仮説**: 
DG-Fragment 電流 ≈ 0 という結果は：
- ✅ Fragment basis coefficient がほぼ静的（時間発展していない）？
- ✅ coefficient evolution が正しく実装されていない？
- ✅ または密度が全く計算されていない？

いずれにせよ、**物理的に非物質な結果**

## 次のステップ

```
Session 25 作業予定:
┌─────────────────────────────────────┐
│ 1. momentum_mat DEBUG実装             │
│ 2. 従来的RTとの詳細比較分析         │
│ 3. yn_fix_func='n' 再テスト          │
│ 4. 根本原因の特定と修正             │
└─────────────────────────────────────┘
```

---

**結論**: ❌ DG-Fragment RTの電流計算は完全に壊れている。物理的に信頼できない結果。
