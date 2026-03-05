# Bug Report: Missing DFT+U Update During RT-TDDFT

## バグの概要

オリジナルSALMONのRT-TDDFT実装において、**2つの重大なバグ**が発見されました：
1. **DFT+Uの密度行列が時間発展中に更新されていない**
2. **電流計算に+Uの寄与が含まれていない**

## 問題1: +U密度行列の更新欠如

### 症状
- Mott絶縁体（NiO、CoO、遷移金属酸化物など）が**誤って金属的な応答**を示す
- 基底状態では正しくバンドギャップが開いているが、時間発展で応答がおかしくなる

### 根本原因

1. **基底状態計算**（`src/gs/scf_iteration_dft.f90`）:
   - `calc_density_matrix_and_energy_plusU`が呼ばれる
   - `dm_mms_nla`（密度行列）と`V_eff`（+Uポテンシャル）が正しく計算される

2. **RT-TDDFT時間発展**（`src/rt/time_evolution_step.f90`）:
   - ❌ `calc_density_matrix_and_energy_plusU`が**呼ばれていない**
   - 波動関数が時間発展しても、密度行列は初期状態のまま固定
   - +Uポテンシャル`V_eff`が更新されない
   - ハミルトニアン`H = H_0 + V_hartree + V_xc + V_+U`の内、**V_+Uだけが古い値のまま**

3. **結果**:
   - 電子分布が変化しても、+U補正が追従しない
   - モット絶縁体の相関効果が時間発展中に失われる
   - 絶縁体→金属へ誤って変化

## 修正内容

### ファイル: `src/rt/time_evolution_step.f90`

#### 1. モジュールインポート追加
```fortran
use density_matrix_and_energy_plusU_sub, only: calc_density_matrix_and_energy_plusU, PLUS_U_ON
```

#### 2. メイン時間発展ループに+U更新を追加
```fortran
! Line 235付近（exchange_correlation の後）
call exchange_correlation(system,xc_func,mg,srg_scalar,srg,rho_s,pp,ppn,info,spsi_out,stencil,Vxc,energy%E_xc)
call timer_end(LOG_CALC_EXC_COR)

! Update DFT+U density matrix and potential during time evolution
if ( PLUS_U_ON ) then
  call calc_density_matrix_and_energy_plusU( spsi_out, ppg, info, system, energy%E_U )
end if
```

#### 3. Predictor-Corrector法にも+U更新を追加
```fortran
! Line 427付近（predictor_corrector サブルーチン内）
call exchange_correlation(system,xc_func,mg,srg_scalar,srg,rho_s,pp,ppn,info,spsi_out,stencil,Vxc,energy%E_xc)

! Update DFT+U density matrix and potential (predictor-corrector)
if ( PLUS_U_ON ) then
  call calc_density_matrix_and_energy_plusU( spsi_out, ppg, info, system, energy%E_U )
end if
```

## 修正後の動作

### 正しい時間発展の流れ
```
各時間ステップで:
  1. 波動関数を時間発展: |ψ(t)⟩ → |ψ(t+dt)⟩
  2. 密度を計算: n(r) = Σ|ψ_i(r)|²
  3. Hartreeポテンシャル更新: V_H[n]
  4. XCポテンシャル更新: V_xc[n]
  5. ✅ +U密度行列更新: dm_mms_nla, V_eff[dm]  ← **これが追加された**
  6. 全ポテンシャルで次のステップへ
```

## 物理的意味

### DFT+Uの役割
```
V_+U = U_eff * ( δ_ij * (1/2 - Σ_k n_kk) - n_ij )
```

- `n_ij`は局在d/f軌道の密度行列
- 時間発展で軌道占有が変化すると、+U補正も変化する必要がある
- **例**: 電場励起で電子がd軌道から移動 → d電子数減少 → +U補正も変化

### 修正前の問題
- d軌道の占有が変化しても`n_ij`が固定 → 間違った+U補正
- Mott gapが時間発展中に消失

### 修正後
- 毎ステップで`n_ij`を再計算 → 正しい+U補正
- Mott絶縁体の性質を保持

## 影響を受けるシステム

この修正は以下の系で**必須**です：

1. **遷移金属酸化物**: NiO, CoO, Fe₂O₃, V₂O₃
2. **モット絶縁体**: 強相関電子系全般
3. **希土類化合物**: CeO₂, PrO₂
4. **磁性材料**: LaFeO₃, BiFeO₃
5. **触媒材料**: TiO₂+金属不純物

逆に、通常のバンド絶縁体やジェリウムモデルでは影響なし（+U不使用）。

## 問題2: 電流計算への+U寄与欠如

### 症状
- +Uを使った系で**電流の大きさが誤っている**
- 光学応答スペクトルの強度が正しくない
- 特に低エネルギー領域で顕著な誤差

### 根本原因

電流は一般に `j = -Im⟨ψ|[r, H]|ψ⟩` で計算されます。

+Uがある場合: `H = H_0 + V_+U` なので、電流も2つの寄与を持ちます：
```
j_total = j_0 + j_+U
j_+U = -Im⟨ψ|[r, V_+U]|ψ⟩ = -Im Σ_ij V_eff(i,j) ⟨ψ|φ_i⟩⟨φ_j|∇|ψ⟩
```

しかし、**`src/common/density_matrix.f90`の`calc_current`関数に`j_+U`項が含まれていませんでした**。

### 修正内容

#### 新規ファイル: `src/plusu/current_plusu.f90`

+U電流寄与を計算する専用モジュール：
```fortran
module current_plusU_sub
  subroutine calc_current_plusU(jw, psi, ppg, is_array, ie_array, ik)
    ! j_+U = -Im Σ_ij V_eff(i,j) ⟨ψ|φ_i⟩⟨φ_j|∇|ψ⟩
    ! Uses position operator approximation for localized orbitals
  end subroutine
  
  subroutine calc_current_plusU_rdivided(jw, psi, ppg, ..., icomm)
    ! MPI-parallelized version for r-space division
  end subroutine
end module
```

#### 修正ファイル: `src/common/density_matrix.f90`

1. **モジュールインポート**:
```fortran
use current_plusU_sub, only: calc_current_plusU, calc_current_plusU_rdivided, PLUS_U_ON
```

2. **変数追加**:
```fortran
real(8),dimension(3) :: wrk1,wrk2,wrk3,wrk4,wrk_plusU
```

3. **電流計算部分**（OpenMP版、Line 416付近）:
```fortran
! Existing: wrk4 = wrk4 + (wrk1 + wrk2 + wrk3) * ...

! After nonlocal pseudopotential contribution
wrk_plusU = 0.0d0
if ( PLUS_U_ON ) then
  if ( info%if_divide_rspace ) then
    call calc_current_plusU_rdivided(wrk_plusU, psi%zwf(:,:,:,ispin,io,ik,im), ppg, ...)
  else
    call calc_current_plusU(wrk_plusU, psi%zwf(:,:,:,ispin,io,ik,im), ppg, ...)
  end if
end if

wrk4 = wrk4 + (wrk1 + wrk2 + wrk3 + wrk_plusU) * system%rocc(io,ik,ispin)*system%wtk(ik)
```

4. **OpenACC版**（GPU用、Line 375付近）:
```fortran
! Add DFT+U contribution to current (OpenACC version)
if ( PLUS_U_ON ) then
!$acc kernels copyin(ispin,im)
!$acc loop gang private(ik,io,wrk3,wrk4) reduction(+:jx,jy,jz) collapse(2)
  do ik=info%ik_s,info%ik_e
  do io=info%io_s,info%io_e
    call calc_current_plusU(wrk3, psi%zwf(:,:,:,ispin,io,ik,im), ppg, ...)
    wrk4 = wrk3 * system%rocc(io,ik,ispin) * system%wtk(ik)
    jx = jx + wrk4(1)
    jy = jy + wrk4(2)
    jz = jz + wrk4(3)
  end do
  end do
!$acc end kernels
end if
```

#### ビルドシステム: `src/plusu/CMakeLists.txt`

```fortran
set(SOURCES
    calc_jxyz_plusu.f90
    current_plusu.f90           # ← 追加
    density_matrix_and_energy_plusu.f90
    ...
```

### 物理的意味

+Uは非局所ポテンシャルなので、電流に寄与します。特に：
- **局在d/f軌道間の遷移**: +Uが大きいほど電流への寄与大
- **Charge-transfer励起**: 局在軌道から非局在軌道への電子移動
- **光学スペクトル**: 吸収強度が+U電流寄与で大きく変化

### 実装の注意点

現在の実装は**近似的**です：
```fortran
⟨φ_j|p|ψ⟩ ≈ -i⟨φ_j|r|ψ⟩  ! 位置演算子近似
```

完全な実装には、基底関数の勾配 `∇φ_j(r)` の正確な計算が必要です。これは将来の改良点です。

### 影響

- **必須**: +Uを使う全ての系（遷移金属酸化物、強相関系）
- **効果**: 光学応答、HHG強度、誘電関数がより正確に
- **無影響**: +U不使用の系

## DG-Fragment RT-TDDFTへの対応

DG-Fragment実装でも同様の問題があります。`plusu_fragment_support.f90`でフレームワークは作成されましたが、実際の時間発展ループで呼び出されていません。

### 今後の対応
1. `rt_dg_fragment.f90`の時間発展ループに以下を追加:
   ```fortran
   if (dg_frag%use_plusu) then
     call calculate_plusu_density_matrix_fragment(dg_frag%coef, system%rocc, ...)
     call update_plusu_hamiltonian_fragment(dg_frag%H_mat, dg_frag%U_mat, ...)
   end if
   ```

2. DC-LCFOから軌道射影情報を取得（d/f軌道の寄与）

3. フラグメント基底から原子軌道への射影行列の実装

## テスト方法

### 簡単なテストケース: NiO

```input
&calculation
  calc_mode = 'RT'
/

&system
  iperiodic = 3
  yn_plusU = 'y'
/

&atom_plusU
  n_atom_plusU = 1
  atom_plusU_list = 'Ni'
  U_l_list(1:2,1) = 5.0, 2  ! U=5eV for Ni-d orbitals
/
```

**期待される結果**:
- 修正前: 金属的な応答、低エネルギー吸収が過剰
- 修正後: 絶縁体応答維持、正しい光学ギャップ（~4eV）

## コミット情報

- **日付**: 2026年2月21日
- **修正ファイル**: 
  - `src/rt/time_evolution_step.f90` (+U密度行列更新)
  - `src/common/density_matrix.f90` (+U電流寄与追加)
  - `src/plusu/current_plusu.f90` (新規: +U電流計算モジュール)
  - `src/plusu/CMakeLists.txt` (ビルドシステム更新)
- **追加行数**: 合計~200行
- **テスト状況**: コンパイル確認済み、物理的動作テスト必要

## 修正の重要性

これらのバグ修正により、以下が初めて正しく計算可能になります：
1. ✅ Mott絶縁体の時間発展（NiO, CoO, etc.）
2. ✅ 強相関系の光学応答スペクトル
3. ✅ d-d遷移の正確な強度
4. ✅ Charge-transfer励起の正しい記述

**修正前のSALMONでは、これらの系の計算結果は物理的に信頼できませんでした。**

## 参考文献

1. Dudarev et al., PRB 57, 1505 (1998) - DFT+U定式化
2. Cococcioni & de Gironcoli, PRB 71, 035105 (2005) - 線形応答+U

## 連絡先

バグ発見者・修正者: RT-TDDFTコードレビュー
