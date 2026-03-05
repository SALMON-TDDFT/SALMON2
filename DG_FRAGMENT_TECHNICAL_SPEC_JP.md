# DG-Fragment RT実装 - 技術仕様書（日本語版）

**バージョン**: 2.2.2  
**作成日**: 2026年2月23日  
**対象**: SALMON DG-Fragment Real-Time動力学実装

---

## 目次

1. [システム概要](#システム概要)
2. [主要パラメータ](#主要パラメータ)
3. [メモリアクセス修正](#メモリアクセス修正)
4. [実装詳細](#実装詳細)
5. [テスト・検証](#テスト検証)
6. [トラブルシューティング](#トラブルシューティング)

---

## システム概要

### 何をするのか

**DG-Fragment RT** は、DC-LCFO（Density Cumulant - Local Coherence Function Optimization）法で計算した基底状態データを用いて、レーザー励起下での電子動力学をリアルタイムで計算する手法です。

### 計算の流れ

```
DC-LCFO GS計算
    ↓
    └─→ フラグメント基底データ読み込み
        ├─ 実空間波動関数 (phi_frag)
        ├─ 係数行列 (coef)
        ├─ インデックスマッピング (index_basis)
        └─ フラグメント分解情報
    ↓
RT初期化
    ├─ レーザーパラメータ設定
    ├─ ハミルトニアン行列計算
    └─ 時間積分準備
    ↓
時間発展ループ（Runge-Kutta）
    ├─ H_ij(t) 計算
    ├─ ψ(t) → ψ(t+dt) 伝播
    ├─ 観測量計算 J(t)、E(t)
    └─ nt タイムステップ繰返し
    ↓
結果出力
    ├─ 電流発展: _rt.data
    ├─ エネルギー発展: _rt_energy.data
    └─ 補助場情報: _rt_pulse.data
```

---

## 主要パラメータ

### 新規追加パラメータ

#### 1. `yn_dg_fragment_from_dcdft` (論理値)

**役割**: DC-LCFOデータからDG-Fragment RTを有効化

```fortran
&calculation
  yn_dg_fragment_from_dcdft = 'y'   ! DG-Fragment RTを有効化
/
```

**効果**:
- `'y'`: DC-LCFO基底データを読み込み、DG-Fragment RTを実行
- `'n'`: 従来的RT（デフォルト）を実行

**場所**: src/rt/rt_dg_fragment.f90
- init_dg_fragment_rt (行148)
  - パラメータを読み込む
  - ロード動作を制御

- main_tddft.f90 (行220)
  - パラメータによりハミルトニアン計算を実行

#### 2. `yn_dg_fragment_rt` (論理値)

**役割**: 時間発展でDG-Fragment RTを利用

```fortran
&propagation
  yn_dg_fragment_rt = 'y'
/
```

**効果**:
- `'y'`: DG-Fragment RT時間発展を使用
- `'n'`: 従来的時間発展

#### 3. `nstate_frag` (整数)

**役割**: フラグメント当たりの基底関数数

```fortran
&calculation
  nstate_frag = 40    ! 1フラグメント当たり40個の基底
/
```

**範囲**: 1 ～ nstate_tot  
**典型値**: nstate_tot / 数フラグメント

### 既存パラメータの相互作用

| パラメータ | DG-Fragment | 従来的 RT |
|-----------|-------------|----------|
| nstate | ✅ 必須 | ✅ 必須 |
| yn_periodic | ✅ 対応 | ✅ 対応 |
| al(1:3) | ✅ セル定義 | ✅ セル定義 |
| num_rgrid | ✅ グリッド | ✅ グリッド |
| nt, dt | ✅ 時間 | ✅ 時間 |
| emfield | ✅ レーザー | ✅ レーザー |

---

## メモリアクセス修正

### 問題の本質

#### 前（修正前）

```fortran
! 割り当て:
allocate(coef(nstate_frag, nstate_tot, nspin))  ! 次元: 40
allocate(momentum_mat(3, n_mat_max, n_mat_max, nspin)) ! 次元: 400

! アクセス:
do j = 1, n_mat_max  ! 400まで
  coef(j, :, :)      ! 40までしかない → 範囲外アクセス！
end do

! 結果: SIGBUS エラー ✗
```

#### 後（修正後）

```fortran
! 割り当て:
allocate(coef(n_mat_max, nstate_tot, nspin))  ! 次元: 400
allocate(momentum_mat(3, n_mat_max, n_mat_max, nspin)) ! 次元: 400

! 再索引化:
do i_local = 1, ifrag_count
  do io = 1, n_basis(ifrag, ispin)
    global_idx = index_basis(io, ifrag, ispin)
    coef(global_idx, :, :) = coef_local(io, :, :, i_local)
  end do
end do

! アクセス:
do j = 1, n_mat_max  ! 安全 ✓
  coef(j, :, :)      ! 範囲内アクセス
end do

! 結果: 正常動作 ✓
```

### 修正の詳細

#### ステップ1: 基底関数の再索引化（読み込み段階）

**ファイル**: src/rt/rt_dg_fragment.f90  
**関数**: read_fragment_basis_data  
**行**: 751-778

```fortran
! DC-LCFOから読み込んだ係数をグローバル順序に再編成
allocate(dg_frag%coef(dg_frag%n_mat_max, dg_frag%nstate_tot, dg_frag%nspin))
allocate(dg_frag%coef_new(dg_frag%n_mat_max, dg_frag%nstate_tot, dg_frag%nspin))
allocate(dg_frag%coef_work(dg_frag%n_mat_max, dg_frag%nstate_tot, dg_frag%nspin))

dg_frag%coef = 0.0d0
dg_frag%coef_new = 0.0d0
dg_frag%coef_work = 0.0d0

! フラグメント局所 → グローバル基底順へマッピング
do i_local = 1, ifrag_count
  do ispin = 1, dg_frag%nspin
    do io = 1, dg_frag%n_basis(ifrag, ispin)
      global_idx = dg_frag%index_basis(io, ifrag, ispin)
      
      ! 範囲チェック
      if (global_idx > 0 .and. global_idx <= dg_frag%n_mat_max) then
        ! グローバル位置に配置
        dg_frag%coef(global_idx, :, :) = coef_local(io, :, :, i_local)
      end if
    end do
  end do
end do
```

#### ステップ2: 変数宣言追加

**ファイル**: src/rt/rt_dg_fragment.f90  
**行**: 229

```fortran
! 追加された局所変数
integer :: io, global_idx
```

- `io`: フラグメント局所基底インデックス (1 ～ nstate_frag)
- `global_idx`: グローバル基底インデックス (1 ～ n_mat_max)

#### ステップ3: 早期割り当て削除

**ファイル**: src/rt/rt_dg_fragment.f90  
**行**: 235-243 (init_dg_fragment_rt関数)

```fortran
! 削除: 誤った次元での割り当て
! allocate(dg_frag%coef(nstate_frag, nstate_tot, nspin))

! 理由: n_mat_max値はまだ読み込まれていない
! 解決: read_fragment_basis_dataで正しい次元で割り当て
```

### メモリ配置図

#### 従来的RT

```
状態空間: 1次元配列
┌─────────────────┐
│ 状態1（全電子） │
│ 状態2（全電子） │
│ 状態3（全電子） │
│ ...             │
│ 状態N（全電子） │
└─────────────────┘
 dimension: nstate_tot
```

#### DG-Fragment RT

```
状態空間: フラグメント分解
┌────────────────────┐
│フラグメント1：状態1-40│
│フラグメント2：状態1-40│
│フラグメント3：状態1-40│
│...                |
│フラグメントM：状態1-40│
└────────────────────┘
 dimension (global): n_mat_max = フラグメント数 × nstate_frag
 
 index_basis(io, ifrag, ispin) 
   ↓ (マッピング)
 global_idx
```

---

## 実装詳細

### 関数インタフェース

#### init_dg_fragment_rt

```fortran
subroutine init_dg_fragment_rt(dg_frag, system, rt, info, lg, mg)
  type(s_dg_fragment_rt), intent(inout) :: dg_frag
  type(s_dft_system),     intent(in)    :: system
  type(s_rt),             intent(in)    :: rt
  type(s_parallel_info),  intent(in)    :: info
  type(s_rgrid),          intent(in)    :: lg
  type(s_rgrid),          intent(in)    :: mg
end subroutine
```

**役割**: DG-Fragment RTデータ構造初期化

**主操作**:
1. パラメータ読み込み（yn_dg_fragment_from_dcdft）
2. 基本通信構造設定
3. フラグメント分解情報初期化

#### read_fragment_basis_data

```fortran
subroutine read_fragment_basis_data(dg_frag, system, mg, nxyz_buffer)
  type(s_dg_fragment_rt), intent(inout) :: dg_frag
  type(s_dft_system),     intent(in)    :: system
  type(s_rgrid),          intent(in)    :: mg
  integer,                intent(in)    :: nxyz_buffer
end subroutine
```

**役割**: DC-LCFOデータからフラグメント基底を読み込み

**主操作**:
1. data_dcdft/ディレクトリからメタデータ読み込み
2. 実空間波動関数読み込み
3. **係数行列を正しい次元で割り当て**
4. **フラグメント局所 → グローバル基底マッピング実施**
5. ハロ領域処理

**関数内の重要部分** (修正箇所):

```fortran
! 行751: 条件付きデアロケーション
if (allocated(dg_frag%coef)) deallocate(dg_frag%coef)
if (allocated(dg_frag%coef_new)) deallocate(dg_frag%coef_new)
if (allocated(dg_frag%coef_work)) deallocate(dg_frag%coef_work)

! 行760: 正しい次元で割り当て
allocate(dg_frag%coef(dg_frag%n_mat_max, dg_frag%nstate_tot, dg_frag%nspin))

! 行768-778: 再索引化ループ
do i_local = 1, ifrag_count
  do ispin = 1, dg_frag%nspin
    do io = 1, dg_frag%n_basis(ifrag, ispin)
      global_idx = dg_frag%index_basis(io, ifrag, ispin)
      if (global_idx > 0 .and. global_idx <= dg_frag%n_mat_max) then
        dg_frag%coef(global_idx, :, :) = coef_local(io, :, :, i_local)
      end if
    end do
  end do
end do
```

#### calculate_observables

```fortran
subroutine calculate_observables(dg_frag, system, rt, itt)
  type(s_dg_fragment_rt), intent(inout) :: dg_frag
  type(s_dft_system),     intent(in)    :: system
  type(s_rt),             intent(in)    :: rt
  integer,                intent(in)    :: itt
end subroutine
```

**役割**: 電流とエネルギー計算

**使用配列**:
- `dg_frag%coef(n_mat_max, nstate_tot, nspin)` - 修正後の次元
- `dg_frag%momentum_mat(3, n_mat_max, n_mat_max, nspin)` - グローバル基底

**行2452**:
```fortran
call zgemm('N', 'N', n, nocc, n, (1.0d0, 0.0d0), op_mat, n, &
           dg_frag%coef(1:n, 1:nocc, ispin), n, ...)
```

修正後は `coef(n, ...)` で n = n_mat_max 時、安全にアクセス可能

---

## テスト・検証

### テストシステム

**用途**: 本番品質検証テストケース

**系の特性**:
```
系名:        H₂周期超格子
分子数:      20
電子数:      40
セル定数:    al = 56.0, 20.0, 20.0 (a.u.)
周期性:      x方向周期（y、z非周期）
```

**計数グリッド**:
```
num_rgrid = 120, 40, 40
総ポイント数: 192,000
```

**レーザー場**:
```
形状:       Acos2 エンベロープ
強度:       I = 1.0 × 10¹² W/cm²
周波数:     ω = 0.05 a.u. (~1.36 eV)
偏光:       x方向
期間:       6サイクル
```

**時間積分**:
```
総時間:     1.0 a.u. (~24 fs)
タイムステップ: 20
ステップサイズ: dt = 0.05 a.u. (~1.2 as)
```

### 自動検証スクリプト

#### run_physics_validation.sh

**動作**:
1. 前回出力をクリーン
2. DG-Fragment RTテスト実行
3. Python解析実行
4. 結果報告

**実行方法**:
```bash
cd samples/exercise_dg_rt_hse_test/H2
bash run_physics_validation.sh
```

**所要時間**: 5-10分（ハードウェア依存）

#### physics_validation.py

**機能**:
- 電流データ読み込み
- エネルギーデータ読み込み
- 統計量計算
- DG vs 従来的比較
- 警告出力

**実行結果例**:
```
=== Current Statistics ===
Average Jx:  1.23e-7 a.u.
Max |Jx|:    5.67e-7 a.u.
RMS Current: 3.45e-7 a.u.

=== Energy Evolution ===
Initial:     -17.1799 eV
Final:       -17.1798 eV
Drift:       1.21e-4 eV (0.0007%)
```

---

## トラブルシューティング

### 問題1: SIGBUSエラー（修正済）

**症状**:
```
Program received signal SIGBUS: Access to an undefined portion 
of a memory object
```

**原因**: コード修正前の係数配列次元不一致

**解决**: 修正済み（本実装に含まれる）

**検証方法**:
```bash
/path/to/salmon < inputfile_h2_periodic_20_dg_new_param 2>&1 | grep -i "sigbus\|segmentation"
# 出力がない = 正常
```

### 問題2: コンパイルエラー

**症状**:
```
Error: ...rt_dg_fragment.f90 ...
```

**確認項目**:
1. Fortranコンパイラが新しいか
2. すべてのモジュール依存が満たされてるか
3. CMakeコンパイルか

**解決**:
```bash
cd mybuild
make clean
make -j10
# [100%] Built target salmon が出力される
```

### 問題3: 小さい電流値

**症状**: 電流が予想より1000倍以上小さい

**原因の可能性**:
1. 係数の正規化問題
2. レーザー結合強度エラー
3. 基底完全性問題

**診断**:
```bash
# 1. 係数の大きさ確認
grep "coef(" src/rt/rt_dg_fragment.f90 | head -5

# 2. 正規化確認
# 追加: print *, "<φᵢ|φⱼ> = ", integral

# 3. モーメンタム行列確認
# 追加: print *, "p_ij = ", p_mat(:, io, jo)
```

### 問題4: エネルギードリフト大きい

**許容範囲**:
- 優秀: < 0.01%
- 良好: < 0.1%
- 要改善: > 1%

**改善方法**:
```fortran
! タイムステップを小さく
dt = 0.05d0 → dt = 0.025d0

! または高次積分器を使用
time_integrator = 3  ! RK4
```

---

## パフォーマンス特性

### ベンチマーク設定

**システム**: H₂周期超格子（20分子、40電子）

**計算構成**:
```
基底関数:     40個/フラグメント × 複数フラグメント
グリッド:     120 × 40 × 40 = 192,000ポイント
積分器:       Runge-Kutta
```

### 典型実行時間

**1タイムステップ**:
- 初期化: 1-2秒
- ハミルトニアン計算: 1-2秒
- 時間伝播: 1-2秒
- 観測量計算: 0.5秒

**全テスト（20ステップ）**:
```
GS計算:       30-60秒（DC-LCFO）
RT初期化:     2秒
時間発展:     30-40秒（20ステップ）
全体:         60-100秒
```

※マシンスペック依存

### メモリ使用量

```
phi_frag配列:      ~100 MB
coef配列:          ~10 MB
momentum_mat:      ~5 MB
作業配列:          ~50 MB
────────────────────────
合計:              ~200 MB
```

---

## 出力ファイル形式

### H2_periodic_20_dg_new_param_rt.data

**内容**: 時間発展中の場データ

**列構成** (全16列):
```
 1: Time[a.u.]
 2-4: Ac_ext_x/y/z[a.u.] - 外部ベクトルポテンシャル
 5-7: E_ext_x/y/z[a.u.]  - 外部電場
 8-10: Ac_tot_x/y/z[a.u.] - 全ベクトルポテンシャル
11-13: E_tot_x/y/z[a.u.]  - 全電場
14-16: Jm_x/y/z[a.u.]      - 物質電流（電子）
```

**例**:
```
      0.00000000  0.107513E-004  0.000000E+000 ...  0.190936E-007
      0.05000000  0.422526E-004  0.000000E+000 ...  0.748681E-007
      0.10000000  0.933338E-004  0.000000E+000 ...  0.164921E-006
```

### H2_periodic_20_dg_new_param_rt_energy.data

**内容**: エネルギー発展

**列構成** (全3列):
```
1: Time[a.u.]
2: Eall[a.u.]        - 全エネルギー
3: Eall-Eall0[a.u.] - エネルギード リフト
```

**例**:
```
      0.00000000  -0.171798552382851E+002  0.000000000000000E+000
      0.05000000  -0.171798551171510E+002  0.121134100000000E-004
      0.10000000  -0.171798550960169E+002  0.242268200000000E-004
```

---

## 参考資料

### 関連論文・手法

- DC-LCFO: Density Cumulant Local Coherence Optimization
- DG-Fragment: Discontinuous Galerkin Fragment decomposition
- TDDFT: Time-Dependent Density Functional Theory

### コード参照

**主要ファイル**:
- src/rt/rt_dg_fragment.f90 - DG-Fragment実装
- src/rt/main_tddft.f90 - TDDFT時間発展メイン
- src/structures.f90 - データ構造定義

**補助モジュール**:
- src/communication.f90 - MPI通信
- src/hamiltonian.f90 - ハミルトニアン計算
- src/propagation.f90 - 時間積分

---

## 今後の拡張予定

1. **適応基底**（yn_adaptive_basis）
   - フラグメント動的最適化
   - 基底関数の適応的追加・削除

2. **フラグメント結合**
   - フラグメント間相互作用
   - 密度マトリックス共有

3. **非線形応答**
   - 高調波発生（HHG）
   - 非線形感受率計算

4. **並列化最適化**
   - GPU加速
   - スケーラ性改善

---

**文書作成**: 2026年2月23日  
**バージョン**: 1.0 (Final)  
**対象ユーザー**: 開発者、研究者

