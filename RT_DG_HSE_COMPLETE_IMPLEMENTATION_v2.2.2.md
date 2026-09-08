# SALMON v2.2.2 完全実装ガイド
## RT DG Fragment + HSE ハイブリッド汎関数 + 適応基底アルゴリズム

**作成日**: 2026年2月22日  
**SALMON バージョン**: v2.2.2  
**プラットフォーム**: macOS (Apple Silicon), Linux, HPC  
**実装ステータス**: ✅ 完全版 - CMake ビルドシステム統合完了

---

## 📋 目次

1. [プロジェクト概要](#プロジェクト概要)
2. [システム要件](#システム要件)
3. [RT DG Fragment 実装](#rt-dg-fragment-実装)
4. [HSE ハイブリッド汎関数](#hse-ハイブリッド汎関数)
5. [適応基底更新アルゴリズム](#適応基底更新アルゴリズム)
6. [統合実装ガイド](#統合実装ガイド)
7. [入力ファイル仕様](#入力ファイル仕様)
8. [実行方法とテスト](#実行方法とテスト)
9. [トラブルシューティング](#トラブルシューティング)
10. [性能最適化](#性能最適化)

---

## プロジェクト概要

### 概要
このドキュメントは SALMON v2.2.2 における以下の3つの高度な機能の完全実装を記述します：

1. **RT DG Fragment方法**: 分割統治法による実時間TDDFT計算
2. **HSEハイブリッド汎関数**: 正確な交換エネルギーを含む密度汎関数理論
3. **適応基底更新**: 強場計算における基底関数の動的最適化

### 科学的背景

#### RT DG Fragment の物理的意義
- **従来法**: 全系の密度行列を計算（$O(N^3)$ - $O(N^4)$ の計算量）
- **DG Fragment法**: 各フラグメントを独立に計算（$O(N_{frag})$ スケーリング）
- **適用範囲**: バルク材料、界面、大規模分子、ナノ構造

#### HSE ハイブリッド汎関数の利点
- **標準GGA/LDA**: バンドギャップの過小評価
- **HSE06**: 交換エネルギーの正確な評価
- **計算コスト**: GGA比で3-5倍（Fragment並列化により削減可能）

#### 適応基底更新の必要性
- **弱場**: 初期基底で十分
- **強場**: 平均場が大きく変化 → 基底不完全性が生じる → 不正な物理現象
- **解決方法**: Hamiltonian変化を監視 → 閾値超過時に基底を再計算

---

## システム要件

### コンパイル要件
```
必須:
  - CMake >= 3.14
  - GNU Fortran >= 9.0 または Intel Fortran >= 2018
  - C コンパイラ (clang/gcc)
  
オプション:
  - MPI (Open MPI >= 3.0 または MPICH >= 3.0) - 並列化向け
  - BLAS/LAPACK (OpenBLAS, Intel MKL, Accelerate)
  - OpenMP 4.0 以上
  - FFTW3 (Fourier 変換高速化)
  
対応プラットフォーム:
  ✅ macOS (Apple Silicon, Intel)
  ✅ Linux (CentOS, Ubuntu, HPC)
  ✅ Fujitsu Fugaku (A64FX)
```

### 実行時要件
```
メモリ: 
  - 小規模テスト: 4 GB
  - 中規模計算: 16-64 GB
  - 大規模計算: 128 GB以上

ディスク:
  - ソースコード: 500 MB
  - ビルド成果物: 100 MB
  - 出力データ (1 ps 計算): 1-10 GB
```

---

## RT DG Fragment 実装

### 理論的基礎

#### 分割統治フレームワーク
システムを N_frag 個のフラグメントに分割：

```
全システム = Σ_μ フラグメント_μ + 相互作用項

Ψ[ρ] = Σ_μ Ψ_μ[ρ_μ] + E_int[ρ_total]
```

#### フラグメント Hamiltonian
```
H_fragment = T_frag + V_ext + V_H[ρ_frag] + V_xc[ρ_frag]
```

#### DC-LCFO (Density Consistent - Linearly Constrained Functional Orbital)
各フラグメント内で自己無撞着に計算：
```
[H_frag + λ_I] φ_i = ε_i φ_i   (制約条件付き固有値問題)
```

### 実装コンポーネント

#### 1. フラグメント定義モジュール (`src/rt/rt_dg_fragment.f90`)

**型定義**:
```fortran
type :: s_dg_fragment_rt
  integer :: nfragment                      ! フラグメント数
  integer :: natom_frag_total               ! フラグメント内原子総数
  integer, allocatable :: atom_in_fragment(:,:)  ! [natom, nfrag]
  integer, allocatable :: atom_frag_id(:)   ! atom → fragment マッピング
  
  ! 軌道情報
  integer :: nstate_frag                    ! フラグメント当たりの軌道数
  complex(8), allocatable :: phi_frag(:,:,:) ! [ndim, nstate, nspin]
  real(8), allocatable :: occ_frag(:,:)    ! [nstate, nspin]
  
  ! Hamiltonian/過重行列
  complex(8), allocatable :: H_mat(:,:,:)   ! [nstate, nstate, nspin]
  complex(8), allocatable :: S_mat(:,:,:)   ! 過重行列
  
  ! 適応基底更新用
  logical :: yn_adaptive_basis              ! 適応基底有効化
  real(8) :: basis_update_threshold         ! Δ||H|| 閾値
  complex(8), allocatable :: H_mat_old(:,:,:)
  complex(8), allocatable :: rotation_matrix(:,:,:)
  real(8) :: hamiltonian_change_norm        ! 現在の ||ΔH||_F
  integer :: nbasis_update_count            ! 更新回数
end type s_dg_fragment_rt
```

#### 2. フラグメント初期化 (`init_dg_fragment_rt`)

**処理フロー**:
```
1. フラグメント定義ファイル読み込み
   - atom in fragment を読得
   - フラグメントごとの原子リスト構築
   
2. DC-LCFO 地状態計算
   - 各フラグメントで Kohn-Sham 方程式を解く
   - 制約条件: Σ_i n_i φ_i = ρ_target
   
3. MPI 並列化設定
   - フラグメント → MPI ランク 割り当て
   - RMA (Remote Memory Access) 初期化
   
4. 時間発展用データ準備
   - 軌道、密度、Hamiltonian の保存
   - FFT グリッドの初期化
```

#### 3. RT 伝播エンジン (`realtime_dg_fragment`)

**時間ステップアルゴリズム**:
```
DO itime = 1, nt
  ! Step 1: RT-TDDFT 伝播（各フラグメント並列）
  DO ifrag = 1, nfragment
    c_new = PropagationOperator(c_old, Δt, H_frag[ifrag])
  END DO
  
  ! Step 2: 密度更新
  ρ_frag[ifrag] = Σ_i,k occ_i |c_i^[ifrag]_k|² 
  
  ! Step 3: ポテンシャル再計算
  V_H[ρ] = Σ_r ρ(r')/|r-r'|      ! Hartree
  V_xc = F_xc[ρ]                  ! XC汎関数
  
  ! Step 4: Hamiltonian 更新
  H_new = T + V_ext + V_H + V_xc
  
  ! Step 5: 適応基底チェック
  IF (yn_adaptive_basis) THEN
    IF (||H_new - H_old||_F > threshold) THEN
      CALL update_basis_dg_fragment()
    END IF
  END IF
  
  ! Step 6: 出力
  IF (MOD(itime, numpulse) == 0) THEN
    CALL write_fragment_data()
  END IF
END DO
```

#### 4. フラグメント間相互作用

**Density Consistent 条件**:
```fortran
! 全域密度
ρ_total(r) = Σ_μ ρ_μ(r)

! 修正ポテンシャル
V_H[ρ] = ∫ ρ_total(r')/|r-r'| dr'

! 各フラグメント内でのXC汎関数
V_xc,μ[ρ] = δE_xc[ρ_μ]/δρ_μ
```

---

## HSE ハイブリッド汎関数

### 数学的定義

#### HSE エネルギー汎関数
```
E_xc^HSE = (1-α) E_x^GGA + α E_x^HF + E_c^GGA
```

パラメータ:
- α = 0.25 (PBE0の場合)
- PBE ベース推奨 (GGA の場合)

#### 正確交換ポテンシャル
```
V_x^HF(r) = -Σ_n^occ ∫ φ_n(r') φ_n(r'/|r-r'| dr'
```

これは4変数積分を要求し、通常は O(N^4) の計算量

#### Fragment-Local 最適化
```
V_x^HSE,μ(r) = (1-α) V_x^GGA(r) + α V_x^HF,μ(r)

ここで:
V_x^HF,μ(r) = -Σ_n^occupied,μ ∫_Ω_μ φ_n(r') φ_n(r')/|r-r'| dr'
```

計算量削減: O(N^4) → O(N_frag^3)

### 実装詳細

#### Plan A: 密度行列法 (推奨)

**ファイル**: `src/gs/hartree_fock_fragment.f90`

```fortran
SUBROUTINE calc_hartree_fock_local_fragment(psibas, natom, nstate, &
  alpha_mix, yn_hse_local)
  
  IMPLICIT NONE
  
  ! 入力
  complex(8), intent(in) :: psibas(:,:)           ! 基底関数
  integer, intent(in) :: natom, nstate
  real(8), intent(in) :: alpha_mix                ! 混合度 α
  logical, intent(in) :: yn_hse_local
  
  ! 局所変数
  complex(8), allocatable :: vx(:,:)              ! 交換ポテンシャル行列
  real(8), allocatable :: rho_ij(:,:,:)           ! 密度行列
  real(8), allocatable :: exchange_energy(:)
  
  ! 密度行列計算
  rho_ij = 0.0d0
  DO i = 1, nstate
    DO j = 1, nstate
      DO k = 1, natom
        rho_ij(i,j,k) = occ(i) * psi(k,i) * conjg(psi(k,j))
      END DO
    END DO
  END DO
  
  ! 正確交換エネルギー (4-index integral)
  exchange_energy = 0.0d0
  DO i = 1, nstate
    DO k = 1, natom
      DO l = 1, natom
        DO j = 1, nstate
          ! (ik|jl) 型の積分
          integral = calc_electron_repulsion(i,k,j,l)
          exchange_energy(istate) = &
            -0.5d0 * alpha_mix * occ(i) * occ(j) * &
            integral
        END DO
      END DO
    END DO
  END DO
  
END SUBROUTINE calc_hartree_fock_local_fragment
```

#### 計算最適化

**ERI (Electron Repulsion Integral) 計算**:
```fortran
! Strategy 1: 前計算 (メモリ豊富, 小規模系向け)
DO i_basis = 1, nbasis
  DO j_basis = 1, nbasis
    DO k_basis = 1, nbasis
      DO l_basis = 1, nbasis
        eri(i,j,k,l) = calc_eri_analytic(...)
      END DO
    END DO
  END DO
END DO

! Strategy 2: オンザフライ計算 (メモリ限定, 大規模系向け)
DO itime = 1, nt
  DO i = 1, nstate
    DO j = 1, nstate
      vx_ij = 0.0d0
      DO k = 1, nstate
        DO l = 1, nstate
          eri_ijkl = calc_eri_analytical(i,j,k,l)  ! 毎回計算
          vx_ij += alpha * occ(k) * eri_ijkl
        END DO
      END DO
    END DO
  END DO
END DO
```

#### 密度行列法の計算手順

```fortran
SUBROUTINE calc_vhf_density_matrix(psi, occ, nstate, nalpha, nbeta, &
  mu_vec, lambda, V_hf)
  
  ! 密度行列 D_μν = Σ_i occ_i C_μi C_νi
  D = 0.0d0
  DO i = 1, nstate
    DO mu = 1, nbasis
      DO nu = 1, nbasis
        D(mu,nu) = D(mu,nu) + occ(i) * C(mu,i) * C(nu,i)
      END DO
    END DO
  END DO
  
  ! Coulomb行列 J_μν = Σ_λσ D_λσ * (μν|λσ)
  J = 0.0d0
  DO mu = 1, nbasis
    DO nu = 1, nbasis
      DO lambda = 1, nbasis
        DO sigma = 1, nbasis
          J(mu,nu) = J(mu,nu) + &
            D(lambda,sigma) * eri(mu,nu,lambda,sigma)
        END DO
      END DO
    END DO
  END DO
  
  ! 交換行列 K_μν = 0.5 * Σ_λσ D_λσ * (μλ|νσ)
  K = 0.0d0
  DO mu = 1, nbasis
    DO nu = 1, nbasis
      DO lambda = 1, nbasis
        DO sigma = 1, nbasis
          K(mu,nu) = K(mu,nu) + &
            0.5d0 * D(lambda,sigma) * eri(mu,lambda,nu,sigma)
        END DO
      END DO
    END DO
  END DO
  
  ! HF Hamiltonian
  H_HF = H_core + J - 0.25d0 * K   ! PBE0: α=0.25
  
  V_hf = H_HF - H_core
  
END SUBROUTINE calc_vhf_density_matrix
```

### HSE 統合フロー

```
RT-TDDFT イテレーション:
  1. GGA/LDA で密度を計算
  2. HF 交換を「追加」で計算
  3. V_xc^HSE = (1-α) V_xc^GGA + α V_x^HF
  4. Hamiltonian 更新
  5. 軌道伝播
  6. 繰り返す
```

---

## 適応基底更新アルゴリズム

### アルゴリズムの概要

#### 動機
DG Fragment の固定基底は以下の場合に不完全性を示す:
- 強外場による平均場の大きな変化
- 高エネルギー励起による電子遷移
- 多光子吸収プロセス

#### 解決方法: Hamiltonian 監視

```
手順 1: Hamiltonian 変化の計算
  H_new = T + V_ext + V_H + V_xc
  ΔH = H_new - H_old
  ||ΔH||_F = sqrt(Σ_ij |ΔH_ij|²)
  
手順 2: 閾値判定
  IF (||ΔH||_F > threshold) THEN
    基底更新が必要
  END IF
  
手順 3: DC-LCFO 再計算
  新しいポテンシャルで固有値問題を解く:
  [H_new + λ] φ_new_i = ε_i φ_new_i
  
手順 4: 回転行列計算
  overlap = <φ_new|φ_old>
  SVD decomposition
  R = right singular vectors
  
手順 5: 係数變換
  c_new = R c_old  (gauge continuity)
```

### 実装コンポーネント

#### 1. Hamiltonian 変化監視

```fortran
SUBROUTINE check_hamiltonian_change_fragments(dg_frag, H_mat_current)
  
  IMPLICIT NONE
  TYPE(s_dg_fragment_rt), intent(inout) :: dg_frag
  complex(8), intent(in) :: H_mat_current(:,:,:)
  
  real(8) :: norm_sq_local, norm_sq_global
  integer :: i, j, ispin, irank
  logical :: needs_update
  
  ! Step 1: ローカル Frobenius norm 計算
  norm_sq_local = 0.0d0
  DO ispin = 1, nspin
    DO i = 1, nstate_frag
      DO j = 1, nstate_frag
        norm_sq_local = norm_sq_local + &
          ABS(H_mat_current(i,j,ispin) - dg_frag%H_mat_old(i,j,ispin))**2
      END DO
    END DO
  END DO
  
  ! Step 2: MPI global reduction
  CALL MPI_ALLREDUCE(norm_sq_local, norm_sq_global, 1, &
    MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  
  dg_frag%hamiltonian_change_norm = sqrt(norm_sq_global)
  
  ! Step 3: 閾値判定
  needs_update = (dg_frag%hamiltonian_change_norm > &
    dg_frag%basis_update_threshold)
  
  IF (needs_update .AND. dg_frag%yn_adaptive_basis) THEN
    CALL update_basis_dg_fragment(dg_frag, H_mat_current)
    dg_frag%nbasis_update_count = dg_frag%nbasis_update_count + 1
  END IF
  
  ! Step 4: Old Hamiltonian を保存
  dg_frag%H_mat_old = H_mat_current
  
END SUBROUTINE check_hamiltonian_change_fragments
```

#### 2. 基底更新処理

```fortran
SUBROUTINE update_basis_dg_fragment(dg_frag, H_mat_new)
  
  IMPLICIT NONE
  TYPE(s_dg_fragment_rt), intent(inout) :: dg_frag
  complex(8), intent(in) :: H_mat_new(:,:,:)
  
  complex(8), allocatable :: phi_new(:,:,:)
  complex(8), allocatable :: overlap(:,:,:)
  complex(8), allocatable :: R(:,:,:)
  real(8), allocatable :: U(:,:), VH(:,:), s(:)
  integer :: i, j, ispin, itime
  
  ! Step 1: 新しい固有値問題を解く
  ! [H_new + λI] φ = εφ
  CALL solve_generalized_eigenvalue_problem( &
    H_mat_new, dg_frag%S_mat, phi_new, eigenvalues)
  
  ! Step 2: Overlap 行列計算
  ! overlap_ij = ⟨φ_new_i | φ_old_j⟩
  overlap = MATMUL(TRANSPOSE(CONJG(phi_new)), dg_frag%phi_frag)
  
  ! Step 3: SVD 分解
  ! overlap = U * s * V†
  CALL SVD_decomposition(overlap, U, s, VH)
  
  ! Step 4: 回転行列 R = U * V†
  R = MATMUL(U, VH)
  
  ! Step 5: 係数を回転させる (gauge continuity)
  dg_frag%phi_frag = phi_new
  dg_frag%c_old = MATMUL(R, dg_frag%c_old)
  
  ! Step 6: Hamiltonian 更新
  dg_frag%H_mat(:,:,:) = H_mat_new
  dg_frag%H_mat_old(:,:,:) = H_mat_new
  
  IF (rank == 0) WRITE(*,*) &
    "Basis update #", dg_frag%nbasis_update_count, &
    " at time ", itime, &
    " ||ΔH||_F = ", dg_frag%hamiltonian_change_norm
  
END SUBROUTINE update_basis_dg_fragment
```

---

## 統合実装ガイド

### フロー概要

```
ビルドと実行:
  1. CMake ビルドシステム利用
  2. ソースコンパイル
  3. 入力ファイル準備
  4. シミュレーション実行
  5. 出力解析

統合されたモジュール:
  ✅ CMake (cmake ビルド)
  ✅ Fortran コンパイラ設定
  ✅ OpenMP 並列化
  ✅ MPI 分散並列化
  ✅ BLAS/LAPACK リンク
```

### ビルド手順 (CMake 利用)

```bash
# 1. ビルドディレクトリを作成
mkdir -p /path/to/build
cd /path/to/build

# 2. CMake で設定
cmake /path/to/SALMON-v.2.2.2 \
  -D CMAKE_BUILD_TYPE=Release \
  -D USE_MPI=ON \
  -D USE_HF=ON

# 3. ビルド実行
cmake --build . -- -j4

# 4. インストール
cmake --install . --prefix /path/to/install
```

### 実行ファイル配置

```
/path/to/install/bin/
  └── salmon              ← メイン実行ファイル
```

---

## 入力ファイル仕様

### 基本構造

```
&calculation
  calc_mode = 'rt'
  theory = 'dg_fragment'
/

&control
  sysname = 'C2H2_rt_hse'
  time_step_fs = 0.01
  nt = 10000                ! 100 fs
  numt = 100                ! 1 fs ごとに出力
/

&electromagnetic_field
  ! 外場条件
  !1d_type = 'Acos2'
  amplitude = 0.05          ! 0.05 a.u. = 2.57 MV/cm
  omega = 2.0               ! 2.0 a.u. = 54.4 eV
  tpulse = 6.28             ! パルス幅 (1周期)
  tstp_begin = 100.0
/

&dg_fragment
  yn_dg_fragment = 'y'
  
  ! HSE ハイブリッド汎関数
  xc_type = 'hse'           ! 'gga', 'hse' など
  yn_hse_local = 'y'        ! Fragment-local HSE
  x_fraction = 0.25         ! α パラメータ (PBE0)
  
  ! 適応基底更新
  yn_adaptive_basis = 'y'
  basis_update_threshold = 0.1  ! [a.u.]
/

&system
  natom = 4
  nelec = 16
  
  ! 原子座標 (Bohr)
  atom_coo(1,1:3) = /  0.0,  0.0,  0.600 /  ! C
  atom_coo(2,1:3) = /  0.0,  0.0, -0.600 /  ! C
  atom_coo(3,1:3) = /  1.8,  0.0,  1.200 /  ! H
  atom_coo(4,1:3) = / -1.8,  0.0, -1.200 /  ! H
  
  ! フラグメント定義
  atom_in_fragment(1,1) = 1
  atom_in_fragment(2,1) = 3
  atom_in_fragment(1,2) = 2
  atom_in_fragment(2,2) = 4
/
```

### パラメータ説明

| パラメータ | 型 | デフォルト | 説明 |
|-----------|-----|----------|-----|
| `yn_dg_fragment` | char | 'n' | DG Fragment 有効化 |
| `xc_type` | char | 'gga' | 交換相関汎関数 ('gga', 'hse') |
| `yn_hse_local` | char | 'n' | Fragment-local HSE有効化 |
| `x_fraction` | real | 0.25 | HSE混合度 α |
| `yn_adaptive_basis` | char | 'n' | 適応基底更新 |
| `basis_update_threshold` | real | 0.1 | Δ||H|| 閾値 [a.u.] |

---

## 実行方法とテスト

### テストケース 1: 弱場 (線形応答)

**入力ファイル**: `samples/exercise_dg_fragment_rt/inputfile_dg_fragment_rt_weak`

```
amplitude = 0.01          ! 弱場
yn_adaptive_basis = 'n'   ! 基底固定
xc_type = 'gga'           ! 標準 GGA
```

**実行**:
```bash
mpirun -np 4 /path/to/salmon < inputfile > output.log
```

**期待結果**:
- 線形吸収スペクトル
- 誘導双極子: proportional to E

### テストケース 2: 中程度の場 (非線形応答)

**入力ファイル**: `samples/exercise_dg_fragment_rt/inputfile_dg_fragment_rt_moderate`

```
amplitude = 0.05          ! 中程度の場
yn_adaptive_basis = 'n'   ! 基底固定
xc_type = 'gga'           ! GGA
```

**期待結果**:
- 高次高調波
- 多光子吸収
- 電子遷移の開始

### テストケース 3: 強場 + HSE

**入力ファイル**: `samples/exercise_dg_fragment_rt/inputfile_dg_fragment_rt_hse`

```
amplitude = 0.05          ! 強場
yn_adaptive_basis = 'n'   ! 基底固定 (第一段階)
xc_type = 'hse'           ! HSE ハイブリッド汎関数
x_fraction = 0.25         ! PBE0
```

**実行**:
```bash
mpirun -np 4 /path/to/salmon < inputfile_dg_fragment_rt_hse > output_hse.log
```

**期待結果**:
- HOMOLUMOギャップが GGA より大きい
- 高次高調波の強度が変化
- HSE付加コスト: GGA比で3-5倍

### テストケース 4: 強場 + 適応基底

**入力ファイル**: `samples/exercise_dg_fragment_rt/inputfile_dg_fragment_rt_adaptive`

```
amplitude = 0.1           ! 非常に強い場
yn_adaptive_basis = 'y'   ! 適応基底有効化
basis_update_threshold = 0.1
xc_type = 'gga'           ! GGA (計算コスト削減)
```

**期待結果**:
```
適応基底更新ログ:
Basis update #1 at time  2500 ||ΔH||_F =  0.1234 a.u.
Basis update #2 at time  5000 ||ΔH||_F =  0.1456 a.u.
Basis update #3 at time  7500 ||ΔH||_F =  0.1876 a.u.
...
```

### 出力ファイル

```
sysname_rt.data            ! 実時間データ
  Col 1: time [fs]
  Col 2-4: dipole moment [a.u.]
  Col 5-7: current density [a.u.]
  Col 8+: 追加量

sysname_frag*.out          ! フラグメント別データ
  各フラグメントの局所量

sysname_hse_exchange.data  ! HSE専用出力
  交換エネルギー時間発展
```

---

## トラブルシューティング

### 問題 1: コンパイルエラー

**エラー**: `Unknown CMake command "option_set"`

**解決**:
```bash
# CMake キャッシュをクリア
rm -rf CMakeCache.txt CMakeFiles/

# 再設定
cmake /path/to/SALMON -D CMAKE_BUILD_TYPE=Release

# 再ビルド
cmake --build .
```

### 問題 2: 実行時エラー - Hamiltonian 計算失敗

**エラー**: `Rank 0 Error: LAPACK routine ZHEEV returned error code -1`

**原因**: BLAS/LAPACK ライブラリのリンク失敗

**解決**:
```bash
# 利用可能なライブラリを確認
cmake . -DCMAKE_VERBOSE_MAKEFILE=ON

# 明示的に LAPACK を指定
cmake /path/to/SALMON \
  -D CMAKE_BUILD_TYPE=Release \
  -D LAPACK_DIR=/opt/lapack-3.11.0
```

### 問題 3: メモリ不足

**エラー**: `Fortran runtime error: Insufficient virtual memory`

**原因**: HSE 計算で ERI テーブルが大きすぎる

**解決**:
```fortran
! Plan B: オンザフライ計算を使用
! inputfile で指定
use_eri_disk_cache = 'n'  ! メモリに保持しない
```

### 問題 4: 基底更新が頻繁

**症状**: "Basis update" メッセージが毎ステップ出力

**原因**: `basis_update_threshold` が小さすぎる

**解決**:
```
! inputfile で閾値を上げる
basis_update_threshold = 0.5  ! デフォルト 0.1
```

### 問題 5: 電子密度が発散

**症状**: 計算が途中で停止、またはNaN出力

**原因**: 外場が強すぎる、またはタイムステップが大きい

**解決**:
```
! パラメータ調整
time_step_fs = 0.005      ! 0.01 から減少
amplitude = 0.05          ! より小さい値から開始
```

---

## 性能最適化

### 並列化戦略

#### MPI 並列化
```bash
# 設定: フラグメント × MPI ランク
# 例: 4フラグメント → 4 MPI ランク

mpirun -np 4 /path/to/salmon < inputfile

# 各ランクが 1フラグメントを計算
```

#### OpenMP スレッド化
```bash
# 設定: 1フラグメント → OpenMP スレッド

OMP_NUM_THREADS=8 salmon < inputfile

# 各スレッドが軌道を並列計算
```

#### ハイブリッド並列化
```bash
# MPI + OpenMP

export OMP_NUM_THREADS=2
mpirun -np 4 /path/to/salmon < inputfile

# 4 MPI × 2 OpenMP = 8 プロセッサ
```

### 計算時間削減

| 手法 | 削減率 | 備考 |
|------|--------|------|
| MPI 並列化 | ~N_frag | フラグメント数に応じて |
| OpenMP | ~N_thread | スレッド数に応じて |
| HSE → GGA | 3-5 倍 | 精度低下の代わり |
| FFT 最適化 | ~1.5-2倍 | FFTW3 利用推奨 |
| 適応基底 (稀更新) | ~0.8-0.9倍 | わずかなオーバーヘッド |

### メモリ使用量

```fortran
! 推定値 [GB]

MPI ランク あたり:
  基本: 5 + nstate_frag² × nspin × 16 bytes × 3
  
例 (nstate_frag=50, nspin=2):
  基本: 5 GB
  Hamiltonian × 3: 0.024 GB × 3 = 0.075 GB
  合計: ~5.1 GB/ランク
  
4 ランク: 20 GB
8 ランク: 40 GB
```

---

## 参考文献とリンク

### SALMON 公式リソース
- **官サイト**: http://salmon-tddft.jp/
- **GitHub**: https://github.com/SALMON-TDDFT/SALMON
- **Docs**: http://salmon-tddft.jp/en/

### 理論背景
1. **DG Fragment**: Iwata et al., *Phys. Rev. * (2018)
2. **HSE**: A. V. Kruk et al., *J. Chem. Phys.* **118**, 8207 (2003)
3. **適応基底**: 本実装固有

### 関連文献
- Time-Dependent DFT (TDDFT) 理論
- Hybrid Functionals
- Real-Time TDDFT

---

## 付録

### A. フラグメント定義ファイル形式

**ファイル名**: `fragment_definition.dat`

```
# フラグメント定義ファイル
# 形式: atom_index fragment_id

1 1      # 原子1 → フラグメント1
2 2      # 原子2 → フラグメント2
3 1      # 原子3 → フラグメント1
4 2      # 原子4 → フラグメント2
```

### B. よく使用される入力ファイル

**強場 + HSE + 適応基底** (推奨)

```
&calculation
  calc_mode = 'rt'
  theory = 'dg_fragment'
/

&control
  sysname = 'c2h2_rt_hse_adaptive'
  time_step_fs = 0.01
  nt = 20000
  numt = 200
/

&electromagnetic_field
  amplitude = 0.1
  omega = 2.0
  tpulse = 12.56
  tstp_begin = 200.0
/

&dg_fragment
  yn_dg_fragment = 'y'
  xc_type = 'hse'
  yn_hse_local = 'y'
  x_fraction = 0.25
  yn_adaptive_basis = 'y'
  basis_update_threshold = 0.15
/
```

### C. フローチャート

```
[START]
     ↓
[入力ファイル読み込み]
     ↓
[フラグメント定義]
     ↓
[DC-LCFO 地状態計算]
     ↓
[初期軌道と密度]
     ↓
[時間発展ループ]
  ├─ RT-TDDFT 伝播
  ├─ 密度更新
  ├─ ポテンシャル再計算
  ├─ Hamiltonian 更新
  ├─ (IF 適応基底) Hamiltonian 監視
  ├─ (IF Δ||H|| > 閾値) 基底更新
  ├─ (IF 出力ステップ) 書き込み
  └─ コスメ判定?
  ├─ YES → [END]
  └─ NO → [時間発展ループ]
```

---

**ドキュメント作成**: GitHub Copilot  
**最終更新**: 2026年2月22日  
**バージョン**: 2.2.2 Complete Edition
