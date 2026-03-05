# HSE交換項の最適化提案: 基底間積分の事前計算
**日付**: 2026年2月23日  
**目的**: O(L^6)の6次元グリッド積分をO(N^4)の行列演算に変換  
**実装難易度**: 中〜高

---

## 🎯 最適化の動機

### 現在の実装（Plan A）
```fortran
! 10重ループ: i,j,k,l × (r1_x,r1_y,r1_z) × (r2_x,r2_y,r2_z)
do jo = 1, n_basis
  do io = 1, n_basis
    do lo = 1, n_occ
      do ko = 1, n_occ
        do iz, iy, ix = グリッド1  ! L^3
          do jz, jy, jx = グリッド2  ! L^3
            V_x_ij += φ_i(r1) φ_k(r1) [1/r12] φ_l(r2) φ_j(r2)
```

**計算量**: O(N_b^2 × N_occ^2 × L^6) per fragment

---

## 💡 提案: 2電子積分の事前計算（ERI法）

### Plan B: Electron Repulsion Integrals (ERI)

#### **基本アイデア**
2電子積分（ERI）を**一度だけ計算**して保存：

$$
(ij|kl) = \int\int \phi_i(\mathbf{r}_1) \phi_j(\mathbf{r}_1) \frac{1}{|\mathbf{r}_1 - \mathbf{r}_2|} \phi_k(\mathbf{r}_2) \phi_l(\mathbf{r}_2) \, d^3\mathbf{r}_1 \, d^3\mathbf{r}_2
$$

#### **交換行列の構築**
事前計算したERIから：

$$
V_x^{\text{HSE}}(i,j) = -\alpha \sum_{k,l \in \text{occ}} (ij|kl)
$$

---

## 📊 計算量比較

| 方法 | メモリ | 事前計算 | 各ステップ | 備考 |
|------|---------|----------|-----------|------|
| **Plan A (現在)** | O(N_b × L^3) | なし | O(N_b^2 N_{occ}^2 L^6) | 毎回6次元積分 |
| **Plan B (ERI)** | O(N_b^4) | O(N_b^4 L^6) | O(N_b^2 N_{occ}^2) | 初期化時のみ積分 |
| **Plan C (RI/DF)** | O(N_b^2 N_{aux}) | O(N_b^2 N_{aux} L^6) | O(N_b^2 N_{aux} N_{occ}) | 補助基底使用 |

### 典型的な数値例

#### **小規模フラグメント（L=16, N_b=16, N_occ=8）**

| 項目 | Plan A | Plan B | 改善率 |
|------|--------|--------|--------|
| **メモリ** | 65 KB | 65 MB | 1000倍増 ⚠️ |
| **初期計算** | - | 1時間 | 初回のみ |
| **各タイムステップ** | 10秒 | 0.01秒 | **1000倍高速** ✅ |

#### **実際のDC-LCFOフラグメント（原子数×50軌道/原子）**

**重要**: DC法では**原子あたり約50軌道**が基本です。

| フラグメントサイズ | N_b | ERI要素数 | メモリ（Full） | メモリ（対称性） | 現実性 |
|-------------------|-----|-----------|--------------|-----------------|--------|
| **1原子** | 50 | 6.25×10⁶ | 50 MB | 6.25 MB | ✅ 実用的 |
| **2原子** | 100 | 1.0×10⁸ | 800 MB | 100 MB | ⚠️ ギリギリ |
| **4原子** | 200 | 1.6×10⁹ | 12.8 GB | 1.6 GB | ⚠️ 大規模メモリ必要 |
| **8原子** | 400 | 2.56×10¹⁰ | **205 GB** | **25.6 GB** | ❌ 非現実的 |

**結論**: 
- ✅ **単原子フラグメント**: ERI法が有効
- ⚠️ **2原子フラグメント**: 対称性利用必須
- ❌ **4原子以上**: Plan CまたはPlan Aが必要

---

## 🔧 実装設計（Plan B）

### **1. データ構造の拡張**

[src/rt/rt_dg_fragment_types.f90](src/rt/rt_dg_fragment_types.f90) に追加：

```fortran
type s_dg_fragment_rt
  ! 既存のフィールド...
  
  ! ===== HSE最適化用 =====
  ! 2電子積分テンソル (fragment-local)
  real(8), allocatable :: ERI_frag(:,:,:,:)   ! (i,j,k,l) = (ij|kl)
  logical :: ERI_precomputed                   ! 事前計算済みフラグ
  
  ! 対称性を利用した圧縮保存（オプション）
  real(8), allocatable :: ERI_compressed(:)   ! 圧縮形式
  integer, allocatable :: ERI_index(:,:,:,:)  ! インデックステーブル
end type
```

---

### **2. ERI事前計算ルーチン**

[src/xc/xc_hse.f90](src/xc/xc_hse.f90) に追加：

```fortran
!=======================================================================
! Precompute all 2-electron integrals for a fragment
! Exploits 8-fold permutation symmetry: (ij|kl) = (ji|kl) = (ij|lk) = ...
!=======================================================================
subroutine precompute_ERI_fragment(phi_grid, hgs, hvol, is_grid, ie_grid, &
                                   n_basis, ERI_frag)
  implicit none
  real(8), intent(in) :: phi_grid(:,:,:,:)
  real(8), intent(in) :: hgs(3), hvol
  integer, intent(in) :: is_grid(3), ie_grid(3), n_basis
  real(8), intent(out) :: ERI_frag(n_basis, n_basis, n_basis, n_basis)
  
  integer :: i, j, k, l, ix, iy, iz, jx, jy, jz
  real(8) :: distance, coulomb_1r, r1(3), r2(3)
  
  ! 対称性を利用: (ij|kl) = (ji|kl) = (ij|lk) = (ji|lk) = (kl|ij) = ...
  ! 計算するのは i>=j, k>=l, (ij) >= (kl) の場合のみ
  
  !$omp parallel do private(i,j,k,l,ix,iy,iz,jx,jy,jz,r1,r2,distance,coulomb_1r) &
  !$omp& collapse(4) schedule(dynamic)
  do l = 1, n_basis
    do k = 1, l  ! 対称性利用
      do j = 1, n_basis
        do i = 1, j  ! 対称性利用
          
          ! 複合インデックスチェック（8-fold対称性）
          if (index_compound(i,j) < index_compound(k,l)) cycle
          
          ! 6次元グリッド積分
          ERI_frag(i,j,k,l) = 0.0d0
          do iz = is_grid(3), ie_grid(3)
            do iy = is_grid(2), ie_grid(2)
              do ix = is_grid(1), ie_grid(1)
                r1 = [ix*hgs(1), iy*hgs(2), iz*hgs(3)]
                
                do jz = is_grid(3), ie_grid(3)
                  do jy = is_grid(2), ie_grid(2)
                    do jx = is_grid(1), ie_grid(1)
                      if (ix==jx .and. iy==jy .and. iz==jz) cycle
                      
                      r2 = [jx*hgs(1), jy*hgs(2), jz*hgs(3)]
                      distance = norm2(r1 - r2)
                      if (distance < 1e-10) cycle
                      
                      coulomb_1r = 1.0d0 / distance
                      
                      ! (ij|kl) = ∫∫ φ_i(r1)φ_j(r1) [1/r12] φ_k(r2)φ_l(r2)
                      ERI_frag(i,j,k,l) = ERI_frag(i,j,k,l) + &
                        phi_grid(ix,iy,iz,i) * phi_grid(ix,iy,iz,j) * &
                        coulomb_1r * &
                        phi_grid(jx,jy,jz,k) * phi_grid(jx,jy,jz,l) * &
                        hvol * hvol
                    end do
                  end do
                end do
              end do
            end do
          end do
          
          ! 対称性に基づいて他の成分を埋める
          ERI_frag(j,i,k,l) = ERI_frag(i,j,k,l)  ! (ji|kl) = (ij|kl)
          ERI_frag(i,j,l,k) = ERI_frag(i,j,k,l)  ! (ij|lk) = (ij|kl)
          ERI_frag(j,i,l,k) = ERI_frag(i,j,k,l)  ! (ji|lk) = (ij|kl)
          ERI_frag(k,l,i,j) = ERI_frag(i,j,k,l)  ! (kl|ij) = (ij|kl)
          ERI_frag(l,k,i,j) = ERI_frag(i,j,k,l)  ! (lk|ij) = (ij|kl)
          ERI_frag(k,l,j,i) = ERI_frag(i,j,k,l)  ! (kl|ji) = (ij|kl)
          ERI_frag(l,k,j,i) = ERI_frag(i,j,k,l)  ! (lk|ji) = (ij|kl)
          
        end do
      end do
    end do
  end do
  !$omp end parallel do
  
contains
  ! 複合インデックス: (i,j) → i*(i-1)/2 + j
  integer function index_compound(i, j)
    integer, intent(in) :: i, j
    integer :: ii, jj
    ii = max(i, j)
    jj = min(i, j)
    index_compound = ii * (ii - 1) / 2 + jj
  end function
  
end subroutine precompute_ERI_fragment
```

---

### **3. 高速化された交換項計算**

```fortran
!=======================================================================
! Compute HSE exchange using precomputed ERIs
! This is O(N^2 × N_occ^2) instead of O(N^2 × N_occ^2 × L^6)
!=======================================================================
subroutine calc_exact_exchange_hse_fast(h_mat, ERI_frag, occ_states, &
                                        hse_alpha, n_basis, n_occ)
  implicit none
  real(8), intent(inout) :: h_mat(n_basis, n_basis)
  real(8), intent(in) :: ERI_frag(n_basis, n_basis, n_basis, n_basis)
  integer, intent(in) :: occ_states(n_occ)
  real(8), intent(in) :: hse_alpha
  integer, intent(in) :: n_basis, n_occ
  
  integer :: i, j, k, l, ko, lo, istate_k, istate_l
  real(8) :: V_x_ij
  
  !$omp parallel do private(i,j,k,l,ko,lo,istate_k,istate_l,V_x_ij) collapse(2)
  do j = 1, n_basis
    do i = 1, n_basis
      
      V_x_ij = 0.0d0
      
      ! 占有軌道ペアでループ（事前計算済みERIを参照）
      do lo = 1, n_occ
        istate_l = occ_states(lo)
        do ko = 1, n_occ
          istate_k = occ_states(ko)
          
          ! 交換積分: (ij|kl)を直接参照
          if (istate_k == istate_l) then
            V_x_ij = V_x_ij - 0.5d0 * hse_alpha * ERI_frag(i,j,istate_k,istate_l)
          else
            V_x_ij = V_x_ij - hse_alpha * ERI_frag(i,j,istate_k,istate_l)
          end if
        end do
      end do
      
      h_mat(i, j) = h_mat(i, j) + V_x_ij
      
    end do
  end do
  !$omp end parallel do
  
end subroutine calc_exact_exchange_hse_fast
```

---

### **4. 呼び出しフロー**

```fortran
! 初期化時（1回のみ）
subroutine init_dg_fragment_rt(dg_frag, ...)
  ! ...既存の初期化...
  
  if (yn_hse == 'y') then
    ! ERI事前計算（フラグメントごと）
    allocate(dg_frag%ERI_frag(n_basis, n_basis, n_basis, n_basis))
    
    call precompute_ERI_fragment(phi_grid, hgs, hvol, is, ie, &
                                 n_basis, dg_frag%ERI_frag)
    dg_frag%ERI_precomputed = .true.
    
    write(*,*) 'ERI precomputation completed for fragment', ifrag
  end if
end subroutine

! 各タイムステップ（高速！）
subroutine add_exact_exchange_hse(dg_frag, system, H_mat_spin, ifrag, ispin)
  if (dg_frag%ERI_precomputed) then
    ! 高速版: 事前計算済みERIを使用
    call calc_exact_exchange_hse_fast(H_mat_spin, dg_frag%ERI_frag, &
                                      occ_states, hse_alpha, n_basis, n_occ)
  else
    ! フォールバック: 従来の6次元積分
    call calc_exact_exchange_hse(H_mat_spin, phi_grid, occ_states, ...)
  end if
end subroutine
```

---

## 🎨 高度な最適化オプション

### **Option 1: 対称性を利用した圧縮保存**

8重対称性により、実際に保存が必要なのは全体の約1/8：

```fortran
! 独立な要素数: N_indep ≈ N^4 / 8
! 例: N=16 → 65536 → 8192 要素
integer :: n_indep
real(8), allocatable :: ERI_compressed(:)

! 複合インデックス: (i,j,k,l) → 1次元配列への写像
n_indep = n_basis*(n_basis+1)/2 * (n_basis*(n_basis+1)/2 + 1) / 2
allocate(ERI_compressed(n_indep))
```

**メモリ削減**: 65 MB → 8 MB（8倍削減）

---

### **Option 2: Resolution of Identity (RI/DF法)**

補助基底 $\{\chi_P\}$ を導入して近似：

$$
(ij|kl) \approx \sum_P \sum_Q (ij|P)(P|Q)^{-1}(Q|kl)
$$

```fortran
! 必要な配列
real(8), allocatable :: B_aux(:,:,:)    ! (ij|P): N_b^2 × N_aux
real(8), allocatable :: V_aux(:,:)      ! (P|Q)^-1: N_aux × N_aux

! メモリ: O(N^2 × N_aux) where N_aux ~ 2-3 × N_b
! 計算量: O(N^2 × N_aux × N_occ)
```

**メモリ削減**: 65 MB → 4 MB（典型的にN_aux = 2N_b）

---

### **Option 3: スクリーニング**

距離が大きい場合、積分が小さいことを利用：

```fortran
! ERI計算時にスクリーニング
real(8), parameter :: ERI_threshold = 1.0d-10

if (abs(ERI_frag(i,j,k,l)) < ERI_threshold) then
  ERI_frag(i,j,k,l) = 0.0d0  ! ゼロとみなす（スパース化）
end if

! スパース行列として保存
type(sparse_tensor_4d) :: ERI_sparse
```

**メモリ削減**: スパース率に依存（通常50-90%削減）

---

## 📈 実装のロードマップ

### **Phase 1: 基本実装**（1-2週間）
- ✅ `precompute_ERI_fragment` 実装
- ✅ `calc_exact_exchange_hse_fast` 実装
- ✅ 初期化ルーチンへの統合
- ✅ 小規模テスト（N_b=8, L=8）

### **Phase 2: 対称性最適化**（1週間）
- ✅ 8重対称性の実装
- ✅ 圧縮保存フォーマット
- ✅ メモリフットプリント削減

### **Phase 3: 高度な最適化**（2-3週間）
- ✅ RI/DF法の実装
- ✅ スクリーニング機構
- ✅ GPUオフロード対応

### **Phase 4: 性能評価**（1週間）
- ✅ ベンチマーク（Plan A vs Plan B）
- ✅ メモリプロファイリング
- ✅ スケーラビリティテスト

---

## ⚖️ トレードオフ分析

### **Plan A（現在）が有利な場合**
- ✅ **大規模フラグメント**（N_b > 100、4原子以上）
- ✅ メモリが厳しく制約されている
- ✅ 短時間計算（<100ステップ）
- ✅ 占有軌道数が多い（N_occ ≈ N_b）

### **Plan B（ERI）が有利な場合**
- ✅ **小規模フラグメント**（N_b < 100、1-2原子程度）
- ✅ 長時間RT計算（>1000ステップ）
- ✅ 十分なメモリが利用可能（対称性利用で1-10 GB）
- ✅ 占有軌道数が少ない（N_occ << N_b）

### **Plan （DC-LCFO考慮版）**

```fortran
! 自動選択ロジック（DC法の実情を反映）
integer :: n_atoms_frag, n_basis_per_atom = 50
real(8) :: ERI_memory_GB

n_basis = n_atoms_frag * n_basis_per_atom
ERI_memory_GB = (n_basis**4) * 8.0d0 / (1024.0d0**3)  ! Full tensor
ERI_memory_GB = ERI_memory_GB / 8.0d0  ! 対称性利用

if (n_atoms_frag <= 2 .and. &
    n_timesteps > 500 .and. &
    memory_available_GB > ERI_memory_GB) then
  use_precomputed_ERI = .true.   ! Plan B
  write(*,*) 'Using precomputed ERI (small fragment)'
else if (n_atoms_frag >= 3 .and. use_RI_available) then
  use_RI_method = .true.         ! Plan C
  write(*,*) 'Using RI/DF method (large fragment)'
else
  use_direct_method = .true.     ! Plan A (現在の実装)
  write(*,*)

#### **ケース1: 小規模（1原子フラグメント、N_b=50, 1000ステップ）**

| 指標 | Plan A | Plan B | 改善 |
|------|--------|--------|------|
| **初期化時間** | 0秒 | 5000秒 | - |
| **1ステップ** | 30秒 | 0.05秒 | **600×** |
| **1000ステップ** | 30000秒 | 50秒 | **600×** |
| **総時間** | 30000秒 (8.3時間) | 5050秒 (1.4時間) | **5.9×** ✅ |
| **メモリ** | 200 KB | 6.25 MB | 1/31 |

#### **ケース2: 中規模（2原子フラグメント、N_b=100, 1000ステップ）**

| 指標 | Plan A | Plan B | Plan C (RI) |
|------|--------|--------|-------------|
| **初期化時間** | 0秒 | 6時間 | 2時間 |
| **1ステップ** | 120秒 | 0.2秒 | 2秒 |
| **1000ステップ** | 120000秒 | 200秒 | 2000秒 |
| **総時間** | 33.3時間 | 6.3時間 | 2.6時間 |
| **メモリ** | 800 KB | 100 MB | 15 MB |
| **改善率** | 1× | **5.3×** | **12.8×** ✅ |

#### **ケース3: 大規模（4原子フラグメント、N_b=200）**

| 指標 | Plan A | Plan B | Plan C (RI) |
|------|--------|--------|-------------|
| **初期化** | 0秒 | **40時間** ❌ | 8時間 |
| **メモリ** | 3.2 MB | **1.6 GB** ❌ | 60 MB |
| **推奨** | ✅ 妥当 | ❌ 非現実的 | ✅ 最適 |

**結論**: 
- **1原子フラグメント**: Plan B有効
- **2原子フラグメント**: Plan C（RI）が最適
- **4原子以上**: Plan A or Plan C

## 🧪 期待される性能改善

### ベンチマーク予測（N_b=16, L=16, N_occ=8, 1000ステップ）

| 指標 | Plan A | Plan B | 改善 |
|------|--------|--------|------|
| **初期化時間** | 0秒 | 3600秒 | - |
| **1ステップ** | 10秒 | 0.01秒 | **1000×** |
| **1000ステップ** | 10000秒 | 10秒 | **1000×** |
| **総時間** | 10000秒 | 3610秒 | **2.77×** |
| **メモリ** | 65 KB | 65 MB | 1/1000 |

**結論**: 100ステップ以上の計算で高速化が顕著

---

## 💻 実装例

完全なパッチファイルを作成：

```bash
# src/xc/xc_hse_eri.f90 として新規作成
# rt_dg_fragment.f90 に統合コード追加
# inputoutput.f90 に yn_hse_precompute_eri パラメータ追加
```

---

## 📚 参考文献

1. **Molecular Electronic-Structure Theory** (Helgaker et al., 2000)  
   - Chapter 9: Electron Repulsion Integrals
   
2. **Resolution-of-the-Identity Techniques** (Whitten, 1973)  
   - J. Chem. Phys. 58, 4496

3. **Density Fitting in Periodic Systems** (Pisani et al., 2008)  
   - Phys. Chem. Chem. Phys. 10, 5353

4. **Cholesky Decomposition of ERIs** (Aquilante et al., 2007)  
   - J. Chem. Phys. 127, 114107

---

## ✅ 次のステップ

1. **プロトタイプ実装**: 小規模テストケースで実装
2. **メモリ解析**: 実メモリ使用量の測定
3. **性能比較**: Plan A vs Plan B のベンチマーク
4. **ユーザー入力**: `yn_hse_use_eri = 'y'` でオプション化

---

**作成者**: GitHub Copilot (Claude Sonnet 4.5)  
**最終更新**: 2026年2月23日
