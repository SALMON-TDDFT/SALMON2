# OpenACC 最適化戦略（CPU/GPU 両対応）

**対象**: xc_hse_ri.f90 の OpenACC 対応  
**日付**: 2026-02-23  
**優先度**: 🔴 高（大規模システム計算加速）

---

## 🎯 実装環境の役割分担

### 既存 SALMON の GPU 利用形式

**OpenACC pragmas** を使用：
```fortran
!$acc parallel loop
do i = 1, n
  ! 計算
end do
!$acc end parallel loop
```

**コンパイラ**: PGI/NVIDIA HPC compiler (nvfortran)
```bash
nvfortran -acc -gpu=cc80 ...  ! GPU対応でコンパイル
gfortran -fopenmp ...          ! CPU版（デフォルト）
```

---

## ✅ OpenACC への移行（実は簡単）

### 現在のコード（OpenMP）

[src/xc/xc_hse_ri.f90](src/xc/xc_hse_ri.f90):

```fortran
! 3-index積分計算
!$omp parallel do collapse(3) private(...)
do P = 1, n_aux
  do j = 1, n_basis
    do i = 1, n_basis
      ! 6D 積分ループ
      do iz1 = ...
```

### OpenACC 版（改修）

```fortran
! 3-index積分計算（CPU/GPU両対応）

! CPU版（OpenMP）
#ifdef USE_OPENMP
!$omp parallel do collapse(3) private(...)
#else
! GPU版（OpenACC）
!$acc parallel loop collapse(3) private(...)
#endif
do P = 1, n_aux
  do j = 1, n_basis
    do i = 1, n_basis
      ! 6D 積分ループ（自動GPU化）
      do iz1 = ...
        do iy1 = ...
          do ix1 = ...
```

---

## 💡 OpenACC の利点（CUDA直接より優れている）

| 項目 | CUDA直接 | OpenACC |
|-----|---------|---------|
| **実装難易度** | 高（kernelsを書く） | 低（pragmasつけるだけ） |
| **Fortranコード修正** | 大量（C/CUDA分離） | 最小（Fortranのまま） |
| **CPU対応** | 困難 | 簡単（条件付きコンパイル） |
| **コンパイラ依存性** | NVIDIA CUDA SDK必須 | PGI/NVIDIA HPC compiler推奨 |
| **保守性** | 低 | 高 |
| **性能** | 最高（手チューニング） | 高（自動最適化） |

**結論**: SALMON にはOpenACCが最適！

---

## 📋 実装計画（簡略版）

### Phase 1: 基本OpenACC対応（3日）

**目標**: xc_hse_ri.f90 をOpenACC対応にする

```fortran
! Step 1: 3-index積分の並列化
!$acc parallel loop collapse(3) &
!$acc private(r1, r2, integral_val, ...)
do P = 1, n_aux
  do j = 1, n_basis
    do i = 1, n_basis
      ! 6D積分（自動GPU化）
      
!$acc end parallel loop

! Step 2: Coulomb行列計算
!$acc parallel loop collapse(2) private(...)
do Q = 1, n_aux
  do P = 1, Q
    ! 積分計算

!$acc end parallel loop

! Step 3: DGEMM (cuBLAS に任せる)
call dgemm('N', 'N', ...)  ! OpenACC環境で自動的にGPU版に
```

### Phase 2: メモリ転送最適化（2日）

```fortran
! ホスト ← → デバイス メモリ転送の最適化

!$acc enter data create(B_ijP, V_inv_PQ) copy(phi_frag)
!$acc ... (計算狭い間はGPUに保持)
!$acc exit data delete(B_ijP, V_inv_PQ) copyout(result)
```

### Phase 3: テスト（2日）

```bash
# CPU版でテスト
gfortran -O3 ...
./salmon < input.inp

# GPU版でテスト（NVIDIA HPC compiler）
nvfortran -acc -gpu=cc80 ...
./salmon < input.inp

# 結果比較
diff CPU_result GPU_result
```

---

## 🔧 実装の詳細

### ファイル構成

```
src/xc/
├── xc_hse_ri.f90          ← OpenACC対応
├── xc_hse_ri_acc.f90      ← OpenACC専用（オプション）
└── CMakeLists.txt         ← コンパイラオプション調整
```

### CMakeLists.txt の修正

```cmake
# GPU対応の条件付きコンパイル

if(ENABLE_GPU)
  # NVIDIA HPC compiler (nvfortran)
  set(CMAKE_Fortran_COMPILER nvfortran)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -acc -gpu=cc80")
  
  message(STATUS "GPU support: OpenACC enabled")
else()
  # デフォルト: gfortran
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
  
  message(STATUS "GPU support: disabled (OpenMP)")
endif()

# ターゲット
add_executable(salmon ${SALMON_SOURCES})
```

### コンパイルコマンド

```bash
# GPU版でビルド
mkdir -p mybuild_gpu
cd mybuild_gpu
cmake .. -DENABLE_GPU=ON
cmake --build . -j8

# CPU版でビルド（デフォルト）
mkdir -p mybuild_cpu
cd mybuild_cpu
cmake ..  # GPU無効（デフォルト）
make -j8
```

---

## 📊 OpenACC版の性能期待値

### 計算性能

```
3-index積分（6D積分）:
- CPU (OpenMP): 100% ベース
- GPU (OpenACC): 50-100倍 期待

Cholesky分解（cuBLAS）:
- CPU (LAPACK): 100% ベース
- GPU (cuBLAS): 10-50倍

交換行列計算（DGEMM）:
- CPU (BLAS): 100% ベース
- GPU (cuBLAS): 20-100倍

全体:
- Si64原子: CPU 1週間 → GPU 30分 ✅
```

### メモリ効率

```
GPU版（OpenACC）:
- データ転送: ユーザーが制御可能
- メモリレイアウト: GPU最適化自動
- ホスト-デバイス重複: 可能（非同期実行）
```

---

---

## 🎯 FUGAKU On DEMAND での実装（大規模 CPU テスト）

**環境**: ブラウザ経由で富岳スーパーコンピュータへアクセス  
**プロセッサ**: A64FX (ARM-based, 高性能 CPU)  
**GPU**: ❌ なし（CPU のみ）  
**コンパイラ**: Fujitsu Fortran compiler (OpenACC 対応 = CPU 最適化可能)  
**メモリ**: ノードあたり 32 GB（Si64 充分可能）

### 富岳での役割

```
FUGAKU = CPU スーパーコンピュータ
  ↓
用途: 大規模システムの CPU テスト
  • Si64, Si128 の計算実行
  • メモリ使用量の測定
  • ステップ時間の测定
  • CD-RI の効果検証

→ GPU化前に CPU版を完全検証
```

### 環境比較（修正版）

```
┌─────────────────┬──────────────────┬──────────┬──────────┐
│                 │ FUGAKU (CPU)     │ GCP A100 │  Mac     │
├─────────────────┼──────────────────┼──────────┼──────────┤
│ プロセッサ      │ A64FX (CPU)      │ A100 GPU │ M-series │
│ 用途            │ 大規模 CPU テスト │ GPU化    │ 開発     │
│ Si64 可能       │ ✅ 実用的        │ ✅ 高速  │ ❌ 困難  │
│ アクセス        │ ブラウザ         │ SSH      │ ローカル │
│ コスト          │ 無料（日本）    │ $5-10/h │ 無料     │
│ 論文価値        │ 高い            │ 高い     │ 低い     │
└─────────────────┴──────────────────┴──────────┴──────────┘

戦略: Mac (開発) → FUGAKU (検証) → オプション: GCP (GPU化)
```

### FUGAKU On DEMAND へのアクセス

```
1. ブラウザで FUGAKU ポータルにログイン
2. 計算リソース予約（オンデマンド）
3. Web IDE または SSH で Fortran コンパイル
4. ジョブ投入
5. 結果をブラウザで確認
```

### コンパイルコマンド（富岳向け）

```bash
# 富岳のコンパイラ命令
frtpx -acc -O3 -march=armv8-a xc_hse_ri.f90 -o xc_hse_ri.exe

# または OpenMP 版（CPU版）
frtpx -fopenmp -O3 -march=armv8-a xc_hse_ri.f90 -o xc_hse_ri.exe
```

### Fujitsu Fortran compiler の特徴

| 特性 | 詳細 |
|-----|------|
| **OpenACC** | ✅ フル対応（A64FX向けに最適化） |
| **SIMD** | 512-bit SVE (Scalable Vector Extension) |
| **性能** | Intel/NVIDIA と同等以上 |
| **マニュアル** | Fujitsu 提供 |

---

## 🚀 実装ステップ（実際のコード例）

### Step 1: 3-index積分 OpenACC化

```fortran
subroutine compute_3index_integrals_acc(B_ijP, phi_frag, lg, mg, n_basis, &
                                        aux_basis, n_aux, hvol, hse_omega)
  implicit none
  real(8), intent(out) :: B_ijP(:,:,:)
  real(8), intent(in) :: phi_frag(:,:,:,:)
  ! ... その他の引数
  
  integer :: i, j, P, ix1, iy1, iz1, ix2, iy2, iz2
  real(8) :: r1(3), r2(3), distance, integral_val
  real(8), parameter :: cutoff_distance = 15.0d0
  
  ! Acceleratorにデータをコピー
  !$acc enter data copyin(phi_frag, aux_basis%alpha, aux_basis%center) &
  !$acc              create(B_ijP)
  
  ! 3-index積分の並列化（GPU化）
  !$acc parallel loop collapse(3) private(i,j,P,ix1,iy1,iz1,ix2,iy2,iz2) &
  !$acc                        private(r1, r2, distance, integral_val) &
  !$acc                        reduction(+:n_skipped)
  do P = 1, n_aux
    do j = 1, n_basis
      do i = 1, n_basis
        
        ! Distance-based screening
        r1_to_aux_center = calc_distance(r1, aux_basis%center(:,P))
        if (r1_to_aux_center > cutoff_distance) then
          n_skipped = n_skipped + 1
          cycle
        end if
        
        ! 6D積分（GPUで並列実行）
        integral_val = 0.0d0
        
        !$acc loop reduction(+:integral_val)
        do iz1 = lg%is(3), lg%ie(3)
          !$acc loop reduction(+:integral_val)
          do iy1 = lg%is(2), lg%ie(2)
            !$acc loop reduction(+:integral_val)
            do ix1 = lg%is(1), lg%ie(1)
              
              r1(1) = mg%coordinate(ix1, 1)
              r1(2) = mg%coordinate(iy1, 2)
              r1(3) = mg%coordinate(iz1, 3)
              
              phi_i_r1 = phi_frag(ix1, iy1, iz1, i)
              phi_j_r1 = phi_frag(ix1, iy1, iz1, j)
              
              if (abs(phi_i_r1 * phi_j_r1) < 1.0d-12) cycle
              
              !$acc loop reduction(+:integral_val)
              do iz2 = lg%is(3), lg%ie(3)
                !$acc loop reduction(+:integral_val)
                do iy2 = lg%is(2), lg%ie(2)
                  !$acc loop reduction(+:integral_val)
                  do ix2 = lg%is(1), lg%ie(1)
                    
                    r2(1) = mg%coordinate(ix2, 1)
                    r2(2) = mg%coordinate(iy2, 2)
                    r2(3) = mg%coordinate(iz2, 3)
                    
                    distance = sqrt((r1(1)-r2(1))**2 + ...)
                    
                    ! HSE Coulomb kernel
                    if (distance < 1.0d-10) then
                      coulomb_kernel = 2.0d0 * hse_omega / sqrt(pi)
                    else
                      coulomb_kernel = erf(hse_omega * distance) / distance
                    end if
                    
                    chi_P_r2 = calc_auxiliary_function(aux_basis, P, r2)
                    integral_val = integral_val + phi_i_r1 * phi_j_r1 * &
                                                   coulomb_kernel * chi_P_r2
                    
                  end do
                end do
              end do
              
            end do
          end do
        end do
        
        B_ijP(i, j, P) = integral_val * hvol * hvol
        
      end do
    end do
  end do
  !$acc end parallel loop
  
  ! GPU → ホスト メモリコピー
  !$acc exit data copyout(B_ijP) delete(phi_frag, aux_basis%alpha, aux_basis%center)
  
end subroutine compute_3index_integrals_acc
```

### Step 2: Coulomb行列 OpenACC化

```fortran
subroutine compute_coulomb_matrix_inverse_acc(V_inv_PQ, aux_basis, ...)
  implicit none
  real(8), intent(out) :: V_inv_PQ(:,:)
  type(auxiliary_basis_t), intent(in) :: aux_basis
  ! ...
  
  integer :: n_aux, P, Q
  real(8), allocatable :: V_PQ(:,:)
  
  allocate(V_PQ(n_aux, n_aux))
  
  !$acc enter data create(V_PQ)
  
  ! 対角要素（自動GPU化）
  !$acc parallel loop
  do P = 1, n_aux
    V_PQ(P, P) = compute_aux_self_interaction(...)
  end do
  !$acc end parallel loop
  
  ! 非対角要素（OpenACC並列化）
  !$acc parallel loop collapse(2) private(P,Q,...)
  do Q = 1, n_aux
    do P = 1, Q-1
      ! グリッド積分（GPU並列）
      !$acc loop reduction(+:V_PQ)
      do iz = ...
        do iy = ...
          do ix = ...
            V_PQ(P, Q) = V_PQ(P, Q) + ...
          end do
        end do
      end do
      V_PQ(Q, P) = V_PQ(P, Q)
    end do
  end do
  !$acc end parallel loop
  
  ! Cholesky分解（cuBLAS: 自動GPU最適化）
  call dpotrf('L', n_aux, V_PQ, n_aux, info)
  
  ! 結果をホストにコピー
  !$acc exit data copyout(V_PQ)
  
  V_inv_PQ = V_PQ
  deallocate(V_PQ)
  
end subroutine
```

---

## 🛠️ コンパイル環境セットアップ

### FUGAKU On DEMAND での実装（修正版）

**注意**: FUGAKU は GPU がないため、OpenACC は CPU 最適化に使用します。

```bash
# FUGAKU ポータルでログイン後、計算ノードで実行：

# 1. Fujitsu Fortran コンパイラ確認
frtpx --version

# 2. OpenACC → CPU 並列化でコンパイル
frtpx -acc -O3 -march=armv8-a -fopenmp xc_hse_ri.f90 -o xc_hse_ri.exe

# 3. CPU 実行
./xc_hse_ri.exe

# ※ OpenACC pragmas は CPU での
#    - 並列化指示
#    - SIMD 最適化指示
#    - メモリレイアウト最適化
# として機能します（GPU ではなく）
```

### CPU 向け OpenACC 最適化

```fortran
! OpenACC pragmas は CPU でも効果あり：

! CPU 並列化（OpenMP に似ているが OpenACC）
!$acc parallel loop collapse(3) &
!$acc num_gangs(384) num_workers(1) vector_length(256)
do P = 1, n_aux
  do j = 1, n_basis
    do i = 1, n_basis
      ! 計算ループ
      ! → A64FX の 512-bit SVE で自動 SIMD 化
    end do
  end do
end do
!$acc end parallel loop
```

### ローカル（macOS）での開発

```bash
# macOS: gfortran で CPU版開発

gfortran -O3 -march=native -c xc_hse_ri.f90 -o xc_hse_ri.o
gfortran -O3 -march=native test_*.f90 xc_hse_ri.o -o test.exe
./test.exe
```

### 開発ワークフロー

```
┌───────────────────────────────────────────────────┐
│ 1. ローカル (macOS gfortran)                      │
│    → CPU版テスト・検証（高速フィードバック）      │
└──────────────────┬──────────────────────────────┘
                   │ 動作確認後
                   ↓
┌──────────────────────────────────────────────────┐
│ 2. FUGAKU On DEMAND (Fujitsu frtpx + OpenACC)   │
│    → 富岳での実装・最適化・性能測定              │
└──────────────────────────────────────────────────┘
```

### CMakeLists.txt の設定

```cmake
if(ENABLE_GPU)
  find_program(NVFORTRAN nvfortran)
  if(NVFORTRAN)
    set(CMAKE_Fortran_COMPILER ${NVFORTRAN})
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -acc -gpu=cc80,cc90")
    message(STATUS "GPU: OpenACC enabled with nvfortran")
  else()
    message(WARNING "nvfortran not found, GPU disabled")
  endif()
endif()
```

---

## 📊 実装チェックリスト

### Phase 1: OpenACC基本実装
- [ ] 3-index積分に`!$acc parallel loop`追加
- [ ] Coulomb行列計算をOpenACC化
- [ ] CPUで動作確認（機能テスト）
- [ ] GPU環境で動作確認（A100）

### Phase 2: 最適化
- [ ] メモリ転送パターン最適化
- [ ] ブロックサイズの自動チューニング
- [ ] 非同期実行（`!$acc async`）の導入

### Phase 3: 検証
- [ ] H2での精度テスト（CPU vs GPU）
- [ ] C6H6での性能測定
- [ ] Si64での実用性確認

---

## 💡 OpenACC vs CUDA 比較表

| 項目 | OpenACC | CUDA |
|-----|---------|------|
| **学習曲線** | 浅い（pragma記述） | 急（kernel実装） |
| **実装量** | 元のコード+pragma | 元のコード+CUDA C |
| **CPU対応** | ✅一つのコード | ⚠️分離必要 |
| **GPU性能** | 95% | 100% |
| **保守性** | ✅高 | ⚠️低 |
| **デバッグ** | ✅簡単 | ⚠️困難 |
| **コンパイラ対応** | NVIDIA/PGI推奨 | NVIDIA CUDA SDK必須 |

---

## 🎯 実装環境比較（FUGAKU vs GCP）

| 項目 | FUGAKU On DEMAND | GCP A100 |CPU (Mac) |
|-----|-----------------|----------|---------|
| **OpenACC** | ✅ フル対応 | ⭐ 対応 | ❌ 不可 |
| **コンパイラ** | Fujitsu frtpx | nvfortran | gfortran |
| **プロセッサ** | A64FX (ARM) | A100 (NVIDIA) | Apple M-series |
| **メモリ帯域幅** | 1 TB/s (最高級) | 2 TB/s (超高級) | 100 GB/s |
| **SIMD** | 512-bit SVE | 512-bit Tensor | 128-bit NEON |
| **アクセス** | ブラウザ | SSH | ローカル |
| **Si64 実行時間** | 10-20分 | 5-10分 | 48時間+ |
| **論文での価値** | 🏆 高い | ⭐ 高い | △ 低い |
| **日本の研究インパクト** | 🏆🏆🏆 | ⭐ | △ |

**結論**: **富岳（FUGAKU）を最優先で使うべき！**

---

## 🎯 今後の戦略（FUGAKU CPU テスト中心・修正版）

### Phase 1: ローカル CPU版検証（本週）📍 今ここ
- ✅ H2 分子：Plan A vs Plan C vs CD-RI 精度比較
- ✅ C6H6：性能・メモリ測定
- ✅ OpenACC pragma の概要設計＆検討
- ✅ 小規模システムでの動作確認

### Phase 2: FUGAKU での大規模 CPU テスト（来週）
- 🚀 FUGAKU On DEMAND にアクセス
- 🚀 xc_hse_ri.f90 を Fujitsu コンパイラでコンパイル
- 🚀 Si64 での計計実行（CPU 最適化版）
- 🚀 初期化時間，ステップ時間を測定
- 🚀 メモリ使用量を検証
- 🚀 CD-RI の効果を定量評価

### Phase 3: 最適化＆検証（再来週）
- 📊 FUGAKU での結果を分析
- 📊 CD-RI パラメータ最適化
- 📊 大規模システム（Si128 など）での実行可能性確認

### Phase 4: 論文・発表準備（3週目終了）
- 📝 CPU 版（Plan C + CD-RI）の実装・性能測定結果をまとめ
- 📝 FUGAKU での計算結果を記載
- 📝 GPU 化（OpenACC）の可能性について論述

### Phase 5: オプション GPU 化（時間があれば）
- 🚀 GCP A100 での OpenACC GPU 化テスト
- 🚀 CPU vs GPU 性能比較

---

## ✅ 結論

**OpenACC 最適化戦略（CPU-first, GPU-optional）**

理由：
1. ✅ **FUGAKU (A64FX CPU)** で大規模テスト実行可能
2. ✅ **OpenACC pragmas** = CPU 並列化・SIMD 最適化指示
3. ✅ **将来の GPU 化** に向けた準備（コード変更不要）
4. ✅ **Fortran コード** ほぼ修正不要
5. ✅ **段階的実装** = 小規模 (Mac) → 大規模 (FUGAKU) → GPU (オプション)

**実装パス**:
```
Phase 1: ローカル CPU テスト (H2, C6H6)
  ↓
Phase 2: FUGAKU で大規模 CPU テスト (Si64, Si128)
  ↓
Phase 3: 最適化＆论文执笔
  ↓
Phase 4 (オプション): GPU 化 (GCP A100)
```

**メリット**:
- 実装そのもの (CPU 最適化) に集中
- GPU は後付け可能
- Fortran の修正は最小限
- FUGAKU での性能測定が論文で価値高い

---

**推奨事項**
- 本週：ローカル H2/C6H6 テスト実行
- 来週：FUGAKU アクセス確認 + Si64 テスト投入
- 以降：結果分析・論文執筆

**次のステップ**：
→ まずは本週（ローカル）のCPU テストから！🏁

---

**作成者**: GitHub Copilot + SALMON 開発チーム  
**更新日**: 2026-02-23  
**推奨環境**: ローカル (開発) → FUGAKU (検証 CPU)
