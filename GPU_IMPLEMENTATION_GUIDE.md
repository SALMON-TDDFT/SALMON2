# GPU実装ガイド: RI/DF HSE交換項の高速化

**対象**: Plan C (RI/DF法) のGPU実装  
**期待効果**: CPU版の40-60倍高速化  
**実装期間**: 2-3週間（CUDA経験者）  
**更新日**: 2026年2月23日

---

## 🎯 GPU実装の目的

### 現状（CPU版）
- **2原子フラグメント**: 1000ステップ → 1時間
- **4原子フラグメント**: 1000ステップ → 4時間
- **8原子フラグメント**: ほぼ不可能（20時間）

### GPU版目標
- **2原子フラグメント**: 1000ステップ → **1分** (60倍高速)
- **4原子フラグメント**: 1000ステップ → **5分** (48倍高速)
- **8原子フラグメント**: 1000ステップ → **20分** (60倍高速)

---

## 📊 性能見積もり

### GPU加速の根拠

#### 1. 3-index積分計算（初期化）
```
計算：B_ijP = ∫∫ φ_i(r1) φ_j(r1) [1/r12] χ_P(r2) dr1 dr2
複雑度：O(N² × N_aux × L⁶)
並列性：N² × N_aux の独立計算
GPU加速：100-200倍（完全並列化可能）
```

**実装**:
- 各 (i, j, P) 組み合わせを独立スレッドで処理
- グリッド積分のループ展開
- 共有メモリでφ_i, φ_j値をキャッシュ

#### 2. Coulomb行列計算
```
計算：V_PQ = (P|Q)
複雑度：O(N²_aux × L⁶)
並列性：N²_aux の独立計算
GPU加速：50-100倍
```

#### 3. BLAS演算（各タイムステップ）
```
計算：C_klQ = Σ_P B_klP × V^(-1)_PQ
使用：DGEMM (行列積)
GPU加速：10-50倍（cuBLAS最適化）
```

### マルチGPU効果

| GPU数 | 理論加速 | 実効加速 | 通信オーバーヘッド |
|-------|---------|---------|------------------|
| 1 GPU | 1.0x | 1.0x | - |
| 2 GPU | 2.0x | 1.8x | ~10% |
| 4 GPU | 4.0x | 3.5x | ~12% |
| 8 GPU | 8.0x | 6.5x | ~18% |

**推奨**: 2-4 GPU（コスト対効果最適）

---

## 💻 実装アーキテクチャ

### GPUメモリ配置

```fortran
! Device memory (GPU上)
real(8), device, allocatable :: phi_frag_d(:,:,:,:)    ! Basis functions
real(8), device, allocatable :: B_ijP_d(:,:,:)         ! 3-index integrals
real(8), device, allocatable :: V_inv_PQ_d(:,:)        ! Coulomb inverse
real(8), device, allocatable :: H_mat_d(:,:)           ! Hamiltonian matrix

! Auxiliary basis (constant)
real(8), device, allocatable :: aux_alpha_d(:)         ! Exponents
real(8), device, allocatable :: aux_center_d(:,:)      ! Centers
integer, device, allocatable :: aux_lm_d(:,:)          ! Quantum numbers
```

### CUDAカーネル設計

#### Kernel 1: 3-index積分
```cuda
__global__ void compute_B_ijP_kernel(
    const double* phi_frag,        // [Lx, Ly, Lz, N_basis]
    const double* aux_alpha,       // [N_aux]
    const double* aux_center,      // [3, N_aux]
    const int* aux_lm,             // [2, N_aux] (l, m)
    double* B_ijP,                 // [N_basis, N_basis, N_aux]
    const int N_basis, const int N_aux,
    const int Lx, const int Ly, const int Lz,
    const double hx, const double hy, const double hz,
    const double omega
) {
    // Thread index: (i, j, P)
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int P = blockIdx.z * blockDim.z + threadIdx.z;
    
    if (i >= N_basis || j >= N_basis || P >= N_aux) return;
    
    double integral = 0.0;
    
    // 6D grid integration (can be further parallelized)
    for (int iz1 = 0; iz1 < Lz; iz1++) {
        for (int iy1 = 0; iy1 < Ly; iy1++) {
            for (int ix1 = 0; ix1 < Lx; ix1++) {
                double r1[3] = {ix1*hx, iy1*hy, iz1*hz};
                double phi_i = phi_frag[ix1 + Lx*(iy1 + Ly*(iz1 + Lz*i))];
                double phi_j = phi_frag[ix1 + Lx*(iy1 + Ly*(iz1 + Lz*j))];
                
                if (fabs(phi_i * phi_j) < 1e-12) continue;
                
                for (int iz2 = 0; iz2 < Lz; iz2++) {
                    for (int iy2 = 0; iy2 < Ly; iy2++) {
                        for (int ix2 = 0; ix2 < Lx; ix2++) {
                            double r2[3] = {ix2*hx, iy2*hy, iz2*hz};
                            double r12 = distance(r1, r2);
                            double coulomb = (r12 < 1e-10) ? 
                                2.0*omega/sqrt(M_PI) : erf(omega*r12)/r12;
                            
                            double chi_P = eval_aux_function(P, r2, aux_alpha, 
                                                             aux_center, aux_lm);
                            
                            integral += phi_i * phi_j * coulomb * chi_P;
                        }
                    }
                }
            }
        }
    }
    
    B_ijP[i + N_basis*(j + N_basis*P)] = integral * hx*hy*hz * hx*hy*hz;
}

// Launch configuration
dim3 block(8, 8, 8);  // 512 threads/block
dim3 grid((N_basis+7)/8, (N_basis+7)/8, (N_aux+7)/8);
compute_B_ijP_kernel<<<grid, block>>>(phi_frag_d, ...);
```

#### Kernel 2: 補助関数評価
```cuda
__device__ double eval_aux_function(
    int P, const double r[3],
    const double* aux_alpha,
    const double* aux_center,
    const int* aux_lm
) {
    double dx = r[0] - aux_center[3*P + 0];
    double dy = r[1] - aux_center[3*P + 1];
    double dz = r[2] - aux_center[3*P + 2];
    double r2 = dx*dx + dy*dy + dz*dz;
    double rlen = sqrt(r2);
    
    int l = aux_lm[2*P + 0];
    int m = aux_lm[2*P + 1];
    double alpha = aux_alpha[P];
    
    // Normalization
    double norm = pow(2.0*alpha/M_PI, 0.75) * pow(4.0*alpha, 0.5*l);
    
    // Radial part
    double radial = (rlen < 1e-10 && l > 0) ? 0.0 : 
                    pow(rlen, l) * exp(-alpha * r2);
    
    // Angular part (spherical harmonic)
    double angular = spherical_harmonic(l, m, dx, dy, dz, rlen);
    
    return norm * radial * angular;
}
```

#### Kernel 3: 交換行列構築（cuBLAS統合）
```cuda
// Use cuBLAS for DGEMM
cublasHandle_t handle;
cublasCreate(&handle);

// C_klQ = B_klP × V^(-1)_PQ
cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
            N_basis*N_basis, N_aux, N_aux,
            &alpha, B_ijP_d, N_basis*N_basis,
                    V_inv_PQ_d, N_aux,
            &beta,  C_klQ_d, N_basis*N_basis);

// H_mat += B_ijQ × C_klQ × D_kl (custom kernel)
__global__ void add_exchange_kernel(...) {
    // Thread per (i, j)
    // Accumulate over Q, k, l
}
```

---

## 🔧 実装ステップ

### Phase 1: 基本CUDAインターフェース（3日）

**目標**: Fortran-CUDA連携確立

```fortran
! src/xc/xc_hse_ri_gpu.f90
module xc_hse_ri_gpu
  use iso_c_binding
  use cudafor
  implicit none
  
  interface
    ! C/C++ wrapper functions
    subroutine compute_B_ijP_cuda(phi_frag, B_ijP, N_basis, N_aux, ...) bind(C)
      use iso_c_binding
      type(c_ptr), value :: phi_frag, B_ijP
      integer(c_int), value :: N_basis, N_aux
    end subroutine
  end interface
  
contains
  
  subroutine init_hse_ri_gpu(ri_data, phi_frag, ...)
    ! Allocate GPU memory
    ! Copy phi_frag to GPU
    ! Launch CUDA kernels
  end subroutine
  
end module
```

**C++/CUDA実装**:
```cpp
// src/xc/cuda_kernels/hse_ri_kernels.cu
extern "C" {
  void compute_B_ijP_cuda(
      const double* phi_frag,
      double* B_ijP,
      int N_basis, int N_aux,
      int Lx, int Ly, int Lz,
      double hx, double hy, double hz,
      double omega
  ) {
      // Launch kernel
      dim3 block(8, 8, 8);
      dim3 grid((N_basis+7)/8, (N_basis+7)/8, (N_aux+7)/8);
      compute_B_ijP_kernel<<<grid, block>>>(
          phi_frag, B_ijP, N_basis, N_aux, Lx, Ly, Lz, hx, hy, hz, omega
      );
      cudaDeviceSynchronize();
  }
}
```

### Phase 2: 3-index積分GPU化（5日）

1. **補助関数評価のGPU実装**
2. **6D積分の並列化**
3. **共有メモリ最適化**
4. **テスト**: CPU版との精度比較

### Phase 3: cuBLAS統合（2日）

1. **Cholesky分解**: `cusolverDnDpotrf`
2. **行列積**: `cublasDgemm`
3. **メモリ転送最適化**

### Phase 4: マルチGPU対応（3日）

```fortran
! NCCL for multi-GPU communication
use nccl
use mpi

! Initialize NCCL
call ncclCommInitRank(nccl_comm, nranks, nccl_id, rank)

! Distribute fragments across GPUs
do ifrag = ifrag_start, ifrag_end
  igpu = mod(ifrag, ngpu)
  call cudaSetDevice(igpu)
  
  ! Compute on GPU igpu
  call compute_B_ijP_cuda(...)
end do

! All-reduce H_mat contributions
call ncclAllReduce(H_mat_d, H_mat_result_d, size_H, ncclDouble, 
                   ncclSum, nccl_comm, stream)
```

### Phase 5: 最適化とテスト（3日）

1. **カーネル占有率の最適化**
2. **メモリコアレッシング**
3. **ストリーム並列実行**
4. **性能測定と比較**

---

## 📈 期待される性能プロファイル

### GPU使用率

| カーネル | GPU占有率 | メモリ帯域 | 計算/メモリ比 |
|---------|----------|-----------|-------------|
| **3-index積分** | 85-95% | 高 | 高（計算律速） |
| **補助関数評価** | 90-98% | 低 | 超高 |
| **cuBLAS DGEMM** | 98-100% | 超高 | 最適 |

### メモリ使用量（GPU VRAM）

| システム | N_basis | N_aux | 必要VRAM | 推奨GPU |
|---------|---------|-------|---------|---------|
| 2原子 | 100 | 300 | 1 GB | V100 (16GB) ✅ |
| 4原子 | 200 | 600 | 4 GB | A100 (40GB) ✅ |
| 8原子 | 400 | 1200 | 15 GB | A100 (40GB) ✅ |
| 16原子 | 800 | 2400 | 60 GB | H100 (80GB) ✅ |

---

## 🎯 実装チェックリスト

### 必須機能
- [ ] CUDAカーネル: 3-index積分
- [ ] CUDAカーネル: 補助関数評価
- [ ] CUDAカーネル: 交換行列構築
- [ ] cuBLAS統合: DGEMM, Cholesky
- [ ] Fortran-CUDA インターフェース
- [ ] メモリ管理: デバイス割り当て/解放

### 最適化機能
- [ ] 共有メモリの活用
- [ ] ストリーム並列実行
- [ ] カーネル融合（B_ijP計算）
- [ ] 非同期メモリ転送
- [ ] マルチGPU対応（NCCL）

### テスト・検証
- [ ] 単体テスト: 各カーネル
- [ ] 精度テスト: CPU版との比較（<1e-10誤差）
- [ ] 性能テスト: 加速率測定
- [ ] スケーリングテスト: 1/2/4 GPU
- [ ] 大規模テスト: 8-16原子フラグメント

---

## 💰 コスト対効果分析

### ハードウェアコスト

| GPU | 価格 | 性能 | コストパフォーマンス |
|-----|------|------|-------------------|
| **V100 (16GB)** | $5000 | 1.0x | 基準 |
| **A100 (40GB)** | $15000 | 2.5x | **最適** ✅ |
| **H100 (80GB)** | $30000 | 4.0x | 大規模系向け |

### 投資回収期間

**前提**: 
- 研究者の時間コスト: $50/時間
- 計算待ち時間削減: 1日4時間 → 5分（約50倍）

**計算**:
```
節約時間 = 4時間/日 × 0.95 = 3.8時間/日
年間節約 = 3.8時間 × 250日 = 950時間
金銭価値 = 950時間 × $50 = $47,500/年

A100投資: $15,000 × 2台 = $30,000
回収期間: 30,000 / 47,500 ≈ 7.6ヶ月 ✅
```

**結論**: **1年以内に投資回収可能**

---

## 🚀 次のステップ

### 即座の実装
1. **Phase 1開始**: Fortran-CUDA連携
2. **プロトタイプ**: 小規模系でのGPU版実装
3. **検証**: CPU版との精度・性能比較

### 将来展望
- **Tensor Core活用**: さらに2-4倍高速化
- **Mixed Precision**: FP16/TF32による高速化
- **動的負荷分散**: 不均一なフラグメントサイズへの最適化
- **クラウドGPU統合**: AWS/GCP/Azureでの大規模計算

---

**文書作成**: GitHub Copilot (Claude Sonnet 4.5)  
**最終更新**: 2026年2月23日  
**実装状態**: 計画段階（実装準備完了）
