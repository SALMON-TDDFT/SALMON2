//
//  Copyright 2026 SALMON developers
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//
// Stencil part of calc_current on a periodic, undivided grid, streamed like hpsi_stencil.cu; block
// partials are summed in a fixed order, so the result does not depend on scheduling.
#include <cuda_runtime.h>
#include <cuda_pipeline.h>
#include "cuda_stencil.h"

namespace {

constexpr int BY = 4;                  // rows per block
constexpr int ND = 4;                  // one-sided stencil reach

struct Args {
    const double2* psi;                // at (is1,is2,is3) of orbital 0
    int nx, ny, nz;
    long long sx2, sx3, sx5;           // element strides: y, z, orbital
    double nab[3][ND];
};

struct alignas(32) Sum4 { double r, x, y, z; };

__device__ __forceinline__ double im_conj_mul(double2 c, double2 s) { return c.x * s.y - c.y * s.x; }

template <int TW, int D>
__global__ void __launch_bounds__(TW * BY)
current_kernel(Args a, Sum4* __restrict__ part)
{
    constexpr int NB = D + 1, TR = BY + ND, RW = TW + ND;
    static_assert(D >= ND + 1, "plane z+4 must be in flight before it is read");
    extern __shared__ double2 smem[];  // [NB][TR][RW]: column c holds x = (x0 + c) mod nx
    const int tx = threadIdx.x, ty = threadIdx.y;
    const int x0 = blockIdx.y * TW;
    const int y0 = blockIdx.x * BY;
    const bool own = x0 + tx < a.nx && y0 + ty < a.ny;
    const bool halo = ty < ND;
    const double2* __restrict__ p = a.psi + (size_t)blockIdx.z * a.sx5;
    const size_t yo = (size_t)((y0 + ty) % a.ny) * a.sx2;
    const size_t yh = (size_t)((y0 + BY + ty) % a.ny) * a.sx2;
    const size_t xc = (size_t)((x0 + tx) % a.nx);            // columns past nx hold wrapped neighbours
    const size_t xr = (size_t)((x0 + TW + tx) % a.nx);       // tx < 4: the 4 columns past the tile
    auto T = [&](int b, int r, int c) -> double2* { return smem + ((size_t)b * TR + r) * RW + c; };
    // Plane z (wrapped) into buffer z % NB; always commits one group so the
    // wait counts stay fixed.
    auto issue = [&](int z) {
        if (z < a.nz + ND) {
            const int b = z % NB;
            const size_t zo = (size_t)(z % a.nz) * a.sx3;
            __pipeline_memcpy_async(T(b, ty, tx), p + zo + yo + xc, 16);
            if (tx < ND) __pipeline_memcpy_async(T(b, ty, TW + tx), p + zo + yo + xr, 16);
            if (halo) {
                __pipeline_memcpy_async(T(b, BY + ty, tx), p + zo + yh + xc, 16);
                if (tx < ND) __pipeline_memcpy_async(T(b, BY + ty, TW + tx), p + zo + yh + xr, 16);
            }
        }
        __pipeline_commit();
    };

    for (int z = 0; z < D; ++z) issue(z);
    __pipeline_wait_prior(D - ND);
    __syncthreads();
    double2 r[ND];                     // plane z+k in slot (z+k) % 4
    #pragma unroll
    for (int k = 1; k < ND; ++k) r[k] = *T(k % NB, ty, tx);

    double rho = 0, gx = 0, gy = 0, gz = 0;
    for (int z0 = 0; z0 < a.nz; z0 += ND) {
        #pragma unroll
        for (int s = 0; s < ND; ++s) {
            const int z = z0 + s;
            if (z >= a.nz) break;
            __pipeline_wait_prior(D - ND - 1);   // plane z+4 landed for this thread
            __syncthreads();                     // ... for all; plane z-1's buffer is free
            issue(z + D);
            const int bz = z % NB;
            r[s] = *T((z + ND) % NB, ty, tx);
            const double2 c = *T(bz, ty, tx);
            double2 sx = make_double2(0, 0), sy = make_double2(0, 0), sz = make_double2(0, 0);
            #pragma unroll
            for (int k = 1; k <= ND; ++k) {
                const double2 nx_ = *T(bz, ty, tx + k);
                const double2 ny_ = *T(bz, ty + k, tx);
                const double2 nz_ = r[(s + k) % ND];
                sx.x += a.nab[0][k - 1] * nx_.x; sx.y += a.nab[0][k - 1] * nx_.y;
                sy.x += a.nab[1][k - 1] * ny_.x; sy.y += a.nab[1][k - 1] * ny_.y;
                sz.x += a.nab[2][k - 1] * nz_.x; sz.y += a.nab[2][k - 1] * nz_.y;
            }
            if (own) {
                rho += c.x * c.x + c.y * c.y;
                gx += im_conj_mul(c, sx);
                gy += im_conj_mul(c, sy);
                gz += im_conj_mul(c, sz);
            }
        }
    }
    __pipeline_wait_prior(0);

    constexpr int NW = TW * BY / 32;
    __shared__ Sum4 wsum[NW];
    #pragma unroll
    for (int o = 16; o > 0; o >>= 1) {
        rho += __shfl_down_sync(0xffffffffu, rho, o);
        gx += __shfl_down_sync(0xffffffffu, gx, o);
        gy += __shfl_down_sync(0xffffffffu, gy, o);
        gz += __shfl_down_sync(0xffffffffu, gz, o);
    }
    const int tid = ty * TW + tx;
    if ((tid & 31) == 0) wsum[tid >> 5] = Sum4{rho, gx, gy, gz};
    __syncthreads();
    if (tid == 0) {
        Sum4 t = wsum[0];
        for (int i = 1; i < NW; ++i) { t.r += wsum[i].r; t.x += wsum[i].x; t.y += wsum[i].y; t.z += wsum[i].z; }
        part[((size_t)blockIdx.z * gridDim.y + blockIdx.y) * gridDim.x + blockIdx.x] = t;
    }
}

// One block: each thread sums its orbitals' partials in block order and applies
// j = (kAc rho + BT 2g) rocc wtk; the block then reduces in a fixed tree.
__global__ void current_finish(const Sum4* __restrict__ part, int nblk, int norb,
                               const double* __restrict__ rocc, double wtk,
                               double k0, double k1, double k2, const double* __restrict__ bt,
                               double* __restrict__ out)
{
    constexpr int NT = 256;
    __shared__ double s[3][NT];
    double j0 = 0, j1 = 0, j2 = 0;
    for (int o = threadIdx.x; o < norb; o += NT) {
        Sum4 t{0, 0, 0, 0};
        for (int b = 0; b < nblk; ++b) {
            const Sum4 v = part[(size_t)o * nblk + b];
            t.r += v.r; t.x += v.x; t.y += v.y; t.z += v.z;
        }
        const double w = rocc[o] * wtk;
        const double g0 = 2 * t.x, g1 = 2 * t.y, g2 = 2 * t.z;     // bt is BT(3,3), column-major
        j0 += (k0 * t.r + bt[0] * g0 + bt[3] * g1 + bt[6] * g2) * w;
        j1 += (k1 * t.r + bt[1] * g0 + bt[4] * g1 + bt[7] * g2) * w;
        j2 += (k2 * t.r + bt[2] * g0 + bt[5] * g1 + bt[8] * g2) * w;
    }
    s[0][threadIdx.x] = j0; s[1][threadIdx.x] = j1; s[2][threadIdx.x] = j2;
    __syncthreads();
    for (int n = NT / 2; n > 0; n >>= 1) {
        if (threadIdx.x < n)
            for (int c = 0; c < 3; ++c) s[c][threadIdx.x] += s[c][threadIdx.x + n];
        __syncthreads();
    }
    if (threadIdx.x < 3) out[threadIdx.x] = s[threadIdx.x][0];
}

Sum4* g_part = nullptr;
size_t g_part_n = 0;
double* g_rocc = nullptr;              // rocc for this spin and k-point, copied per call
size_t g_rocc_n = 0;
double* g_small = nullptr;             // BT(3,3), then the 3 results

template <int TW, int D>
size_t smem_bytes()
{
    return (size_t)(D + 1) * (BY + ND) * (TW + ND) * sizeof(double2);
}

template <int TW, int D>
cuda_stencil::Fit fit()
{
    return cuda_stencil::probe(current_kernel<TW, D>, smem_bytes<TW, D>());
}

// Returns the number of block partials per orbital.
template <int TW, int D>
int run(const Args& a, int norb, cudaStream_t s)
{
    const size_t smem = smem_bytes<TW, D>();
    auto k = current_kernel<TW, D>;
    cuda_stencil::check(cudaFuncSetAttribute(k, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)smem),
                        "current_stencil", "cudaFuncSetAttribute");
    dim3 block(TW, BY), grid((a.ny + BY - 1) / BY, (a.nx + TW - 1) / TW, norb);
    k<<<grid, block, smem, s>>>(a, g_part);
    cuda_stencil::check(cudaGetLastError(), "current_stencil", "launch");
    return grid.x * grid.y;
}

struct Variant { int tw, d; cuda_stencil::Fit (*fit)(); int (*run)(const Args&, int, cudaStream_t); };

#define CUR_VARIANT(W, DD) {W, DD, fit<W, DD>, run<W, DD>}
const Variant kVariants[] = {
    CUR_VARIANT(32, 7), CUR_VARIANT(32, 6), CUR_VARIANT(32, 5),
    CUR_VARIANT(64, 7), CUR_VARIANT(64, 6), CUR_VARIANT(64, 5),
    CUR_VARIANT(96, 7), CUR_VARIANT(96, 6), CUR_VARIANT(96, 5),
};
#undef CUR_VARIANT

// 32-wide tiles. CURRENT_STENCIL_TW pins another width at compile time, so widths
// can be compared on the same shape.
int tile_width(int nx)
{
    (void)nx;
#ifdef CURRENT_STENCIL_TW
    return CURRENT_STENCIL_TW;
#else
    return 32;
#endif
}

// Deepest lookahead that fits without spilling, else the deepest that fits.
// Chosen once per tile width.
const Variant* choose(int tw)
{
    static const Variant* chosen[4] = {};
    static bool decided[4] = {};
    const int slot = tw / 32;
    if (decided[slot]) return chosen[slot];
    const Variant* fallback = nullptr;
    for (const Variant& v : kVariants) {             // listed deepest first
        if (v.tw != tw) continue;
        const cuda_stencil::Fit f = v.fit();
        if (!f.ok) continue;
        if (!f.spills && !chosen[slot]) chosen[slot] = &v;
        if (!fallback) fallback = &v;
    }
    if (!chosen[slot]) chosen[slot] = fallback;
    decided[slot] = true;
    return chosen[slot];
}

template <class T>
void reserve(T** buf, size_t* have, size_t need)
{
    if (need <= *have) return;
    if (*buf) cuda_stencil::check(cudaFree(*buf), "current_stencil", "cudaFree");
    cuda_stencil::check(cudaMalloc(buf, need * sizeof(T)), "current_stencil", "cudaMalloc");
    *have = need;
}

} // namespace

// Adds the orbital-summed stencil current of one spin and k-point to out3 and
// returns cuda_stencil::kOk, or returns another status (the caller then runs the
// OpenACC kernel). nab is coef_nab(4,3); bt is BT(3,3); rocc is rocc(io_s:io_e)
// for this k-point and spin; cfg receives the tile width and lookahead used.
extern "C" int current_stencil(const void* psi, int nx, int ny, int nz, int norb,
                               long long sx2, long long sx3, long long sx5,
                               const double* nab, const double* kac, const double* bt,
                               const double* rocc, double wtk, double* out3,
                               unsigned long long stream, int* cfg)
{
    if (nx < ND || ny < ND || nz < ND + 1 || norb < 1 || norb > 65535) return cuda_stencil::kShape;
    static int device = -1, memory = -1;
    if (device < 0) device = cuda_stencil::device_ok();
    if (!device) return cuda_stencil::kDevice;
    if (memory < 0) memory = cuda_stencil::address_ok(psi);
    if (!memory) return cuda_stencil::kMemory;
    const int tw = tile_width(nx);
    const Variant* v = choose(tw);
    if (!v) return cuda_stencil::kNoFit;

    Args a;
    a.psi = (const double2*)psi; a.nx = nx; a.ny = ny; a.nz = nz;
    a.sx2 = sx2; a.sx3 = sx3; a.sx5 = sx5;
    for (int d = 0; d < 3; ++d)
        for (int k = 0; k < ND; ++k) a.nab[d][k] = nab[k + ND * d];
    cudaStream_t s = (cudaStream_t)stream;

    size_t one = 0;
    reserve(&g_part, &g_part_n, (size_t)norb * ((ny + BY - 1) / BY) * ((nx + tw - 1) / tw));
    reserve(&g_rocc, &g_rocc_n, (size_t)norb);
    if (!g_small) reserve(&g_small, &one, 12);

    const int nblk = v->run(a, norb, s);
    cuda_stencil::check(cudaMemcpyAsync(g_small, bt, 9 * sizeof(double), cudaMemcpyDefault, s),
                        "current_stencil", "cudaMemcpyAsync");
    cuda_stencil::check(cudaMemcpyAsync(g_rocc, rocc, norb * sizeof(double), cudaMemcpyDefault, s),
                        "current_stencil", "cudaMemcpyAsync");
    current_finish<<<1, 256, 0, s>>>(g_part, nblk, norb, g_rocc, wtk, kac[0], kac[1], kac[2], g_small, g_small + 9);
    cuda_stencil::check(cudaGetLastError(), "current_stencil", "launch");
    double o[3];
    cuda_stencil::check(cudaMemcpyAsync(o, g_small + 9, 3 * sizeof(double), cudaMemcpyDefault, s),
                        "current_stencil", "cudaMemcpyAsync");
    cuda_stencil::check(cudaStreamSynchronize(s), "current_stencil", "cudaStreamSynchronize");
    out3[0] += o[0]; out3[1] += o[1]; out3[2] += o[2];
    cfg[0] = v->tw;
    cfg[1] = v->d;
    return cuda_stencil::kOk;
}
