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
// Fused orthogonal stencil plus V for one spin of a block of orbitals on a periodic, undivided grid:
// each block streams z over an x tile of 4 rows, staging planes D ahead and keeping z neighbours in registers.
#include <cuda_runtime.h>
#include <cuda_pipeline.h>
#include <cuda/barrier>
#include "cuda_stencil.h"

namespace {

constexpr int BY = 4;                  // output rows per block

struct Args {
    const double2* tr; double2* thr; const double* V;   // at (is1,is2,is3) of orbital 0
    int nxi, nyi, nzi;
    long long sx2, sx3, sx5, vx2, vx3;                  // element strides: psi y, z, orbital; V y, z
    double lap0, lx[4], mx[4], ly[4], my[4], lz[4], mz[4];
};

#define HPSI_TAP(L, M, A, B) { const double2 t0 = (A), t1 = (B); \
    v.x = (L) * (t0.x + t1.x) + v.x;  v.y = (L) * (t0.y + t1.y) + v.y; \
    u.x = (M) * (t0.x - t1.x) + u.x;  u.y = (M) * (t0.y - t1.y) + u.y; }

template <int TW, int D>
__global__ void __launch_bounds__(TW * BY, 1)
hpsi_stencil_kernel(const __grid_constant__ Args a)
{
    constexpr int RW = TW + 8;                         // staged row: 4 + TW + 4 columns
    constexpr int NHALO = 8 * BY + 8 * TW;             // x halo of the centre rows, then the y-halo rows
    constexpr int NE = (NHALO + TW * BY - 1) / (TW * BY);
    constexpr int NB = D + 2;                          // a thread past wait(n) knows all finished n-2
    static_assert(D >= 5 && D <= 8, "the ring refill waits on plane n+4");
    constexpr int PSTR = (BY + 8) * RW;
    constexpr int VSTR = BY * TW;
    extern __shared__ double2 planes[];
    double* vbuf = reinterpret_cast<double*>(planes + NB * PSTR);
#pragma nv_diag_suppress static_var_with_dynamic_init
    __shared__ cuda::barrier<cuda::thread_scope_block> bar;   // initialised by init() below
    const int tx = threadIdx.x, ty = threadIdx.y;
    const int tid = ty * TW + tx;
    const int nxi = a.nxi, nyi = a.nyi, nzi = a.nzi;
    const int x0 = blockIdx.y * TW;
    const int ybase = blockIdx.x * BY;
    const int b0 = ybase + ty;
    const int gx = (x0 + tx) % nxi;                    // columns past nx hold wrapped neighbours
    const bool live = (x0 + tx < nxi) && (b0 < nyi);
    const long long sx2 = a.sx2, sx3 = a.sx3, span = a.sx3 * nzi, vspan = a.vx3 * nzi;
    const double2 zero = make_double2(0.0, 0.0);
    if (tid == 0) init(&bar, TW * BY);
    __syncthreads();

    const int bw = b0 % nyi;                           // rows past ny only feed neighbours
    const long long korb = a.sx5 * (long long)blockIdx.z;
    const long long kxy  = (long long)gx + sx2 * bw;
    double2* po = a.thr + kxy + korb;

    int zc = nzi - 4;                                  // ring prologue: planes -4 .. +3
    const double2* pc = a.tr + kxy + korb + sx3 * (long long)zc;
    double2 r[11];
    #pragma unroll
    for (int j = 0; j < 8; ++j) {
        r[j] = live ? *pc : zero;
        pc += sx3; if (++zc == nzi) { zc = 0; pc -= span; }
    }

    const double2* ph[NE]; int hslot[NE]; bool hon[NE];
    #pragma unroll
    for (int j = 0; j < NE; ++j) {
        const int e = tid + j * TW * BY;
        hon[j] = e < NHALO;
        int row, c, gy;
        if (e < 8 * BY) {                              // x halo of centre row e/8
            const int k = e % 8;
            row = 4 + e / 8; c = k < 4 ? k : TW + k;
            gy = (ybase + e / 8) % nyi;
        } else {                                       // y-halo rows, tile columns
            const int f = hon[j] ? e - 8 * BY : 0, rr = f / TW;
            row = rr < 4 ? rr : BY + rr; c = 4 + f % TW;
            gy = (rr < 4 ? ybase - 4 + rr : ybase + BY + (rr - 4)) % nyi;
            if (gy < 0) gy += nyi;
        }
        int x = (x0 - 4 + c) % nxi; if (x < 0) x += nxi;
        hslot[j] = row * RW + c;
        ph[j] = a.tr + (long long)x + sx2 * gy + korb;
    }
    const int vofs = ty * TW + tx;
    const int cofs = (ty + 4) * RW + tx + 4;
    const double* pv = a.V + (long long)gx + a.vx2 * bw;
    const double2* pk = a.tr + kxy + korb;             // own centre: the ring refill reads it back
    int zh = 0;

    // halo, own centre and own V of plane zh into buffer B, then advance one plane
#define HPSI_ISSUE(B) do { \
        double2* plf = planes + (B) * PSTR; \
        _Pragma("unroll") for (int j = 0; j < NE; ++j) \
            if (hon[j]) __pipeline_memcpy_async(plf + hslot[j], ph[j], 16); \
        __pipeline_memcpy_async(plf + cofs, pk, 16); \
        if (live) __pipeline_memcpy_async(vbuf + (B) * VSTR + vofs, pv, 8); \
        __pipeline_commit(); \
        _Pragma("unroll") for (int j = 0; j < NE; ++j) ph[j] += sx3; \
        pk += sx3; pv += a.vx3; \
        if (++zh == nzi) { \
            zh = 0; \
            _Pragma("unroll") for (int j = 0; j < NE; ++j) ph[j] -= span; \
            pk -= span; pv -= vspan; \
        } \
    } while (0)
    #pragma unroll
    for (int p = 0; p < D; ++p) HPSI_ISSUE(p);

    __pipeline_wait_prior(D - 1);                      // own copies of plane 0
    cuda::barrier<cuda::thread_scope_block>::arrival_token tok = bar.arrive();

    int n = 0;
    // ring phase K: slot (K+j)%11 holds plane n-4+j
#define HPSI_STEP(K) do { \
        const int b = n % NB; \
        bar.wait(cuda::std::move(tok));                /* plane n visible, plane n-2 free */ \
        HPSI_ISSUE((n + D) % NB); \
        __pipeline_wait_prior(D - 4);                  /* own copies up to plane n+4 */ \
        tok = bar.arrive(); \
        r[((K) + 8) % 11] = planes[((n + 4) % NB) * PSTR + cofs]; \
        if (live) { \
            const double2* pl = planes + b * PSTR; \
            const double2 c = r[((K) + 4) % 11]; \
            const double vq = vbuf[b * VSTR + vofs]; \
            double2 v = zero, u = zero; \
            _Pragma("unroll") for (int i = 0; i < 4; ++i) \
                HPSI_TAP(a.lx[i], a.mx[i], pl[cofs + i + 1], pl[cofs - i - 1]); \
            double2 h = make_double2(vq * c.x + a.lap0 * c.x - 0.5 * v.x + u.y, \
                                     vq * c.y + a.lap0 * c.y - 0.5 * v.y - u.x); \
            v = zero; u = zero; \
            _Pragma("unroll") for (int i = 0; i < 4; ++i) \
                HPSI_TAP(a.ly[i], a.my[i], pl[cofs + (i+1) * RW], pl[cofs - (i+1) * RW]); \
            HPSI_TAP(a.lz[0], a.mz[0], r[((K) + 5) % 11], r[((K) + 3) % 11]); \
            HPSI_TAP(a.lz[1], a.mz[1], r[((K) + 6) % 11], r[((K) + 2) % 11]); \
            HPSI_TAP(a.lz[2], a.mz[2], r[((K) + 7) % 11], r[((K) + 1) % 11]); \
            HPSI_TAP(a.lz[3], a.mz[3], r[((K) + 8) % 11], r[((K) + 0) % 11]); \
            h.x = h.x - 0.5 * v.x + u.y; \
            h.y = h.y - 0.5 * v.y - u.x; \
            *po = h; \
        } \
        po += sx3; ++n; \
    } while (0)

    for (int zb = 0; zb < nzi; zb += 11) {
        HPSI_STEP(0);
        if (zb + 1  < nzi) HPSI_STEP(1);  else break;
        if (zb + 2  < nzi) HPSI_STEP(2);  else break;
        if (zb + 3  < nzi) HPSI_STEP(3);  else break;
        if (zb + 4  < nzi) HPSI_STEP(4);  else break;
        if (zb + 5  < nzi) HPSI_STEP(5);  else break;
        if (zb + 6  < nzi) HPSI_STEP(6);  else break;
        if (zb + 7  < nzi) HPSI_STEP(7);  else break;
        if (zb + 8  < nzi) HPSI_STEP(8);  else break;
        if (zb + 9  < nzi) HPSI_STEP(9);  else break;
        if (zb + 10 < nzi) HPSI_STEP(10); else break;
    }
    __pipeline_wait_prior(0);                          // copies issued past the last plane
    bar.wait(cuda::std::move(tok));
#undef HPSI_STEP
#undef HPSI_ISSUE
}
#undef HPSI_TAP

template <int TW, int D>
size_t smem_bytes()
{
    return (size_t)(D + 2) * ((size_t)(BY + 8) * (TW + 8) * sizeof(double2) + (size_t)BY * TW * sizeof(double));
}

template <int TW, int D>
cuda_stencil::Fit fit()
{
    return cuda_stencil::probe(hpsi_stencil_kernel<TW, D>, smem_bytes<TW, D>());
}

template <int TW, int D>
void run(const Args& a, int norb, cudaStream_t s)
{
    const size_t smem = smem_bytes<TW, D>();
    auto k = hpsi_stencil_kernel<TW, D>;
    cuda_stencil::check(cudaFuncSetAttribute(k, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)smem),
                        "hpsi_stencil", "cudaFuncSetAttribute");
    dim3 block(TW, BY), grid((a.nyi + BY - 1) / BY, (a.nxi + TW - 1) / TW, norb);
    k<<<grid, block, smem, s>>>(a);
    cuda_stencil::check(cudaGetLastError(), "hpsi_stencil", "launch");
    cuda_stencil::check(cudaStreamSynchronize(s), "hpsi_stencil", "cudaStreamSynchronize");
}

struct Variant { int tw, d; cuda_stencil::Fit (*fit)(); void (*run)(const Args&, int, cudaStream_t); };

#define HPSI_VARIANT(W, DD) {W, DD, fit<W, DD>, run<W, DD>}
const Variant kVariants[] = {
    HPSI_VARIANT(32, 7), HPSI_VARIANT(32, 6), HPSI_VARIANT(32, 5),
    HPSI_VARIANT(64, 7), HPSI_VARIANT(64, 6), HPSI_VARIANT(64, 5),
    HPSI_VARIANT(96, 7), HPSI_VARIANT(96, 6), HPSI_VARIANT(96, 5),
};
#undef HPSI_VARIANT

// One tile when x fits in 96, else 32-wide tiles.
int tile_width(int nxi)
{
    return nxi <= 96 ? 32 * ((nxi + 31) / 32) : 32;
}

// Deepest lookahead that fits without spilling, else the deepest that fits; Hopper
// tries D=6 first. Chosen once per tile width.
const Variant* choose(int tw)
{
    static const Variant* chosen[4] = {};
    static bool decided[4] = {};
    const int slot = tw / 32;
    if (decided[slot]) return chosen[slot];
    const bool hopper = cuda_stencil::device_attr(cudaDevAttrComputeCapabilityMajor) == 9;
    const int order[3] = {hopper ? 6 : 7, hopper ? 7 : 6, 5};
    const Variant* fallback = nullptr;
    for (int d : order)
        for (const Variant& v : kVariants) {
            if (v.tw != tw || v.d != d) continue;
            const cuda_stencil::Fit f = v.fit();
            if (!f.ok) continue;
            if (!f.spills && !chosen[slot]) chosen[slot] = &v;
            if (!fallback) fallback = &v;
        }
    if (!chosen[slot]) chosen[slot] = fallback;
    decided[slot] = true;
    return chosen[slot];
}

} // namespace

// Applies the stencil and returns cuda_stencil::kOk, or returns another status (the
// caller then runs the OpenACC kernels). lapt/nabt are (4,3), x then y then z; cfg
// receives the tile width and lookahead used.
extern "C" int hpsi_stencil(const void* tpsi, void* htpsi, const void* vloc,
                            int nxi, int nyi, int nzi, int norb,
                            long long sx2, long long sx3, long long sx5,
                            long long vx2, long long vx3,
                            double lap0, const double* lapt, const double* nabt,
                            unsigned long long stream, int* cfg)
{
    if (nxi < 5 || nyi < 1 || nzi < 5 || norb < 1 || norb > 65535) return cuda_stencil::kShape;
    static int device = -1, memory = -1;
    if (device < 0) device = cuda_stencil::device_ok();
    if (!device) return cuda_stencil::kDevice;
    if (memory < 0)
        memory = cuda_stencil::address_ok(tpsi) && cuda_stencil::address_ok(htpsi) && cuda_stencil::address_ok(vloc);
    if (!memory) return cuda_stencil::kMemory;
    const Variant* v = choose(tile_width(nxi));
    if (!v) return cuda_stencil::kNoFit;

    Args a;
    a.tr = (const double2*)tpsi; a.thr = (double2*)htpsi; a.V = (const double*)vloc;
    a.nxi = nxi; a.nyi = nyi; a.nzi = nzi;
    a.sx2 = sx2; a.sx3 = sx3; a.sx5 = sx5; a.vx2 = vx2; a.vx3 = vx3; a.lap0 = lap0;
    for (int i = 0; i < 4; ++i) {
        a.lx[i] = lapt[i]; a.ly[i] = lapt[4 + i]; a.lz[i] = lapt[8 + i];
        a.mx[i] = nabt[i]; a.my[i] = nabt[4 + i]; a.mz[i] = nabt[8 + i];
    }
    v->run(a, norb, (cudaStream_t)stream);
    cfg[0] = v->tw;
    cfg[1] = v->d;
    return cuda_stencil::kOk;
}
