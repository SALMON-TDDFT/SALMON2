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
// Nonlocal back-projection, htpsi(v,io) += sum_k zekr_uV(j_k,ilma_k) uVpsibox(ilma_k,io), on tiles
// of consecutive support points. A tile's uVpsibox rows are staged in shared memory a chunk of
// orbitals ahead with cp.async, so each row is fetched once per tile rather than once per point.
#include "cuda_stencil.h"
#include <cuda_pipeline.h>
#include <algorithm>
#include <vector>

namespace {

constexpr int TP = 64;     // support points per tile
constexpr int G = 4;       // orbital groups per tile
constexpr int OB = 2;      // orbitals per thread per chunk
constexpr int RCAP = 96;   // uVpsibox rows per tile
constexpr int KMAX = 32;   // projector entries a point can hold
constexpr int NT = TP * G, OC = G * OB;
constexpr size_t kSmem = sizeof(double2) * (2 * OC * RCAP + KMAX * TP) + sizeof(int) * TP
                       + sizeof(unsigned short) * KMAX * TP;

// Built from ppg's fixed index tables; the caller rebuilds it when they change.
struct Plan {
    int max_vi = 0, ntile = 0;
    int* n = nullptr;        // entries per point
    int2* ent = nullptr;     // [k][vi]: zekr_uV offset, row slot within the tile
    int* rows = nullptr;     // each tile's ilma rows (0-based), concatenated
    int* rowoff = nullptr;   // ntile+1 offsets into rows
    int* tv0 = nullptr;      // ntile+1 first points of the tiles
};
Plan g_plan;

struct Args {
    const int* __restrict__ n; const int2* __restrict__ ent; const int* __restrict__ rows;
    const int* __restrict__ rowoff; const int* __restrict__ tv0; const int* __restrict__ v2j;
    const double2* __restrict__ zekr; const double2* __restrict__ ubox; double2* __restrict__ htpsi;
    int max_vi, nspin, nio, nik, lb1, lb2, lb3;
    long long ub2, ub3, ub4, ub5, z3, s2, s3, s4, s5, s6, s7, zoff0, hoff0;
};

// A block owns one tile of at most TP points and RCAP rows. Thread (p, g) accumulates point p
// for OB orbitals of each chunk of OC, adding the projector entries in their original order.
__global__ void __launch_bounds__(NT) back_projection_kernel(Args a)
{
    extern __shared__ double2 smem[];
    double2* U = smem;                                     // [2][OC][RCAP]
    double2* s_zk = U + 2 * OC * RCAP;                     // [KMAX][TP]
    int* s_n = reinterpret_cast<int*>(s_zk + KMAX * TP);   // [TP]
    unsigned short* s_slot = reinterpret_cast<unsigned short*>(s_n + TP); // [KMAX][TP]

    const int tile = blockIdx.x, t = threadIdx.x, p = t % TP, g = t / TP;
    int q = blockIdx.y;
    const int isp = q % a.nspin; q /= a.nspin;
    const int ik = q % a.nik, im = q / a.nik;
    const int v0 = a.tv0[tile], np = a.tv0[tile + 1] - v0;
    const int r0 = a.rowoff[tile], R = a.rowoff[tile + 1] - r0;
    const long long uo = isp * a.ub2 + ik * a.ub4 + im * a.ub5;

    auto stage = [&](int io0, int buf) {
        const int oc = min(OC, a.nio - io0);
        double2* Ub = U + buf * OC * RCAP;
        for (int o = 0; o < oc; ++o)
            for (int r = t; r < R; r += NT)
                __pipeline_memcpy_async(Ub + o * RCAP + r, a.ubox + uo + (io0 + o) * a.ub3 + a.rows[r0 + r], 16);
        __pipeline_commit();
    };
    stage(0, 0);

    const long long zoff = a.zoff0 + ik * a.z3;
    for (int i = t; i < TP; i += NT) s_n[i] = i < np ? a.n[v0 + i] : 0;
    __syncthreads();
    for (int i = t; i < KMAX * TP; i += NT) {
        const int k = i / TP, pp = i % TP;
        if (k < s_n[pp]) {
            const int2 e = a.ent[(long long)k * a.max_vi + v0 + pp];
            s_zk[i] = a.zekr[zoff + e.x];
            s_slot[i] = (unsigned short)e.y;
        }
    }
    const int n = s_n[p];
    double2* h = a.htpsi;
    if (n > 0) {
        const int* c = a.v2j + 3 * (v0 + p);
        h += (c[0] - a.lb1) + (long long)(c[1] - a.lb2) * a.s2 + (long long)(c[2] - a.lb3) * a.s3
           + isp * a.s4 + a.hoff0 + ik * a.s6 + im * a.s7;
    }

    int buf = 0;
    for (int io0 = 0; io0 < a.nio; io0 += OC, buf ^= 1) {
        const int oc = min(OC, a.nio - io0);
        if (io0 + OC < a.nio) {
            stage(io0 + OC, buf ^ 1);
            __pipeline_wait_prior(1);
        } else {
            __pipeline_wait_prior(0);
        }
        __syncthreads();
        if (n > 0) {
            const double2* Ub = U + buf * OC * RCAP + g * OB * RCAP;
            double2 acc[OB], hv[OB];
#pragma unroll
            for (int b = 0; b < OB; ++b) {
                acc[b] = make_double2(0.0, 0.0);
                if (g * OB + b < oc) hv[b] = h[(io0 + g * OB + b) * a.s5];
            }
            for (int k = 0; k < n; ++k) {
                const double2 z = s_zk[k * TP + p];
                const double2* Uk = Ub + s_slot[k * TP + p];
#pragma unroll
                for (int b = 0; b < OB; ++b) {
                    const double2 u = Uk[b * RCAP];
                    acc[b].x += u.x * z.x - u.y * z.y;
                    acc[b].y += u.x * z.y + u.y * z.x;
                }
            }
#pragma unroll
            for (int b = 0; b < OB; ++b) {
                if (g * OB + b < oc) {
                    hv[b].x += acc[b].x;
                    hv[b].y += acc[b].y;
                    h[(io0 + g * OB + b) * a.s5] = hv[b];
                }
            }
        }
        __syncthreads();
    }
}

void release(Plan& p)
{
    cudaFree(p.n); cudaFree(p.ent); cudaFree(p.rows); cudaFree(p.rowoff); cudaFree(p.tv0);
    p = Plan();
}

template <class T>
bool upload(T*& d, const std::vector<T>& h)
{
    if (h.empty()) return true;
    if (cudaMalloc(&d, sizeof(T) * h.size()) != cudaSuccess) {
        cudaGetLastError();
        return false;
    }
    cuda_stencil::check(cudaMemcpy(d, h.data(), sizeof(T) * h.size(), cudaMemcpyHostToDevice),
                        "back_projection", "cudaMemcpy");
    return true;
}

} // namespace

extern "C" {

// Builds the tile plan from ppg's tables, whose Fortran layouts are v2nlma(0:max_vi-1) and
// k2ilma, k2j(0:max_vi-1, kmax), 1-based values. Returns a cuda_stencil::Status.
int back_projection_plan(int max_vi, int kmax, const int* v2nlma, const int* k2ilma, const int* k2j, long long z2)
{
    using namespace cuda_stencil;
    release(g_plan);
    if (!device_ok()) return kDevice;
    if (kmax > KMAX) return kShape;
    if (kSmem > (size_t)device_attr(cudaDevAttrMaxSharedMemoryPerBlockOptin)) return kNoFit;
    std::vector<int> n(max_vi);
    int nl = 0;
    for (int v = 0; v < max_vi; ++v) {
        n[v] = std::max(0, v2nlma[v]);
        for (int k = 0; k < n[v]; ++k) nl = std::max(nl, k2ilma[(size_t)k * max_vi + v]);
    }
    if ((long long)nl * z2 > 2147483647LL) return kShape;  // zekr_uV offsets are stored as int
    // Consecutive points, closed at TP points or before a tile would exceed RCAP distinct rows.
    std::vector<int> tv0{0}, rows, rowoff{0}, tr, mark(nl, -1);
    std::vector<int2> ent((size_t)KMAX * max_vi, make_int2(0, 0));
    auto close_tile = [&](int v1) {
        std::sort(tr.begin(), tr.end());
        for (int v = tv0.back(); v < v1; ++v)
            for (int k = 0; k < n[v]; ++k) {
                const int il = k2ilma[(size_t)k * max_vi + v] - 1, j = k2j[(size_t)k * max_vi + v];
                const int s = (int)(std::lower_bound(tr.begin(), tr.end(), il) - tr.begin());
                ent[(size_t)k * max_vi + v] = make_int2((int)(il * z2 + j - 1), s);
            }
        rows.insert(rows.end(), tr.begin(), tr.end());
        rowoff.push_back((int)rows.size());
        tv0.push_back(v1);
        tr.clear();
    };
    for (int v = 0; v < max_vi; ++v) {
        const int tag = (int)tv0.size();
        int fresh = 0;
        for (int k = 0; k < n[v]; ++k) fresh += mark[k2ilma[(size_t)k * max_vi + v] - 1] != tag;
        if (v - tv0.back() == TP || (int)tr.size() + fresh > RCAP) close_tile(v);
        for (int k = 0; k < n[v]; ++k) {
            const int il = k2ilma[(size_t)k * max_vi + v] - 1;
            if (mark[il] != (int)tv0.size()) {
                mark[il] = (int)tv0.size();
                tr.push_back(il);
            }
        }
    }
    if (max_vi > 0) close_tile(max_vi);

    Plan& p = g_plan;
    p.max_vi = max_vi;
    p.ntile = (int)tv0.size() - 1;
    if (!(upload(p.n, n) && upload(p.ent, ent) && upload(p.rows, rows) && upload(p.rowoff, rowoff)
          && upload(p.tv0, tv0))) {
        release(p);
        return kMemory;
    }
    check(cudaFuncSetAttribute(back_projection_kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)kSmem),
          "back_projection", "cudaFuncSetAttribute");
    return kOk;
}

// Adds the back-projection to htpsi on the stream and waits for it. zoff0 and hoff0 are the
// 0-based offsets of the first k-point (and band group) in zekr_uV and htpsi.
int back_projection(const void* zekr, const void* ubox, void* htpsi, const int* v2j,
                    int nspin, int nio, int nik, int nim,
                    long long ub2, long long ub3, long long ub4, long long ub5, long long z3,
                    int lb1, int lb2, int lb3, long long s2, long long s3, long long s4, long long s5,
                    long long s6, long long s7, long long zoff0, long long hoff0, unsigned long long stream)
{
    using namespace cuda_stencil;
    static int memory = -1;
    if (memory < 0)
        memory = address_ok(zekr) && address_ok(ubox) && address_ok(htpsi) && address_ok(v2j);
    if (!memory) return kMemory;
    const Plan& p = g_plan;
    if (p.ntile == 0 || nio == 0) return kOk;
    Args a{p.n, p.ent, p.rows, p.rowoff, p.tv0, v2j, (const double2*)zekr, (const double2*)ubox,
           (double2*)htpsi, p.max_vi, nspin, nio, nik, lb1, lb2, lb3,
           ub2, ub3, ub4, ub5, z3, s2, s3, s4, s5, s6, s7, zoff0, hoff0};
    cudaStream_t s = (cudaStream_t)stream;
    back_projection_kernel<<<dim3(p.ntile, nspin * nik * nim), NT, kSmem, s>>>(a);
    check(cudaGetLastError(), "back_projection", "launch");
    check(cudaStreamSynchronize(s), "back_projection", "cudaStreamSynchronize");
    return kOk;
}

} // extern "C"
