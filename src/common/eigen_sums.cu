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
// calc_eigen_energy's sums in one pass, per orbital: A = sum Re(conj(t) h), B = sum Re(conj(t) s),
// C = sum Re(conj(t) (h - (s + V t))). Blocks reduce z planes; planes are added in z order.
#include "cuda_stencil.h"

namespace {

constexpr int NT = 256;          // threads per block
constexpr int NW = NT / 32;
constexpr int U = 4;             // points per thread in flight

struct Args {
    const double2* t; const double2* h; const double2* s;
    const double* v0; const double* v1;
    int npl, nz, nspin, nio;
    long long pz, ps, po, pk;    // orbital strides in points: plane, spin, orbital, k
    long long vz;                // V plane stride
    double* part;                // [slab][z][3]
};

struct Acc { double a, b, c; };

__device__ __forceinline__ void add(Acc& r, double2 t, double2 h, double2 s, double v)
{
    r.a = fma(t.x, h.x, fma(t.y, h.y, r.a));
    r.b = fma(t.x, s.x, fma(t.y, s.y, r.b));
    const double ux = h.x - (s.x + v * t.x), uy = h.y - (s.y + v * t.y);
    r.c = fma(t.x, ux, fma(t.y, uy, r.c));
}

// The orbitals are read once, so their loads are evict-first; V is shared by every slab.
__global__ void __launch_bounds__(NT) plane_sums(Args a)
{
    const int o = blockIdx.y, z = blockIdx.x;
    const int is = o % a.nspin, io = (o / a.nspin) % a.nio, ik = o / (a.nspin * a.nio);
    const long long off = z * a.pz + is * a.ps + io * a.po + ik * a.pk;
    const double2 *t = a.t + off, *h = a.h + off, *s = a.s + off;
    const double* v = (is == 0 ? a.v0 : a.v1) + z * a.vz;
    Acc r{0, 0, 0};
    int i = threadIdx.x;
    for (; i + (U - 1) * NT < a.npl; i += U * NT) {
        double2 tt[U], hh[U], ss[U]; double vv[U];
#pragma unroll
        for (int u = 0; u < U; ++u) {
            tt[u] = __ldcs(t + i + u * NT); hh[u] = __ldcs(h + i + u * NT); ss[u] = __ldcs(s + i + u * NT);
            vv[u] = __ldg(v + i + u * NT);
        }
#pragma unroll
        for (int u = 0; u < U; ++u) add(r, tt[u], hh[u], ss[u], vv[u]);
    }
    for (; i < a.npl; i += NT) add(r, __ldcs(t + i), __ldcs(h + i), __ldcs(s + i), __ldg(v + i));

    __shared__ double sh[3][NW];
    for (int d = 16; d > 0; d >>= 1) {
        r.a += __shfl_down_sync(0xffffffffu, r.a, d);
        r.b += __shfl_down_sync(0xffffffffu, r.b, d);
        r.c += __shfl_down_sync(0xffffffffu, r.c, d);
    }
    if (threadIdx.x % 32 == 0) { sh[0][threadIdx.x / 32] = r.a; sh[1][threadIdx.x / 32] = r.b; sh[2][threadIdx.x / 32] = r.c; }
    __syncthreads();
    if (threadIdx.x < 3) {
        double x = 0.0;
        for (int w = 0; w < NW; ++w) x += sh[threadIdx.x][w];
        a.part[((long long)o * a.nz + z) * 3 + threadIdx.x] = x;
    }
}

__global__ void slab_sums(const double* part, int nz, int n, double* out)
{
    const int o = blockIdx.x * blockDim.x + threadIdx.x;
    if (o >= n) return;
    double a = 0, b = 0, c = 0;
    for (int z = 0; z < nz; ++z) {
        const double* p = part + ((long long)o * nz + z) * 3;
        a += p[0]; b += p[1]; c += p[2];
    }
    out[3 * o] = a; out[3 * o + 1] = b; out[3 * o + 2] = c;
}

double* g_part = nullptr;
double* g_out = nullptr;
size_t g_cap = 0;

} // namespace

extern "C" {

// sums(3, nspin*nio*nik) receives (A, B, C) per slab, spin fastest, then orbital, then k.
// Returns a cuda_stencil::Status; on anything but kOk, sums is untouched.
int eigen_sums(const double* t, const double* h, const double* s, const double* v0, const double* v1,
               int npl, int nz, int nspin, int nio, int nik,
               long long pz, long long ps, long long po, long long pk, long long vz,
               double* sums, unsigned long long stream)
{
    using namespace cuda_stencil;
    static int memory = -1;
    if (!device_ok()) return kDevice;
    if (memory < 0)
        memory = address_ok(t) && address_ok(h) && address_ok(s) && address_ok(v0) && address_ok(v1);
    if (!memory) return kMemory;
    const int n = nspin * nio * nik;
    if (n == 0) return kOk;
    const size_t need = (size_t)n * nz * 3;
    if (need > g_cap) {
        if (g_part) {
            check(cudaFree(g_part), "eigen_sums", "cudaFree");
            check(cudaFree(g_out), "eigen_sums", "cudaFree");
        }
        check(cudaMalloc(&g_part, need * sizeof(double)), "eigen_sums", "cudaMalloc");
        check(cudaMalloc(&g_out, (size_t)n * 3 * sizeof(double)), "eigen_sums", "cudaMalloc");
        g_cap = need;
    }
    Args a{(const double2*)t, (const double2*)h, (const double2*)s, v0, v1, npl, nz, nspin, nio,
           pz, ps, po, pk, vz, g_part};
    cudaStream_t st = (cudaStream_t)stream;
    plane_sums<<<dim3(nz, n), NT, 0, st>>>(a);
    check(cudaGetLastError(), "eigen_sums", "plane_sums");
    slab_sums<<<(n + 127) / 128, 128, 0, st>>>(g_part, nz, n, g_out);
    check(cudaGetLastError(), "eigen_sums", "slab_sums");
    check(cudaMemcpyAsync(sums, g_out, (size_t)n * 3 * sizeof(double), cudaMemcpyDeviceToHost, st),
          "eigen_sums", "cudaMemcpyAsync");
    check(cudaStreamSynchronize(st), "eigen_sums", "cudaStreamSynchronize");
    return kOk;
}

} // extern "C"
