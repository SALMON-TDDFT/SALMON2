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
// Periodic Hartree potential with cuFFT on a rank that holds the whole grid: rhoG = FFT(rho)/N,
// Vh = Re IFFT(coef rhoG), with poisson_ft's sign and normalization.
#include "../common/cuda_stencil.h"
#include <cufft.h>

namespace {

constexpr int NT = 256;

__global__ void load_rho(const double* __restrict__ rho, double2* __restrict__ w, long long n)
{
    const long long i = (long long)blockIdx.x * NT + threadIdx.x;
    if (i < n) w[i] = make_double2(rho[i], 0.0);
}

__global__ void apply_coef(const double* __restrict__ coef, double2* __restrict__ zrhog,
                           double2* __restrict__ w, double inv_n, long long n)
{
    const long long i = (long long)blockIdx.x * NT + threadIdx.x;
    if (i < n) {
        const double2 g = make_double2(w[i].x * inv_n, w[i].y * inv_n);
        zrhog[i] = g;
        w[i] = make_double2(coef[i] * g.x, coef[i] * g.y);
    }
}

__global__ void store_vh(const double2* __restrict__ w, double* __restrict__ vh, long long n)
{
    const long long i = (long long)blockIdx.x * NT + threadIdx.x;
    if (i < n) vh[i] = w[i].x;
}

struct Plan { cufftHandle h; int nx, ny, nz; bool made; };
Plan g_plan{0, 0, 0, 0, false};

void fft_check(cufftResult r, const char* what)
{
    if (r != CUFFT_SUCCESS) {
        fprintf(stderr, "hartree_fft: %s: cuFFT error %d\n", what, (int)r);
        abort();
    }
}

} // namespace

extern "C" {

// rho, vh, coef, zrhog and work are contiguous nx*ny*nz arrays, x fastest.
// Returns a cuda_stencil::Status; on anything but kOk the caller's arrays are untouched.
int hartree_fft(const double* rho, double* vh, const double* coef, double* zrhog, double* work,
                int nx, int ny, int nz, unsigned long long stream)
{
    using namespace cuda_stencil;
    static int memory = -1;
    if (!device_ok()) return kDevice;
    if (memory < 0)
        memory = address_ok(rho) && address_ok(vh) && address_ok(coef) && address_ok(zrhog) && address_ok(work);
    if (!memory) return kMemory;
    if (!g_plan.made || g_plan.nx != nx || g_plan.ny != ny || g_plan.nz != nz) {
        if (g_plan.made) fft_check(cufftDestroy(g_plan.h), "cufftDestroy");
        g_plan.made = false;
        if (cufftPlan3d(&g_plan.h, nz, ny, nx, CUFFT_Z2Z) != CUFFT_SUCCESS) return kNoFit;
        g_plan = {g_plan.h, nx, ny, nz, true};
    }
    const long long n = (long long)nx * ny * nz;
    const unsigned nb = (unsigned)((n + NT - 1) / NT);
    cudaStream_t s = (cudaStream_t)stream;
    cufftDoubleComplex* w = (cufftDoubleComplex*)work;
    fft_check(cufftSetStream(g_plan.h, s), "cufftSetStream");
    load_rho<<<nb, NT, 0, s>>>(rho, (double2*)w, n);
    check(cudaGetLastError(), "hartree_fft", "load_rho");
    fft_check(cufftExecZ2Z(g_plan.h, w, w, CUFFT_FORWARD), "forward transform");
    apply_coef<<<nb, NT, 0, s>>>(coef, (double2*)zrhog, (double2*)w, 1.0 / (double)n, n);
    check(cudaGetLastError(), "hartree_fft", "apply_coef");
    fft_check(cufftExecZ2Z(g_plan.h, w, w, CUFFT_INVERSE), "inverse transform");
    store_vh<<<nb, NT, 0, s>>>((const double2*)w, vh, n);
    check(cudaGetLastError(), "hartree_fft", "store_vh");
    check(cudaStreamSynchronize(s), "hartree_fft", "cudaStreamSynchronize");
    return kOk;
}

} // extern "C"
