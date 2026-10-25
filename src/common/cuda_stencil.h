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
// Host-side checks shared by the CUDA stencil kernels.
// A kernel that cannot run here reports a status and the caller uses OpenACC.
#pragma once
#include <cuda_runtime.h>
#include <cstddef>
#include <cstdio>
#include <cstdlib>

namespace cuda_stencil {

// Keep in sync with stencil_cuda.f90's report().
enum Status { kOk = 0, kShape = 1, kDevice = 2, kMemory = 3, kNoFit = 4 };

inline void check(cudaError_t e, const char* who, const char* what)
{
    if (e != cudaSuccess) {
        fprintf(stderr, "%s: %s: %s\n", who, what, cudaGetErrorString(e));
        abort();
    }
}

inline int device_attr(cudaDeviceAttr attr)
{
    int dev = 0, v = 0;
    check(cudaGetDevice(&dev), "cuda_stencil", "cudaGetDevice");
    check(cudaDeviceGetAttribute(&v, attr, dev), "cuda_stencil", "cudaDeviceGetAttribute");
    return v;
}

// cp.async needs sm_80 or newer.
inline bool device_ok()
{
    return device_attr(cudaDevAttrComputeCapabilityMajor) >= 8;
}

// The kernels dereference the caller's addresses: managed or device memory, or host
// memory the GPU reaches coherently (pageableMemoryAccess, e.g. Grace).
inline bool address_ok(const void* p)
{
    cudaPointerAttributes at;
    if (cudaPointerGetAttributes(&at, p) != cudaSuccess) {
        cudaGetLastError();
        return false;
    }
    if (at.type == cudaMemoryTypeManaged || at.type == cudaMemoryTypeDevice) return true;
    return device_attr(cudaDevAttrPageableMemoryAccess) != 0;
}

// Whether one kernel instance can run here: an image for this device, its shared memory
// within the opt-in limit. Spilling instances run but rank after those that do not.
struct Fit { bool ok; bool spills; };

template <class Kernel>
Fit probe(Kernel k, size_t dyn_smem)
{
    cudaFuncAttributes fa;
    if (cudaFuncGetAttributes(&fa, k) != cudaSuccess) {
        cudaGetLastError();
        return {false, false};
    }
    const size_t optin = (size_t)device_attr(cudaDevAttrMaxSharedMemoryPerBlockOptin);
    if (dyn_smem + fa.sharedSizeBytes > optin) return {false, false};
    return {true, fa.localSizeBytes > 0};
}

} // namespace cuda_stencil
