#include <cstdio>
#include <cstdlib>
#include <complex>
#include <cuComplex.h>
#include "array_index.h"
#include <thrust/device_vector.h>

extern "C" {
    static constexpr unsigned block_size = 128;
    static constexpr unsigned Nd = 4;
    static_assert(block_size != 0 && (block_size & (block_size-1)) == 0, "block_size must be 2^n");

    __device__ void thread_reduction(const int threadId, double* const d_res, double* const smem, const int blockDim) {
        for(int i = 1; i < blockDim; i *= 2) { // blockDim must be 2^n
            if ((threadId % (2*i)) == 0) {
                smem[threadId] += smem[threadId + i];
            }
            __syncthreads();
        }
        if(threadId == 0) {
            atomicAdd(d_res, smem[0]);
        }
    }

    __global__ void stencil_current_kernel(double* const d_res, const cuDoubleComplex* const psi_data, const int xlen,
                                           const int ylen, const int xoffset, const int yoffset, const int zoffset, const int xsize, const int ysize, const int zsize, const int maxlen,
                                           const double* const nabt, const int* const idx, const int* const idy, const int* const idz, const double nabt0, const double nabt1,
                                           const double nabt2, const double nabt3, const double nabt4, const double nabt5, const double nabt6, const double nabt7, const double nabt8,
                                           const double nabt9, const double nabt10, const double nabt11) {
        const int tid = blockIdx.x * blockDim.x + threadIdx.x;
        const int threadId = threadIdx.x;
        if (tid >= maxlen) {
            return;
        }
        extern __shared__ double smem[];
        const int ix = (tid % xsize) + xoffset;
        const int iy = ((tid/xsize) % ysize) + yoffset;
        const int iz = ((tid/(xsize*ysize)) % zsize) + zoffset;
        cuDoubleComplex tmp;
        const cuDoubleComplex cpsi = cuConj(psi_data[ARRAY_INDEX_3D(ix, iy, iz, 1, xlen, 1, ylen, 1)]);
        tmp = psi_data[ARRAY_INDEX_3D(ix, iy, iz, 1, xlen, 1, ylen, 1)]; // psi_data[ix][iy][iz]
        const double psi_abs = cuCabs(tmp);
        smem[threadId] = psi_abs*psi_abs;
        __syncthreads();
        thread_reduction(threadId, &d_res[0], smem, blockDim.x);

        tmp = psi_data[ARRAY_INDEX_3D(idx[ARRAY_INDEX_1D(ix+1, 1)], iy, iz, 1, xlen, 1, ylen, 1)];
        smem[threadId] = nabt0*(cpsi.y*tmp.x + cpsi.x*tmp.y);
        tmp = psi_data[ARRAY_INDEX_3D(idx[ARRAY_INDEX_1D(ix+2, 1)], iy, iz, 1, xlen, 1, ylen, 1)];
        smem[threadId] += nabt1*(cpsi.y*tmp.x + cpsi.x*tmp.y);
        tmp = psi_data[ARRAY_INDEX_3D(idx[ARRAY_INDEX_1D(ix+3, 1)], iy, iz, 1, xlen, 1, ylen, 1)];
        smem[threadId] += nabt2*(cpsi.y*tmp.x + cpsi.x*tmp.y);
        tmp = psi_data[ARRAY_INDEX_3D(idx[ARRAY_INDEX_1D(ix+4, 1)], iy, iz, 1, xlen, 1, ylen, 1)];
        smem[threadId] += nabt3*(cpsi.y*tmp.x + cpsi.x*tmp.y);
        __syncthreads();
        thread_reduction(threadId, &d_res[1], smem, blockDim.x);

        tmp = psi_data[ARRAY_INDEX_3D(ix, idy[ARRAY_INDEX_1D(iy+1, 1)], iz, 1, xlen, 1, ylen, 1)];
        smem[threadId] = nabt4*(cpsi.y*tmp.x + cpsi.x*tmp.y);
        tmp = psi_data[ARRAY_INDEX_3D(ix, idy[ARRAY_INDEX_1D(iy+2, 1)], iz, 1, xlen, 1, ylen, 1)];
        smem[threadId] += nabt5*(cpsi.y*tmp.x + cpsi.x*tmp.y);
        tmp = psi_data[ARRAY_INDEX_3D(ix, idy[ARRAY_INDEX_1D(iy+3, 1)], iz, 1, xlen, 1, ylen, 1)];
        smem[threadId] += nabt6*(cpsi.y*tmp.x+cpsi.x*tmp.y);
        tmp = psi_data[ARRAY_INDEX_3D(ix, idy[ARRAY_INDEX_1D(iy+4, 1)], iz, 1, xlen, 1, ylen, 1)];
        smem[threadId] += nabt7*(cpsi.y*tmp.x + cpsi.x*tmp.y);
        __syncthreads();
        thread_reduction(threadId, &d_res[2], smem, blockDim.x);

        tmp = psi_data[ARRAY_INDEX_3D(ix, iy, idz[ARRAY_INDEX_1D(iz+1, 1)], 1, xlen, 1, ylen, 1)];
        smem[threadId] = nabt8*(cpsi.y*tmp.x + cpsi.x*tmp.y);
        tmp = psi_data[ARRAY_INDEX_3D(ix, iy, idz[ARRAY_INDEX_1D(iz+2, 1)], 1, xlen, 1, ylen, 1)];
        smem[threadId] += nabt9*(cpsi.y*tmp.x + cpsi.x*tmp.y);
        tmp = psi_data[ARRAY_INDEX_3D(ix, iy, idz[ARRAY_INDEX_1D(iz+3, 1)], 1, xlen, 1, ylen, 1)];
        smem[threadId] += nabt10*(cpsi.y*tmp.x + cpsi.x*tmp.y);
        tmp = psi_data[ARRAY_INDEX_3D(ix, iy, idz[ARRAY_INDEX_1D(iz+4, 1)], 1, xlen, 1, ylen, 1)];
        smem[threadId] += nabt11*(cpsi.y*tmp.x + cpsi.x*tmp.y);
        __syncthreads();
        thread_reduction(threadId, &d_res[3], smem, blockDim.x);
    }

    void stencil_current_core_gpu(const int ik_s, const int ik_e, const int io_s, const int io_e, const double* const vec_k, const double* const vec_Ac,
                                  const int* const is_array, const int* const ie_array, const int* const is, const int* const ie, const int* const idx, const int* const idy,
                                  const int* const idz, const double* const nabt, const int ispin, const int im, const int spin_len, const cuDoubleComplex* const psi, 
                                  const double* const BT, const double* const rocc, const double* const wtk, double* const jx, double* const  jy, double* const jz) {
        const int xlen = ie_array[0] - is_array[0] + 1;
        const int ylen = ie_array[1] - is_array[1] + 1;
        const int zlen = ie_array[2] - is_array[2] + 1;
        const int xsize = ie[0] - is[0] + 1;
        const int ysize = ie[1] - is[1] + 1;
        const int zsize = ie[2] - is[2] + 1;
        const int maxlen = xsize*ysize*zsize;
        int psi_index = xlen*ylen*zlen*(ispin - 1) + xlen*ylen*zlen*spin_len*(ik_e - ik_s + 1)*(io_e - io_s + 1)*(im - 1);
        double res[4];
        double *d_nabt, *d_res;
        int *d_idx, *d_idy, *d_idz;

        cudaMalloc(reinterpret_cast<void**>(&d_nabt), sizeof(double)*12);
        cudaMalloc(reinterpret_cast<void**>(&d_idx), sizeof(int)*(ie[0] - is[0] + 1 + Nd));
        cudaMalloc(reinterpret_cast<void**>(&d_idy), sizeof(int)*(ie[1] - is[1] + 1 + Nd));
        cudaMalloc(reinterpret_cast<void**>(&d_idz), sizeof(int)*(ie[2] - is[2] + 1 + Nd));
        cudaMalloc(reinterpret_cast<void**>(&d_res), sizeof(double)*4);

        cudaMemcpy(d_nabt, nabt, sizeof(double)*12, cudaMemcpyHostToDevice);
        cudaMemcpy(d_idx, idx, sizeof(int)*(ie[0] - is[0] + 1 + Nd), cudaMemcpyHostToDevice);
        cudaMemcpy(d_idy, idy, sizeof(int)*(ie[1] - is[1] + 1 + Nd), cudaMemcpyHostToDevice);
        cudaMemcpy(d_idz, idz, sizeof(int)*(ie[2] - is[2] + 1 + Nd), cudaMemcpyHostToDevice);
        for(int ik = ik_s; ik <= ik_e; ik++) {
            double kAc[3];
            for(int i = 1; i <= 3; i++) {
                kAc[ARRAY_INDEX_1D(i, 1)] = vec_k[ARRAY_INDEX_2D(i, ik, 1, 3, 1)] + vec_Ac[ARRAY_INDEX_1D(i, 1)];
            }
            for(int io = io_s; io <= io_e; io++) {
                const int grid_size = ((maxlen+block_size - 1)/block_size);
                // gpu kernel
                cudaMemset(reinterpret_cast<void*>(d_res), static_cast<double>(0), sizeof(double)*4);
                stencil_current_kernel<<<grid_size, block_size, sizeof(double)*block_size>>>(
                        d_res, &psi[psi_index], xlen, ylen, is[0], is[1], is[2], xsize, ysize, zsize, maxlen, d_nabt, d_idx, d_idy, d_idz,
                        nabt[0], nabt[1], nabt[2], nabt[3], nabt[4], nabt[5], nabt[6], nabt[7], nabt[8], nabt[9], nabt[10], nabt[11]);
                cudaMemcpy(res, d_res, sizeof(double)*4, cudaMemcpyDeviceToHost);
                psi_index += xlen*ylen*zlen*spin_len;
                double wrk1[3], wrk2[3], wrk3[3], wrk4[3];
                for(int i = 0; i < 3; i++) {
                    wrk1[i] = kAc[i]*res[0];
                }
                wrk2[0] = res[1]*static_cast<double>(2);
                wrk2[1] = res[2]*static_cast<double>(2);
                wrk2[2] = res[3]*static_cast<double>(2);
                for(int i = 0; i < 3; i++) {
                    wrk3[i] = BT[i]*wrk2[0] + BT[i + 3]*wrk2[1] + BT[i + 6]*wrk2[2];
                }
                for(int i = 0; i < 3; i++) {
                    wrk4[i] = (wrk1[i] + wrk3[i])*rocc[ARRAY_INDEX_2D(io, ik, io_s, io_e, 1)]*wtk[ARRAY_INDEX_1D(ik, 1)];
                }
                *jx += wrk4[0];
                *jy += wrk4[1];
                *jz += wrk4[2];
            }
        }
    }
}
