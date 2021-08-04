#include <cstdio>
#include <cstdlib>
#include <complex>
#include <cuComplex.h>
#include <cub/cub.cuh>

extern "C" {
    static constexpr unsigned block_size = 128;
    static constexpr unsigned Nd = 4;
    static_assert(block_size != 0 && (block_size & (block_size-1)) == 0, "block_size must be 2^n");

    void stencil_current_core_cpu(const int ik_s, const int ik_e, const int io_s, const int io_e, const double* const vec_k, const double* const vec_Ac,
                                  const int* const is_array, const int* const ie_array, const int* const is, const int* const ie, const int* const idx, const int* const idy,
                                  const int* const idz, const double* const nabt, const int ispin, const int im, const int spin_len, const std::complex<double>* const psi, const double* const BT,
                                  const double* const rocc, const double* const wtk, double* const jx, double* const jy, double* const jz) {
        const int xlen = ie_array[0] - is_array[0] + 1;
        const int ylen = ie_array[1] - is_array[1] + 1;
        const int zlen = ie_array[2] - is_array[2] + 1;
        int psi_index = xlen*ylen*zlen*(ispin - 1) + xlen*ylen*zlen*spin_len*(ik_e - ik_s + 1)*(io_e - io_s + 1)*(im - 1);
        for(int ik = ik_s-1; ik < ik_e; ik++) {
            double kAc[3];
            for(int i = 0; i < 3; i++) {
                kAc[i] = vec_k[i+3*ik] + vec_Ac[i];
            }
            for(int io = io_s-1; io < io_e; io++) {
                double rtmp = 0.0;
                std::complex<double> xtmp(0.0, 0.0);
                std::complex<double> ytmp(0.0, 0.0);
                std::complex<double> ztmp(0.0, 0.0);
                std::complex<double> cpsi;
                for(int iz = is[2]-1; iz < ie[2]; iz++) {
                    for(int iy = is[1]-1; iy < ie[1]; iy++) {
                        for(int ix = is[0]-1; ix < ie[0]; ix++) {
                            rtmp = rtmp + std::abs(psi[psi_index + ix + iy*xlen + iz*xlen*ylen])*std::abs(psi[psi_index + ix + iy*xlen + iz*xlen*ylen]);
                            cpsi = std::conj(psi[psi_index + ix + iy*xlen + iz*xlen*ylen]);
                            xtmp = xtmp + nabt[0] * cpsi * (psi[psi_index + (idx[ix + 1] - 1) + iy*xlen + iz*xlen*ylen])
                                        + nabt[1] * cpsi * (psi[psi_index + (idx[ix + 2] - 1) + iy*xlen + iz*xlen*ylen])
                                        + nabt[2] * cpsi * (psi[psi_index + (idx[ix + 3] - 1) + iy*xlen + iz*xlen*ylen])
                                        + nabt[3] * cpsi * (psi[psi_index + (idx[ix + 4] - 1) + iy*xlen + iz*xlen*ylen]);

                            ytmp = ytmp + nabt[4] * cpsi * (psi[psi_index + ix + (idy[iy + 1] - 1)*xlen + iz*xlen*ylen])
                                        + nabt[5] * cpsi * (psi[psi_index + ix + (idy[iy + 2] - 1)*xlen + iz*xlen*ylen])
                                        + nabt[6] * cpsi * (psi[psi_index + ix + (idy[iy + 3] - 1)*xlen + iz*xlen*ylen])
                                        + nabt[7] * cpsi * (psi[psi_index + ix + (idy[iy + 4] - 1)*xlen + iz*xlen*ylen]);

                            ztmp = ztmp + nabt[8]  * cpsi * (psi[psi_index + ix + iy*xlen + (idz[iz + 1] - 1)*xlen*ylen])
                                        + nabt[9]  * cpsi * (psi[psi_index + ix + iy*xlen + (idz[iz + 2] - 1)*xlen*ylen])
                                        + nabt[10] * cpsi * (psi[psi_index + ix + iy*xlen + (idz[iz + 3] - 1)*xlen*ylen])
                                        + nabt[11] * cpsi * (psi[psi_index + ix + iy*xlen + (idz[iz + 4] - 1)*xlen*ylen]);
                        }
                    }
                }
                psi_index += xlen*ylen*zlen*spin_len;
                double wrk1[3], wrk2[3], wrk3[3], wrk4[3];
                for(int i = 0; i < 3; i++) {
                    wrk1[i] = kAc[i]*rtmp;
                }
                wrk2[0] = xtmp.imag()*static_cast<double>(2);
                wrk2[1] = ytmp.imag()*static_cast<double>(2);
                wrk2[2] = ztmp.imag()*static_cast<double>(2);

                for(int i = 0; i < 3; i++) {
                    wrk3[i] = BT[i]*wrk2[0]+BT[i + 3]*wrk2[1]+BT[i + 6]*wrk2[2];
                }
                for(int i = 0; i < 3; i++) {
                    wrk4[i] = (wrk1[i] + wrk3[i])*rocc[io + ik*(io_e - io_s + 1)]*wtk[ik];
                }
                *jx += wrk4[0];
                *jy += wrk4[1];
                *jz += wrk4[2];
            }
        }
    }

    __global__ void stencil_current_kernel(double* const d_res, const cuDoubleComplex* const psi_data, const int xlen,
                                           const int ylen, const int xsize, const int ysize, const int zsize, const int maxlen, const double* const nabt, const int* const idx, const int* const idy,
                                           const int* const idz, const double nabt0, const double nabt1, const double nabt2, const double nabt3, const double nabt4, const double nabt5, const double nabt6,
                                           const double nabt7, const double nabt8, const double nabt9, const double nabt10, const double nabt11) {
        const int tid = blockIdx.x * blockDim.x + threadIdx.x;
        const int threadId = threadIdx.x;
        if (tid >= maxlen) {
            return;
        }
        typedef cub::BlockReduce<double, block_size> BlockReduce;
        __shared__ typename BlockReduce::TempStorage temp_storage;
        const int ix = tid % xsize;
        const int iy = (tid/xsize) % ysize;
        const int iz = (tid/(xsize*ysize)) % zsize;
        cuDoubleComplex tmp;
        double val, block_sum;
        const cuDoubleComplex cpsi = cuConj(psi_data[ix + iy*xlen + iz*xlen*ylen]);
        tmp = psi_data[ix + iy*xlen + iz*xlen*ylen]; // psi_data[ix][iy][iz]
        const double psi_abs = cuCabs(tmp);
        val = psi_abs*psi_abs;
        __syncthreads();
        block_sum = BlockReduce(temp_storage).Sum(val);
        if(threadId == 0) {
             atomicAdd(&d_res[0], block_sum);
        }

        tmp = psi_data[(idx[ix + 1] - 1) + iy*xlen + iz*xlen*ylen];
        val = nabt0*(cpsi.y*tmp.x + cpsi.x*tmp.y);
        tmp = psi_data[(idx[ix + 2] - 1) + iy*xlen + iz*xlen*ylen];
        val += nabt1*(cpsi.y*tmp.x + cpsi.x*tmp.y);
        tmp = psi_data[(idx[ix + 3] - 1) + iy*xlen + iz*xlen*ylen];
        val += nabt2*(cpsi.y*tmp.x + cpsi.x*tmp.y);
        tmp = psi_data[(idx[ix + 4] - 1) + iy*xlen + iz*xlen*ylen];
        val += nabt3*(cpsi.y*tmp.x + cpsi.x*tmp.y);
        __syncthreads();
        block_sum = BlockReduce(temp_storage).Sum(val);
        if(threadId == 0) {
             atomicAdd(&d_res[1], block_sum);
        }
        
        tmp = psi_data[ix+(idy[iy + 1] - 1)*xlen + iz*xlen*ylen];
        val = nabt4*(cpsi.y*tmp.x + cpsi.x*tmp.y);
        tmp = psi_data[ix + (idy[iy + 2] - 1)*xlen + iz*xlen*ylen];
        val += nabt5*(cpsi.y*tmp.x + cpsi.x*tmp.y);
        tmp = psi_data[ix + (idy[iy + 3] - 1)*xlen + iz*xlen*ylen];
        val += nabt6*(cpsi.y*tmp.x+cpsi.x*tmp.y);
        tmp = psi_data[ix + (idy[iy + 4] - 1)*xlen + iz*xlen*ylen];
        val += nabt7*(cpsi.y*tmp.x + cpsi.x*tmp.y);
        __syncthreads();
        block_sum = BlockReduce(temp_storage).Sum(val);
        if(threadId == 0) {
             atomicAdd(&d_res[2], block_sum);
        }
        
        tmp = psi_data[ix + iy*xlen + (idz[iz + 1] - 1)*xlen*ylen];
        val = nabt8*(cpsi.y*tmp.x + cpsi.x*tmp.y);
        tmp = psi_data[ix + iy*xlen + (idz[iz + 2] - 1)*xlen*ylen];
        val += nabt9*(cpsi.y*tmp.x + cpsi.x*tmp.y);
        tmp = psi_data[ix + iy*xlen + (idz[iz + 3] - 1)*xlen*ylen];
        val += nabt10*(cpsi.y*tmp.x + cpsi.x*tmp.y);
        tmp = psi_data[ix + iy*xlen + (idz[iz + 4] - 1)*xlen*ylen];
        val += nabt11*(cpsi.y*tmp.x + cpsi.x*tmp.y);
        __syncthreads();
        block_sum = BlockReduce(temp_storage).Sum(val);
        if(threadId == 0) {
             atomicAdd(&d_res[3], block_sum);
        }
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
        for(int ik = ik_s - 1; ik < ik_e; ik++) {
            double kAc[3];
            for(int i = 0; i < 3; i++) {
                kAc[i] = vec_k[i + 3*ik] + vec_Ac[i];
            }
            for(int io = io_s - 1; io < io_e; io++) {
                const int grid_size = ((maxlen+block_size - 1)/block_size);
                // gpu kernel
                cudaMemset(reinterpret_cast<void*>(d_res), static_cast<double>(0), sizeof(double)*4);
                stencil_current_kernel<<<grid_size, block_size, sizeof(double)*block_size>>>(
                        d_res, &psi[psi_index], xlen, ylen, xsize, ysize, zsize, maxlen, d_nabt, d_idx, d_idy, d_idz,
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
                    wrk4[i] = (wrk1[i] + wrk3[i])*rocc[io + ik*(io_e - io_s + 1)]*wtk[ik];
                }
                *jx += wrk4[0];
                *jy += wrk4[1];
                *jz += wrk4[2];
            }
        }
    }
}
