#include <stdio.h>
#include <complex>
#include <math.h>
#include <cuComplex.h>

extern "C" {
static const unsigned block_size = 256;
__global__ void update_kvector_nonlocalpt_kernel(
		// Output
		cuDoubleComplex* ppg_zekr_uV,
		// Input (size)
		const int ik_s,
		const int ik_e,
		const int natom,
		const int ppg_nlma,
		const int ppg_nps,
		const int ia_tbl_size,
		// Input (ptr)
		const int* const ppg_ia_tbl,
		const int* const ppg_mps,
		const double* const ppg_rxyz,
		const double* ppg_uv,
		const double* kAc
		) {
	const unsigned ik = (blockIdx.x / ppg_nlma) + ik_s - 1;
	if (ik >= ik_e) return;
	const unsigned ilma = blockIdx.x % ppg_nlma;

	const int iatom = ppg_ia_tbl[ilma] - 1;
	const int j_end = ppg_mps[iatom];

	const double kAc_x = kAc[0 + 3 * (ik - (ik_s - 1))];
	const double kAc_y = kAc[1 + 3 * (ik - (ik_s - 1))];
	const double kAc_z = kAc[2 + 3 * (ik - (ik_s - 1))];

	for (int j_offset = 0; j_offset < j_end; j_offset += block_size) {
		const int j = j_offset + threadIdx.x;
		if (j >= j_end) return;
		const double x = ppg_rxyz[0 + j * (3) + iatom * (3 * ppg_nps)];
		const double y = ppg_rxyz[1 + j * (3) + iatom * (3 * ppg_nps)];
		const double z = ppg_rxyz[2 + j * (3) + iatom * (3 * ppg_nps)];

		const double theta = kAc_x * x +kAc_y * y + kAc_z * z;
		const cuDoubleComplex ekr_conj = {cos(theta), -sin(theta)};
		const double tmp = ppg_uv[j + ilma * (ppg_nps)];
		const cuDoubleComplex res = {ekr_conj.x * tmp, ekr_conj.y * tmp};
		ppg_zekr_uV[j + ilma * (ppg_nps) + ik * (ppg_nps * ppg_nlma)] = res;
	}
}
void update_kvector_nonlocalpt_core(
		// Output
		cuDoubleComplex* ppg_zekr_uV,
		// Input (size)
		const int ik_s,
		const int ik_e,
		const int natom,
		const int ppg_nlma,
		const int ppg_nps,
		const int ia_tbl_size,
		// Input (ptr)
		const int* const ppg_ia_tbl,
		const int* const ppg_mps,
		const double* const ppg_rxyz,
		const double* ppg_uv,
		const double* kAc
		) {
	dim3 grid_dim((ik_e - (ik_s - 1) + 1) * ppg_nlma);
	dim3 block_dim(block_size);
	update_kvector_nonlocalpt_kernel<<<grid_dim, block_dim>>>(
			ppg_zekr_uV,
			ik_s, ik_e,
			natom,
			ppg_nlma, ppg_nps,
			ia_tbl_size,
			ppg_ia_tbl, ppg_mps,
			ppg_rxyz, ppg_uv,
			kAc
			);
	cudaDeviceSynchronize();
}
}
