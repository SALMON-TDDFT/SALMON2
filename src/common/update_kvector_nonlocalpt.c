#include <stdio.h>
#include <complex.h>
#include <math.h>

void update_kvector_nonlocalpt_kernel(
		// Output
		double complex* ppg_zekr_uV,
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
	for (int ik = ik_s - 1; ik < ik_e; ik++) {
		for (int ilma = 0; ilma < ppg_nlma; ilma++) {
			const int iatom = ppg_ia_tbl[ilma] - 1;
			for (int j = 0; j < ppg_mps[iatom]; j++) {
				const double x = ppg_rxyz[0 + j * (3) + iatom * (3 * ppg_nps)];
				const double y = ppg_rxyz[1 + j * (3) + iatom * (3 * ppg_nps)];
				const double z = ppg_rxyz[2 + j * (3) + iatom * (3 * ppg_nps)];

				const double complex ekr = cexp(
						I * (
							kAc[0 + 3 * (ik - (ik_s - 1))] * x +
							kAc[1 + 3 * (ik - (ik_s - 1))] * y +
							kAc[2 + 3 * (ik - (ik_s - 1))] * z)
						);
				ppg_zekr_uV[j + ilma * (ppg_nps) + ik * (ppg_nps * ppg_nlma)] = conj(ekr) * ppg_uv[j + ilma * (ppg_nps)];
			}
		}
	}
}
