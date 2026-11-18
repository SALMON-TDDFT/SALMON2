/* SALMONpy retained standard GS protocol1. Host MPI, collective main-thread calls. */
#ifndef SALMONPY_STANDARD_GS_H
#define SALMONPY_STANDARD_GS_H
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
int sp_gs_protocol(void);
/* Immutable relative input and dependencies must be visible on every subgroup rank.
 * Create computes SCF. Unconverged SCF still succeeds and retains the last iterate.
 * Handles retain numerical/XC/MPI owners, with unique tokens invalidated at close. */
int sp_gs_create(const char *filename,int64_t *token);
/* stats: iterations, native convergence flag, atoms, bands.
 * values: total Hartree energy, residual, tolerance. */
int sp_gs_info(int64_t token,int stats[4],double values[3]);
/* q:1 density,2 occupations,3 eigenvalues,4 force,5 orbital. Fortran-order arrays;
 * shape/lower use7 slots; ndim selects meaningful slots. Local density/orbitals
 * may include halos; eigen/occupation/force are native replicated arrays. */
int sp_gs_array_info(int64_t token,int quantity,int *ndim,int64_t shape[7],int64_t lower[7]);
int sp_gs_read(int64_t token,int quantity,double *real_part,double *imag_part,int64_t length);
/* close destroys only this owner. Native stale/mixed tokens return10; Python
 * balanced close is idempotent. Caller closes every owner before detaching MPI. */
int sp_gs_close(int64_t token);
/* states,live child communicators,management slot high-water. */
int sp_gs_resources(int counts[3]);
#ifdef __cplusplus
}
#endif
#endif
