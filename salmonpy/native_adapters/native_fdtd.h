/* Copyright 2026 SALMONpy contributors. Apache-2.0. Current experimental FDTD ABI. */
#ifndef SALMONPY_NATIVE_FDTD_H
#define SALMONPY_NATIVE_FDTD_H
#include "native_context.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Caller-owned buffers, not retained. All state operations are sequential on
 * the main thread. Collective on field group except current_mask/plane_read;
 * protocol queries are local and require no field/MPI initialization.
 * indices: C [points][3] int, ONE-BASED global grid triples.
 * Fields: C [points][3] double = Fortran (3,points); follow COUPLED current,
 * normalization/sign and time rules. J is physical charge-current density j
 * in atomic units at previous half-step, no extra cell-volume factor.
 * Buffers must have the declared capacity.
 * Return0 success; otherwise operation-specific status, outputs not usable.
 * Protected rollback is a coordinator contract, not every nonzero return.
 */
int sp_eh_init(const char *path,int points,const int *indices);
int sp_eh_step(int points,const double *J,double *E);
int sp_eh_step_checked(int points,const double *J,const double *A,double *E);
/* J contains ascending mask-selected local rows; A/E contain all mapped rows. */
int sp_eh_step_owned(int points,int local_points,const double *J,const double *A,double *E);
int sp_eh_current_mask(int points,int *mask);
/* Accepted step>=1; raw macro-centered H at (n+1/2)dt; no advancement. */
int sp_eh_sample_h(int points,double *H);
/* Local only. info int[21]; values double[count]. count0 queries metadata.
 * axis1..3, local-grid offset; committed sampled field, no pending restore. */
int sp_eh_plane_read(int axis,int offset,int count,int *info,double *values);
int sp_eh_checkpoint_protocol(void);
/* Each rank has a distinct path; caller validates input/runtime/topology.
 * Modes0 save,1 load,2 stage,3 commit,4 cancel. No multi-file atomic publication. */
int sp_eh_checkpoint(const char *path,int mode);
int sp_eh_step_transaction_protocol(void);
/* Modes0 capture,1 rollback,2 release,3 no-op. Supported numerical field/LD
 * profiles only; unsupported output modes return15. No general MPI/OS recovery. */
int sp_eh_step_transaction(int mode);
int sp_eh_close(void);
#ifdef __cplusplus
}
#endif
#endif
