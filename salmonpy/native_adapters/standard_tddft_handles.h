/* Copyright 2026 SALMONpy contributors. Apache-2.0. */
#ifndef SALMONPY_STANDARD_TDDFT_HANDLES_H
#define SALMONPY_STANDARD_TDDFT_HANDLES_H
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
/* See standard_tddft_models.h for the additive explicit heterogeneous registry.
 * Experimental version1: one immutable periodic, non-symmetric, unpolarized,
 * fixed-ion semilocal model (PZ/TB-mBJ/Libxc-PZ/PBE), middlepoint/transverse.
 * dt, capacity, decomposition and initial A are common. New states read original
 * GS, not a sibling excitation. Only RT/current and energy streams are owned;
 * optional global-name analysis outputs are rejected. Legacy ABI unchanged.
 * Every call except protocol queries: main thread, same order on all subgroup
 * ranks, within one attachment lifetime. Host initializes/finalizes MPI. Close
 * every state before detach; don't copy owning native objects.
 * A is Cartesian A/c (atomic units). All handle j values are CHARGE current
 * density = minus SALMON matter current; legacy sp_tddft_* still return matter J.
 * Checkpoint is same-binary/input/topology, caller-authenticated, not portable.
 * 0 success,11 mixed-rank arguments,12 stale/invalid token or lifetime,
 * 13 small label buffer or exhausted propagation capacity,14 transaction conflict,
 * 15 unsupported profile,16 registry allocation failure.
 * Other status meanings follow existing native TDDFT, including safe numerical
 * rejection. No process/MPI/unreturned STOP/initialization-failure recovery.
 */
int sp_attach(int64_t fortran_communicator);
int sp_detach(void);
int sp_tddft_handle_protocol(void);
int sp_tddft_model_protocol(void);
int sp_tddft_handle_init(const char *input,double dt,int capacity,
                        const double initial_A[3],const int decomposition[5],
                        int64_t *token,double j[3]);
int sp_tddft_model_new(int64_t *token,double j[3]);
int sp_tddft_handle_step(int64_t token,const double Aold[3],const double Anew[3],double j[3]);
/* Additive read-only collective density integral; includes occupation/k weights and volume. */
int sp_tddft_electron_number(double *value);
/* Optional protocol1: rank-local completed total/init/density/Hartree seconds.
 * Handle-local diagnostics; not checkpoint state; queries do not switch banks. */
int sp_tddft_legacy_timer_protocol(void);
int sp_tddft_legacy_timers(double seconds[4]);
int sp_tddft_timer_protocol(void);
int sp_tddft_handle_timers(int64_t token,double seconds[4]);
int sp_tddft_handle_electron_number(int64_t token,double *value);
int sp_tddft_handle_info(int64_t token,int *step,double A[3],double j[3]);
/* Modes0 save,1 fresh load,2 stage,3 commit,4 cancel; paths are rank-local. */
int sp_tddft_handle_checkpoint(int64_t token,const char *path,int mode);
/* 0 capture,1 numeric rollback/close streams,3 reopen,2 release.
 * Host restores only this owner's stream bytes between1 and3. */
int sp_tddft_handle_transaction(int64_t token,int mode);
/* Prefix only; combine with immutable base_directory and _rt.data/_rt_energy.data. */
int sp_tddft_handle_label(int64_t token,char *label,int capacity);
int sp_tddft_handle_close(int64_t token);
#ifdef __cplusplus
}
#endif
#endif
