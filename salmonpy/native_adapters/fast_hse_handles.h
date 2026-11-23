/* Copyright 2026 SALMONpy contributors. Apache-2.0. */
#ifndef SALMONPY_FAST_HSE_HANDLES_V1_H
#define SALMONPY_FAST_HSE_HANDLES_V1_H
#include <stdint.h>
#include "native_context.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Explicit host-MPI attachment is required. All calls except protocol are
 * collective, serial, main-thread operations on the same attached communicator.
 * Model inputs remain immutable; host supplies the model I/O directory per call. Positive tokens are opaque, not pointers.
 * handle initialization and single-state initialization are exclusive.
 * new starts from the original GS, not the current state of its live model token.
 * Single/shared calls select one immutable model; registry selects independent models.
 * Each context owns dt, capacity, numerical state, output label and transactions.
 * dt is atomic time; A is Cartesian A/c in atomic units at step endpoints.
 * Aold must exactly match the accepted endpoint (query info after init/restore).
 * j is charge-current density in atomic units, minus SALMON matter current.
 * energy is total energy in Hartree. info returns the accepted step/A/j.
 * Close consumes a token without detaching/finalizing MPI; a stale token stays
 * invalid even if its slot is reused. Heterogeneous models require registry; no threading.
 * 0=success; 11=mixed rank arguments; 12=invalid/stale handle or attachment;
 * 13=label buffer too small; 16=registry allocation failure. Other codes retain
 * existing fast-PR status meanings, including recoverable numerical failures.
 */
/* Host initializes MPI and passes MPI_Comm_c2f(comm); keep comm alive until
 * every handle is closed, then detach before host MPI finalization. */
int sp_attach(int64_t fortran_communicator);
int sp_detach(void);
/* Optional protocol1, rank-local completed total/init/density/Hartree seconds. */
int sp_fast_hse_legacy_timer_protocol(void);
int sp_fast_hse_legacy_timers(double seconds[4]);
int sp_fast_hse_timer_protocol(void);
int sp_fast_hse_handle_timers(int64_t token,double seconds[4]);
/* Additive read-only queries; no attachment or MPI initialization is required.
 * Missing required symbols reject; rebuild/use a matching current runtime. */
int sp_fast_hse_checkpoint_protocol(void);
int sp_fast_hse_recovery_protocol(void);
int sp_fast_hse_handle_protocol(void); /* version 1 */
int sp_fast_hse_handle_init(const char *path, double dt, int capacity,
                          int allow_symmetry, int64_t *handle, double j[3]);
/* Model allocator remains usable after all states close, before detach.
 * Requires a successful handle_init and model protocol1. */
/* Optional capability: immutable same-symmetry model with independent states. */
int sp_fast_hse_shared_symmetry_protocol(void);
/* Context initialization children retire at close; diagnostic live/slot counts. */
int sp_fast_hse_context_mpi_protocol(void);
int sp_fast_hse_mpi_resources(int values[2]);
int sp_fast_hse_model_protocol(void);
int sp_fast_hse_model_new(double dt, int capacity, int64_t *handle, double j[3]);
/* Explicit heterogeneous registry: negative model tokens, positive state tokens.
 * Same supported attached MPI group/main thread; host supplies each model I/O directory.
 * Model close stops new states while live children retain their inputs/resources.
 * Choose registry registration or single/shared initialization for an attachment lifetime. */
int sp_fast_hse_registry_protocol(void);
int sp_fast_hse_model_register(const char *path,double dt,int capacity,int symmetry,
                              int64_t *model,int64_t *state,double j[3]);
int sp_fast_hse_model_spawn(int64_t model,double dt,int capacity,int64_t *state,double j[3]);
int sp_fast_hse_model_close(int64_t model);
/* Optional retained-GS bridge: copies a value seed; later source GS close is safe.
 * path describes the same GS material; no GS restart file is read or generated.
 * GS boundary 0 parks RT settings, 1 restores them around a subsequent GS import. */
int sp_fast_hse_memory_gs_protocol(void);
int sp_fast_hse_model_from_gs(int64_t gs,const char *path,double dt,int capacity,int symmetry,
                            int64_t *model,int64_t *state,double j[3]);
int sp_fast_hse_gs_boundary(int mode);


int sp_fast_hse_handle_new(int64_t model, double dt, int capacity,
                         int64_t *handle, double j[3]);
/* Single-state unchecked update; symmetry admissibility remains. No rollback promise. */
int sp_fast_hse_step_performance(const double Aold[3],const double Anew[3],double j[3],double *energy);
/* Unchecked handle update; local generation bounds and symmetry are retained. */
int sp_fast_hse_handle_step_performance(int64_t token,const double Aold[3],const double Anew[3],double j[3],double *energy);
int sp_fast_hse_handle_step(int64_t handle, const double a_old[3],
                          const double a_new[3], double j[3], double *energy);
int sp_fast_hse_handle_info(int64_t handle, int *step, double a[3], double j[3]);
int sp_fast_hse_handle_history(int64_t handle, int values[8]);
/* Caller owns all buffers; no pointer is retained; consume outputs only on
 * success. Nonzero status is not universal rollback or MPI-failure recovery. */
/* transaction modes: 0 capture, 1 numeric rollback/close logs,
 * 2 reopen/flush/release (retained on failure), 3 reopen. Host restores log bytes.
 * checkpoint modes: 0 save, 1 fresh direct load, 2 stage, 3 commit, 4 cancel.
 * Checkpoint paths may differ by rank. Authenticity/identity remain host-owned.
 */
int sp_fast_hse_handle_transaction(int64_t handle, int mode);
int sp_fast_hse_handle_checkpoint(int64_t handle, const char *path, int mode);
int sp_fast_hse_handle_label(int64_t handle, char *buffer, int capacity);
int sp_fast_hse_handle_close(int64_t handle);
#ifdef __cplusplus
}
#endif
#endif
