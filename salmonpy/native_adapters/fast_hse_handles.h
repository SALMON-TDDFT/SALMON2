/* Copyright 2026 SALMONpy contributors. Apache-2.0. */
#ifndef SALMONPY_FAST_HSE_HANDLES_V1_H
#define SALMONPY_FAST_HSE_HANDLES_V1_H
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
/* Explicit host-MPI attachment is required. All calls except protocol are
 * collective, serial, main-thread operations on the same attached communicator.
 * Model inputs and cwd remain immutable. Positive tokens are opaque, not pointers.
 * init and the legacy singleton init are exclusive over a process lifetime.
 * new starts from the original GS, not the current state of its live model token.
 * One model initialization/process; multi-context symmetry is not yet supported.
 * Each context owns dt, capacity, numerical state, output label and transactions.
 * dt is atomic time; A is Cartesian A/c in atomic units at step endpoints.
 * Aold must exactly match the accepted endpoint (query info after init/restore).
 * j is charge-current density in atomic units, minus SALMON matter current.
 * energy is total energy in Hartree. info returns the accepted step/A/j.
 * Close consumes a token without detaching/finalizing MPI; a stale token stays
 * invalid even if its slot is reused. No heterogeneous/threaded guarantee.
 * 0=success; 11=mixed rank arguments; 12=invalid/stale handle or attachment;
 * 13=label buffer too small; 16=registry allocation failure. Other codes retain
 * existing fast-PR status meanings, including recoverable numerical failures.
 */
/* Host initializes MPI and passes MPI_Comm_c2f(comm); keep comm alive until
 * every handle is closed, then detach before host MPI finalization. */
int sp_attach(int64_t fortran_communicator);
int sp_detach(void);
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
int sp_fast_hse_handle_new(int64_t model, double dt, int capacity,
                         int64_t *handle, double j[3]);
int sp_fast_hse_handle_step(int64_t handle, const double a_old[3],
                          const double a_new[3], double j[3], double *energy);
int sp_fast_hse_handle_info(int64_t handle, int *step, double a[3], double j[3]);
int sp_fast_hse_handle_history(int64_t handle, int values[8]);
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
