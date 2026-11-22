/* Copyright 2026 SALMONpy contributors. Apache-2.0. Experimental native ABI. */
#ifndef SALMONPY_NATIVE_CONTEXT_H
#define SALMONPY_NATIVE_CONTEXT_H
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
/* Host descriptor version1, returned locally without MPI calls/allocation.
 * Supply at least10 int64_t words. Success0 writes only words0..9; short
 * capacity returns13 without writing. No prior-version rescue is provided.
 * Words: version,C int bits,C double bits,int64 bits,Fortran INTEGER bits,
 * real(kind=8) bits,communicator encoding,host MPI ownership,thread policy,
 * directory protocol required for I/O operations (not an I/O capability claim).
 * Current values:1,32,64,64,32,64,1,1,1,3.
 * Encoding1 means MPI_Comm_c2f(comm) widened to int64_t, never raw MPI_Comm.
 * Ownership1: host initializes/finalizes MPI. Thread1: main thread, sequential.
 * A match establishes host layout/rules only; query operation protocols too.
 */
enum { SALMONPY_HOST_ABI_VERSION=1, SALMONPY_HOST_ABI_WORDS=10 };
int sp_host_abi(int capacity,int64_t *values);
/* Current directory protocol3 is required. Query does not initialize/attach MPI.
 * Enter an existing absolute POSIX work directory on the main thread, then
 * leave after that call. Nested enter is rejected; never leave another owner's
 * scope after a failed enter. It does not change the host's current directory. */
int sp_io_directory_protocol(void);
int sp_io_directory_enter(const char *directory);
void sp_io_directory_leave(void);
#ifdef __cplusplus
}
#endif
#endif
