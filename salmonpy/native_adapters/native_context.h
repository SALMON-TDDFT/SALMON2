/* Copyright 2026 SALMONpy contributors. Apache-2.0. Experimental native ABI. */
#ifndef SALMONPY_NATIVE_CONTEXT_H
#define SALMONPY_NATIVE_CONTEXT_H
#ifdef __cplusplus
extern "C" {
#endif
/* Query before use. Version2/3 support cwd-free standard/HSE operations;
 * native E/H requires3. This does not initialize or attach MPI.
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
