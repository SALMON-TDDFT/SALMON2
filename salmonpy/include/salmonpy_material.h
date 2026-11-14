#ifndef SALMONPY_MATERIAL_H
#define SALMONPY_MATERIAL_H
/* Preview material ABI1. Atomic units, Ac=A/c -> charge-current density j.
 * Serial calls, no MPI collectives or process termination. The library owns
 * positive opaque handles; 0 is invalid. Configuration bytes are plugin-defined.
 * Destroy states before models. Errors leave accepted state unchanged; failed
 * create sets its output handle to0. Error text is optional, truncated and NUL
 * terminated when buffer/capacity permit. No C++ exceptions cross this ABI. */
#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
#define SALMONPY_MATERIAL_ABI 1

typedef int64_t salmonpy_handle;
enum salmonpy_status {
    SALMONPY_OK=0, SALMONPY_ARGUMENT=1, SALMONPY_HANDLE=2,
    SALMONPY_MEMORY=3, SALMONPY_NUMERICAL=4, SALMONPY_BUSY=5,
    SALMONPY_ABI=6, SALMONPY_INTERNAL=7
};

typedef struct salmonpy_material_api_v1 {
    int32_t abi_version;
    int32_t struct_size;
    const char *backend; /* Immutable NUL-terminated algorithm/version name. */
    int32_t input_lag;  /* Accepted-state steps, as in salmonpy.time/1. */
    int32_t uses_initial_A; /*0 or1; physical initial current is always required. */
    int32_t (*model_create)(const void *,size_t,salmonpy_handle *,char *,size_t);
    int32_t (*model_destroy)(salmonpy_handle,char *,size_t);
    int32_t (*state_create)(salmonpy_handle,double,double,const double[3],salmonpy_handle *,char *,size_t);
    int32_t (*state_observe)(salmonpy_handle,double *,double[3],char *,size_t);
    int32_t (*state_step)(salmonpy_handle,const double[3],double *,double[3],char *,size_t);
    int32_t (*state_destroy)(salmonpy_handle,char *,size_t);
} salmonpy_material_api_v1;

/* Query must check requested ABI and output capacity before writing the table.
 * ABI1 fields and callbacks must all be present. On success struct_size equals
 * sizeof(salmonpy_material_api_v1) for this platform. Host and plugin use the
 * platform C ABI; this is not a wire/file format. Callback metadata lives until
 * library unload. State/config/output buffers belong to the caller and must not
 * be retained (except the static backend/table metadata). All numerical double
 * values use IEEE754 binary64. Input/output buffers must not overlap.
 * Initial/returned time is t0+n*dt; dt is fixed and positive.
 * Output clock/current may be ignored on nonzero status. */
int32_t salmonpy_material_query(int32_t,size_t,salmonpy_material_api_v1 *);
#ifdef __cplusplus
}
#endif
#endif
