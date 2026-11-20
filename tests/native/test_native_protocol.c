/* Query-only C host: no solver state, MPI_Init or MPI_Finalize. */
#include <stdio.h>
#include <stdint.h>
#include <mpi.h>
#ifdef SALMONPY_ABI_FAST_HSE
#include "fast_hse_handles.h"
#define QUERY1 sp_fast_hse_handle_protocol
#define QUERY2 sp_fast_hse_model_protocol
#define QUERY3 sp_fast_hse_registry_protocol
#define QUERY4 sp_fast_hse_checkpoint_protocol
#else
#include "standard_tddft_handles.h"
#include "standard_tddft_models.h"
int sp_eh_checkpoint_protocol(void);
#define QUERY1 sp_tddft_handle_protocol
#define QUERY2 sp_tddft_model_protocol
#define QUERY3 sp_tddft_registry_protocol
#define QUERY4 sp_eh_checkpoint_protocol
#endif
int main(void) {
    int initialized, values[4];
    if (sizeof(int64_t) != 8) return 40;
    if (MPI_Initialized(&initialized) != MPI_SUCCESS || initialized) return 41;
    values[0] = QUERY1(); values[1] = QUERY2();
    values[2] = QUERY3(); values[3] = QUERY4();
    for (int i=0; i<4; ++i) if (values[i] != 1) return 42;
    if (MPI_Initialized(&initialized) != MPI_SUCCESS || initialized) return 43;
    printf("%d %d %d %d\n", values[0], values[1], values[2], values[3]);
    return 0;
}
