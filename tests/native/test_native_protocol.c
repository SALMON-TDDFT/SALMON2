/* Query-only C host: no solver state, MPI_Init or MPI_Finalize. */
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <inttypes.h>
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
    int initialized, values[4],i;
    int64_t host[12];
    const int64_t expected[10]={1,sizeof(int)*CHAR_BIT,sizeof(double)*CHAR_BIT,64,
                               sizeof(MPI_Fint)*CHAR_BIT,64,1,1,1,3};
    for(i=0;i<12;i++)host[i]=-91;
    if(sp_host_abi(9,host)!=13 || sp_host_abi(0,NULL)!=13) return 44;
    for(i=0;i<12;i++)if(host[i]!=-91)return 45;
    if(sp_host_abi(12,host)!=0)return 46;
    for(i=0;i<10;i++)if(host[i]!=expected[i])return 47;
    if(host[10]!=-91 || host[11]!=-91)return 48;
    if (sizeof(int64_t) != 8) return 40;
    if (MPI_Initialized(&initialized) != MPI_SUCCESS || initialized) return 41;
    values[0] = QUERY1(); values[1] = QUERY2();
    values[2] = QUERY3(); values[3] = QUERY4();
    for (int i=0; i<4; ++i) if (values[i] != 1) return 42;
    if (MPI_Initialized(&initialized) != MPI_SUCCESS || initialized) return 43;
    for(i=0;i<10;i++)printf("%" PRId64 " ",host[i]);
    printf("%d %d %d %d\n", values[0], values[1], values[2], values[3]);
    return 0;
}
