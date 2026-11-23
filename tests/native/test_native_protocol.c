/* Query-only C host: no solver state, MPI_Init or MPI_Finalize. */
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <inttypes.h>
#include <string.h>
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
static int host_operations(void) {
    int rank,size,nrank,nsize,active,finalized;
    char version[8192],short_version[1];
    int64_t communicator;
    if(MPI_Init(NULL,NULL)!=MPI_SUCCESS)return 50;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
    sp_library(version,sizeof(version),&active);
    if(!active || !version[0])return 51;
    short_version[0]='X';sp_library(short_version,1,&active);
    if(short_version[0]!='\0' || !active)return 52;
    /* A range rejection is local and must not consume the attachment. */
    if(sp_attach(INT64_MAX)!=-1)return 53;
    communicator=(int64_t)MPI_Comm_c2f(MPI_COMM_WORLD);
    if(sp_attach(communicator)!=0)return 54;
    if(sp_attach(communicator)==0)return 55;
    sp_info(&nrank,&nsize);
    if(nrank!=rank || nsize!=size)return 56;
    if(sp_sum((double)(rank+1))!=(double)size*(size+1)/2)return 57;
    sp_barrier();
    if(sp_detach()!=0)return 58;
    MPI_Finalized(&finalized);if(finalized)return 59;
    if(sp_attach(communicator)!=0 || sp_detach()!=0)return 60;
    if(MPI_Barrier(MPI_COMM_WORLD)!=MPI_SUCCESS)return 61;
    if(rank==0)puts("common host operations PASS");
    return MPI_Finalize()==MPI_SUCCESS?0:62;
}
int main(int argc,char **argv) {
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
    if(argc==2 && strcmp(argv[1],"--host-operations")==0)return host_operations();
    return 0;
}
