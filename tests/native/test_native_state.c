/* Functional TDDFT/fast-PR-HSE host. Existing GS only; no disk checkpoints. */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef SALMONPY_STATE_SYMMETRY
#define SALMONPY_STATE_SYMMETRY 0
#endif
#ifdef SALMONPY_FAST_HSE
#include "fast_hse_handles.h"
#define HOST_DT .001
#define REJECTED_OLD_A 1
/* Match the standard fixture locally while exercising the exact HSE C ABI. */
static int sp_tddft_handle_init(const char *path,double dt,int capacity,
                              const double a[3],const int dec[5],int64_t *h,double j[3]) {
    (void)a;(void)dec;
    return sp_fast_hse_handle_init(path,dt,capacity,SALMONPY_STATE_SYMMETRY,h,j);
}
static int sp_tddft_model_new(int64_t *h,double j[3]) {
    return sp_fast_hse_model_new(HOST_DT,16,h,j);
}
static int sp_tddft_handle_step(int64_t h,const double old[3],const double a[3],double j[3]) {
    double energy;
    return sp_fast_hse_handle_step(h,old,a,j,&energy);
}
#define sp_tddft_handle_info sp_fast_hse_handle_info
#define sp_tddft_handle_close sp_fast_hse_handle_close
#else
#include "standard_tddft_handles.h"
#define HOST_DT .003
#define REJECTED_OLD_A 10
#endif
static const char *work;
static void require(int got, int wanted) {
    if (got != wanted) { fprintf(stderr,"native status %d, expected %d\n",got,wanted); MPI_Abort(MPI_COMM_WORLD,90); }
}
#define CALL(expr,wanted) do { require(sp_io_directory_enter(work),0); int rc_=(expr); sp_io_directory_leave(); require(rc_,wanted); } while(0)
static void same(const double a[3],const double b[3]) {
    for(int k=0;k<3;k++) if(a[k]!=b[k]) MPI_Abort(MPI_COMM_WORLD,91);
}
static void drive(int i,int point,double a[3]) {
    double factor=point==0?1.0:-0.5;
    a[0]=factor*(i%5)/1048576.0;
    a[1]=-factor*(i%7)/2097152.0;
    a[2]=factor*(i%3)/4194304.0;
#if SALMONPY_STATE_SYMMETRY
    /* The retained fixture operations preserve z, so drive that direction. */
    a[0]=0.;a[1]=0.;
#endif
}
int main(int argc,char **argv) {
    int provided,rank,size,finalized,point,first,last,step,before_step,after_step;
    int64_t h[2]={0,0},fresh;
    double zero[3]={0,0,0},a[3],j[2][3],old[2][3]={{0}},initial[3],before_a[3],before_j[3],after_a[3],after_j[3],history[4][3];
    FILE *trace=NULL;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
    if(provided<MPI_THREAD_FUNNELED || argc!=3) MPI_Abort(MPI_COMM_WORLD,92);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
    work=argv[1];point=atoi(argv[2]);first=point?point-1:0;last=point?first:1;
    if(point<0 || point>2) MPI_Abort(MPI_COMM_WORLD,93);
    int dec[5]={1,size,1,1,1};
    require(sp_attach((int64_t)MPI_Comm_c2f(MPI_COMM_WORLD)),0);
    if(sp_io_directory_protocol()<2) MPI_Abort(MPI_COMM_WORLD,94);
    if(rank==0) { char path[8192];snprintf(path,sizeof(path),"%s/trace.tsv",work);trace=fopen(path,"wx");if(!trace)MPI_Abort(MPI_COMM_WORLD,95); }
    CALL(sp_tddft_handle_init("rt.inp",HOST_DT,16,zero,dec,&h[first],j[first]),0);
    for(int k=0;k<3;k++)initial[k]=j[first][k];
    if(!point) { CALL(sp_tddft_model_new(&h[1],j[1]),0);same(j[1],initial); }
    for(int i=1;i<=12;i++) for(int k=first;k<=last;k++) {
        int sibling=1-k;
        if(!point) CALL(sp_tddft_handle_info(h[sibling],&before_step,before_a,before_j),0);
        drive(i,k,a);
        if(i==3 && k==first) {
            double bad[3],rejected[3];for(int d=0;d<3;d++)bad[d]=old[k][d];bad[0]+=1.;
            CALL(sp_tddft_handle_step(h[k],bad,a,rejected),REJECTED_OLD_A);
            CALL(sp_tddft_handle_info(h[k],&step,after_a,after_j),0);
            require(step,i-1);same(after_a,old[k]);same(after_j,j[k]);
        }
        CALL(sp_tddft_handle_step(h[k],old[k],a,j[k]),0);
        for(int d=0;d<3;d++)old[k][d]=a[d];
        if(!point) {
            CALL(sp_tddft_handle_info(h[sibling],&after_step,after_a,after_j),0);
            require(after_step,before_step);same(after_a,before_a);same(after_j,before_j);
            if(k==0 && i<=4)for(int d=0;d<3;d++)history[i-1][d]=j[k][d];
        }
        if(trace)fprintf(trace,"%d %d %.17g %.17g %.17g %.17g %.17g %.17g\n",i,k+1,a[0],a[1],a[2],j[k][0],j[k][1],j[k][2]);
    }
    if(!point) {
        CALL(sp_tddft_handle_info(h[1],&before_step,before_a,before_j),0);
        CALL(sp_tddft_handle_close(h[0]),0);
        CALL(sp_tddft_handle_info(h[0],&step,after_a,after_j),12);
        CALL(sp_tddft_model_new(&fresh,after_j),0);same(after_j,initial);
        if(fresh==h[0])MPI_Abort(MPI_COMM_WORLD,96);
        CALL(sp_tddft_handle_info(h[0],&step,after_a,after_j),12);
        double previous[3]={0,0,0};
        for(int i=1;i<=4;i++) {
            drive(i,0,a);CALL(sp_tddft_handle_step(fresh,previous,a,after_j),0);
            for(int d=0;d<3;d++){if(fabs(after_j[d]-history[i-1][d])>1e-12)MPI_Abort(MPI_COMM_WORLD,97);previous[d]=a[d];}
        }
        CALL(sp_tddft_handle_close(fresh),0);
        CALL(sp_tddft_handle_info(h[1],&after_step,after_a,after_j),0);
        require(after_step,before_step);same(after_a,before_a);same(after_j,before_j);
        CALL(sp_tddft_handle_close(h[1]),0);
    } else {
        CALL(sp_tddft_handle_close(h[first]),0);
        CALL(sp_tddft_handle_info(h[first],&step,after_a,after_j),12);
    }
    require(sp_detach(),0);
    MPI_Finalized(&finalized);if(finalized)MPI_Abort(MPI_COMM_WORLD,98);
    MPI_Barrier(MPI_COMM_WORLD);
    if(trace)fclose(trace);
    MPI_Finalize();return 0;
}
