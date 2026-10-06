#include "config.h"
program test_metric
  use omp_lib, only: omp_get_max_threads
  use exx_blas_threads, only: scope_entries,scope_active
  use mpi
  use exx_sparse_orbitals, only: s_sparse_orbitals,sparse_pack,sparse_clear
  use exx_ace, only: s_exx_ace,exx_ace_clear,exx_ace_ready
  use exx_orbitals, only: orbital_ace_build,orbital_ace_apply
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  type(s_exx_ace) :: reference,distributed,copy
  type(s_sparse_orbitals) :: sparse_u
  complex(8),allocatable :: u(:,:,:),w(:,:,:),t(:,:,:),a(:,:,:),b(:,:,:)
  integer :: ierr,rank,np,nr,ro,rr,cr,co,n,g,lo,hi,gs,ge,j,mode,status,packed,trial,tl,th,bytes,maxbytes
  real(8) :: error,total,tolerance
  character(32) :: arg
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
  n=7
  if(command_argument_count()>0)then
    call get_command_argument(1,arg);read(arg,*)n
  endif
  do mode=1,3
    nr=1
    if(mode==2)nr=np
    if(mode==3.and.mod(np,2)==0)nr=2
    rr=mod(rank,nr);ro=rank/nr
    call MPI_Comm_split(MPI_COMM_WORLD,ro,rr,cr,ierr)
    call MPI_Comm_split(MPI_COMM_WORLD,rr,ro,co,ierr)
    lo=ro*n/(np/nr)+1;hi=(ro+1)*n/(np/nr)
    gs=rr*(n+2)/nr+1;ge=(rr+1)*(n+2)/nr
    ! A different target count tests streamed action independently of training.
    tl=ro*3/(np/nr)+1;th=(ro+1)*3/(np/nr)
    allocate(u(ge-gs+1,hi-lo+1,1),w(ge-gs+1,hi-lo+1,1))
    allocate(t(ge-gs+1,th-tl+1,1),a(ge-gs+1,th-tl+1,1),b(ge-gs+1,th-tl+1,1))
    do j=tl,th
      do g=gs,ge
        t(g-gs+1,j-tl+1,1)=cmplx(cos(dble(g*j)),sin(dble(g+j)),8)
      enddo
    enddo
    do trial=1,6
      u=0d0
      do j=lo,hi
        do g=gs,ge
          if(trial/=3)u(g-gs+1,j-lo+1,1)=cmplx(sin(dble(g*j)),cos(dble(g+2*j)),8)*0.01d0
          if(g==j)u(g-gs+1,j-lo+1,1)=u(g-gs+1,j-lo+1,1)+1d0
          if(trial==3.and.g==n)u(g-gs+1,j-lo+1,1)=u(g-gs+1,j-lo+1,1)*sqrt(2d-12)
          w(g-gs+1,j-lo+1,1)=-u(g-gs+1,j-lo+1,1)
          if(trial/=3)w(g-gs+1,j-lo+1,1)=dble(g)*w(g-gs+1,j-lo+1,1)
        enddo
      enddo
      if(trial==2)w=0d0
      if(trial==4.and.hi==n.and.hi>=lo)w(:,hi-lo+1,1)=0d0
      if(trial==5.and.lo<=1.and.hi>=1)w(:,1-lo+1,1)=cmplx(0d0,1d0,8)*w(:,1-lo+1,1)
      if(trial==6.and.rank==np-1.and.size(w)>0)w(1,1,1)=ieee_value(0d0,ieee_quiet_nan)
      do packed=0,1
        call orbital_ace_build(reference,u,w,0.5d0,cr,co,status,packed=packed==1)
        if(trial<=3)call check(status==0,'reference build')
        if(trial>3)call check(status/=0,'reference reject')
        call orbital_ace_build(distributed,u,w,0.5d0,cr,co,status,packed=packed==1,comm_matrix=MPI_COMM_WORLD)
        if(trial>3)then
          call check(status/=0,'distributed reject')
          call check(.not.exx_ace_ready(distributed),'failed state cleared')
          cycle
        endif
        call check(status==0,'distributed build')
#ifdef USE_SCALAPACK
        if(np>1.and.packed==1)then
          call check(distributed%metric_distributed,'distributed path exercised')
          bytes=16*size(distributed%metric_factor)
          call MPI_Allreduce(bytes,maxbytes,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
          call check(maxbytes<16*n*n,'factor not replicated')
          if(rank==0.and.mode==1.and.trial==1)write(*,*) 'TILE bytes N P replicated/max',n,np,16*n*n,maxbytes
        endif
#endif
        call orbital_ace_apply(reference,t,a,cr,co,status)
        call check(status==0,'reference apply')
        copy=distributed
        call exx_ace_clear(distributed)
        call orbital_ace_apply(copy,t,b,cr,co,status)
        call check(status==0,'copied distributed apply')
        error=0d0
        if(size(a)>0)error=maxval(abs(a-b))
        call MPI_Allreduce(error,total,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
        tolerance=1d-10
        if(trial==3)tolerance=1d-9
        call check(total<tolerance,'action equality')
        if(rank==0)write(*,*) 'PASS N/layout/trial/packed/error',n,mode,trial,packed,total
        call exx_ace_clear(copy)
        call sparse_pack(sparse_u,u(:,:,1))
        call orbital_ace_build(distributed,w,w,0.5d0,cr,co,status,packed=packed==1, &
          comm_matrix=MPI_COMM_WORLD,sparse_u=sparse_u)
        call check(status==0,'sparse training build')
        call orbital_ace_apply(distributed,t,b,cr,co,status)
        call check(status==0,'sparse training apply')
        error=0d0
        if(size(a)>0)error=maxval(abs(a-b))
        call MPI_Allreduce(error,total,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
        call check(total<tolerance,'sparse training equality')
        call exx_ace_clear(distributed)
        call sparse_clear(sparse_u)
      enddo
    enddo
    deallocate(u,w,t,a,b)
    call MPI_Comm_free(cr,ierr);call MPI_Comm_free(co,ierr)
  enddo
  call MPI_Finalize(ierr)
contains
  subroutine check(ok,message)
    implicit none
    logical,intent(in) :: ok
    character(*),intent(in) :: message
    if(scope_active)error stop 'unrestored ACE thread scope'
    if(omp_get_max_threads()>1.and.scope_entries==0)error stop 'ACE thread scope not entered'
    if(ok)return
    write(*,*) 'FAIL ',rank,message
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  end subroutine
end program
