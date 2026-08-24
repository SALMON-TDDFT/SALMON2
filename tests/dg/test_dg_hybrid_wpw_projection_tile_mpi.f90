#include "config.h"
program test_dg_hybrid_wpw_projection_tile_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_wannier_complement,only:compute_dg_hybrid_wannier_projection_tile
  implicit none
  integer,parameter::nglobal=6,nw=2,np=3
  integer::comm,rank,nproc,ierr,nlocal,p,j
  integer(int64),allocatable::row_ids(:)
  real(real64),allocatable::weights(:)
  complex(real64),allocatable::wannier(:,:),pw(:,:),coefficients(:,:)
  complex(real64)::wref(nw,nglobal),pref(np,nglobal),expected(nw,np)
  integer(int64)::workspace,fingerprint
  logical::ok
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  wref=(0d0,0d0);wref(1,1)=1d0;wref(2,2)=1d0
  pref=(0d0,0d0);pref(1,1)=0.5d0;pref(1,3)=1d0
  pref(2,2)=cmplx(0.25d0,0.1d0,real64);pref(2,4)=1d0
  pref(3,1)=0.2d0;pref(3,2)=-0.3d0;pref(3,5)=1d0
  expected=matmul(conjg(wref),transpose(pref))
  nlocal=count([(mod(p-1,nproc)==rank,p=1,nglobal)])
  allocate(row_ids(nlocal),weights(nlocal),wannier(nw,nlocal),pw(np,nlocal));j=0
  do p=1,nglobal
    if(mod(p-1,nproc)/=rank)cycle
    j=j+1;row_ids(j)=p;weights(j)=1d0;wannier(:,j)=wref(:,p);pw(:,j)=pref(:,p)
  enddo
  call compute_dg_hybrid_wannier_projection_tile(comm,nglobal,row_ids,weights,wannier,pw,&
    701_int64,709_int64,11,1d-12,coefficients,workspace,fingerprint,ok,message)
  call require(ok,trim(message));call require(all(shape(coefficients)==[nw,np]),'projection tile shape mismatch')
  call require(maxval(abs(coefficients-expected))<1d-14,'projection tile coefficient mismatch')
  call require(workspace>0_int64.and.fingerprint/=0_int64,'projection tile receipts are empty')
  if(rank==0)write(*,'(a,i0,a,i0)')'HYBRID_WPW_PROJECTION_TILE ranks=',nproc,' fingerprint=',fingerprint
  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid WPW projection tile on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,text)
    logical,intent(in)::condition;character(*),intent(in)::text;logical::global_condition
    call MPI_Allreduce(condition,global_condition,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(.not.global_condition)then;if(rank==0)write(0,'(a)')trim(text);call MPI_Abort(comm,1,ierr);endif
  end subroutine require
end program test_dg_hybrid_wpw_projection_tile_mpi
