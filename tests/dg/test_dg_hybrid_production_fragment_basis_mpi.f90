#include "config.h"
program test_dg_hybrid_production_fragment_basis_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
  use dg_hybrid_production_fragment_basis,only:build_dg_hybrid_production_fragment_basis
  implicit none
  integer,parameter::nglobal=4,nw=2,np=2
  integer::comm,rank,nproc,ierr,nlocal,p,j
  integer(int64),allocatable::point_ids(:),physical_point_ids(:)
  real(real64),allocatable::weights(:)
  complex(real64),allocatable::wannier(:,:),pw(:,:)
  type(s_dg_hybrid_fragment_basis)::basis
  integer::packet_ids(np),near_offsets(np+1),near_ids(2)
  integer(int64)::workspace,fingerprint
  logical::ok
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  nlocal=count([(mod(p-1,nproc)==rank,p=1,nglobal)])
  allocate(point_ids(nlocal),physical_point_ids(nlocal),weights(nlocal),wannier(nw,nlocal),pw(np,nlocal));j=0
  do p=1,nglobal
    if(mod(p-1,nproc)/=rank)cycle
    j=j+1;point_ids(j)=p;physical_point_ids(j)=10+p;weights(j)=1d0;wannier(:,j)=(0d0,0d0);pw(:,j)=(0d0,0d0)
    if(p==1)then;wannier(1,j)=1d0;pw(1,j)=1d0;endif
    if(p==2)then;wannier(2,j)=1d0;pw(2,j)=1d0;endif
    if(p==3)pw(1,j)=1d0
    if(p==4)pw(2,j)=1d0
  enddo
  packet_ids=[1,2];near_offsets=[1,2,3];near_ids=[1,2]
  call build_dg_hybrid_production_fragment_basis(comm,nglobal,point_ids,weights,wannier,pw,&
    packet_ids,near_offsets,near_ids,701_int64,709_int64,100_int64,1,1d-12,basis,workspace,fingerprint,ok,message,&
    physical_point_ids)
  call require(ok,trim(message));call require(basis%generation==1.and.basis%fragment_id==1,'fragment receipt mismatch')
  call require(size(basis%global_ids)==count([(mod(p-1,nproc)==rank,p=1,nw+np)]),&
    'production basis is not column-owner distributed')
  call require(all(basis%buffer_point_ids==[11_int64,12_int64,13_int64,14_int64]),&
    'production physical support IDs were not preserved')
  do j=1,size(basis%global_ids)
    p=int(basis%global_ids(j)-100_int64)
    call require(basis%sector(j)==merge(1,2,p<=nw),'production WF/PW sector mismatch')
    if(p==1)call require(abs(basis%buffer_values(1,j)-1d0)<1d-14,'WF1 value mismatch')
    if(p==2)call require(abs(basis%buffer_values(2,j)-1d0)<1d-14,'WF2 value mismatch')
    if(p==3)call require(abs(basis%buffer_values(3,j)-1d0)<1d-14.and.&
      sum(abs(basis%buffer_values(:,j)))<1d0+1d-14,'PW1 complement mismatch')
    if(p==4)call require(abs(basis%buffer_values(4,j)-1d0)<1d-14.and.&
      sum(abs(basis%buffer_values(:,j)))<1d0+1d-14,'PW2 complement mismatch')
  enddo
  call require(workspace>0_int64.and.fingerprint/=0_int64,'production basis receipts are empty')
  if(rank==0)write(*,'(a,i0,a,i0)')'HYBRID_PRODUCTION_FRAGMENT_BASIS ranks=',nproc,' fingerprint=',fingerprint
  if(rank==0)write(*,'(a,i0,a)')'PASS production fragment basis on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,text)
    logical,intent(in)::condition;character(*),intent(in)::text;logical::global_condition
    call MPI_Allreduce(condition,global_condition,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(.not.global_condition)then;if(rank==0)write(0,'(a)')trim(text);call MPI_Abort(comm,1,ierr);endif
  end subroutine require
end program test_dg_hybrid_production_fragment_basis_mpi
