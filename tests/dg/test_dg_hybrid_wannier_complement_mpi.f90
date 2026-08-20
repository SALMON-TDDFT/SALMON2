#include "config.h"
program test_dg_hybrid_wannier_complement_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_wannier_complement,only:project_dg_hybrid_wannier_complement
  implicit none
  integer,parameter::nrow=11,nw=3,np=2
  integer::comm,rank,nproc,ierr,i,p,nlocal,pos
  integer(int64),allocatable::row_ids(:),duplicate_ids(:)
  real(real64),allocatable::weights(:),duplicate_weights(:)
  complex(real64),allocatable::wannier(:,:),pw(:,:),projected(:,:),dense(:,:),duplicate_w(:,:),duplicate_p(:,:)
  integer::packet_ids(np),offsets(np+1),near_ids(3)
  real(real64)::x,pi,tail,defect,epsilon_tail
  integer(int64)::workspace,fingerprint,reference_fingerprint
  logical::ok
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  pi=acos(-1d0);nlocal=count([(mod(i-1,nproc)==rank,i=1,nrow)])
  allocate(row_ids(nlocal),weights(nlocal),wannier(nw,nlocal),pw(np,nlocal),dense(np,nlocal))
  pos=0;epsilon_tail=1d-10
  do i=nrow,1,-1
    if(mod(i-1,nproc)/=rank)cycle
    pos=pos+1;row_ids(pos)=i;weights(pos)=1d0
    x=2d0*pi*real(i-1,real64)/real(nrow,real64)
    wannier(1,pos)=exp(cmplx(0d0,x,real64))/sqrt(real(nrow,real64))
    wannier(2,pos)=exp(cmplx(0d0,2d0*x,real64))/sqrt(real(nrow,real64))
    wannier(3,pos)=exp(cmplx(0d0,3d0*x,real64))/sqrt(real(nrow,real64))
    pw(1,pos)=0.3d0*wannier(1,pos)+(0.4d0,0.1d0)*wannier(2,pos)+&
      epsilon_tail*wannier(3,pos)+exp(cmplx(0d0,4d0*x,real64))/sqrt(real(nrow,real64))
    pw(2,pos)=(-0.2d0,0.3d0)*wannier(3,pos)+&
      exp(cmplx(0d0,5d0*x,real64))/sqrt(real(nrow,real64))
    dense(1,pos)=exp(cmplx(0d0,4d0*x,real64))/sqrt(real(nrow,real64))
    dense(2,pos)=exp(cmplx(0d0,5d0*x,real64))/sqrt(real(nrow,real64))
  enddo
  packet_ids=[1,1];offsets=[1,3,4];near_ids=[1,2,3]
  call project_dg_hybrid_wannier_complement(comm,nrow,row_ids,weights,wannier,pw,101_int64,202_int64,packet_ids,&
    offsets,near_ids,.true.,1d-8,projected,tail,workspace,fingerprint,ok,message)
  call require(ok,trim(message));reference_fingerprint=fingerprint
  defect=maxval(abs(projected-dense));call require(defect<2d-9,'local complement differs from dense oracle')
  call require(tail>0d0.and.tail<1d-8,'omitted-tail receipt is incorrect')
  call require(workspace>0_int64.and.fingerprint/=0_int64,'complement receipts are invalid')

  pw(1,:)=pw(1,:)+(1d-4-epsilon_tail)*wannier(3,:)
  call project_dg_hybrid_wannier_complement(comm,nrow,row_ids,weights,wannier,pw,101_int64,202_int64,packet_ids,&
    offsets,near_ids,.true.,1d-8,projected,tail,workspace,fingerprint,ok,message)
  call require(.not.ok,'large omitted Wannier tail was accepted')
  pw(1,:)=pw(1,:)-(1d-4-epsilon_tail)*wannier(3,:)

  wannier(1,:)=2d0*wannier(1,:)
  call project_dg_hybrid_wannier_complement(comm,nrow,row_ids,weights,wannier,pw,101_int64,202_int64,packet_ids,&
    offsets,near_ids,.true.,1d-8,projected,tail,workspace,fingerprint,ok,message)
  call require(.not.ok,'nonorthonormal retained Wannier frame was accepted')
  wannier(1,:)=0.5d0*wannier(1,:)

  if(rank==0)then
    allocate(duplicate_ids(nlocal+1),duplicate_weights(nlocal+1),duplicate_w(nw,nlocal+1),duplicate_p(np,nlocal+1))
    duplicate_ids(1:nlocal)=row_ids;duplicate_ids(nlocal+1)=row_ids(1)
    duplicate_weights(1:nlocal)=weights;duplicate_weights(nlocal+1)=weights(1)
    duplicate_w(:,1:nlocal)=wannier;duplicate_w(:,nlocal+1)=wannier(:,1)
    duplicate_p(:,1:nlocal)=pw;duplicate_p(:,nlocal+1)=pw(:,1)
  else
    allocate(duplicate_ids(nlocal),duplicate_weights(nlocal),duplicate_w(nw,nlocal),duplicate_p(np,nlocal))
    duplicate_ids=row_ids;duplicate_weights=weights;duplicate_w=wannier;duplicate_p=pw
  endif
  call project_dg_hybrid_wannier_complement(comm,nrow,duplicate_ids,duplicate_weights,duplicate_w,duplicate_p,&
    101_int64,202_int64,&
    packet_ids,offsets,near_ids,.true.,1d-8,projected,tail,workspace,fingerprint,ok,message)
  call require(.not.ok,'same-rank duplicate spatial row was accepted')

  if(rank==0)then
    write(*,'(a,i0,a,i0)')'HYBRID_COMPLEMENT ranks=',nproc,' fingerprint=',reference_fingerprint
    write(*,'(a,i0,a)')'PASS hybrid Wannier complement on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label
    integer::lf,gf
    lf=merge(0,1,condition);call MPI_Allreduce(lf,gf,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gf/=0)error stop label
  end subroutine require
end program test_dg_hybrid_wannier_complement_mpi
