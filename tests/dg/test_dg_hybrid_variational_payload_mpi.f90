#include "config.h"
program test_dg_hybrid_variational_payload_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_variational_payload,only:s_dg_hybrid_fixed_payload,s_dg_hybrid_variational_iterate,&
    freeze_dg_hybrid_variational_payload,compose_dg_hybrid_variational_hamiltonian
  implicit none
  integer,parameter::n=4
  integer::comm,rank,nproc,ierr,i,j,nrow,p
  integer(int64),allocatable::rows(:)
  complex(real64),allocatable::s(:,:),t(:,:),vnl(:,:),sipg(:,:),vlocal(:,:)
  complex(real64)::dense_s(n,n),dense_t(n,n),dense_nl(n,n),dense_i(n,n),dense_v(n,n)
  type(s_dg_hybrid_fixed_payload)::fixed
  type(s_dg_hybrid_variational_iterate)::zero,full
  logical::ok,checks_ok
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  nrow=count([(mod(i-1,nproc)==rank,i=1,n)]);allocate(rows(nrow),s(nrow,n),t(nrow,n),&
    vnl(nrow,n),sipg(nrow,n),vlocal(nrow,n))
  dense_s=(0d0,0d0);dense_t=(0d0,0d0);dense_nl=(0d0,0d0);dense_i=(0d0,0d0);dense_v=(0d0,0d0)
  do i=1,n
    dense_s(i,i)=1d0+0.1d0*i;dense_t(i,i)=0.4d0+0.05d0*i
    dense_v(i,i)=-0.2d0+0.03d0*i
  enddo
  dense_t(1,2)=cmplx(0.04d0,0.01d0,real64);dense_t(2,1)=conjg(dense_t(1,2))
  dense_t(3,4)=cmplx(-0.02d0,0.015d0,real64);dense_t(4,3)=conjg(dense_t(3,4))
  dense_nl(1,3)=cmplx(0.025d0,-0.01d0,real64);dense_nl(3,1)=conjg(dense_nl(1,3))
  dense_i(2,3)=cmplx(-0.08d0,0.03d0,real64);dense_i(3,2)=conjg(dense_i(2,3))
  p=0
  do i=1,n
    if(mod(i-1,nproc)/=rank)cycle
    p=p+1;rows(p)=int(i,int64);s(p,:)=dense_s(i,:);t(p,:)=dense_t(i,:)
    vnl(p,:)=dense_nl(i,:);sipg(p,:)=dense_i(i,:);vlocal(p,:)=dense_v(i,:)
  enddo
  call freeze_dg_hybrid_variational_payload(comm,n,rows,s,t,vnl,sipg,91_int64,92_int64,93_int64,&
    fixed,ok,message)
  call require(ok,trim(message))
  call compose_dg_hybrid_variational_hamiltonian(comm,fixed,vlocal,0d0,1,zero,ok,message)
  call require(ok,trim(message))
  call compose_dg_hybrid_variational_hamiltonian(comm,fixed,vlocal,1d0,2,full,ok,message)
  call require(ok,trim(message))
  checks_ok=.true.
  do p=1,nrow
    i=int(rows(p))
    checks_ok=checks_ok.and.maxval(abs(zero%hamiltonian_rows(p,:)-&
      (dense_t(i,:)+dense_nl(i,:)+dense_v(i,:))))<1d-14
    checks_ok=checks_ok.and.maxval(abs(full%hamiltonian_rows(p,:)-&
      (dense_t(i,:)+dense_nl(i,:)+dense_v(i,:)+dense_i(i,:))))<1d-14
  enddo
  checks_ok=checks_ok.and.maxval(abs(zero%hamiltonian_rows-full%hamiltonian_rows+sipg))<1d-14
  call require(checks_ok,'variational Hamiltonian composition is incorrect')
  call require(fixed%frozen.and.fixed%fingerprint/=0_int64.and.zero%epoch==1.and.full%epoch==2,&
    'variational payload provenance is incomplete')
  if(nrow>0)fixed%kinetic_rows(1,1)=fixed%kinetic_rows(1,1)+1d-3
  call compose_dg_hybrid_variational_hamiltonian(comm,fixed,vlocal,0.5d0,3,full,ok,message)
  call require(.not.ok.and.index(message,'fingerprint')>0,'mutated fixed variational payload was accepted')
  if(rank==0)write(*,'(a,i0,a)')'PASS variational DG payload on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label
    integer::bad,global_bad
    bad=merge(0,1,condition);call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;if(rank==0)write(0,'(a)')trim(label);error stop 1;endif
  end subroutine require
end program test_dg_hybrid_variational_payload_mpi
