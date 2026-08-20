#include "config.h"
program test_dg_hybrid_generalized_eigensystem_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_get_halting_mode,ieee_set_halting_mode,ieee_set_flag,&
    ieee_invalid,ieee_divide_by_zero,ieee_overflow
  use dg_hybrid_generalized_eigensystem,only:solve_dg_hybrid_generalized_scalapack
  implicit none
  integer,parameter::n=4,nstate=2
  integer::comm,rank,nproc,ierr,nowned,row,i,j,position
  integer(int64),allocatable::row_ids(:)
  complex(real64),allocatable::hrows(:,:),srows(:,:),coefficients(:,:)
  complex(real64)::s(n,n),h(n,n),u(n,n),phase
  real(real64)::expected(n),eigenvalues(nstate),residual,orthogonality,projector_defect
  integer(int64)::workspace,fingerprint,reference_fingerprint
  logical::ok
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  s=(0d0,0d0);h=(0d0,0d0);u=(0d0,0d0)
  do i=1,n;s(i,i)=1.4d0+0.1d0*i;u(i,i)=1d0;enddo
  s(1,2)=cmplx(0.08d0,0.03d0,real64);s(2,1)=conjg(s(1,2))
  s(3,4)=cmplx(-0.05d0,0.02d0,real64);s(4,3)=conjg(s(3,4))
  h(1,1)=-0.7d0;h(2,2)=-0.2d0;h(3,3)=0.4d0;h(4,4)=0.9d0
  h(1,3)=cmplx(0.06d0,-0.04d0,real64);h(3,1)=conjg(h(1,3))
  h(2,4)=cmplx(-0.03d0,0.05d0,real64);h(4,2)=conjg(h(2,4))
  expected=0d0;call dense_oracle(h,s,expected)
  nowned=count([(mod(row-1,nproc)==rank,row=1,n)])
  allocate(row_ids(nowned),hrows(nowned,n),srows(nowned,n));position=0
  do row=n,1,-1
    if(mod(row-1,nproc)/=rank)cycle
    position=position+1;row_ids(position)=row;hrows(position,:)=h(row,:);srows(position,:)=s(row,:)
  enddo
  call solve_dg_hybrid_generalized_scalapack(comm,n,nstate,row_ids,hrows,srows,1d-11,coefficients,eigenvalues,&
    residual,orthogonality,projector_defect,workspace,fingerprint,ok,message)
  call require(ok,trim(message));reference_fingerprint=fingerprint
  call require(maxval(abs(eigenvalues-expected(1:nstate)))<2d-11,'ScaLAPACK eigenvalues differ from ZHEGV')
  call require(residual<2d-10.and.orthogonality<2d-10.and.projector_defect<2d-10,&
    'generalized eigensystem receipts are invalid')
  if(nproc>1)then
    call solve_dg_hybrid_generalized_scalapack(comm,n,nstate,row_ids,hrows,srows,merge(2d-11,1d-11,rank==0),&
      coefficients,eigenvalues,residual,orthogonality,projector_defect,workspace,fingerprint,ok,message)
    call require(.not.ok,'rank-disagreeing generalized tolerance was accepted')
  endif
  if(nowned>0)hrows(1,mod(int(row_ids(1)),n)+1)=hrows(1,mod(int(row_ids(1)),n)+1)+(0.2d0,0.1d0)
  call solve_dg_hybrid_generalized_scalapack(comm,n,nstate,row_ids,hrows,srows,1d-11,coefficients,eigenvalues,&
    residual,orthogonality,projector_defect,workspace,fingerprint,ok,message)
  call require(.not.ok,'non-Hermitian hybrid Hamiltonian was accepted')
  do i=1,nowned;hrows(i,:)=h(int(row_ids(i)),:);enddo
  s(4,:)=(0d0,0d0);s(:,4)=(0d0,0d0)
  do i=1,nowned;srows(i,:)=s(int(row_ids(i)),:);enddo
  call solve_dg_hybrid_generalized_scalapack(comm,n,nstate,row_ids,hrows,srows,1d-11,coefficients,eigenvalues,&
    residual,orthogonality,projector_defect,workspace,fingerprint,ok,message)
  call require(.not.ok,'rank-deficient reference metric was accepted')
  if(rank==0)then
    write(*,'(a,i0,a,i0)')'HYBRID_GENERALIZED_EIGENSYSTEM ranks=',nproc,' fingerprint=',reference_fingerprint
    write(*,'(a,i0,a)')'PASS hybrid generalized eigensystem on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine dense_oracle(a,b,w)
    complex(real64),intent(in)::a(:,:),b(:,:)
    real(real64),intent(out)::w(:)
    complex(real64)::acopy(n,n),bcopy(n,n),query(1)
    complex(real64),allocatable::work(:)
    real(real64)::rwork(max(1,3*n-2))
    integer::info,lwork
    logical::halt_invalid,halt_zero,halt_overflow
    external::zhegv
    acopy=a;bcopy=b;lwork=-1
    call ieee_get_halting_mode(ieee_invalid,halt_invalid);call ieee_get_halting_mode(ieee_divide_by_zero,halt_zero)
    call ieee_get_halting_mode(ieee_overflow,halt_overflow)
    call ieee_set_halting_mode(ieee_invalid,.false.);call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
    call ieee_set_halting_mode(ieee_overflow,.false.)
    call zhegv(1,'V','U',n,acopy,n,bcopy,n,w,query,lwork,rwork,info);call require(info==0,'ZHEGV query failed')
    lwork=max(1,int(real(query(1))));allocate(work(lwork));acopy=a;bcopy=b
    call zhegv(1,'V','U',n,acopy,n,bcopy,n,w,work,lwork,rwork,info);call require(info==0,'ZHEGV oracle failed')
    call ieee_set_flag(ieee_invalid,.false.);call ieee_set_flag(ieee_divide_by_zero,.false.)
    call ieee_set_flag(ieee_overflow,.false.);call ieee_set_halting_mode(ieee_invalid,halt_invalid)
    call ieee_set_halting_mode(ieee_divide_by_zero,halt_zero);call ieee_set_halting_mode(ieee_overflow,halt_overflow)
  end subroutine dense_oracle
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label
    integer::local_bad,global_bad
    local_bad=merge(0,1,condition);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)error stop label
  end subroutine require
end program test_dg_hybrid_generalized_eigensystem_mpi
