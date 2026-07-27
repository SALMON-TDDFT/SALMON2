#include "config.h"
program test_dg_overlapping_wannier_solver_mpi
  use mpi
  use dg_overlapping_wannier_solver,only:solve_dg_overlapping_wannier_coefficients
  use dg_overlapping_wannier_density,only:reconstruct_dg_overlapping_wannier_density
  implicit none
  integer::comm,rank,nproc,ierr,nlocal,i,j,p,index,lapack_info
  integer(8),allocatable::row_ids(:),point_ids(:),bad_ids(:)
  integer,allocatable::tail_generation(:,:)
  complex(8),allocatable::hrows(:,:),srows(:,:),coeff(:,:),values(:,:),smetric(:,:),hmetric(:,:),&
    reference_h(:,:),reference_s(:,:),work(:)
  complex(8),allocatable::degenerate_projector(:,:)
  complex(8)::gauge(4,4)
  real(8),allocatable::eval(:),weights(:),rho(:),reference_eval(:),rwork(:)
  real(8)::residual,orthogonality,condition,charge,trace_charge,pi
  integer(8)::signature
  logical::ok
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  pi=acos(-1d0)

  nlocal=count([(mod(i-1,nproc)==rank,i=1,4)])
  allocate(row_ids(nlocal),hrows(nlocal,4),srows(nlocal,4))
  smetric=reshape([cmplx(2d0,0d0,8),cmplx(0.2d0,0.1d0,8),cmplx(0d0,0d0,8),cmplx(0d0,0d0,8),&
    cmplx(0.2d0,-0.1d0,8),cmplx(1.5d0,0d0,8),cmplx(0.1d0,0d0,8),cmplx(0d0,0d0,8),&
    cmplx(0d0,0d0,8),cmplx(0.1d0,0d0,8),cmplx(1.2d0,0d0,8),cmplx(0.05d0,0.02d0,8),&
    cmplx(0d0,0d0,8),cmplx(0d0,0d0,8),cmplx(0.05d0,-0.02d0,8),cmplx(0.9d0,0d0,8)],[4,4])
  hmetric=reshape([cmplx(0.4d0,0d0,8),cmplx(0.03d0,0.02d0,8),cmplx(0d0,0d0,8),cmplx(0d0,0d0,8),&
    cmplx(0.03d0,-0.02d0,8),cmplx(0.9d0,0d0,8),cmplx(0.04d0,0d0,8),cmplx(0d0,0d0,8),&
    cmplx(0d0,0d0,8),cmplx(0.04d0,0d0,8),cmplx(1.6d0,0d0,8),cmplx(0.02d0,0.01d0,8),&
    cmplx(0d0,0d0,8),cmplx(0d0,0d0,8),cmplx(0.02d0,-0.01d0,8),cmplx(2.3d0,0d0,8)],[4,4])
  index=0
  do i=1,4
    if(mod(i-1,nproc)/=rank)cycle
    index=index+1;row_ids(index)=i;hrows(index,:)=hmetric(i,:);srows(index,:)=smetric(i,:)
  enddo
  allocate(coeff(4,3),eval(3))
  call solve_dg_overlapping_wannier_coefficients(comm,row_ids,hrows,srows,3,300,1d-9,1d-10,&
    coeff,eval,residual,orthogonality,condition,ok,message)
  call require(ok,trim(message))
  call require(residual<1d-9,'generalized residual')
  call require(orthogonality<1d-9,'S orthonormality')
  call require(all(eval(2:3)>=eval(1:2)),'ordered occupied and empty targets')
  allocate(reference_h,source=hmetric);allocate(reference_s,source=smetric)
  allocate(reference_eval(4),work(8),rwork(10))
  call dense_generalized_reference(reference_h,reference_s,reference_eval,work,rwork,lapack_info)
  call require(lapack_info==0.and.maxval(abs(eval-reference_eval(1:3)))<1d-8,&
    'distributed solver agrees with dense LAPACK fixture reference')

  gauge=(0d0,0d0);gauge(1,1)=sqrt(0.5d0);gauge(1,2)=sqrt(0.5d0)
  gauge(2,1)=-sqrt(0.5d0);gauge(2,2)=sqrt(0.5d0);gauge(3,3)=1d0;gauge(4,4)=1d0
  hmetric=matmul(conjg(transpose(gauge)),matmul(hmetric,gauge))
  smetric=matmul(conjg(transpose(gauge)),matmul(smetric,gauge))
  hrows=hmetric(row_ids,:);srows=smetric(row_ids,:)
  call solve_dg_overlapping_wannier_coefficients(comm,row_ids,hrows,srows,3,300,1d-9,1d-10,&
    coeff,eval,residual,orthogonality,condition,ok,message)
  call require(ok.and.maxval(abs(eval-reference_eval(1:3)))<1d-8,&
    'candidate-basis rotation invariant generalized spectrum')
  hmetric=matmul(gauge,matmul(hmetric,conjg(transpose(gauge))))
  smetric=matmul(gauge,matmul(smetric,conjg(transpose(gauge))))
  hrows=hmetric(row_ids,:);srows=smetric(row_ids,:)

  smetric=(0d0,0d0);hmetric=(0d0,0d0)
  do i=1,4;smetric(i,i)=1d0;enddo
  hmetric(1,1)=1d0;hmetric(2,2)=1d0;hmetric(3,3)=3d0;hmetric(4,4)=4d0
  hrows=hmetric(row_ids,:);srows=smetric(row_ids,:)
  call solve_dg_overlapping_wannier_coefficients(comm,row_ids,hrows,srows,2,300,1d-10,1d-10,&
    coeff(:,1:2),eval(1:2),residual,orthogonality,condition,ok,message)
  call require(ok.and.maxval(abs(eval(1:2)-1d0))<1d-10,'degenerate target convergence')
  degenerate_projector=matmul(coeff(:,1:2),conjg(transpose(coeff(:,1:2))))
  call require(maxval(abs(degenerate_projector-reshape([cmplx(1d0,0d0,8),cmplx(0d0,0d0,8),&
    cmplx(0d0,0d0,8),cmplx(0d0,0d0,8),cmplx(0d0,0d0,8),cmplx(1d0,0d0,8),&
    cmplx(0d0,0d0,8),cmplx(0d0,0d0,8),cmplx(0d0,0d0,8),cmplx(0d0,0d0,8),&
    cmplx(0d0,0d0,8),cmplx(0d0,0d0,8),cmplx(0d0,0d0,8),cmplx(0d0,0d0,8),&
    cmplx(0d0,0d0,8),cmplx(0d0,0d0,8)],[4,4])))<1d-9,'degenerate S-projector')
  signature=nint(sum([(real(i,8)*real(degenerate_projector(i,i),8),i=1,4)])*1d10,8)

  smetric=(0d0,0d0)
  do i=1,4;smetric(i,i)=1d0;enddo
  smetric(1,2)=1d0-1d-12;smetric(2,1)=1d0-1d-12
  srows=smetric(row_ids,:);hrows=hmetric(row_ids,:)
  call solve_dg_overlapping_wannier_coefficients(comm,row_ids,hrows,srows,3,50,1d-10,1d-8,&
    coeff,eval,residual,orthogonality,condition,ok,message)
  call require(.not.ok,'ill-conditioned metric rejection')
  srows=smetric(row_ids,:)

  nlocal=count([(mod(p-1,nproc)==rank,p=1,6)])
  allocate(point_ids(nlocal),bad_ids(nlocal),weights(nlocal),values(4,nlocal),tail_generation(4,nlocal),rho(nlocal))
  index=0
  do p=1,6
    if(mod(p-1,nproc)/=rank)cycle
    index=index+1;point_ids(index)=p;bad_ids(index)=p;weights(index)=0.5d0
    values(:,index)=[cmplx(1d0,0d0,8),cmplx(cos(2*pi*p/6d0),sin(2*pi*p/6d0),8),&
      cmplx(cos(4*pi*p/6d0),sin(4*pi*p/6d0),8),cmplx((-1d0)**p,0d0,8)]/sqrt(3d0)
  enddo
  tail_generation=17
  smetric=(0d0,0d0)
  do p=1,nlocal
    do j=1,4;do i=1,4
      smetric(i,j)=smetric(i,j)+weights(p)*conjg(values(i,p))*values(j,p)
    enddo;enddo
  enddo
  call MPI_Allreduce(MPI_IN_PLACE,smetric,16,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  hmetric=smetric
  hrows=hmetric(row_ids,:);srows=smetric(row_ids,:)
  call solve_dg_overlapping_wannier_coefficients(comm,row_ids,hrows,srows,2,300,1d-9,1d-10,&
    coeff(:,1:2),eval(1:2),residual,orthogonality,condition,ok,message)
  call require(ok,trim(message))
  call reconstruct_dg_overlapping_wannier_density(comm,point_ids,weights,values,tail_generation,17,&
    coeff(:,1:2),[2d0,1d0],6_8,smetric,1d-10,rho,charge,trace_charge,ok,message)
  call require(ok,trim(message))
  call require(maxval(abs(rho-sum(abs(matmul(transpose(values),coeff(:,1:2)))**2*&
    spread([2d0,1d0],1,nlocal),dim=2)))<1d-10,'all covering tails in density')
  call require(abs(charge-3d0)<1d-9.and.abs(trace_charge-3d0)<1d-9,'electron count identities')
  if(nproc>1)then
    if(rank==0)coeff(1,1)=coeff(1,1)+1d-3
    call reconstruct_dg_overlapping_wannier_density(comm,point_ids,weights,values,tail_generation,17,&
      coeff(:,1:2),[2d0,1d0],6_8,smetric,1d-10,rho,charge,trace_charge,ok,message)
    call require(.not.ok,'rank-corrupted replicated density payload rejection')
    if(rank==0)coeff(1,1)=coeff(1,1)-1d-3
  endif

  if(rank==0.and.nlocal>0)bad_ids(1)=2_8
  call reconstruct_dg_overlapping_wannier_density(comm,bad_ids,weights,values,tail_generation,17,&
    coeff(:,1:2),[2d0,1d0],6_8,smetric,1d-10,rho,charge,trace_charge,ok,message)
  call require(.not.ok,'duplicate or missing core write rejection')
  if(nlocal>0)tail_generation(4,1)=0
  call reconstruct_dg_overlapping_wannier_density(comm,point_ids,weights,values,tail_generation,17,&
    coeff(:,1:2),[2d0,1d0],6_8,smetric,1d-10,rho,charge,trace_charge,ok,message)
  call require(.not.ok,'missing tail pair rejection')
  tail_generation=17
  if(nproc>1)then
    call reconstruct_dg_overlapping_wannier_density(comm,point_ids,weights,values,tail_generation,17,&
      coeff(:,1:2),[2d0,1d0],merge(7_8,6_8,rank==0),smetric,1d-10,rho,charge,trace_charge,ok,message)
    call require(.not.ok,'rank-mismatched density count rejection')
  endif

  if(rank==0)then
    write(*,'(a,i0,a,i0)')'SOLVER ranks=',nproc,' signature=',signature
    write(*,'(a,i0,a)')'PASS overlapping-Wannier solver and density on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine dense_generalized_reference(h,s,e,w,rw,info)
    complex(8),intent(inout)::h(4,4),s(4,4),w(:)
    real(8),intent(out)::e(4),rw(:)
    integer,intent(out)::info
    complex(8)::l(4,4),linv(4,4),rhs(4),x(4)
    real(8)::pivot
    integer::ii,jj,kk
    l=(0d0,0d0);info=0
    do jj=1,4
      pivot=real(s(jj,jj)-sum(l(jj,1:jj-1)*conjg(l(jj,1:jj-1))),8)
      if(pivot<=0d0)then;info=1;return;endif
      l(jj,jj)=sqrt(pivot)
      do ii=jj+1,4
        l(ii,jj)=(s(ii,jj)-sum(l(ii,1:jj-1)*conjg(l(jj,1:jj-1))))/l(jj,jj)
      enddo
    enddo
    linv=(0d0,0d0)
    do jj=1,4
      rhs=(0d0,0d0);rhs(jj)=1d0;x=(0d0,0d0)
      do ii=1,4
        x(ii)=(rhs(ii)-sum(l(ii,1:ii-1)*x(1:ii-1)))/l(ii,ii)
      enddo
      linv(:,jj)=x
    enddo
    h=matmul(linv,matmul(h,conjg(transpose(linv))))
    h=0.5d0*(h+conjg(transpose(h)))
    call zheev('V','U',4,h,4,e,w,size(w),rw,info)
  end subroutine

  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::lf,gf
    lf=merge(0,1,condition);call MPI_Allreduce(lf,gf,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(gf/=0)error stop label
  end subroutine
end program
