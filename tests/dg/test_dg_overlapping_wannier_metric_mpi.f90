#include "config.h"
program test_dg_overlapping_wannier_metric_mpi
  use mpi
  use dg_overlapping_wannier_metric,only:assemble_dg_overlapping_wannier_metric
  implicit none
  integer::comm,rank,nproc,ierr,i,nlocal,owned,rejected,reference_owned
  integer(8),allocatable::ids(:)
  real(8),allocatable::weights(:)
  complex(8),allocatable::values(:,:),metric(:,:),vectors(:,:)
  complex(8),allocatable::base_values(:,:),reference_metric(:,:)
  real(8),allocatable::spectrum(:),reference_spectrum(:)
  complex(8)::rotation(3,3)
  logical,allocatable::pairs(:,:)
  real(8)::minimum,condition
  logical::ok
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  nlocal=count([(mod(i-1,nproc)==rank,i=1,4)])
  allocate(ids(nlocal),weights(nlocal),values(3,nlocal),pairs(3,nlocal))
  nlocal=0
  do i=1,4
    if(mod(i-1,nproc)/=rank)cycle
    nlocal=nlocal+1;ids(nlocal)=i;weights(nlocal)=0.5d0+0.1d0*i
    values(:,nlocal)=[cmplx(1d0+0.1d0*i,0.05d0*i,8),&
      cmplx((-1d0)**i*0.3d0,0.02d0*i,8),cmplx(0.2d0*i,-0.04d0*i,8)]
  enddo
  pairs=.true.
  call assemble_dg_overlapping_wannier_metric(comm,3,ids,weights,values,pairs,4_8,1d-10,&
    metric,vectors,spectrum,minimum,condition,rejected,owned,ok,message)
  call require(ok,trim(message));call require(owned==4,'unique ownership')
  call require(maxval(abs(metric-conjg(transpose(metric))))<1d-13,'Hermitian metric')
  call require(abs(metric(1,2))>1d-8,'off-fragment block')
  call require(minimum>0d0.and.condition>=1d0.and.rejected==0,'positive rank revelation')
  call require_same_matrix(metric)
  reference_spectrum=spectrum;base_values=values;reference_metric=metric;reference_owned=owned

  values=base_values;values(2,:)=-values(2,:)
  call check_invariant('sign invariance')
  values=base_values([2,1,3],:)
  call check_invariant('permutation invariance')
  rotation=(0d0,0d0);rotation(1,1)=sqrt(0.5d0);rotation(1,2)=sqrt(0.5d0)
  rotation(2,1)=-sqrt(0.5d0);rotation(2,2)=sqrt(0.5d0);rotation(3,3)=1d0
  values=matmul(rotation,base_values)
  call check_invariant('unitary candidate-window invariance')

  values=base_values;values(3,:)=values(1,:)
  call assemble_dg_overlapping_wannier_metric(comm,3,ids,weights,values,pairs,4_8,1d-10,&
    metric,vectors,spectrum,minimum,condition,rejected,owned,ok,message)
  call require(ok.and.rejected==1,'positive-metric null-space rank revelation')

  values=base_values
  if(rank==0.and.size(ids)>0)pairs(1,1)=.false.
  call assemble_dg_overlapping_wannier_metric(comm,3,ids,weights,values,pairs,4_8,1d-10,&
    metric,vectors,spectrum,minimum,condition,rejected,owned,ok,message)
  call require(.not.ok,'missing owner-pair rejection')
  pairs=.true.

  if(nproc>1)then
    call assemble_dg_overlapping_wannier_metric(comm,merge(4,3,rank==0),ids,weights,values,pairs,&
      4_8,1d-10,metric,vectors,spectrum,minimum,condition,rejected,owned,ok,message)
    call require(.not.ok,'rank-inconsistent metric contract rejection')
  endif

  if(rank==0.and.size(ids)>0)then
    ids=[ids,ids(1)];weights=[weights,weights(1)]
    values=reshape([values,values(:,1)],[3,size(ids)])
    pairs=reshape([pairs,pairs(:,1)],[3,size(ids)])
  endif
  call assemble_dg_overlapping_wannier_metric(comm,3,ids,weights,values,pairs,4_8,1d-10,&
    metric,vectors,spectrum,minimum,condition,rejected,owned,ok,message)
  call require(.not.ok,'duplicate core quadrature rejection')
  if(rank==0)then
    write(*,'(a,i0,a,i0)')'METRIC ranks=',nproc,' ownership=',reference_owned
    write(*,'(*(es24.16,1x))')[(real(reference_metric(i,1)),aimag(reference_metric(i,1)),&
      i=1,size(reference_metric,1))]
    write(*,'(a,i0,a)')'PASS overlapping-Wannier metric on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine check_invariant(label)
    character(*),intent(in)::label
    call assemble_dg_overlapping_wannier_metric(comm,3,ids,weights,values,pairs,4_8,1d-10,&
      metric,vectors,spectrum,minimum,condition,rejected,owned,ok,message)
    call require(ok,label)
    call require(size(spectrum)==size(reference_spectrum),label)
    call require(maxval(abs(spectrum-reference_spectrum))<1d-12,label)
  end subroutine
  subroutine require(c,label)
    logical,intent(in)::c;character(*),intent(in)::label
    integer::l,g
    l=merge(0,1,c);call MPI_Allreduce(l,g,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(g/=0)error stop label
  end subroutine
  subroutine require_same_matrix(a)
    complex(8),intent(in)::a(:,:)
    complex(8)::reference(size(a,1),size(a,2))
    reference=a
    call MPI_Bcast(reference,size(a),MPI_DOUBLE_COMPLEX,0,comm,ierr)
    call require(maxval(abs(reference-a))<1d-14,'deterministic matrix agreement')
  end subroutine
end program
