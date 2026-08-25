#include "config.h"
program test_dg_hybrid_lcfo_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis,build_dg_hybrid_fragment_basis
  use dg_hybrid_lcfo,only:assemble_dg_hybrid_lcfo_rows
  implicit none
  integer,parameter::npoint=6,nbasis=8
  integer::comm,rank,nproc,ierr,i,j,nowned,position,maximum_callback_columns
  integer(int64),allocatable::wf_ids(:),pw_ids(:),row_ids(:)
  complex(real64),allocatable::wf_values(:,:),pw_values(:,:),hrows(:,:),srows(:,:)
  complex(real64)::vectors(npoint,nbasis),hop(npoint,npoint),sop(npoint,npoint),href(nbasis,nbasis),sref(nbasis,nbasis)
  type(s_dg_hybrid_fragment_basis)::basis
  integer(int64)::peak_elements,fingerprint
  logical::ok
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  maximum_callback_columns=0
  vectors=(0d0,0d0)
  do j=1,nbasis
    vectors(mod(j-1,npoint)+1,j)=cmplx(1d0+0.03d0*j,0.02d0*j,real64)
    vectors(mod(j,npoint)+1,j)=cmplx(-0.18d0,0.01d0*j,real64)
  enddo
  hop=(0d0,0d0);sop=(0d0,0d0)
  do i=1,npoint
    hop(i,i)=1.5d0+0.2d0*i ! local potential
    sop(i,i)=1d0
  enddo
  do i=1,npoint-1
    hop(i,i+1)=cmplx(-0.35d0,0.04d0,real64) ! kinetic/boundary coupling
    hop(i+1,i)=conjg(hop(i,i+1))
  enddo
  hop(2,5)=cmplx(0.11d0,-0.07d0,real64) ! neighboring-core nonlocal projector
  hop(5,2)=conjg(hop(2,5))
  sop(1,4)=cmplx(0.025d0,0.01d0,real64);sop(4,1)=conjg(sop(1,4))
  href=matmul(conjg(transpose(vectors)),matmul(hop,vectors))
  sref=matmul(conjg(transpose(vectors)),matmul(sop,vectors))

  nowned=count([(mod(j-1,nproc)==rank,j=1,nbasis)])
  allocate(wf_ids(count([(mod(j-1,nproc)==rank.and.j<=4,j=1,nbasis)])),&
    pw_ids(count([(mod(j-1,nproc)==rank.and.j>4,j=1,nbasis)])))
  allocate(wf_values(npoint,size(wf_ids)),pw_values(npoint,size(pw_ids)))
  position=0
  do j=1,4
    if(mod(j-1,nproc)/=rank)cycle
    position=position+1;wf_ids(position)=j;wf_values(:,position)=vectors(:,j)
  enddo
  position=0
  do j=5,nbasis
    if(mod(j-1,nproc)/=rank)cycle
    position=position+1;pw_ids(position)=j;pw_values(:,position)=vectors(:,j)
  enddo
  call build_dg_hybrid_fragment_basis(comm,1,wf_ids,wf_values,pw_ids,pw_values,2,3,basis,ok,message)
  call require(ok,trim(message))
  do i=1,npoint
    basis%buffer_point_ids(i)=int(mod(i+rank-1,npoint)+1,int64)
    do j=1,nowned
      basis%buffer_values(i,j)=vectors(int(basis%buffer_point_ids(i)),int(basis%global_ids(j)))
    enddo
  enddo
  allocate(row_ids(nowned));row_ids=basis%global_ids
  call assemble_dg_hybrid_lcfo_rows(comm,basis,row_ids,apply_h,apply_s,hrows,srows,peak_elements,fingerprint,ok,message)
  call require(ok,trim(message))
  call require(size(hrows,1)==nowned.and.size(hrows,2)==nbasis,'LCFO Hamiltonian is not row distributed')
  call require(all(shape(srows)==shape(hrows)),'LCFO metric row shape mismatch')
  do i=1,nowned
    call require(maxval(abs(hrows(i,:)-href(int(row_ids(i)),:)))<2d-13,'LCFO Hamiltonian contribution mismatch')
    call require(maxval(abs(srows(i,:)-sref(int(row_ids(i)),:)))<2d-13,'LCFO metric contribution mismatch')
  enddo
  call require(peak_elements==int(2*nowned*nbasis,int64),'LCFO persistent matrix receipt is not row local')
  call require(maximum_callback_columns==1,'LCFO assembly did not stream operator columns')
  call require(fingerprint/=0_int64,'LCFO operator fingerprint is empty')
  if(rank==0)write(*,'(a,i0,a,i0)')'HYBRID_LCFO ranks=',nproc,' fingerprint=',fingerprint
  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid LCFO assembly on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine apply_h(input,output,callback_ok)
    complex(real64),intent(in)::input(:,:);complex(real64),intent(out)::output(:,:)
    logical,intent(out)::callback_ok
    complex(real64)::global_input(npoint,size(input,2)),global_output(npoint,size(input,2))
    integer::p
    maximum_callback_columns=max(maximum_callback_columns,size(input,2))
    global_input=(0d0,0d0)
    do p=1,size(input,1);global_input(int(basis%buffer_point_ids(p)),:)=input(p,:);enddo
    global_output=matmul(hop,global_input)
    do p=1,size(output,1);output(p,:)=global_output(int(basis%buffer_point_ids(p)),:);enddo
    callback_ok=.true.
  end subroutine apply_h
  subroutine apply_s(input,output,callback_ok)
    complex(real64),intent(in)::input(:,:);complex(real64),intent(out)::output(:,:)
    logical,intent(out)::callback_ok
    complex(real64)::global_input(npoint,size(input,2)),global_output(npoint,size(input,2))
    integer::p
    maximum_callback_columns=max(maximum_callback_columns,size(input,2))
    global_input=(0d0,0d0)
    do p=1,size(input,1);global_input(int(basis%buffer_point_ids(p)),:)=input(p,:);enddo
    global_output=matmul(sop,global_input)
    do p=1,size(output,1);output(p,:)=global_output(int(basis%buffer_point_ids(p)),:);enddo
    callback_ok=.true.
  end subroutine apply_s
  subroutine require(condition,text)
    logical,intent(in)::condition;character(*),intent(in)::text
    logical::global_condition
    call MPI_Allreduce(condition,global_condition,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(.not.global_condition)then
      if(rank==0)write(0,'(a)')trim(text)
      call MPI_Abort(comm,1,ierr)
    endif
  end subroutine require
end program test_dg_hybrid_lcfo_mpi
