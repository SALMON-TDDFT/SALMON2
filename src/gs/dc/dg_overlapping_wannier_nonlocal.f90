#include "config.h"
module dg_overlapping_wannier_nonlocal
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::assemble_dg_overlapping_wannier_nonlocal
contains
  subroutine assemble_dg_overlapping_wannier_nonlocal(comm,nwann,projector_ids,strength,overlap,&
      complete_tail_overlap,expected_projector_count,matrix,ownership_count,ok,message)
    integer,intent(in)::comm,nwann
    integer(int64),intent(in)::projector_ids(:),expected_projector_count
    real(real64),intent(in)::strength(:)
    complex(real64),intent(in)::overlap(:,:)
    logical,intent(in)::complete_tail_overlap(:,:)
    complex(real64),allocatable,intent(out)::matrix(:,:)
    integer,intent(out)::ownership_count
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nproc,ierr,local_bad,global_bad,total_count,i,j,p,matrix_count
    integer,allocatable::counts(:),displs(:)
    integer(int64),allocatable::all_ids(:)
    integer(int64)::matrix_count64
    complex(real64),allocatable::local_matrix(:,:)
    logical::shape_ok,finite_overlap
    ok=.false.;message='';ownership_count=0;local_bad=0
    shape_ok=size(strength)==size(projector_ids).and.size(overlap,1)==nwann.and.&
      size(overlap,2)==size(projector_ids).and.size(complete_tail_overlap,1)==nwann.and.&
      size(complete_tail_overlap,2)==size(projector_ids)
    if(nwann<=0.or.expected_projector_count<=0_int64.or..not.shape_ok)local_bad=1
    if(any(projector_ids<=0_int64).or.any(.not.ieee_is_finite(strength)))local_bad=1
    finite_overlap=.true.
    if(shape_ok)then
      do p=1,size(projector_ids);do i=1,nwann
        finite_overlap=finite_overlap.and.ieee_is_finite(real(overlap(i,p))).and.&
          ieee_is_finite(aimag(overlap(i,p)))
      enddo;enddo
      if(.not.all(complete_tail_overlap))local_bad=1
    endif
    if(.not.finite_overlap)local_bad=1
    if(nwann>0.and.int(nwann,int64)<=huge(1_int64)/int(nwann,int64))then
      matrix_count64=int(nwann,int64)*int(nwann,int64)
      if(matrix_count64>int(huge(matrix_count),int64))local_bad=1
    else
      matrix_count64=0_int64;local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid or incomplete tail-projector overlaps';return;endif
    matrix_count=int(matrix_count64)
    call MPI_Comm_size(comm,nproc,ierr);allocate(counts(nproc),displs(nproc))
    call MPI_Allgather(size(projector_ids),1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,ierr)
    total_count=0;displs(1)=0
    do i=1,nproc
      if(counts(i)<0.or.total_count>huge(total_count)-counts(i))then;local_bad=1;exit;endif
      if(i>1)displs(i)=total_count
      total_count=total_count+counts(i)
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.int(total_count,int64)/=expected_projector_count)then
      message='missing or extra atom/projector owner';return
    endif
    allocate(all_ids(total_count))
    call MPI_Allgatherv(projector_ids,size(projector_ids),MPI_INTEGER8,all_ids,counts,displs,&
      MPI_INTEGER8,comm,ierr)
    call sort_ids(all_ids)
    do i=1,total_count
      if(all_ids(i)/=int(i,int64))then;message='duplicate or missing atom/projector owner';return;endif
    enddo
    allocate(local_matrix(nwann,nwann),matrix(nwann,nwann));local_matrix=(0d0,0d0)
    do p=1,size(projector_ids);do j=1,nwann;do i=1,nwann
      local_matrix(i,j)=local_matrix(i,j)+strength(p)*conjg(overlap(i,p))*overlap(j,p)
    enddo;enddo;enddo
    call MPI_Allreduce(local_matrix,matrix,matrix_count,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    matrix=0.5d0*(matrix+conjg(transpose(matrix)))
    ownership_count=total_count;ok=.true.
#else
    ok=.false.;message='overlapping-Wannier nonlocal assembly requires MPI';ownership_count=0
#endif
  end subroutine
  subroutine sort_ids(ids)
    integer(int64),intent(inout)::ids(:)
    integer::i,j
    integer(int64)::key
    do i=2,size(ids)
      key=ids(i);j=i-1
      do while(j>=1)
        if(ids(j)<=key)exit
        ids(j+1)=ids(j);j=j-1
      enddo
      ids(j+1)=key
    enddo
  end subroutine
end module dg_overlapping_wannier_nonlocal
