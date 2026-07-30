#include "config.h"
module dg_overlapping_wannier_nonlocal
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::assemble_dg_overlapping_wannier_nonlocal,assemble_dg_overlapping_wannier_nonlocal_rows
contains
  subroutine assemble_dg_overlapping_wannier_nonlocal_rows(comm,nwann,row_ids,projector_ids,strength,&
      overlap,complete_tail_overlap,expected_projector_count,matrix_rows,ownership_count,ok,message)
    integer,intent(in)::comm,nwann
    integer(int64),intent(in)::row_ids(:),projector_ids(:),expected_projector_count
    real(real64),intent(in)::strength(:)
    complex(real64),intent(in)::overlap(:,:)
    logical,intent(in)::complete_tail_overlap(:,:)
    complex(real64),allocatable,intent(out)::matrix_rows(:,:)
    integer,intent(out)::ownership_count
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer,parameter::row_batch_size=32
    integer::rank,nproc,ierr,local_bad,global_bad,total_projectors,total_rows,r,i,j,p,nrows,&
      batch_first,batch_count,nwann_min,nwann_max
    integer,allocatable::projector_counts(:),projector_displs(:),row_counts(:),row_displs(:)
    integer(int64),allocatable::all_projector_ids(:),all_row_ids(:),validation_ids(:)
    integer(int64)::expected_min,expected_max
    complex(real64),allocatable::partial(:,:),reduced(:,:)
    logical::shape_ok

    ok=.false.;message='';ownership_count=0;local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(nwann,nwann_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(nwann,nwann_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(expected_projector_count,expected_min,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(expected_projector_count,expected_max,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    if(nwann_min/=nwann_max.or.expected_min/=expected_max)local_bad=1
    shape_ok=size(strength)==size(projector_ids).and.size(overlap,1)==nwann.and.&
      size(overlap,2)==size(projector_ids).and.size(complete_tail_overlap,1)==nwann.and.&
      size(complete_tail_overlap,2)==size(projector_ids)
    if(nwann<=0.or.expected_projector_count<0_int64.or..not.shape_ok)local_bad=1
    if(any(row_ids<1_int64).or.any(row_ids>int(nwann,int64)).or.any(projector_ids<=0_int64).or.&
        any(.not.ieee_is_finite(strength)))local_bad=1
    if(shape_ok)then
      if(.not.all(complete_tail_overlap).or..not.all(ieee_is_finite(real(overlap))).or.&
          .not.all(ieee_is_finite(aimag(overlap))))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid row-owned nonlocal payload';return;endif
    allocate(projector_counts(nproc),projector_displs(nproc),row_counts(nproc),row_displs(nproc))
    call MPI_Allgather(size(projector_ids),1,MPI_INTEGER,projector_counts,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allgather(size(row_ids),1,MPI_INTEGER,row_counts,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0)then;message='nonlocal ownership metadata collective failed';return;endif
    total_projectors=0;total_rows=0
    do r=1,nproc
      projector_displs(r)=total_projectors;row_displs(r)=total_rows
      if(projector_counts(r)<0.or.row_counts(r)<0.or.&
          total_projectors>huge(total_projectors)-projector_counts(r).or.&
          total_rows>huge(total_rows)-row_counts(r))local_bad=1
      if(local_bad==0)then
        total_projectors=total_projectors+projector_counts(r);total_rows=total_rows+row_counts(r)
      endif
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.int(total_projectors,int64)/=expected_projector_count.or.total_rows/=nwann)then
      message='missing or extra row/projector owner in nonlocal assembly';return
    endif
    allocate(all_projector_ids(total_projectors),all_row_ids(total_rows))
    call MPI_Allgatherv(projector_ids,size(projector_ids),MPI_INTEGER8,all_projector_ids,&
      projector_counts,projector_displs,MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allgatherv(row_ids,size(row_ids),MPI_INTEGER8,all_row_ids,row_counts,row_displs,&
      MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0)then;message='nonlocal ownership payload collective failed';return;endif
    allocate(validation_ids(max(total_projectors,total_rows)))
    validation_ids(1:total_projectors)=all_projector_ids
    call sort_ids(validation_ids(1:total_projectors))
    do i=1,total_projectors
      if(validation_ids(i)/=int(i,int64))local_bad=1
    enddo
    validation_ids(1:total_rows)=all_row_ids
    call sort_ids(validation_ids(1:total_rows))
    do i=1,total_rows
      if(validation_ids(i)/=int(i,int64))local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='duplicate or missing row/projector owner in nonlocal assembly';return;endif
    allocate(matrix_rows(size(row_ids),nwann));matrix_rows=(0d0,0d0)
    do r=0,nproc-1
      nrows=row_counts(r+1)
      do batch_first=1,nrows,row_batch_size
        batch_count=min(row_batch_size,nrows-batch_first+1)
        allocate(partial(batch_count,nwann),reduced(batch_count,nwann));partial=(0d0,0d0)
        do p=1,size(projector_ids);do j=1,nwann;do i=1,batch_count
          partial(i,j)=partial(i,j)+strength(p)*conjg(overlap(&
            int(all_row_ids(row_displs(r+1)+batch_first+i-1)),p))*overlap(j,p)
        enddo;enddo;enddo
        call MPI_Reduce(partial,reduced,batch_count*nwann,MPI_DOUBLE_COMPLEX,MPI_SUM,r,comm,ierr)
        if(ierr/=MPI_SUCCESS)local_bad=1
        if(rank==r)matrix_rows(batch_first:batch_first+batch_count-1,:)=reduced
        deallocate(partial,reduced)
      enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0)then;message='nonlocal row reduction failed';return;endif
    ownership_count=total_projectors;ok=.true.
#else
    ok=.false.;message='row-owned overlapping-Wannier nonlocal assembly requires MPI'
    ownership_count=0
#endif
  end subroutine

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
    integer::nproc,ierr,local_bad,global_bad,total_count,i,j,p,matrix_count,nwann_min,nwann_max
    integer,allocatable::counts(:),displs(:)
    integer(int64),allocatable::all_ids(:)
    integer(int64)::matrix_count64,expected_min,expected_max
    complex(real64),allocatable::local_matrix(:,:)
    logical::shape_ok,finite_overlap
    real(real64)::hermiticity_defect,scale
    ok=.false.;message='';ownership_count=0;local_bad=0
    call MPI_Allreduce(nwann,nwann_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(nwann,nwann_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(expected_projector_count,expected_min,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(expected_projector_count,expected_max,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(nwann_min/=nwann_max.or.expected_min/=expected_max)then
      message='inconsistent nonlocal assembly contract across ranks';return
    endif
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
    scale=max(1d0,maxval(abs(matrix)))
    hermiticity_defect=maxval(abs(matrix-conjg(transpose(matrix))))
    if(hermiticity_defect>1d-12*scale)then
      message='nonlocal Hermiticity defect exceeds tolerance';return
    endif
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
