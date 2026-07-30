#include "config.h"
module dg_overlapping_wannier_operators
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::assemble_dg_overlapping_wannier_weak_operators,&
    assemble_dg_overlapping_wannier_weak_operator_rows
contains
  subroutine assemble_dg_overlapping_wannier_weak_operator_rows(comm,nwann,row_ids,core_ids,weights,&
      values,gradients,local_potential,expected_core_count,kinetic_rows,potential_rows,&
      ownership_count,ok,message)
    integer,intent(in)::comm,nwann
    integer(int64),intent(in)::row_ids(:),core_ids(:),expected_core_count
    real(real64),intent(in)::weights(:),local_potential(:)
    complex(real64),intent(in)::values(:,:),gradients(:,:,:)
    complex(real64),allocatable,intent(out)::kinetic_rows(:,:),potential_rows(:,:)
    integer,intent(out)::ownership_count
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer,parameter::row_batch_size=32
    integer::rank,nproc,ierr,local_bad,global_bad,total_count,total_rows,p,i,j,r,nrows,&
      batch_first,batch_count,nwann_min,nwann_max
    integer,allocatable::core_counts(:),core_displs(:),row_counts(:),row_displs(:)
    integer(int64),allocatable::all_core_ids(:),all_row_ids(:),validation_ids(:)
    integer(int64)::expected_min,expected_max
    complex(real64),allocatable::partial_t(:,:),partial_v(:,:),reduced_t(:,:),reduced_v(:,:)
    logical::finite_payload,shapes_valid

    ok=.false.;message='';ownership_count=0;local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(nwann,nwann_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(nwann,nwann_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(expected_core_count,expected_min,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(expected_core_count,expected_max,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    if(nwann_min/=nwann_max.or.expected_min/=expected_max)local_bad=1
    shapes_valid=size(values,1)==nwann.and.size(values,2)==size(core_ids).and.&
      size(gradients,1)==3.and.size(gradients,2)==nwann.and.size(gradients,3)==size(core_ids)
    finite_payload=all(ieee_is_finite(weights)).and.all(ieee_is_finite(local_potential))
    if(nwann<=0.or.expected_core_count<=0_int64.or.size(weights)/=size(core_ids).or.&
        size(local_potential)/=size(core_ids).or..not.shapes_valid) local_bad=1
    if(any(row_ids<1_int64).or.any(row_ids>int(nwann,int64)).or.any(core_ids<=0_int64).or.&
        any(weights<=0d0))local_bad=1
    if(shapes_valid)then
      finite_payload=finite_payload.and.all(ieee_is_finite(real(values))).and.&
        all(ieee_is_finite(aimag(values))).and.all(ieee_is_finite(real(gradients))).and.&
        all(ieee_is_finite(aimag(gradients)))
    endif
    if(.not.finite_payload)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid row-owned weak-operator payload';return;endif

    allocate(core_counts(nproc),core_displs(nproc),row_counts(nproc),row_displs(nproc))
    call MPI_Allgather(size(core_ids),1,MPI_INTEGER,core_counts,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allgather(size(row_ids),1,MPI_INTEGER,row_counts,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0)then;message='weak-operator ownership metadata collective failed';return;endif
    total_count=0;total_rows=0
    do r=1,nproc
      core_displs(r)=total_count;row_displs(r)=total_rows
      if(core_counts(r)<0.or.row_counts(r)<0.or.total_count>huge(total_count)-core_counts(r).or.&
          total_rows>huge(total_rows)-row_counts(r))local_bad=1
      if(local_bad==0)then
        total_count=total_count+core_counts(r);total_rows=total_rows+row_counts(r)
      endif
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.int(total_count,int64)/=expected_core_count.or.total_rows/=nwann)then
      message='missing or extra row/core owner in weak-operator assembly';return
    endif
    allocate(all_core_ids(total_count),all_row_ids(total_rows))
    call MPI_Allgatherv(core_ids,size(core_ids),MPI_INTEGER8,all_core_ids,core_counts,core_displs,&
      MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allgatherv(row_ids,size(row_ids),MPI_INTEGER8,all_row_ids,row_counts,row_displs,&
      MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0)then;message='weak-operator ownership payload collective failed';return;endif
    allocate(validation_ids(max(total_count,total_rows)))
    validation_ids(1:total_count)=all_core_ids;call sort_ids(validation_ids(1:total_count))
    do i=1,total_count
      if(validation_ids(i)/=int(i,int64))local_bad=1
    enddo
    validation_ids(1:total_rows)=all_row_ids;call sort_ids(validation_ids(1:total_rows))
    do i=1,total_rows
      if(validation_ids(i)/=int(i,int64))local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='duplicate or missing row/core owner in weak-operator assembly';return;endif

    allocate(kinetic_rows(size(row_ids),nwann),potential_rows(size(row_ids),nwann))
    kinetic_rows=(0d0,0d0);potential_rows=(0d0,0d0)
    do r=0,nproc-1
      nrows=row_counts(r+1)
      do batch_first=1,nrows,row_batch_size
        batch_count=min(row_batch_size,nrows-batch_first+1)
        allocate(partial_t(batch_count,nwann),partial_v(batch_count,nwann),&
          reduced_t(batch_count,nwann),reduced_v(batch_count,nwann))
        partial_t=(0d0,0d0);partial_v=(0d0,0d0)
        do p=1,size(core_ids)
          do j=1,nwann;do i=1,batch_count
            partial_t(i,j)=partial_t(i,j)+0.5d0*weights(p)*sum(conjg(gradients(:,&
              int(all_row_ids(row_displs(r+1)+batch_first+i-1)),p))*gradients(:,j,p))
            partial_v(i,j)=partial_v(i,j)+weights(p)*local_potential(p)*conjg(values(&
              int(all_row_ids(row_displs(r+1)+batch_first+i-1)),p))*values(j,p)
          enddo;enddo
        enddo
        call MPI_Reduce(partial_t,reduced_t,batch_count*nwann,MPI_DOUBLE_COMPLEX,MPI_SUM,r,comm,ierr)
        if(ierr/=MPI_SUCCESS)local_bad=1
        call MPI_Reduce(partial_v,reduced_v,batch_count*nwann,MPI_DOUBLE_COMPLEX,MPI_SUM,r,comm,ierr)
        if(ierr/=MPI_SUCCESS)local_bad=1
        if(rank==r)then
          kinetic_rows(batch_first:batch_first+batch_count-1,:)=reduced_t
          potential_rows(batch_first:batch_first+batch_count-1,:)=reduced_v
        endif
        deallocate(partial_t,partial_v,reduced_t,reduced_v)
      enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0)then;message='weak-operator row reduction failed';return;endif
    ownership_count=total_count;ok=.true.
#else
    ok=.false.;message='row-owned weak overlapping-Wannier operators require MPI'
    ownership_count=0
#endif
  end subroutine

  subroutine assemble_dg_overlapping_wannier_weak_operators(comm,nwann,core_ids,weights,values,&
      gradients,local_potential,expected_core_count,kinetic,potential,ownership_count,ok,message)
    integer,intent(in)::comm,nwann
    integer(int64),intent(in)::core_ids(:),expected_core_count
    real(real64),intent(in)::weights(:),local_potential(:)
    complex(real64),intent(in)::values(:,:),gradients(:,:,:)
    complex(real64),allocatable,intent(out)::kinetic(:,:),potential(:,:)
    integer,intent(out)::ownership_count
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nproc,ierr,local_bad,global_bad,total_count,i,j,p,matrix_count,nwann_min,nwann_max
    integer,allocatable::counts(:),displs(:)
    integer(int64),allocatable::all_ids(:)
    integer(int64)::matrix_count64,expected_min,expected_max
    complex(real64),allocatable::local_kinetic(:,:),local_potential_matrix(:,:)
    logical::finite_payload,shapes_valid
    real(real64)::hermiticity_defect,scale
    ok=.false.;message='';ownership_count=0;local_bad=0
    call MPI_Allreduce(nwann,nwann_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(nwann,nwann_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(expected_core_count,expected_min,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(expected_core_count,expected_max,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(nwann_min/=nwann_max.or.expected_min/=expected_max)then
      message='inconsistent weak-operator assembly contract across ranks';return
    endif
    finite_payload=all(ieee_is_finite(weights)).and.all(ieee_is_finite(local_potential))
    if(nwann<=0.or.expected_core_count<=0_int64)local_bad=1
    if(size(weights)/=size(core_ids).or.size(local_potential)/=size(core_ids))local_bad=1
    shapes_valid=size(values,1)==nwann.and.size(values,2)==size(core_ids).and.&
      size(gradients,1)==3.and.size(gradients,2)==nwann.and.size(gradients,3)==size(core_ids)
    if(.not.shapes_valid)local_bad=1
    if(shapes_valid)then
      do p=1,size(values,2);do i=1,size(values,1)
        finite_payload=finite_payload.and.ieee_is_finite(real(values(i,p))).and.&
          ieee_is_finite(aimag(values(i,p)))
        finite_payload=finite_payload.and.all(ieee_is_finite(real(gradients(:,i,p)))).and.&
          all(ieee_is_finite(aimag(gradients(:,i,p))))
      enddo;enddo
    endif
    if(any(core_ids<=0_int64).or.any(weights<=0d0).or..not.finite_payload)local_bad=1
    if(nwann>0.and.int(nwann,int64)<=huge(1_int64)/int(nwann,int64))then
      matrix_count64=int(nwann,int64)*int(nwann,int64)
      if(matrix_count64>int(huge(matrix_count),int64))local_bad=1
    else
      matrix_count64=0_int64;local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid weak-operator unique-core payload';return;endif
    matrix_count=int(matrix_count64)

    call MPI_Comm_size(comm,nproc,ierr);allocate(counts(nproc),displs(nproc))
    call MPI_Allgather(size(core_ids),1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,ierr)
    total_count=0;displs(1)=0
    do i=1,nproc
      if(counts(i)<0.or.total_count>huge(total_count)-counts(i))then
        local_bad=1;exit
      endif
      if(i>1)displs(i)=total_count
      total_count=total_count+counts(i)
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.int(total_count,int64)/=expected_core_count)then
      message='missing or extra weak-operator core owner';return
    endif
    allocate(all_ids(total_count))
    call MPI_Allgatherv(core_ids,size(core_ids),MPI_INTEGER8,all_ids,counts,displs,MPI_INTEGER8,comm,ierr)
    call sort_ids(all_ids)
    do i=1,total_count
      if(all_ids(i)/=int(i,int64))then
        message='duplicate or missing weak-operator core owner';return
      endif
    enddo

    allocate(local_kinetic(nwann,nwann),local_potential_matrix(nwann,nwann))
    local_kinetic=(0d0,0d0);local_potential_matrix=(0d0,0d0)
    do p=1,size(core_ids)
      do j=1,nwann;do i=1,nwann
        local_kinetic(i,j)=local_kinetic(i,j)+0.5d0*weights(p)*&
          sum(conjg(gradients(:,i,p))*gradients(:,j,p))
        local_potential_matrix(i,j)=local_potential_matrix(i,j)+weights(p)*local_potential(p)*&
          conjg(values(i,p))*values(j,p)
      enddo;enddo
    enddo
    allocate(kinetic(nwann,nwann),potential(nwann,nwann))
    call MPI_Allreduce(local_kinetic,kinetic,matrix_count,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_potential_matrix,potential,matrix_count,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    scale=max(1d0,max(maxval(abs(kinetic)),maxval(abs(potential))))
    hermiticity_defect=max(maxval(abs(kinetic-conjg(transpose(kinetic)))),&
      maxval(abs(potential-conjg(transpose(potential)))))
    if(hermiticity_defect>1d-12*scale)then
      message='weak-operator Hermiticity defect exceeds tolerance';return
    endif
    kinetic=0.5d0*(kinetic+conjg(transpose(kinetic)))
    potential=0.5d0*(potential+conjg(transpose(potential)))
    ownership_count=total_count;ok=.true.
#else
    ok=.false.;message='weak overlapping-Wannier operators require MPI';ownership_count=0
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
end module dg_overlapping_wannier_operators
