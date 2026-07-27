#include "config.h"
module dg_overlapping_wannier_types
  use iso_fortran_env, only: int64,real64
  use,intrinsic :: ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private

  type,public :: s_dg_wannier_tail
    integer :: center_fragment=0
    integer :: center_owner_rank=-1
    integer :: local_index=0
    integer :: generation=0
    integer(int64) :: geometry_fingerprint=0_int64
    integer(int64),allocatable :: physical_grid_ids(:)
    real(real64),allocatable :: value(:)
    real(real64),allocatable :: gradient(:,:)
  end type

  type,public :: s_dg_overlapping_wannier_basis
    integer :: candidate_rank=0
    integer :: target_rank=0
    integer :: retained_rank=0
    integer :: generation=0
    integer(int64) :: geometry_fingerprint=0_int64
    integer(int64) :: basis_fingerprint=0_int64
    integer(int64),allocatable :: owned_core_physical_ids(:)
    type(s_dg_wannier_tail),allocatable :: tail(:)
  end type

  public :: initialize_dg_wannier_tail
  public :: validate_dg_overlapping_wannier_basis
  public :: checked_dg_wannier_extent_product
  public :: release_dg_overlapping_wannier_basis

contains

  subroutine checked_dg_wannier_extent_product(extents,product_value,ok)
    integer(int64),intent(in) :: extents(:)
    integer(int64),intent(out) :: product_value
    logical,intent(out) :: ok
    integer :: i
    product_value=1_int64
    ok=.false.
    do i=1,size(extents)
      if(extents(i)<0_int64)return
      if(extents(i)==0_int64)then
        product_value=0_int64
        ok=.true.
        return
      endif
      if(product_value>huge(product_value)/extents(i))return
      product_value=product_value*extents(i)
    enddo
    ok=.true.
  end subroutine checked_dg_wannier_extent_product

  subroutine initialize_dg_wannier_tail(tail,center_fragment,center_owner_rank,local_index, &
      generation,geometry_fingerprint,physical_grid_ids,values,gradients,ok,message)
    type(s_dg_wannier_tail),intent(out) :: tail
    integer,intent(in) :: center_fragment,center_owner_rank,local_index,generation
    integer(int64),intent(in) :: geometry_fingerprint,physical_grid_ids(:)
    real(real64),intent(in) :: values(:),gradients(:,:)
    logical,intent(out) :: ok
    character(*),intent(out) :: message
    integer(int64) :: allocation_size
    logical :: extent_ok
    ok=.false.;message=''
    if(center_fragment<=0.or.center_owner_rank<0.or.local_index<=0.or.generation<=0)then
      message='invalid tail center metadata';return
    endif
    if(geometry_fingerprint==0_int64)then
      message='missing tail geometry fingerprint';return
    endif
    if(size(physical_grid_ids)/=size(values).or.size(gradients,1)/=3.or. &
        size(gradients,2)/=size(values))then
      message='incomplete tail value/gradient map';return
    endif
    if(any(physical_grid_ids<=0_int64))then
      message='invalid physical grid id';return
    endif
    if(.not.all(ieee_is_finite(values)).or..not.all(ieee_is_finite(gradients)))then
      message='non-finite tail value or gradient';return
    endif
    call checked_dg_wannier_extent_product([3_int64,int(size(values),int64)],allocation_size,extent_ok)
    if(.not.extent_ok)then
      message='tail allocation extent overflow';return
    endif
    tail%center_fragment=center_fragment
    tail%center_owner_rank=center_owner_rank
    tail%local_index=local_index
    tail%generation=generation
    tail%geometry_fingerprint=geometry_fingerprint
    allocate(tail%physical_grid_ids(size(physical_grid_ids)),tail%value(size(values)))
    allocate(tail%gradient(3,size(values)))
    tail%physical_grid_ids=physical_grid_ids
    tail%value=values
    tail%gradient=gradients
    ok=.true.
  end subroutine initialize_dg_wannier_tail

  subroutine validate_dg_overlapping_wannier_basis(basis,comm,expected_core_count,ok,message, &
      ownership_count,fingerprint)
    type(s_dg_overlapping_wannier_basis),intent(in) :: basis
    integer,intent(in) :: comm
    integer(int64),intent(in) :: expected_core_count
    logical,intent(out) :: ok
    character(*),intent(out) :: message
    integer,intent(out),optional :: ownership_count
    integer(int64),intent(out),optional :: fingerprint
#ifdef USE_MPI
    integer :: rank,nproc,ierr,local_bad,global_bad,i,j,total_core,total_tail
    integer :: scalar_min(4),scalar_max(4),local_scalar(4)
    integer,allocatable :: core_counts(:),tail_counts(:),core_displs(:),tail_displs(:)
    integer,allocatable :: local_centers(:),global_centers(:)
    integer(int64),allocatable :: global_core(:)
    integer(int64) :: fp_min,fp_max,geo_min,geo_max,total64,tail_total64,expected_min,expected_max
    logical :: extent_ok

    ok=.false.;message=''
    if(present(ownership_count))ownership_count=0
    if(present(fingerprint))fingerprint=0_int64
    call MPI_Comm_rank(comm,rank,ierr)
    call MPI_Comm_size(comm,nproc,ierr)
    call MPI_Allreduce(expected_core_count,expected_min,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(expected_core_count,expected_max,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(expected_min/=expected_max)then
      message='inconsistent expected core ownership count across ranks'
      return
    endif

    local_bad=merge(0,1,locally_valid(basis,rank))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then
      message='invalid or stale overlapping-Wannier metadata before payload collective'
      return
    endif

    local_scalar=[basis%candidate_rank,basis%target_rank,basis%retained_rank,basis%generation]
    call MPI_Allreduce(local_scalar,scalar_min,4,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(local_scalar,scalar_max,4,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(basis%geometry_fingerprint,geo_min,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(basis%geometry_fingerprint,geo_max,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    call MPI_Allreduce(basis%basis_fingerprint,fp_min,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(basis%basis_fingerprint,fp_max,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(any(scalar_min/=scalar_max).or.geo_min/=geo_max.or.fp_min/=fp_max)then
      message='inconsistent overlapping-Wannier generation or fingerprints'
      return
    endif

    allocate(core_counts(nproc),tail_counts(nproc),core_displs(nproc),tail_displs(nproc))
    call MPI_Allgather(size(basis%owned_core_physical_ids),1,MPI_INTEGER,core_counts,1,MPI_INTEGER,comm,ierr)
    call MPI_Allgather(size(basis%tail),1,MPI_INTEGER,tail_counts,1,MPI_INTEGER,comm,ierr)
    call checked_count_sum(core_counts,total64,extent_ok)
    if(.not.extent_ok.or.total64>int(huge(total_core),int64))then
      message='core ownership allocation extent overflow';return
    endif
    total_core=int(total64)
    call checked_count_sum(tail_counts,tail_total64,extent_ok)
    if(extent_ok)call checked_dg_wannier_extent_product([2_int64,tail_total64],total64,extent_ok)
    if(.not.extent_ok.or.tail_total64>int(huge(total_tail),int64))then
      message='tail center allocation extent overflow';return
    endif
    total_tail=int(tail_total64)
    core_displs(1)=0;tail_displs(1)=0
    do i=2,nproc
      core_displs(i)=core_displs(i-1)+core_counts(i-1)
      tail_displs(i)=tail_displs(i-1)+tail_counts(i-1)
    enddo
    allocate(global_core(total_core),local_centers(2*size(basis%tail)),global_centers(2*total_tail))
    do i=1,size(basis%tail)
      local_centers(2*i-1)=basis%tail(i)%center_fragment
      local_centers(2*i)=basis%tail(i)%local_index
    enddo
    call MPI_Allgatherv(basis%owned_core_physical_ids,size(basis%owned_core_physical_ids),MPI_INTEGER8, &
      global_core,core_counts,core_displs,MPI_INTEGER8,comm,ierr)
    tail_counts=2*tail_counts;tail_displs=2*tail_displs
    call MPI_Allgatherv(local_centers,size(local_centers),MPI_INTEGER,global_centers,tail_counts, &
      tail_displs,MPI_INTEGER,comm,ierr)

    call sort_int64(global_core)
    if(expected_core_count<=0_int64.or.total_core/=int(expected_core_count))then
      message='missing or extra physical core owner';return
    endif
    do i=1,total_core
      if(global_core(i)/=int(i,int64))then
        message='duplicate or missing physical core owner';return
      endif
    enddo
    call sort_center_pairs(global_centers)
    do i=2,total_tail
      if(all(global_centers(2*i-1:2*i)==global_centers(2*i-3:2*i-2)))then
        message='duplicate authoritative Wannier center owner';return
      endif
    enddo

    ok=.true.;message=''
    if(present(ownership_count))ownership_count=total_core
    if(present(fingerprint))then
      fingerprint=basis%basis_fingerprint
      fingerprint=ieor(ishftc(fingerprint,7),basis%geometry_fingerprint)
      fingerprint=ieor(ishftc(fingerprint,7),int(basis%generation,int64))
      do i=1,total_core
        fingerprint=ieor(ishftc(fingerprint,7),global_core(i))
      enddo
      do i=1,2*total_tail
        fingerprint=ieor(ishftc(fingerprint,7),int(global_centers(i),int64))
      enddo
      if(fingerprint==0_int64)fingerprint=1_int64
    endif
#else
    ok=.false.;message='overlapping-Wannier metadata validation requires MPI'
    if(present(ownership_count))ownership_count=0
    if(present(fingerprint))fingerprint=0_int64
#endif
  end subroutine validate_dg_overlapping_wannier_basis

  subroutine checked_count_sum(counts,total,ok)
    integer,intent(in) :: counts(:)
    integer(int64),intent(out) :: total
    logical,intent(out) :: ok
    integer :: i
    total=0_int64;ok=.false.
    do i=1,size(counts)
      if(counts(i)<0)return
      if(total>huge(total)-int(counts(i),int64))return
      total=total+int(counts(i),int64)
    enddo
    ok=.true.
  end subroutine checked_count_sum

  logical function locally_valid(basis,rank)
    type(s_dg_overlapping_wannier_basis),intent(in) :: basis
    integer,intent(in) :: rank
    integer :: i,j
    locally_valid=.false.
    if(basis%candidate_rank<=0.or.basis%target_rank<=0.or.basis%retained_rank<=0)return
    if(basis%candidate_rank<basis%target_rank.or.basis%retained_rank>basis%target_rank)return
    if(basis%generation<=0.or.basis%geometry_fingerprint==0_int64.or.basis%basis_fingerprint==0_int64)return
    if(.not.allocated(basis%owned_core_physical_ids).or..not.allocated(basis%tail))return
    if(any(basis%owned_core_physical_ids<=0_int64))return
    do i=1,size(basis%owned_core_physical_ids)
      if(any(basis%owned_core_physical_ids(i)==basis%owned_core_physical_ids(i+1:)))return
    enddo
    do i=1,size(basis%tail)
      if(basis%tail(i)%center_fragment<=0.or.basis%tail(i)%local_index<=0)return
      if(basis%tail(i)%center_owner_rank/=rank)return
      if(basis%tail(i)%generation/=basis%generation)return
      if(basis%tail(i)%geometry_fingerprint/=basis%geometry_fingerprint)return
      if(.not.allocated(basis%tail(i)%physical_grid_ids).or..not.allocated(basis%tail(i)%value).or. &
          .not.allocated(basis%tail(i)%gradient))return
      if(size(basis%tail(i)%physical_grid_ids)/=size(basis%tail(i)%value))return
      if(size(basis%tail(i)%gradient,1)/=3.or. &
          size(basis%tail(i)%gradient,2)/=size(basis%tail(i)%value))return
      if(any(basis%tail(i)%physical_grid_ids<=0_int64))return
      if(.not.all(ieee_is_finite(basis%tail(i)%value)).or.&
          .not.all(ieee_is_finite(basis%tail(i)%gradient)))return
      do j=1,size(basis%tail(i)%physical_grid_ids)
        if(any(basis%tail(i)%physical_grid_ids(j)==basis%tail(i)%physical_grid_ids(j+1:)))return
      enddo
    enddo
    locally_valid=.true.
  end function locally_valid

  subroutine sort_int64(values)
    integer(int64),intent(inout) :: values(:)
    integer :: i,j
    integer(int64) :: key
    do i=2,size(values)
      key=values(i);j=i-1
      do while(j>=1)
        if(values(j)<=key)exit
        values(j+1)=values(j);j=j-1
      enddo
      values(j+1)=key
    enddo
  end subroutine sort_int64

  subroutine sort_center_pairs(values)
    integer,intent(inout) :: values(:)
    integer :: i,j,key1,key2
    do i=2,size(values)/2
      key1=values(2*i-1);key2=values(2*i);j=i-1
      do while(j>=1)
        if(values(2*j-1)<key1.or.(values(2*j-1)==key1.and.values(2*j)<=key2))exit
        values(2*j+1:2*j+2)=values(2*j-1:2*j);j=j-1
      enddo
      values(2*j+1)=key1;values(2*j+2)=key2
    enddo
  end subroutine sort_center_pairs

  subroutine release_dg_overlapping_wannier_basis(basis)
    type(s_dg_overlapping_wannier_basis),intent(inout) :: basis
    integer :: i
    if(allocated(basis%tail))then
      do i=1,size(basis%tail)
        if(allocated(basis%tail(i)%physical_grid_ids))deallocate(basis%tail(i)%physical_grid_ids)
        if(allocated(basis%tail(i)%value))deallocate(basis%tail(i)%value)
        if(allocated(basis%tail(i)%gradient))deallocate(basis%tail(i)%gradient)
      enddo
      deallocate(basis%tail)
    endif
    if(allocated(basis%owned_core_physical_ids))deallocate(basis%owned_core_physical_ids)
  end subroutine release_dg_overlapping_wannier_basis
end module dg_overlapping_wannier_types
