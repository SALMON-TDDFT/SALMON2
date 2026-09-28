! Column-streamed linear algebra on a Cartesian spatial/orbital process grid.
! Only band matrices are replicated; mesh arrays contain local orbital columns.
module exx_orbitals
  use communication, only: comm_get_groupinfo,comm_get_max,comm_summation,comm_bcast
  use hse_ace, only: hse_ace_state
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: orbital_layout,orbital_check,orbital_overlap,orbital_rotate
  public :: orbital_ace_build,orbital_ace_apply
contains
  subroutine orbital_check(bad,comm_r,comm_o)
    integer,intent(inout) :: bad
    integer,intent(in) :: comm_r,comm_o
    call comm_get_max(bad,comm_r)
    call comm_get_max(bad,comm_o)
  end subroutine

  subroutine orbital_layout(nlocal,comm_r,comm_o,counts,first,status)
    integer,intent(in) :: nlocal,comm_r,comm_o
    integer,allocatable,intent(out) :: counts(:)
    integer,intent(out) :: first,status
    integer :: rank,np,owner,largest
    call comm_get_groupinfo(comm_o,rank,np)
    allocate(counts(0:np-1))
    do owner=0,np-1
      counts(owner)=nlocal
      call comm_bcast(counts(owner),comm_o,owner)
    enddo
    status=0
    do owner=0,np-1
      largest=counts(owner)
      call comm_get_max(largest,comm_r)
      if(largest/=counts(owner).or.counts(owner)<0)status=1
    enddo
    call orbital_check(status,comm_r,comm_o)
    first=1+sum(counts(:rank-1))
  end subroutine

  subroutine orbital_overlap(left,right,dv,comm_r,comm_o,counts,first,matrix,phase)
    ! left/right have the same contiguous column ownership. Optional phase is
    ! applied to right grid rows. Global small matrix is identical on all peers.
    complex(8),intent(in) :: left(:,:),right(:,:)
    real(8),intent(in) :: dv
    integer,intent(in) :: comm_r,comm_o,counts(0:),first
    complex(8),intent(out) :: matrix(:,:)
    complex(8),intent(in),optional :: phase(:)
    complex(8),allocatable :: column(:),total(:,:)
    integer :: owner,rank,np,j,k
    call comm_get_groupinfo(comm_o,rank,np)
    allocate(column(size(right,1)),total(size(matrix,1),size(matrix,2)))
    matrix=0d0;k=0
    do owner=0,np-1
      do j=1,counts(owner)
        k=k+1
        if(rank==owner)column=right(:,j)
        call comm_bcast(column,comm_o,owner)
        if(present(phase))column=column*phase
        matrix(first:first+size(left,2)-1,k)=matmul(conjg(transpose(left)),column)*dv
      enddo
    enddo
    call comm_summation(matrix,total,size(matrix),comm_r)
    call comm_summation(total,matrix,size(matrix),comm_o)
  end subroutine

  subroutine orbital_rotate(input,matrix,comm_o,counts,first,output,weights)
    complex(8),intent(in) :: input(:,:),matrix(:,:)
    integer,intent(in) :: comm_o,counts(0:),first
    complex(8),intent(out) :: output(:,:)
    real(8),intent(in),optional :: weights(:)
    complex(8),allocatable :: column(:)
    integer :: rank,np,owner,i,j,k
    call comm_get_groupinfo(comm_o,rank,np)
    allocate(column(size(input,1)))
    output=0d0;k=0
    do owner=0,np-1
      do i=1,counts(owner)
        k=k+1
        if(rank==owner)then
          column=input(:,i)
          if(present(weights))column=column*weights(i)
        endif
        call comm_bcast(column,comm_o,owner)
        do j=1,size(output,2)
          output(:,j)=output(:,j)+column*matrix(k,first+j-1)
        enddo
      enddo
    enddo
  end subroutine

  subroutine orbital_ace_build(ace,u,w,dv,comm_r,comm_o,status)
    type(hse_ace_state),intent(inout) :: ace
    complex(8),intent(in) :: u(:,:,:),w(:,:,:)
    real(8),intent(in) :: dv
    integer,intent(in) :: comm_r,comm_o
    integer,intent(out) :: status
    integer,allocatable :: counts(:)
    complex(8),allocatable :: metric(:,:),work(:)
    real(8),allocatable :: e(:),rwork(:)
    integer :: first,n,j,bad,nonzero
    real(8) :: scale
    status=1;bad=0
    if(allocated(ace%factors))deallocate(ace%factors)
    if(any(shape(u)/=shape(w)).or.size(u,3)/=1.or.dv<=0d0.or..not.ieee_is_finite(dv))bad=1
    if(.not.all(ieee_is_finite(real(u))).or..not.all(ieee_is_finite(aimag(u))))bad=1
    if(.not.all(ieee_is_finite(real(w))).or..not.all(ieee_is_finite(aimag(w))))bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    call orbital_layout(size(u,2),comm_r,comm_o,counts,first,bad)
    if(bad/=0.or.sum(counts)<1)return
    n=sum(counts)
    allocate(metric(n,n),work(2*n),e(n),rwork(max(1,3*n-2)))
    nonzero=0
    if(any(w/=(0d0,0d0)))nonzero=1
    call orbital_check(nonzero,comm_r,comm_o)
    ace%dv=dv;ace%condition=0d0
    if(nonzero==0)then
      allocate(ace%factors(size(u,1),size(u,2),1));ace%factors=0d0;status=0;return
    endif
    call orbital_overlap(u(:,:,1),w(:,:,1),-dv,comm_r,comm_o,counts,first,metric)
    if(.not.all(ieee_is_finite(real(metric))).or..not.all(ieee_is_finite(aimag(metric))))bad=1
    scale=sqrt(sum(abs(metric)**2))
    if(scale==0d0.or.sqrt(sum(abs(metric-transpose(conjg(metric)))**2))>1d-10*scale)bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    metric=.5d0*(metric+transpose(conjg(metric)))
    call zheev('V','U',n,metric,n,e,work,size(work),rwork,bad)
    if(bad/=0)bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    if(e(n)<=0d0.or.e(1)<=1d-12*e(n))bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    ace%condition=e(n)/e(1)
    do j=1,n
      metric(:,j)=metric(:,j)/sqrt(e(j))
    enddo
    allocate(ace%factors(size(u,1),size(u,2),1))
    call orbital_rotate(w(:,:,1),metric,comm_o,counts,first,ace%factors(:,:,1))
    status=0
  end subroutine

  subroutine orbital_ace_apply(ace,target,action,comm_r,comm_o,status)
    type(hse_ace_state),intent(in) :: ace
    complex(8),intent(in) :: target(:,:,:)
    complex(8),intent(out) :: action(:,:,:)
    integer,intent(in) :: comm_r,comm_o
    integer,intent(out) :: status
    integer,allocatable :: counts(:),target_counts(:)
    complex(8),allocatable :: column(:),overlap(:),total(:)
    integer :: first,rank,np,owner,i,j,bad
    status=1;action=0d0;bad=0
    if(.not.allocated(ace%factors))bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    if(size(ace%factors,1)/=size(target,1).or.size(ace%factors,3)/=1.or.size(target,3)/=1)bad=1
    if(any(shape(target)/=shape(action)).or.ace%dv<=0d0.or..not.ieee_is_finite(ace%dv))bad=1
    if(.not.all(ieee_is_finite(real(target))).or..not.all(ieee_is_finite(aimag(target))))bad=1
    if(.not.all(ieee_is_finite(real(ace%factors))).or..not.all(ieee_is_finite(aimag(ace%factors))))bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    call orbital_layout(size(target,2),comm_r,comm_o,target_counts,first,bad)
    if(bad/=0)return
    call orbital_layout(size(ace%factors,2),comm_r,comm_o,counts,first,bad)
    if(bad/=0)return
    allocate(column(size(target,1)),overlap(size(target,2)),total(size(target,2)))
    call comm_get_groupinfo(comm_o,rank,np)
    do owner=0,np-1
      do i=1,counts(owner)
        if(rank==owner)column=ace%factors(:,i,1)
        call comm_bcast(column,comm_o,owner)
        overlap=matmul(conjg(column),target(:,:,1))*ace%dv
        call comm_summation(overlap,total,size(overlap),comm_r)
        do j=1,size(target,2)
          action(:,j,1)=action(:,j,1)-column*total(j)
        enddo
      enddo
    enddo
    status=0
  end subroutine
end module
