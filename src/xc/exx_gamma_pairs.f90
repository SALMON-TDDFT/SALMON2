! Exact Gamma orbital pairs. Each unordered pair shares one real-kernel convolution.
! Occupations are source weights (rocc/2), independent of target columns; no division
! by a weight is used, including for empty states and empty orbital partitions.
module exx_gamma_pairs
  use iso_fortran_env,only:int64
  use ieee_arithmetic,only:ieee_is_finite
  use communication,only:comm_get_groupinfo,comm_get_max,comm_bcast,comm_summation
  use fftw_blocks,only:mesh_transform,block_layout
  implicit none
  private
  public::gamma_pair_action
contains
  subroutine gamma_pair_action(n,dims,coords,groups,comm_r,comm_o,target,weight,multiplier,action, &
                               status,pairs,scratch)
    implicit none
    integer,intent(in)::n(3),dims(:),coords(:),groups(:),comm_r,comm_o
    complex(8),intent(in)::target(:,:,:)
    real(8),intent(in)::weight(:),multiplier(:)
    complex(8),intent(out)::action(:,:,:)
    integer,intent(out)::status
    integer(int64),intent(out)::pairs,scratch
    integer,parameter::batch=4
    integer::rank,np,bad,m(3),lo(3),ng,nt,owner,count,i,j,g,nb,first,global_i,total,local_first
    integer,allocatable::counts(:),selected(:)
    complex(8),allocatable::source(:),reverse(:),received(:),density(:,:),spectrum(:,:)
    real(8)::source_weight,weight_hi,weight_neg_hi
    status=1;action=0d0;pairs=0;scratch=0;bad=0
    call comm_get_groupinfo(comm_o,rank,np)
    call block_layout(n,dims,coords,m,lo,bad)
    ng=product(m);nt=size(target,2)
    if(size(target,1)/=ng.or.size(target,3)/=1.or.any(shape(action)/=shape(target)))bad=1
    if(size(weight)/=nt.or.size(multiplier)/=ng)bad=1
    call reject_collectively()
    if(bad/=0)return
    do i=1,size(weight)
      if(.not.ieee_is_finite(weight(i)))then
        bad=1
      else if(weight(i)<0d0.or.weight(i)>1d0)then
        bad=1
      endif
    enddo
    do i=1,size(multiplier)
      if(.not.ieee_is_finite(multiplier(i)))bad=1
    enddo
    do j=1,nt;do g=1,size(target,1)
      if(.not.ieee_is_finite(real(target(g,j,1))).or..not.ieee_is_finite(aimag(target(g,j,1))))bad=1
    enddo;enddo
    call reject_collectively()
    if(bad/=0)return
    allocate(counts(np))
    do owner=0,np-1
      count=nt
      call comm_bcast(count,comm_o,owner)
      counts(owner+1)=count
    enddo
    ! Spatial peers must own the same columns and participate in matching FFT calls.
    do i=1,np
      count=counts(i);call comm_get_max(count,comm_r)
      if(count/=counts(i))bad=1
    enddo
    total=0
    do i=1,np
      if(counts(i)>huge(total)-total)bad=1
      if(bad==0)total=total+counts(i)
    enddo
    call reject_collectively()
    if(bad/=0)return
    do i=1,nt
      weight_hi=weight(i);weight_neg_hi=-weight(i)
      call max_weight(weight_hi);call max_weight(weight_neg_hi)
      if(weight_hi/=-weight_neg_hi)bad=1
    enddo
    call reject_collectively()
    if(bad/=0)return
    local_first=sum(counts(:rank))+1
    allocate(source(ng),reverse(ng),received(ng),density(ng,batch),spectrum(ng,batch),selected(nt))
    scratch=16_int64*ng*(3+2*batch)+int(storage_size(1)/8,int64)*(np+nt)
    do owner=0,np-1
      do i=1,counts(owner+1)
        global_i=sum(counts(:owner))+i
        if(rank==owner)then
          source=target(:,i,1);source_weight=weight(i)
        endif
        call comm_bcast(source,comm_o,owner);call comm_bcast(source_weight,comm_o,owner)
        count=0
        do j=1,nt
          if(local_first+j-1<global_i)cycle
          count=count+1;selected(count)=j
        enddo
        reverse=0d0
        do first=1,count,batch
          nb=min(batch,count-first+1);density=0d0
          do j=1,nb
            density(:,j)=conjg(source)*target(:,selected(first+j-1),1)
          enddo
          call mesh_transform(n,dims,coords,groups,density(:,:nb),spectrum(:,:nb),-1,status)
          bad=status;call comm_get_max(bad,comm_r)
          ! Only spatial peers perform these FFTs; orbital peers can own other counts.
          if(bad/=0)exit
          do j=1,nb
            spectrum(:,j)=spectrum(:,j)*multiplier
          enddo
          call mesh_transform(n,dims,coords,groups,spectrum(:,:nb),density(:,:nb),1,status)
          bad=status;call comm_get_max(bad,comm_r)
          if(bad/=0)exit
          do j=1,nb
            g=selected(first+j-1)
            action(:,g,1)=action(:,g,1)-source_weight*source*density(:,j)
            if(local_first+g-1/=global_i)reverse=reverse-weight(g)*target(:,g,1)*conjg(density(:,j))
          enddo
          pairs=pairs+nb
        enddo
        ! Reduce status before reverse traffic, even when another orbital group failed an FFT.
        call reject_collectively()
        if(bad/=0)return
        call comm_summation(reverse,received,ng,comm_o,owner)
        if(rank==owner)action(:,i,1)=action(:,i,1)+received
      enddo
    enddo
    do j=1,nt;do g=1,ng
      if(.not.ieee_is_finite(real(action(g,j,1))).or..not.ieee_is_finite(aimag(action(g,j,1))))bad=1
    enddo;enddo
    call reject_collectively()
    if(bad==0)status=0
  contains
    subroutine max_weight(w)
      implicit none
      real(8),intent(inout)::w
      real(8)::a(1),b(1)
      a(1)=w;call comm_get_max(a,b,1,comm_r);w=b(1)
    end subroutine
    subroutine reject_collectively()
      implicit none
      call comm_get_max(bad,comm_r);call comm_get_max(bad,comm_o)
    end subroutine
  end subroutine
end module
