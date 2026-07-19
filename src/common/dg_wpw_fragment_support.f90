!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!
module dg_wpw_fragment_support
  implicit none
  private

  type, public :: s_dg_wpw_support_record
    integer :: fragment_id=0
    integer :: image_shift(3)=0
    integer :: periodic_shift(3)=0
  end type s_dg_wpw_support_record

  public :: build_dg_wpw_fragment_support
  public :: dg_wpw_support_overlap_box

contains

  subroutine dg_wpw_support_overlap_box(local_core_lo,local_core_hi,neighbor_displacement,buffer,&
      overlap_lo,overlap_hi,info)
    integer,intent(in)::local_core_lo(3),local_core_hi(3),neighbor_displacement(3),buffer(3)
    integer,intent(out)::overlap_lo(3),overlap_hi(3),info
    integer::extent(3),neighbor_lo(3),neighbor_hi(3)
    overlap_lo=0;overlap_hi=-1;info=1
    extent=local_core_hi-local_core_lo+1
    if(any(extent<=0).or.any(buffer<0).or.any(buffer>extent).or.&
       any(abs(neighbor_displacement)>1).or.all(neighbor_displacement==0))return
    neighbor_lo=local_core_lo+neighbor_displacement*extent
    neighbor_hi=neighbor_lo+extent-1
    overlap_lo=max(local_core_lo,neighbor_lo-buffer)
    overlap_hi=min(local_core_hi,neighbor_hi+buffer)
    if(any(overlap_hi<overlap_lo))then;info=2;return;endif
    info=0
  end subroutine dg_wpw_support_overlap_box

  subroutine build_dg_wpw_fragment_support(fragment_id,nfrag_axis,support_extent,records, &
      unique_fragments,info)
    integer,intent(in)::fragment_id,nfrag_axis(3),support_extent(3)
    type(s_dg_wpw_support_record),allocatable,intent(out)::records(:)
    integer,allocatable,intent(out)::unique_fragments(:)
    integer,intent(out)::info
    integer::coord(3),target(3),raw_target(3),offset(3),lo(3),hi(3)
    integer::dx,dy,dz,nrecord,capacity,destination,nunique,i,j,key

    info=0
    if(any(nfrag_axis<1).or.any(support_extent<0).or.any(support_extent>1).or.&
       fragment_id<1.or.fragment_id>product(nfrag_axis))then
      allocate(records(0),unique_fragments(0));info=1;return
    endif
    coord(1)=(fragment_id-1)/(nfrag_axis(2)*nfrag_axis(3))
    coord(2)=modulo(fragment_id-1,nfrag_axis(2)*nfrag_axis(3))/nfrag_axis(3)
    coord(3)=modulo(fragment_id-1,nfrag_axis(3))
    do i=1,3
      if(nfrag_axis(i)>1)then
        lo(i)=-support_extent(i);hi(i)=support_extent(i)
      else
        lo(i)=0;hi(i)=0
      endif
    enddo
    capacity=product(hi-lo+1)-1
    allocate(records(max(0,capacity)))
    nrecord=0
    do dz=lo(3),hi(3)
      do dy=lo(2),hi(2)
        do dx=lo(1),hi(1)
          offset=[dx,dy,dz]
          if(all(offset==0))cycle
          raw_target=coord+offset
          target=modulo(raw_target,nfrag_axis)
          destination=target(1)*nfrag_axis(2)*nfrag_axis(3)+target(2)*nfrag_axis(3)+target(3)+1
          nrecord=nrecord+1
          records(nrecord)%fragment_id=destination
          records(nrecord)%image_shift=offset
          records(nrecord)%periodic_shift=(raw_target-target)/nfrag_axis
        enddo
      enddo
    enddo

    allocate(unique_fragments(nrecord+1))
    unique_fragments(1)=fragment_id
    do i=1,nrecord;unique_fragments(i+1)=records(i)%fragment_id;enddo
    do i=2,size(unique_fragments)
      key=unique_fragments(i);j=i-1
      do while(j>=1)
        if(unique_fragments(j)<=key)exit
        unique_fragments(j+1)=unique_fragments(j);j=j-1
      enddo
      unique_fragments(j+1)=key
    enddo
    nunique=0
    do i=1,size(unique_fragments)
      if(i==1.or.unique_fragments(i)/=unique_fragments(max(1,i-1)))then
        nunique=nunique+1;unique_fragments(nunique)=unique_fragments(i)
      endif
    enddo
    if(nunique<size(unique_fragments))unique_fragments=unique_fragments(1:nunique)
  end subroutine build_dg_wpw_fragment_support

end module dg_wpw_fragment_support
