module dg_wpw_canonical_face_schedule
  implicit none
  private
  type,public::s_dg_wpw_canonical_face_record
    integer::neighbor_fragment=0,axis=0,side_from_local=0
    integer::neighbor_displacement(3)=0,periodic_shift(3)=0
    logical::canonical_owner=.false.
  end type
  public::build_dg_wpw_canonical_face_schedule
contains
  subroutine build_dg_wpw_canonical_face_schedule(fragment_id,nfrag_axis,faces,info)
    integer,intent(in)::fragment_id,nfrag_axis(3)
    type(s_dg_wpw_canonical_face_record),allocatable,intent(out)::faces(:)
    integer,intent(out)::info
    integer::coord(3),raw(3),target(3),offset(3),axis,side,nface,neighbor
    info=1
    if(any(nfrag_axis<1).or.fragment_id<1.or.fragment_id>product(nfrag_axis))then
      allocate(faces(0));return
    endif
    coord(1)=(fragment_id-1)/(nfrag_axis(2)*nfrag_axis(3))
    coord(2)=modulo(fragment_id-1,nfrag_axis(2)*nfrag_axis(3))/nfrag_axis(3)
    coord(3)=modulo(fragment_id-1,nfrag_axis(3))
    allocate(faces(2*count(nfrag_axis>1)));nface=0
    do axis=1,3
      if(nfrag_axis(axis)<=1)cycle
      do side=-1,1,2
        offset=0;offset(axis)=side;raw=coord+offset;target=modulo(raw,nfrag_axis)
        neighbor=target(1)*nfrag_axis(2)*nfrag_axis(3)+target(2)*nfrag_axis(3)+target(3)+1
        nface=nface+1;faces(nface)%neighbor_fragment=neighbor;faces(nface)%axis=axis
        faces(nface)%side_from_local=side;faces(nface)%neighbor_displacement=offset
        faces(nface)%periodic_shift=(raw-target)/nfrag_axis
        faces(nface)%canonical_owner=fragment_id<neighbor
      enddo
    enddo
    info=0
  end subroutine
end module dg_wpw_canonical_face_schedule
