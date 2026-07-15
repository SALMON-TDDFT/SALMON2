module rt_dg_wpw_canonical_faces
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use rt_dg_wpw_window, only: wpw_face_neighbor_fragment
  implicit none
  private
  type,public::s_wpw_canonical_face_list
    integer::nface=0
    integer,allocatable::k_minus(:),k_plus(:),axis(:),side_from_k_minus(:)
  end type s_wpw_canonical_face_list
  public::build_wpw_canonical_face_list

contains
  subroutine build_wpw_canonical_face_list(dg_frag,faces,info)
    type(s_dg_fragment_rt),intent(in)::dg_frag
    type(s_wpw_canonical_face_list),intent(out)::faces
    integer,intent(out)::info
    integer::ifrag,jfrag,axis,side,n
    info=0; faces%nface=0
    if(dg_frag%n_frag<=0 .or. dg_frag%ifrag_start<1 .or. dg_frag%ifrag_end<dg_frag%ifrag_start .or. &
       dg_frag%ifrag_end>dg_frag%n_frag .or. any(dg_frag%num_fragment<=0) .or. &
       .not.allocated(dg_frag%ixyz_frag) .or. .not.allocated(dg_frag%nxyz_domain)) then
      info=1
      return
    end if
    n=0
    do ifrag=dg_frag%ifrag_start,dg_frag%ifrag_end
      do axis=1,3
        if(dg_frag%num_fragment(axis)<=1) cycle
        do side=-1,1,2
          jfrag=wpw_face_neighbor_fragment(dg_frag,ifrag,axis,side)
          if(jfrag<=0 .or. jfrag>dg_frag%n_frag) then; info=2; return; end if
          if (ifrag >= jfrag) cycle
          n=n+1
        end do
      end do
    end do
    allocate(faces%k_minus(n),faces%k_plus(n),faces%axis(n),faces%side_from_k_minus(n))
    n=0
    do ifrag=dg_frag%ifrag_start,dg_frag%ifrag_end
      do axis=1,3
        if(dg_frag%num_fragment(axis)<=1) cycle
        do side=-1,1,2
          jfrag=wpw_face_neighbor_fragment(dg_frag,ifrag,axis,side)
          if (ifrag >= jfrag) cycle
          n=n+1
          faces%k_minus(n)=ifrag
          faces%k_plus(n)=jfrag
          faces%axis(n)=axis
          faces%side_from_k_minus(n)=side
        end do
      end do
    end do
    faces%nface=n
  end subroutine build_wpw_canonical_face_list
end module rt_dg_wpw_canonical_faces
