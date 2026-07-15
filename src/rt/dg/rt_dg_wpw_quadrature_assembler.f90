module rt_dg_wpw_quadrature_assembler
  use rt_dg_wpw_weak_form_evaluator, only: wpw_volume_weak_pair, wpw_sipg_face_pair
  use rt_dg_wpw_sparse_builder, only: s_dg_wpw_sparse_candidates, wpw_candidate_volume, &
    wpw_candidate_volume_face
  implicit none
  private
  public :: assemble_wpw_volume_point
  public :: assemble_wpw_canonical_face_point
  public :: pack_wpw_point_candidates

contains

  subroutine assemble_wpw_volume_point(w, grad_w, p_owned, grad_p_owned, p_support, grad_p_support, &
                                        potential, weight, wp_h, wp_s, pp_h, pp_s, info)
    complex(8), intent(in) :: w(:), grad_w(:,:), p_owned(:), grad_p_owned(:,:)
    complex(8), intent(in) :: p_support(:), grad_p_support(:,:)
    real(8), intent(in) :: potential, weight
    complex(8), intent(out) :: wp_h(:,:), wp_s(:,:), pp_h(:,:), pp_s(:,:)
    integer, intent(out) :: info
    integer :: iw, ip, jp, pair_info

    wp_h=(0d0,0d0); wp_s=(0d0,0d0); pp_h=(0d0,0d0); pp_s=(0d0,0d0)
    info=0
    if(size(grad_w,1)/=3 .or. size(grad_w,2)/=size(w) .or. &
       size(grad_p_owned,1)/=3 .or. size(grad_p_owned,2)/=size(p_owned) .or. &
       size(grad_p_support,1)/=3 .or. size(grad_p_support,2)/=size(p_support) .or. &
       any(shape(wp_h)/=[size(w),size(p_owned)]) .or. any(shape(wp_s)/=shape(wp_h)) .or. &
       any(shape(pp_h)/=[size(p_owned),size(p_support)]) .or. any(shape(pp_s)/=shape(pp_h))) then
      info=1
      return
    end if
    do ip=1,size(p_owned)
      do iw=1,size(w)
        call wpw_volume_weak_pair(w(iw),grad_w(:,iw),p_owned(ip),grad_p_owned(:,ip), &
          potential,weight,wp_s(iw,ip),wp_h(iw,ip),pair_info)
        if(pair_info/=0) then; info=2; return; end if
      end do
      do jp=1,size(p_support)
        call wpw_volume_weak_pair(p_owned(ip),grad_p_owned(:,ip),p_support(jp),grad_p_support(:,jp), &
          potential,weight,pp_s(ip,jp),pp_h(ip,jp),pair_info)
        if(pair_info/=0) then; info=2; return; end if
      end do
    end do
  end subroutine assemble_wpw_volume_point

  subroutine assemble_wpw_canonical_face_point(k_minus,k_plus,w_minus,w_plus,grad_w_minus,grad_w_plus, &
                                                p_minus,p_plus,grad_p_minus,grad_p_plus,normal,h_normal, &
                                                face_weight,wp_h,info)
    integer, intent(in) :: k_minus,k_plus
    complex(8), intent(in) :: w_minus(:),w_plus(:),grad_w_minus(:,:),grad_w_plus(:,:)
    complex(8), intent(in) :: p_minus(:),p_plus(:),grad_p_minus(:,:),grad_p_plus(:,:)
    real(8), intent(in) :: normal(3),h_normal,face_weight
    complex(8), intent(out) :: wp_h(:,:)
    integer, intent(out) :: info
    integer :: iw,ip,pair_info
    real(8), parameter :: h1_trace_tolerance=1d-12

    wp_h=(0d0,0d0); info=0
    if(k_minus<=0 .or. k_minus >= k_plus) then; info=1; return; end if
    if(size(w_minus)/=size(w_plus) .or. size(p_minus)/=size(p_plus) .or. &
       size(grad_w_minus,1)/=3 .or. size(grad_w_plus,1)/=3 .or. &
       size(grad_p_minus,1)/=3 .or. size(grad_p_plus,1)/=3 .or. &
       size(grad_w_minus,2)/=size(w_minus) .or. size(grad_w_plus,2)/=size(w_plus) .or. &
       size(grad_p_minus,2)/=size(p_minus) .or. size(grad_p_plus,2)/=size(p_plus) .or. &
       any(shape(wp_h)/=[size(w_minus),size(p_minus)])) then; info=2; return; end if
    if(any(abs(p_minus-p_plus)>h1_trace_tolerance)) then; info=4; return; end if
    do ip=1,size(p_minus)
      do iw=1,size(w_minus)
        call wpw_sipg_face_pair(w_minus(iw),w_plus(iw),grad_w_minus(:,iw),grad_w_plus(:,iw), &
          p_minus(ip),p_plus(ip),grad_p_minus(:,ip),grad_p_plus(:,ip),normal,h_normal,face_weight, &
          wp_h(iw,ip),pair_info)
        if(pair_info/=0) then; info=3; return; end if
      end do
    end do
  end subroutine assemble_wpw_canonical_face_point

  subroutine pack_wpw_point_candidates(w_row_ids,p_owned_ids,p_support_ids,wp_h,wp_s,pp_h,pp_s, &
                                       wp_face_h,candidates,info)
    integer, intent(in) :: w_row_ids(:),p_owned_ids(:),p_support_ids(:)
    complex(8), intent(in) :: wp_h(:,:),wp_s(:,:),pp_h(:,:),pp_s(:,:),wp_face_h(:,:)
    type(s_dg_wpw_sparse_candidates), intent(out) :: candidates
    integer, intent(out) :: info
    integer :: iw,ip,jp,ientry,nwp,npp

    info=0
    if(.not. strictly_increasing(w_row_ids) .or. .not. strictly_increasing(p_owned_ids) .or. &
       .not. strictly_increasing(p_support_ids) .or. &
       any(shape(wp_h)/=[size(w_row_ids),size(p_owned_ids)]) .or. any(shape(wp_s)/=shape(wp_h)) .or. &
       any(shape(wp_face_h)/=shape(wp_h)) .or. &
       any(shape(pp_h)/=[size(p_owned_ids),size(p_support_ids)]) .or. any(shape(pp_s)/=shape(pp_h))) then
      info=1
      return
    end if
    nwp=size(w_row_ids)*size(p_owned_ids); npp=size(p_owned_ids)*size(p_support_ids)
    allocate(candidates%wp_w_row_ids(nwp),candidates%wp_pw_col_ids(nwp),candidates%wp_origin(nwp))
    allocate(candidates%wp_h_values(nwp,1),candidates%wp_s_values(nwp,1))
    allocate(candidates%pp_pw_row_ids(npp),candidates%pp_pw_col_ids(npp),candidates%pp_origin(npp))
    allocate(candidates%pp_h_values(npp,1),candidates%pp_s_values(npp,1))
    ientry=0
    do ip=1,size(p_owned_ids)
      do iw=1,size(w_row_ids)
        ientry=ientry+1
        candidates%wp_w_row_ids(ientry)=w_row_ids(iw)
        candidates%wp_pw_col_ids(ientry)=p_owned_ids(ip)
        candidates%wp_origin(ientry)=wpw_candidate_volume_face
        candidates%wp_h_values(ientry,1)=wp_h(iw,ip)+wp_face_h(iw,ip)
        candidates%wp_s_values(ientry,1)=wp_s(iw,ip)
      end do
    end do
    ientry=0
    do ip=1,size(p_owned_ids)
      do jp=1,size(p_support_ids)
        ientry=ientry+1
        candidates%pp_pw_row_ids(ientry)=p_owned_ids(ip)
        candidates%pp_pw_col_ids(ientry)=p_support_ids(jp)
        candidates%pp_origin(ientry)=wpw_candidate_volume
        candidates%pp_h_values(ientry,1)=pp_h(ip,jp)
        candidates%pp_s_values(ientry,1)=pp_s(ip,jp)
      end do
    end do
  end subroutine pack_wpw_point_candidates

  logical function strictly_increasing(ids) result(ok)
    integer,intent(in)::ids(:)
    integer::i
    ok=.true.
    do i=2,size(ids)
      if(ids(i)<=ids(i-1)) then; ok=.false.; return; end if
    end do
  end function strictly_increasing
end module rt_dg_wpw_quadrature_assembler
