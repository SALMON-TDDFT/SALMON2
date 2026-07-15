module rt_dg_wpw_face_trace_scanner
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use rt_dg_wpw_face_trace_provider, only: s_wpw_face_trace_provider, evaluate_wpw_face_traces
  use rt_dg_wpw_quadrature_assembler, only: assemble_wpw_canonical_face_point
  implicit none
  private
  public :: assemble_wpw_canonical_face_grid

contains

  subroutine assemble_wpw_canonical_face_grid(provider,k_minus,k_plus,axis,side_from_k_minus, &
      minus_box_lo,minus_box_hi,plus_box_lo,plus_box_hi,hgs,w_row_ids,p_column_ids,wp_face_h,info)
    type(s_wpw_face_trace_provider), intent(inout) :: provider
    integer, intent(in) :: k_minus,k_plus,axis,side_from_k_minus
    integer, intent(in) :: minus_box_lo(3),minus_box_hi(3),plus_box_lo(3),plus_box_hi(3)
    real(8), intent(in) :: hgs(3)
    integer, intent(in) :: w_row_ids(:),p_column_ids(:)
    complex(8), intent(out) :: wp_face_h(:,:)
    integer, intent(out) :: info
    integer :: tangent(2),face_lo(3),face_hi(3),unwrapped_grid(3)
    integer :: it,jt,trace_info,point_info
    real(8) :: normal(3),h_normal,face_weight
    complex(8), allocatable :: w_minus(:),w_plus(:),grad_w_minus(:,:),grad_w_plus(:,:)
    complex(8), allocatable :: p_minus(:),p_plus(:),grad_p_minus(:,:),grad_p_plus(:,:)
    complex(8), allocatable :: point_block(:,:),temporary_block(:,:)

    wp_face_h=(0d0,0d0); info=0
    if(k_minus<=0 .or. k_minus>=k_plus .or. axis<1 .or. axis>3 .or. &
       abs(side_from_k_minus)/=1 .or. any(hgs<=0d0) .or. &
       any(minus_box_hi<minus_box_lo) .or. any(plus_box_hi<plus_box_lo) .or. &
       .not.strictly_increasing(w_row_ids) .or. .not.strictly_increasing(p_column_ids) .or. &
       any(shape(wp_face_h)/=[size(w_row_ids),size(p_column_ids)])) then
      info=1
      return
    end if

    select case(axis)
    case(1); tangent=[2,3]
    case(2); tangent=[1,3]
    case(3); tangent=[1,2]
    end select
    face_lo=max(minus_box_lo,plus_box_lo)
    face_hi=min(minus_box_hi,plus_box_hi)
    if(any(face_hi(tangent)<face_lo(tangent))) then
      info=2
      return
    end if
    unwrapped_grid=minus_box_lo
    if(side_from_k_minus==1) unwrapped_grid(axis)=minus_box_hi(axis)
    normal=0d0
    normal(axis)=dble(side_from_k_minus)
    h_normal=hgs(axis)
    face_weight=hgs(tangent(1))*hgs(tangent(2))

    allocate(w_minus(size(w_row_ids)),w_plus(size(w_row_ids)))
    allocate(grad_w_minus(3,size(w_row_ids)),grad_w_plus(3,size(w_row_ids)))
    allocate(p_minus(size(p_column_ids)),p_plus(size(p_column_ids)))
    allocate(grad_p_minus(3,size(p_column_ids)),grad_p_plus(3,size(p_column_ids)))
    allocate(point_block(size(w_row_ids),size(p_column_ids)))
    allocate(temporary_block(size(w_row_ids),size(p_column_ids)))
    temporary_block=(0d0,0d0)

    do jt=face_lo(tangent(2)),face_hi(tangent(2))
      unwrapped_grid(tangent(2))=jt
      do it=face_lo(tangent(1)),face_hi(tangent(1))
        unwrapped_grid(tangent(1))=it
        call evaluate_wpw_face_traces(provider,k_minus,k_plus,axis,side_from_k_minus, &
          unwrapped_grid,w_minus,w_plus,grad_w_minus,grad_w_plus,p_minus,p_plus, &
          grad_p_minus,grad_p_plus,trace_info)
        if(trace_info/=0) then
          info=3
          return
        end if
        if(.not.all(finite_complex(w_minus)) .or. .not.all(finite_complex(w_plus)) .or. &
           .not.all(finite_complex(grad_w_minus)) .or. .not.all(finite_complex(grad_w_plus)) .or. &
           .not.all(finite_complex(p_minus)) .or. .not.all(finite_complex(p_plus)) .or. &
           .not.all(finite_complex(grad_p_minus)) .or. .not.all(finite_complex(grad_p_plus))) then
          info=3
          return
        end if
        call assemble_wpw_canonical_face_point(k_minus,k_plus,w_minus,w_plus,grad_w_minus, &
          grad_w_plus,p_minus,p_plus,grad_p_minus,grad_p_plus,normal,h_normal,face_weight, &
          point_block,point_info)
        if(point_info/=0) then
          info=4
          return
        end if
        temporary_block=temporary_block+point_block
      end do
    end do
    wp_face_h=temporary_block
  end subroutine assemble_wpw_canonical_face_grid

  logical function strictly_increasing(ids) result(ok)
    integer, intent(in) :: ids(:)
    integer :: i
    ok=size(ids)>0
    do i=2,size(ids)
      if(ids(i)<=ids(i-1)) then
        ok=.false.
        return
      end if
    end do
  end function strictly_increasing

  elemental logical function finite_complex(value) result(ok)
    complex(8), intent(in) :: value
    ok=ieee_is_finite(real(value,8)) .and. ieee_is_finite(aimag(value))
  end function finite_complex

end module rt_dg_wpw_face_trace_scanner
