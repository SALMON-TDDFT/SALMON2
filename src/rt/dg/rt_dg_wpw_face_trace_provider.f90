module rt_dg_wpw_face_trace_provider
  implicit none
  private

  public :: s_wpw_face_trace_provider
  public :: wpw_face_trace_callback
  public :: bind_wpw_face_trace_provider
  public :: unbind_wpw_face_trace_provider
  public :: evaluate_wpw_face_traces

  abstract interface
    subroutine wpw_face_trace_callback(user_context,k_minus,k_plus,axis,side_from_k_minus, &
                                        unwrapped_grid,w_minus,w_plus,grad_w_minus,grad_w_plus, &
                                        p_minus,p_plus,grad_p_minus,grad_p_plus,info)
      class(*), pointer, intent(inout) :: user_context
      integer, intent(in) :: k_minus,k_plus,axis,side_from_k_minus,unwrapped_grid(3)
      complex(8), intent(out) :: w_minus(:),w_plus(:),grad_w_minus(:,:),grad_w_plus(:,:)
      complex(8), intent(out) :: p_minus(:),p_plus(:),grad_p_minus(:,:),grad_p_plus(:,:)
      integer, intent(out) :: info
    end subroutine wpw_face_trace_callback
  end interface

  type, public :: s_wpw_face_trace_provider
    class(*), pointer :: user_context=>null()
    procedure(wpw_face_trace_callback), pointer, nopass :: callback=>null()
  end type s_wpw_face_trace_provider

contains

  subroutine bind_wpw_face_trace_provider(provider,user_context,callback,info)
    type(s_wpw_face_trace_provider), intent(inout) :: provider
    class(*), pointer, intent(inout) :: user_context
    procedure(wpw_face_trace_callback) :: callback
    integer, intent(out) :: info

    call unbind_wpw_face_trace_provider(provider)
    provider%user_context=>user_context
    ! CALLBACK and USER_CONTEXT must both outlive PROVIDER or be explicitly
    ! detached with unbind_wpw_face_trace_provider before either lifetime ends.
    provider%callback=>callback
    info=0
  end subroutine bind_wpw_face_trace_provider

  subroutine unbind_wpw_face_trace_provider(provider)
    type(s_wpw_face_trace_provider), intent(inout) :: provider

    nullify(provider%callback)
    nullify(provider%user_context)
  end subroutine unbind_wpw_face_trace_provider

  subroutine evaluate_wpw_face_traces(provider,k_minus,k_plus,axis,side_from_k_minus, &
                                       unwrapped_grid,w_minus,w_plus,grad_w_minus,grad_w_plus, &
                                       p_minus,p_plus,grad_p_minus,grad_p_plus,info)
    type(s_wpw_face_trace_provider), intent(inout) :: provider
    integer, intent(in) :: k_minus,k_plus,axis,side_from_k_minus,unwrapped_grid(3)
    complex(8), intent(out) :: w_minus(:),w_plus(:),grad_w_minus(:,:),grad_w_plus(:,:)
    complex(8), intent(out) :: p_minus(:),p_plus(:),grad_p_minus(:,:),grad_p_plus(:,:)
    integer, intent(out) :: info

    w_minus=(0d0,0d0); w_plus=(0d0,0d0)
    grad_w_minus=(0d0,0d0); grad_w_plus=(0d0,0d0)
    p_minus=(0d0,0d0); p_plus=(0d0,0d0)
    grad_p_minus=(0d0,0d0); grad_p_plus=(0d0,0d0)
    if(.not.associated(provider%callback) .or. .not.associated(provider%user_context)) then
      info=1
      return
    end if
    call provider%callback(provider%user_context,k_minus,k_plus,axis,side_from_k_minus, &
      unwrapped_grid,w_minus,w_plus,grad_w_minus,grad_w_plus,p_minus,p_plus, &
      grad_p_minus,grad_p_plus,info)
  end subroutine evaluate_wpw_face_traces

end module rt_dg_wpw_face_trace_provider
