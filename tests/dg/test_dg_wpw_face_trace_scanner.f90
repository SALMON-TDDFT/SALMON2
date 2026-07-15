program test_dg_wpw_face_trace_scanner
  use rt_dg_wpw_face_trace_provider, only: s_wpw_face_trace_provider, &
    bind_wpw_face_trace_provider, unbind_wpw_face_trace_provider
  use rt_dg_wpw_face_trace_scanner, only: assemble_wpw_canonical_face_grid
  use rt_dg_wpw_quadrature_assembler, only: assemble_wpw_canonical_face_point
  implicit none

  type :: s_mock_context
    integer :: point_count=0
    integer :: visited(3,4)=0
    logical :: identity_ok=.true.
  end type s_mock_context

  type(s_mock_context), target :: mock
  type(s_wpw_face_trace_provider) :: provider
  complex(8) :: actual(1,1),expected(1,1),point_value(1,1)
  complex(8) :: wm(1),wp(1),gwm(3,1),gwp(3,1)
  complex(8) :: pm(1),pp(1),gpm(3,1),gpp(3,1)
  integer :: info,point_info,iy,iz,grid(3),expected_grids(3,4)
  integer, parameter :: w_ids(1)=[7],p_ids(1)=[11]
  real(8), parameter :: hgs(3)=[0.5d0,0.25d0,0.4d0]

  expected_grids=reshape([2,1,1, 2,2,1, 2,1,2, 2,2,2],[3,4])
  call bind_wpw_face_trace_provider(provider,mock,mock_face_traces,info)
  if(info/=0) error stop 1
  call assemble_wpw_canonical_face_grid(provider,1,2,1,1,[1,1,1],[2,2,2], &
    [2,1,1],[3,2,2],hgs,w_ids,p_ids,actual,info)
  if(info/=0) error stop 2
  if(mock%point_count/=4 .or. .not.mock%identity_ok) error stop 3
  if(any(mock%visited/=expected_grids)) error stop 4

  expected=(0d0,0d0)
  do iz=1,2
    do iy=1,2
      grid=[2,iy,iz]
      call deterministic_traces(grid,wm,wp,gwm,gwp,pm,pp,gpm,gpp)
      call assemble_wpw_canonical_face_point(1,2,wm,wp,gwm,gwp,pm,pp,gpm,gpp, &
        [1d0,0d0,0d0],hgs(1),hgs(2)*hgs(3),point_value,point_info)
      if(point_info/=0) error stop 5
      expected=expected+point_value
    end do
  end do
  if(maxval(abs(actual-expected))>1d-13) error stop 6
  call unbind_wpw_face_trace_provider(provider)
  print '(a)', 'PASS canonical WPW face trace numerical fixture'

contains

  subroutine mock_face_traces(user_context,k_minus,k_plus,axis,side_from_k_minus, &
      unwrapped_grid,w_minus,w_plus,grad_w_minus,grad_w_plus,p_minus,p_plus, &
      grad_p_minus,grad_p_plus,callback_info)
    class(*), pointer, intent(inout) :: user_context
    integer, intent(in) :: k_minus,k_plus,axis,side_from_k_minus,unwrapped_grid(3)
    complex(8), intent(out) :: w_minus(:),w_plus(:),grad_w_minus(:,:),grad_w_plus(:,:)
    complex(8), intent(out) :: p_minus(:),p_plus(:),grad_p_minus(:,:),grad_p_plus(:,:)
    integer, intent(out) :: callback_info
    select type(context=>user_context)
    type is(s_mock_context)
      context%point_count=context%point_count+1
      if(context%point_count<=size(context%visited,2)) &
        context%visited(:,context%point_count)=unwrapped_grid
      context%identity_ok=context%identity_ok .and. k_minus==1 .and. k_plus==2 .and. &
        axis==1 .and. side_from_k_minus==1
      call deterministic_traces(unwrapped_grid,w_minus,w_plus,grad_w_minus,grad_w_plus, &
        p_minus,p_plus,grad_p_minus,grad_p_plus)
      callback_info=0
    class default
      callback_info=1
    end select
  end subroutine mock_face_traces

  subroutine deterministic_traces(grid,w_minus,w_plus,grad_w_minus,grad_w_plus, &
      p_minus,p_plus,grad_p_minus,grad_p_plus)
    integer, intent(in) :: grid(3)
    complex(8), intent(out) :: w_minus(:),w_plus(:),grad_w_minus(:,:),grad_w_plus(:,:)
    complex(8), intent(out) :: p_minus(:),p_plus(:),grad_p_minus(:,:),grad_p_plus(:,:)
    real(8) :: y,z
    y=dble(grid(2)); z=dble(grid(3))
    w_minus=cmplx(1d0+0.1d0*y+0.01d0*z,0.02d0*y,kind=8)
    w_plus=w_minus+cmplx(0.5d0,-0.1d0,kind=8)
    grad_w_minus(:,1)=[cmplx(0.2d0,0.01d0,8),cmplx(0.3d0,0d0,8),cmplx(0.4d0,0d0,8)]
    grad_w_plus(:,1)=[cmplx(-0.1d0,0.02d0,8),cmplx(0.2d0,0d0,8),cmplx(0.3d0,0d0,8)]
    p_minus=cmplx(2d0+0.05d0*y+0.02d0*z,-0.03d0*z,kind=8)
    p_plus=p_minus
    grad_p_minus(:,1)=[cmplx(0.1d0,0d0,8),cmplx(0.05d0,0d0,8),cmplx(0.02d0,-0.03d0,8)]
    grad_p_plus=grad_p_minus
  end subroutine deterministic_traces

end program test_dg_wpw_face_trace_scanner
