program test_dg_wpw_face_trace_scanner
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use rt_dg_wpw_face_trace_provider, only: s_wpw_face_trace_provider, &
    bind_wpw_face_trace_provider, unbind_wpw_face_trace_provider
  use rt_dg_wpw_face_trace_scanner, only: assemble_wpw_canonical_face_grid
  use rt_dg_wpw_quadrature_assembler, only: assemble_wpw_canonical_face_point
  implicit none

  type :: s_mock_context
    integer :: point_count=0
    integer :: visited(3,4)=0
    logical :: identity_ok=.true.
    integer :: failure_after=huge(1)
    integer :: nan_point=0
    integer :: p_mismatch_point=0
    integer :: expected_side=1
  end type s_mock_context

  type(s_mock_context), target :: mock
  type(s_wpw_face_trace_provider) :: provider
  complex(8) :: actual(1,1),expected(1,1),point_value(1,1),wrapped_minus(1,1),wrapped_plus(1,1)
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

  mock%point_count=0; mock%failure_after=2; actual=cmplx(9d0,0d0,8)
  call assemble_wpw_canonical_face_grid(provider,1,2,1,1,[1,1,1],[2,2,2], &
    [2,1,1],[3,2,2],hgs,w_ids,p_ids,actual,info)
  if(info==0 .or. mock%point_count/=3 .or. any(actual/=(0d0,0d0))) error stop 7
  mock%failure_after=huge(1)

  mock%point_count=0; mock%nan_point=2; actual=cmplx(9d0,0d0,8)
  call assemble_wpw_canonical_face_grid(provider,1,2,1,1,[1,1,1],[2,2,2], &
    [2,1,1],[3,2,2],hgs,w_ids,p_ids,actual,info)
  if(info/=3 .or. any(actual/=(0d0,0d0))) error stop 8
  mock%nan_point=0

  mock%point_count=0; mock%p_mismatch_point=2; actual=cmplx(9d0,0d0,8)
  call assemble_wpw_canonical_face_grid(provider,1,2,1,1,[1,1,1],[2,2,2], &
    [2,1,1],[3,2,2],hgs,w_ids,p_ids,actual,info)
  if(info/=4 .or. any(actual/=(0d0,0d0))) error stop 9
  mock%p_mismatch_point=0

  actual=cmplx(9d0,0d0,8)
  call assemble_wpw_canonical_face_grid(provider,2,1,1,1,[1,1,1],[2,2,2], &
    [2,1,1],[3,2,2],hgs,w_ids,p_ids,actual,info)
  if(info==0 .or. any(actual/=(0d0,0d0))) error stop 10

  actual=cmplx(9d0,0d0,8)
  call assemble_wpw_canonical_face_grid(provider,1,2,4,1,[1,1,1],[2,2,2], &
    [2,1,1],[3,2,2],hgs,w_ids,p_ids,actual,info)
  if(info==0 .or. any(actual/=(0d0,0d0))) error stop 11

  actual=cmplx(9d0,0d0,8)
  call assemble_wpw_canonical_face_grid(provider,1,2,1,0,[1,1,1],[2,2,2], &
    [2,1,1],[3,2,2],hgs,w_ids,p_ids,actual,info)
  if(info==0 .or. any(actual/=(0d0,0d0))) error stop 12

  actual=cmplx(9d0,0d0,8)
  call assemble_wpw_canonical_face_grid(provider,1,2,1,1,[1,1,1],[2,2,2], &
    [2,4,1],[3,5,2],hgs,w_ids,p_ids,actual,info)
  if(info==0 .or. any(actual/=(0d0,0d0))) error stop 13

  mock%point_count=0; mock%identity_ok=.true.; mock%expected_side=-1
  call assemble_wpw_canonical_face_grid(provider,1,2,1,-1,[1,1,1],[2,2,2], &
    [1,1,1],[2,2,2],hgs,w_ids,p_ids,wrapped_minus,info)
  if(info/=0 .or. .not.mock%identity_ok .or. mock%point_count/=4 .or. &
    any(mock%visited(1,:)/=1)) error stop 14
  mock%point_count=0; mock%visited=0; mock%identity_ok=.true.; mock%expected_side=1
  call assemble_wpw_canonical_face_grid(provider,1,2,1,1,[1,1,1],[2,2,2], &
    [1,1,1],[2,2,2],hgs,w_ids,p_ids,wrapped_plus,info)
  if(info/=0 .or. .not.mock%identity_ok .or. mock%point_count/=4 .or. &
    any(mock%visited(1,:)/=2)) error stop 15
  if(abs(wrapped_minus(1,1)-wrapped_plus(1,1))<1d-13) error stop 16
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
      if(context%point_count>context%failure_after) then
        callback_info=7
        return
      end if
      if(context%point_count<=size(context%visited,2)) &
        context%visited(:,context%point_count)=unwrapped_grid
      context%identity_ok=context%identity_ok .and. k_minus==1 .and. k_plus==2 .and. &
        axis==1 .and. side_from_k_minus==context%expected_side
      call deterministic_traces(unwrapped_grid,w_minus,w_plus,grad_w_minus,grad_w_plus, &
        p_minus,p_plus,grad_p_minus,grad_p_plus)
      if(context%point_count==context%nan_point) &
        w_minus(1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,kind=8)
      if(context%point_count==context%p_mismatch_point) p_plus(1)=p_minus(1)+1d-6
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
    real(8) :: x,y,z
    x=dble(grid(1)); y=dble(grid(2)); z=dble(grid(3))
    w_minus=cmplx(1d0+0.03d0*x+0.1d0*y+0.01d0*z,0.02d0*y,kind=8)
    w_plus=w_minus+cmplx(0.5d0,-0.1d0,kind=8)
    grad_w_minus(:,1)=[cmplx(0.2d0,0.01d0,8),cmplx(0.3d0,0d0,8),cmplx(0.4d0,0d0,8)]
    grad_w_plus(:,1)=[cmplx(-0.1d0,0.02d0,8),cmplx(0.2d0,0d0,8),cmplx(0.3d0,0d0,8)]
    p_minus=cmplx(2d0+0.04d0*x+0.05d0*y+0.02d0*z,-0.03d0*z,kind=8)
    p_plus=p_minus
    grad_p_minus(:,1)=[cmplx(0.1d0,0d0,8),cmplx(0.05d0,0d0,8),cmplx(0.02d0,-0.03d0,8)]
    grad_p_plus=grad_p_minus
  end subroutine deterministic_traces

end program test_dg_wpw_face_trace_scanner
