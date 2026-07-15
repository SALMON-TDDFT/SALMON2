program test_dg_wpw_weak_form_evaluator
  use rt_dg_wpw_weak_form_evaluator, only: wpw_volume_weak_pair, wpw_sipg_face_pair
  implicit none
  complex(8) :: u, up, v, vp, gu(3), gup(3), gv(3), gvp(3)
  complex(8) :: suv, huv, svu, hvu, face_uv, face_vu, face_flipped, face_zero
  real(8) :: n(3)
  integer :: info
  real(8), parameter :: tol = 2.0d-13

  u = (1.0d0,0.2d0); up = (-0.1d0,0.4d0)
  v = (0.3d0,-0.5d0); vp = (0.2d0,0.1d0)
  gu = [(0.2d0,0.1d0),(0.3d0,-0.2d0),(-0.1d0,0.4d0)]
  gup = [(0.1d0,-0.3d0),(0.2d0,0.5d0),(0.4d0,0.1d0)]
  gv = [(-0.2d0,0.4d0),(0.6d0,0.1d0),(0.3d0,-0.2d0)]
  gvp = [(0.5d0,0.2d0),(-0.1d0,0.3d0),(0.2d0,0.4d0)]
  n = [1.0d0,0.0d0,0.0d0]

  call wpw_volume_weak_pair(u, gu, v, gv, 0.7d0, 0.25d0, suv, huv, info)
  if (info /= 0) error stop 1
  call wpw_volume_weak_pair(v, gv, u, gu, 0.7d0, 0.25d0, svu, hvu, info)
  if (info /= 0 .or. abs(suv-conjg(svu)) > tol .or. abs(huv-conjg(hvu)) > tol) error stop 2

  call wpw_sipg_face_pair(u, up, gu, gup, v, vp, gv, gvp, n, 0.4d0, 0.4d0, face_uv, info)
  if (info /= 0) error stop 3
  call wpw_sipg_face_pair(v, vp, gv, gvp, u, up, gu, gup, n, 0.4d0, 0.4d0, face_vu, info)
  if (info /= 0 .or. abs(face_uv-conjg(face_vu)) > tol) error stop 4
  call wpw_sipg_face_pair(up, u, gup, gu, vp, v, gvp, gv, -n, 0.4d0, 0.4d0, face_flipped, info)
  if (info /= 0 .or. abs(face_uv-face_flipped) > tol) error stop 5
  call wpw_sipg_face_pair(u, u, gu, gu, v, v, gv, gv, n, 0.4d0, 0.4d0, face_zero, info)
  if (info /= 0 .or. abs(face_zero) > 0.0d0) error stop 6
  write(*,'(a)') 'PASS WPW SIPG weak-form numerical fixture'
end program test_dg_wpw_weak_form_evaluator
