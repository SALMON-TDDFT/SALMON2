program test_calc_nstate_active
  use nstate_active_util, only: calc_nstate_active
  implicit none
  real(8) :: rocc(5,1,1), wtk(1), dropped
  integer :: nact
  ! 占有: [2,2,1.5,1d-5,0] (1 k, 1 spin). threshold 1d-3 -> count(>1d-3)=3
  rocc(:,1,1) = (/2d0, 2d0, 1.5d0, 1d-5, 0d0/); wtk(1)=1d0
  call calc_nstate_active(rocc,5,1,1,wtk, 0, 1d-3, nact, dropped)
  if (nact /= 3) stop 1
  if (abs(dropped - 1d-5) > 1d-12) stop 2
  ! 明示優先: nstate_active_in=4 -> nact=4（threshold無視）
  call calc_nstate_active(rocc,5,1,1,wtk, 4, 1d-3, nact, dropped)
  if (nact /= 4) stop 3
  ! 既定（両0）-> 全 nstate
  call calc_nstate_active(rocc,5,1,1,wtk, 0, 0d0, nact, dropped)
  if (nact /= 5) stop 4
  print *, 'PASS'
end program
