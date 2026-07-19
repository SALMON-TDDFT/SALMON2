program test_dg_wpw_lcfo_volume_operator
  use dg_wpw_lcfo_volume_operator,only:rebuild_dg_wpw_ww_volume,rebuild_dg_wpw_wp_pp_volume
  implicit none
  complex(8)::kin(1,1),nonlocal(1,1),face(1,1),potential(1,1),hww(1,1)
  complex(8)::vkin(2),vnonlocal(2),vpotential(2),hwp(2)
  integer::info
  kin=(1d0,0.2d0);nonlocal=(2d0,-0.1d0);face=(3d0,0.4d0);potential=(4d0,-0.5d0)
  call rebuild_dg_wpw_ww_volume(kin,nonlocal,face,potential,hww,info)
  if(info/=0.or.maxval(abs(hww-(kin+nonlocal+face+potential)))>1d-14)error stop 1
  vkin=[(1d0,0.2d0),(2d0,-0.1d0)];vnonlocal=[(3d0,0.3d0),(4d0,-0.2d0)]
  vpotential=[(5d0,0.4d0),(6d0,-0.3d0)]
  call rebuild_dg_wpw_wp_pp_volume(vkin,vnonlocal,vpotential,hwp,info)
  if(info/=0.or.maxval(abs(hwp-(vkin+vnonlocal+vpotential)))>1d-14)error stop 2
  print '(a)','PASS LCFO WW/WP/PP potential-volume replacement preserves fixed components'
end program
