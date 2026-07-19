program test_dg_wpw_volume_header_bounds
  use,intrinsic::iso_fortran_env,only:int64
  use rt_dg_wpw_volume_halo_provider,only:dg_wpw_volume_header_is_bounded
  implicit none
  if(.not.dg_wpw_volume_header_is_bounded([3,4,1,1,1,2,2,2],4,[2,2,2],8_int64))error stop 1
  if(dg_wpw_volume_header_is_bounded([3,5,1,1,1,2,2,2],4,[2,2,2],8_int64))error stop 2
  if(dg_wpw_volume_header_is_bounded([3,4,1,1,1,3,2,2],4,[2,2,2],8_int64))error stop 3
  if(dg_wpw_volume_header_is_bounded([3,4,1,1,1,50000,50000,2],4,[50000,50000,2],&
    1000_int64))error stop 4
  print '(a)','PASS bounded volume-halo header validation'
end program
