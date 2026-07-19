program test_dg_wpw_volume_halo_box_mapping
  use rt_dg_wpw_volume_halo_provider,only:map_dg_wpw_volume_halo_box_to_canonical
  implicit none
  integer::lo(3),hi(3),info
  call map_dg_wpw_volume_halo_box_to_canonical([17,1,1],[18,2,2],[16,8,8],[1,0,0],lo,hi,info)
  if(info/=0.or.any(lo/=[1,1,1]).or.any(hi/=[2,2,2]))error stop 1
  call map_dg_wpw_volume_halo_box_to_canonical([1,1,1],[2,2,2],[16,8,8],[-1,0,0],lo,hi,info)
  if(info/=0.or.any(lo/=[17,1,1]).or.any(hi/=[18,2,2]))error stop 2
  call map_dg_wpw_volume_halo_box_to_canonical([1,1,1],[2,2,2],[16,8,8],[2,0,0],lo,hi,info)
  if(info==0)error stop 3
  print '(a)','PASS periodic WPW volume-halo canonical box mapping'
end program
