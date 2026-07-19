program test_dg_wpw_volume_halo_packer
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use rt_dg_wpw_volume_halo_provider,only:s_dg_wpw_volume_halo_send,pack_dg_wpw_volume_halo_send
  implicit none
  type(s_dg_wpw_volume_halo_send)::send
  complex(8)::values(2,4),gradients(3,2,4)
  integer::info,i
  do i=1,4
    values(:,i)=[cmplx(10+i,0d0,8),cmplx(20+i,0d0,8)]
    gradients(:,:,i)=cmplx(i,0d0,8)
  enddo
  call pack_dg_wpw_volume_halo_send(3,[1,0,0],[1,0,0],[17,1,1],[18,1,1],[16,8,8],&
    [16,1,1],[19,1,1],[5,9],values,gradients,send,info)
  if(info/=0.or.send%peer/=3.or.any(send%image_shift/=[1,0,0]))error stop 1
  if(any(send%box_lo/=[1,1,1]).or.any(send%box_hi/=[2,1,1]))error stop 2
  if(any(send%w_ids/=[5,9]).or.any(abs(send%values(:,1)-values(:,2))>1d-13).or.&
    any(abs(send%values(:,2)-values(:,3))>1d-13))error stop 3
  if(any(abs(send%gradients(:,:,1)-gradients(:,:,2))>1d-13))error stop 4
  values(1,2)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  call pack_dg_wpw_volume_halo_send(3,[1,0,0],[1,0,0],[17,1,1],[18,1,1],[16,8,8],&
    [16,1,1],[19,1,1],[5,9],values,gradients,send,info)
  if(info==0.or.allocated(send%w_ids).or.allocated(send%values).or.allocated(send%gradients))error stop 5
  values(1,2)=(12d0,0d0)
  call pack_dg_wpw_volume_halo_send(3,[1,0,0],[0,0,0],[5,1,1],[6,1,1],[16,8,8],&
    [4,1,1],[7,1,1],[5,9],values,gradients,send,info)
  if(info/=0.or.any(send%box_lo/=[5,1,1]).or.any(send%box_hi/=[6,1,1]))error stop 6
  print '(a)','PASS bounded WPW volume-halo payload packer'
end program
