program test_dg_wpw_g_modes
  use dg_wpw_g_modes,only:select_dg_wpw_g_modes
  implicit none
  integer,allocatable::indices(:,:)
  real(8),allocatable::gvec(:,:)
  integer::info
  integer,parameter::expected(3,6)=reshape([-1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1],[3,6])
  call select_dg_wpw_g_modes([4d0,4d0,4d0],1.3d0,4,indices,gvec,info)
  if(info/=0)error stop 1
  if(size(indices,2)/=6.or.any(indices/=expected))error stop 2
  if(maxval(abs(gvec-0.5d0*acos(-1d0)*dble(indices)))>1d-13)error stop 3
  call select_dg_wpw_g_modes([4d0,4d0,4d0],1.3d0,1,indices,gvec,info)
  if(info/=0.or.size(indices,2)/=6)error stop 4
  call select_dg_wpw_g_modes([0d0,4d0,4d0],1.3d0,4,indices,gvec,info)
  if(info==0.or.allocated(indices).or.allocated(gvec))error stop 5
  print '(a)','PASS GS-neutral complete-shell WPW G-mode selection'
end program
