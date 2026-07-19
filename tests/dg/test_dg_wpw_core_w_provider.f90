program test_dg_wpw_core_w_provider
  use rt_dg_wpw_volume_halo_provider,only:s_dg_wpw_volume_halo_state
  use dg_wpw_core_w_provider,only:evaluate_dg_wpw_core_w_support
  implicit none
  type(s_dg_wpw_volume_halo_state)::halo(1),duplicate(2)
  integer::owned_ids(1),support_ids(2),info
  complex(8)::owned_value(1),owned_gradient(3,1),values(2),gradients(3,2)
  owned_ids=[4];support_ids=[4,9];owned_value=[(2d0,0d0)];owned_gradient=(1d0,0d0)
  halo(1)%epoch=7;halo(1)%valid=.true.;halo(1)%box_lo=[3,2,1];halo(1)%box_hi=[3,2,1]
  allocate(halo(1)%w_ids(1),halo(1)%values(1,1),halo(1)%gradients(3,1,1))
  halo(1)%w_ids=[9];halo(1)%values(1,1)=(5d0,-1d0);halo(1)%gradients(:,1,1)=[(2d0,0d0),(3d0,0d0),(4d0,0d0)]
  call evaluate_dg_wpw_core_w_support(owned_ids,owned_value,owned_gradient,support_ids,halo,&
    [3,2,1],7,values,gradients,info)
  if(info/=0.or.abs(values(1)-(2d0,0d0))>1d-13.or.abs(values(2)-(5d0,-1d0))>1d-13)error stop 1
  if(any(abs(gradients(:,2)-[ (2d0,0d0),(3d0,0d0),(4d0,0d0) ])>1d-13))error stop 2
  call evaluate_dg_wpw_core_w_support(owned_ids,owned_value,owned_gradient,support_ids,halo,&
    [3,2,1],8,values,gradients,info)
  if(info==0.or.any(abs(values)>0d0).or.any(abs(gradients)>0d0))error stop 3
  duplicate=[halo(1),halo(1)]
  call evaluate_dg_wpw_core_w_support(owned_ids,owned_value,owned_gradient,support_ids,duplicate,&
    [3,2,1],7,values,gradients,info)
  if(info==0.or.any(abs(values)>0d0))error stop 4
  call evaluate_dg_wpw_core_w_support(owned_ids,owned_value,owned_gradient,support_ids,halo,&
    [1,1,1],7,values,gradients,info,zero_outside_halo=.true.)
  if(info/=0.or.abs(values(1)-(2d0,0d0))>1d-13.or.abs(values(2))>1d-13.or.&
    any(abs(gradients(:,2))>1d-13))error stop 5
  halo(1)%w_ids=[8]
  call evaluate_dg_wpw_core_w_support(owned_ids,owned_value,owned_gradient,support_ids,halo,&
    [3,2,1],7,values,gradients,info,zero_outside_halo=.true.)
  if(info/=0.or.abs(values(1)-(2d0,0d0))>1d-13.or.abs(values(2))>1d-13.or.&
    any(abs(gradients(:,2))>1d-13))error stop 6
  print '(a)','PASS owned/halo WPW core W provider is epoch-safe'
end program
