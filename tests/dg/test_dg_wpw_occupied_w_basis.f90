program test_dg_wpw_occupied_w_basis
  use dg_wpw_occupied_w_basis,only:s_dg_wpw_occupied_w_basis,initialize_dg_wpw_occupied_w_basis
  implicit none
  type(s_dg_wpw_occupied_w_basis)::basis
  integer::global_keys(5,4),owners(4),local_keys(5,2),grid_ids(3),info
  integer(8)::saved_fingerprint
  complex(8)::core_values(3,2),buffer_values(5,2),buffer_gradients(3,5,2)

  global_keys(:,1)=[8,9,0,0,0]
  global_keys(:,2)=[1,2,0,0,0]
  global_keys(:,3)=[5,6,-1,0,0]
  global_keys(:,4)=[3,4,0,1,0]
  owners=[2,1,2,1]
  ! Fragment 1 keys in the lexicographic global order: (1,2,...) then (3,4,...).
  local_keys(:,1)=global_keys(:,2);local_keys(:,2)=global_keys(:,4)
  grid_ids=[11,12,13]
  core_values=reshape([(cmplx(dble(info),-0.1d0*dble(info),8),info=1,6)],[3,2])
  buffer_values=reshape([(cmplx(0.2d0*dble(info),0.05d0*dble(info),8),info=1,10)],[5,2])
  buffer_gradients=(0.1d0,0.02d0)
  call initialize_dg_wpw_occupied_w_basis(basis,1,global_keys,owners,local_keys,grid_ids,core_values,&
    [-1,-1,-1],[3,-1,-1],buffer_values,buffer_gradients,2.5d0,7,info)
  if(info/=0.or..not.basis%valid)error stop 1
  if(any(basis%owned_ids/=[1,2]).or.basis%global_count/=4.or.&
    basis%local_count/=2.or.any(basis%stable_keys/=local_keys).or.basis%epoch/=7.or.&
    basis%fingerprint==0_8)error stop 2
  saved_fingerprint=basis%fingerprint

  global_keys(:,4)=global_keys(:,2)
  call initialize_dg_wpw_occupied_w_basis(basis,1,global_keys,owners,local_keys,grid_ids,core_values,&
    [-1,-1,-1],[3,-1,-1],buffer_values,buffer_gradients,2.5d0,8,info)
  if(info==0.or..not.basis%valid.or.basis%fingerprint/=saved_fingerprint.or.basis%epoch/=7)error stop 3

  global_keys(:,4)=[3,4,0,1,0]
  call initialize_dg_wpw_occupied_w_basis(basis,1,global_keys,owners,local_keys,grid_ids,core_values,&
    [-1,-1,-1],[3,-1,-1],buffer_values(1:4,:),buffer_gradients(:,1:4,:),2.5d0,8,info)
  if(info==0.or.basis%fingerprint/=saved_fingerprint.or.basis%epoch/=7)error stop 4
  call initialize_dg_wpw_occupied_w_basis(basis,1,global_keys,owners,local_keys,grid_ids,core_values,&
    [-1,-1,-1],[3,-1,-1],buffer_values,buffer_gradients,-1d0,8,info)
  if(info==0.or.basis%fingerprint/=saved_fingerprint.or.basis%epoch/=7)error stop 5

  print '(a)','PASS transactional occupied-W stable ID descriptor'
end program test_dg_wpw_occupied_w_basis
