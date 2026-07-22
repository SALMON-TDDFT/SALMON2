program test_dg_wpw_occupied_w_basis
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_positive_inf
  use dg_wpw_occupied_w_basis,only:s_dg_wpw_occupied_w_basis,t_dg_wpw_periodic_image_mismatch,&
    initialize_dg_wpw_occupied_w_basis,&
    evaluate_dg_wpw_occupied_w_point,dg_wpw_unwrapped_to_storage_index,&
    reorder_dg_wpw_fragment_buffer,extract_dg_wpw_canonical_cell
  implicit none
  type(s_dg_wpw_occupied_w_basis)::basis,changed_basis
  integer::global_keys(5,4),owners(4),local_keys(5,2),grid_ids(3),info
  integer(8)::saved_fingerprint
  complex(8)::core_values(3,2),buffer_values(5,2),buffer_gradients(3,5,2)
  complex(8)::periodized_values(2),periodized_gradients(3,2)
  complex(8)::storage3d(4,4,4,1),unwrapped3d(4,4,4,1)
  integer::sx,sy,sz,ix,iy,iz
  complex(8),allocatable::p6(:,:,:,:),p10(:,:,:,:),p5(:,:,:,:),canonical6(:,:,:,:),canonical10(:,:,:,:)
  complex(8),allocatable::canonical_bad(:,:,:,:),pzero(:,:,:,:),canonical_zero(:,:,:,:)
  integer::canonical_grid(3),origin(3)
  type(t_dg_wpw_periodic_image_mismatch)::mismatch

  ! Fragment storage is core/positive-buffer followed by negative-buffer.
  ! Logical P must instead be contiguous from 1-B through n+B.
  if(any([(dg_wpw_unwrapped_to_storage_index(info,12,6),info=-5,18)]/=&
    [19,20,21,22,23,24,(info,info=1,18)]))error stop 200
  if(any([(dg_wpw_unwrapped_to_storage_index(info,12,10),info=-9,22)]/=&
    [23,24,25,26,27,28,29,30,31,32,(info,info=1,22)]))error stop 201
  do sz=1,4;do sy=1,4;do sx=1,4
    storage3d(sx,sy,sz,1)=cmplx(sx+10*sy+100*sz,0d0,8)
  enddo;enddo;enddo
  call reorder_dg_wpw_fragment_buffer(storage3d,[2,2,2],[1,1,1],unwrapped3d,info)
  if(info/=0)error stop 202
  do iz=1,4;do iy=1,4;do ix=1,4
    sx=dg_wpw_unwrapped_to_storage_index(ix-1,2,1)
    sy=dg_wpw_unwrapped_to_storage_index(iy-1,2,1)
    sz=dg_wpw_unwrapped_to_storage_index(iz-1,2,1)
    if(abs(unwrapped3d(ix,iy,iz,1)-storage3d(sx,sy,sz,1))>0d0)error stop 203
  enddo;enddo;enddo

  origin=[12,0,12]
  allocate(p6(24,24,24,1),p10(32,32,32,1),p5(22,22,22,1),&
    canonical6(24,24,24,1),canonical10(24,24,24,1))
  do iz=1,24;do iy=1,24;do ix=1,24
    canonical_grid=modulo(origin+[ix-6,iy-6,iz-6]-1,[24,24,24])+1
    p6(ix,iy,iz,1)=cmplx(canonical_grid(1)+100*canonical_grid(2)+10000*canonical_grid(3),0d0,8)
  enddo;enddo;enddo
  do iz=1,32;do iy=1,32;do ix=1,32
    canonical_grid=modulo(origin+[ix-10,iy-10,iz-10]-1,[24,24,24])+1
    p10(ix,iy,iz,1)=cmplx(canonical_grid(1)+100*canonical_grid(2)+10000*canonical_grid(3),0d0,8)
  enddo;enddo;enddo
  call extract_dg_wpw_canonical_cell(p6,[12,12,12],[6,6,6],[24,24,24],origin,canonical6,info)
  if(info/=0.or.abs(canonical6(13,1,13,1)-cmplx(130113d0,0d0,8))>0d0.or.&
    abs(canonical6(24,12,24,1)-cmplx(241224d0,0d0,8))>0d0)error stop 204
  call extract_dg_wpw_canonical_cell(p10,[12,12,12],[10,10,10],[24,24,24],origin,canonical10,info)
  if(info/=0.or.maxval(abs(canonical10-canonical6))>0d0)error stop 205
  p5=(0d0,0d0)
  call extract_dg_wpw_canonical_cell(p5,[12,12,12],[5,5,5],[24,24,24],origin,canonical6,info)
  if(info==0)error stop 206
  call extract_dg_wpw_canonical_cell(p6,[12,12,12],[6,6,6],[24,24,24],[-1,0,0],canonical6,info)
  if(info==0)error stop 207
  allocate(canonical_bad(23,24,24,1))
  call extract_dg_wpw_canonical_cell(p6,[12,12,12],[6,6,6],[24,24,24],origin,canonical_bad,info)
  if(info==0)error stop 208
  p6(1,1,1,1)=cmplx(ieee_value(0d0,ieee_positive_inf),0d0,8)
  call extract_dg_wpw_canonical_cell(p6,[12,12,12],[6,6,6],[24,24,24],origin,canonical6,info)
  if(info==0)error stop 209
  p6(1,1,1,1)=p10(1,1,1,1)
  p10(25,1,1,1)=p10(1,1,1,1)+50d0*epsilon(1d0)*max(1d0,abs(p10(1,1,1,1)))
  call extract_dg_wpw_canonical_cell(p10,[12,12,12],[10,10,10],[24,24,24],origin,canonical10,info)
  if(info/=0)error stop 210
  p10(25,1,1,1)=p10(1,1,1,1)+200d0*epsilon(1d0)*max(1d0,abs(p10(1,1,1,1)))
  call extract_dg_wpw_canonical_cell(p10,[12,12,12],[10,10,10],[24,24,24],origin,canonical10,info,mismatch)
  if(info==0.or..not.mismatch%valid.or.any(mismatch%canonical/=[3,15,3]).or.&
    any(mismatch%first_p/=[1,1,1]).or.any(mismatch%second_p/=[25,1,1]).or.&
    mismatch%w_row/=1.or.mismatch%abs_diff<=0d0)error stop 211
  allocate(pzero(24,24,24,0),canonical_zero(24,24,24,0))
  call extract_dg_wpw_canonical_cell(pzero,[12,12,12],[6,6,6],[24,24,24],origin,canonical_zero,info)
  if(info==0)error stop 212

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
  if(any(basis%gradient_lo/=basis%buffer_lo).or.any(basis%gradient_hi/=basis%buffer_hi))error stop 201
  saved_fingerprint=basis%fingerprint
  call evaluate_dg_wpw_occupied_w_point(basis,[3,-1,-1],&
    periodized_values,periodized_gradients,info)
  if(info/=0.or.maxval(abs(periodized_values-buffer_values(5,:)))>1d-14.or.&
    maxval(abs(periodized_gradients-buffer_gradients(:,5,:)))>1d-14)error stop 202

  buffer_values(3,1)=buffer_values(3,1)+cmplx(0.125d0,-0.25d0,8)
  call initialize_dg_wpw_occupied_w_basis(changed_basis,1,global_keys,owners,local_keys,grid_ids,&
    core_values,[-1,-1,-1],[3,-1,-1],buffer_values,buffer_gradients,2.5d0,7,info)
  if(info/=0.or.changed_basis%fingerprint==saved_fingerprint)error stop 21
  buffer_values(3,1)=buffer_values(3,1)-cmplx(0.125d0,-0.25d0,8)
  call initialize_dg_wpw_occupied_w_basis(changed_basis,1,global_keys,owners,local_keys,grid_ids+1,&
    core_values,[-1,-1,-1],[3,-1,-1],buffer_values,buffer_gradients,2.5d0,7,info)
  if(info/=0.or.changed_basis%fingerprint==saved_fingerprint)error stop 22

  call initialize_dg_wpw_occupied_w_basis(changed_basis,1,global_keys,owners,local_keys(:,[2,1]),&
    grid_ids,core_values(:,[2,1]),[-1,-1,-1],[3,-1,-1],buffer_values(:,[2,1]),&
    buffer_gradients(:,:,[2,1]),2.5d0,7,info)
  if(info/=0.or.any(changed_basis%stable_keys/=local_keys).or.&
    maxval(abs(changed_basis%core_values-core_values))>0d0.or.&
    maxval(abs(changed_basis%buffer_values-buffer_values))>0d0.or.&
    maxval(abs(changed_basis%buffer_gradients-buffer_gradients))>0d0)error stop 23

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
