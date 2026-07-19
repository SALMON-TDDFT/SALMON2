program test_dg_wpw_lcfo_ww_import
  use,intrinsic::iso_fortran_env,only:int64
  use dg_wpw_lcfo_operator_adapter,only:s_dg_wpw_lcfo_ww_components,import_dg_wpw_lcfo_ww_components,&
    validate_dg_wpw_surface_self_policy
  implicit none
  type(s_dg_wpw_lcfo_ww_components)::ww,bad,ww_meta
  integer::info
  integer::n_basis(2),index_basis(2,2)
  real(8)::kinetic(2,2),potential(2,2),nonlocal(2,2),face_self(2,2)
  integer::face_id(2),face_row(2),face_col(2),face_axis(2),face_side(2),face_image(3,2)
  real(8)::face_value(2)

  n_basis=[2,1];index_basis=reshape([11,12,21,0],[2,2])
  call validate_dg_wpw_surface_self_policy(.true.,info);if(info/=0)error stop 9
  call validate_dg_wpw_surface_self_policy(.false.,info);if(info==0)error stop 10
  kinetic=reshape([1d0,2d0,3d0,4d0],[2,2])
  potential=10d0+kinetic;nonlocal=20d0+kinetic;face_self=30d0+kinetic
  face_id=[101,102];face_row=[11,11];face_col=[21,21];face_axis=[1,1];face_side=[-1,1]
  face_image=0;face_image(1,:)=[-1,1];face_value=[-0.25d0,-0.75d0]
  call import_dg_wpw_lcfo_ww_components(1,n_basis,index_basis,kinetic,potential,nonlocal,face_self,&
    face_id,face_row,face_col,face_axis,face_side,face_image,face_value,'orthonormal_ww',ww,info)
  if(info/=0.or..not.ww%valid)error stop 1
  if(any(ww%owned_w_ids/=[11,12]))error stop 2
  if(maxval(abs(ww%kinetic-kinetic))>1d-14.or.maxval(abs(ww%potential-potential))>1d-14)error stop 3
  if(maxval(abs(ww%nonlocal-nonlocal))>1d-14.or.maxval(abs(ww%face_self-face_self))>1d-14)error stop 4
  if(any(ww%cross_face_id/=face_id).or.any(ww%cross_side/=[-1,1]))error stop 5
  if(any(ww%cross_image(1,:)/=[-1,1]).or.any(ww%cross_value/=face_value))error stop 6
  if(trim(ww%metric_convention)/='orthonormal_ww'.or.ww%provenance_fingerprint==0_int64)error stop 7
  call import_dg_wpw_lcfo_ww_components(1,n_basis,index_basis,kinetic,potential,nonlocal,face_self,&
    face_id,[12,11],face_col,face_axis,face_side,face_image,face_value,'orthonormal_ww',ww_meta,info)
  if(info/=0.or.ww_meta%provenance_fingerprint==ww%provenance_fingerprint)error stop 71
  call import_dg_wpw_lcfo_ww_components(1,n_basis,index_basis,kinetic,potential,nonlocal,face_self,&
    [101,101],face_row,face_col,face_axis,face_side,face_image,face_value,'orthonormal_ww',bad,info)
  if(info==0.or.bad%valid)error stop 8
  print '(a)','PASS LCFO WW components preserve canonical periodic faces and provenance'
end program
