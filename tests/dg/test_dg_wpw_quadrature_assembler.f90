program test_dg_wpw_quadrature_assembler
  use rt_dg_wpw_quadrature_assembler, only: assemble_wpw_volume_point, assemble_wpw_canonical_face_point, &
    pack_wpw_point_candidates
  use rt_dg_wpw_sparse_builder, only: s_dg_wpw_sparse_candidates, wpw_candidate_volume_face, &
    build_windowed_sparse_wpw_operators
  use rt_dg_wpw_column_layout, only: s_dg_wpw_column_layout, initialize_wpw_column_layout
  use rt_dg_wpw_sparse_blocks, only: s_dg_wpw_sparse_blocks
  implicit none
  complex(8) :: w(1), gw(3,1), po(1), gpo(3,1), ps(1), gps(3,1)
  complex(8) :: wp_h(1,1), wp_s(1,1), pp_h(1,1), pp_s(1,1)
  complex(8) :: wm(1), wx(1), gwm(3,1), gwx(3,1), pm(1), px(1), gpm(3,1), gpx(3,1)
  complex(8) :: face_h(1,1), reverse_h(1,1), zero_h(1,1)
  integer :: info
  type(s_dg_wpw_sparse_candidates) :: candidates
  type(s_dg_wpw_column_layout) :: layout
  type(s_dg_wpw_sparse_blocks) :: h_blocks,s_blocks
  real(8), parameter :: tol=1d-13

  w(1)=(1d0,0.2d0); gw(:,1)=[(0.3d0,0.1d0),(0d0,-0.2d0),(0.4d0,0d0)]
  po(1)=(0.7d0,-0.3d0); gpo(:,1)=[(0d0,0.1d0),(0.2d0,0d0),(0d0,-0.3d0)]
  ps=po; gps=gpo
  call assemble_wpw_volume_point(w,gw,po,gpo,ps,gps,1.4d0,0.25d0,wp_h,wp_s,pp_h,pp_s,info)
  if(info/=0) error stop 1
  if(abs(wp_s(1,1)-0.25d0*conjg(w(1))*po(1))>tol) error stop 2
  if(abs(pp_s(1,1)-0.25d0*conjg(po(1))*ps(1))>tol) error stop 3

  wm(1)=(1d0,0.1d0); wx(1)=(0.2d0,-0.3d0)
  gwm(:,1)=[(0.4d0,0d0),(0d0,0.1d0),(0d0,0d0)]
  gwx(:,1)=[(-0.2d0,0d0),(0d0,0.1d0),(0d0,0d0)]
  pm(1)=(0.8d0,0.2d0); px=pm
  gpm(:,1)=[(0d0,0.3d0),(0.2d0,0d0),(0d0,0d0)]; gpx=gpm
  call assemble_wpw_canonical_face_point(1,2,wm,wx,gwm,gwx,pm,px,gpm,gpx,[1d0,0d0,0d0],0.5d0,0.4d0,face_h,info)
  if(info/=0) error stop 4
  call assemble_wpw_canonical_face_point(1,2,wx,wm,gwx,gwm,px,pm,gpx,gpm,[-1d0,0d0,0d0],0.5d0,0.4d0,reverse_h,info)
  if(info/=0 .or. maxval(abs(face_h-reverse_h))>tol) error stop 5
  call assemble_wpw_canonical_face_point(1,2,wm,wm,gwm,gwm,pm,px,gpm,gpx,[1d0,0d0,0d0],0.5d0,0.4d0,zero_h,info)
  if(info/=0 .or. maxval(abs(zero_h))>0d0) error stop 6
  px(1)=px(1)+(1d-5,0d0)
  call assemble_wpw_canonical_face_point(1,2,wm,wx,gwm,gwx,pm,px,gpm,gpx,[1d0,0d0,0d0],0.5d0,0.4d0,face_h,info)
  if(info/=4) error stop 7
  px=pm
  call pack_wpw_point_candidates([11],[1],[1],wp_h,wp_s,pp_h,pp_s,face_h,candidates,info)
  if(info/=0) error stop 8
  if(size(candidates%wp_h_values,1)/=1 .or. abs(candidates%wp_h_values(1,1)-wp_h(1,1)-face_h(1,1))>tol) error stop 9
  if(candidates%wp_origin(1)/=wpw_candidate_volume_face .or. candidates%pp_origin(1)/=1) error stop 10
  call initialize_wpw_column_layout(layout,1,1,0,1,info)
  if(info/=0) error stop 11
  call build_windowed_sparse_wpw_operators(layout,0,1,candidates,h_blocks,s_blocks,info)
  if(info/=0 .or. size(h_blocks%wp_values,1)/=1 .or. size(h_blocks%pp_values,1)/=1) error stop 12
  write(*,'(a)') 'PASS one-point support-local WP/PP quadrature assembler'
end program test_dg_wpw_quadrature_assembler
