program test_dg_wpw_grid_point_adapter
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use rt_dg_wpw_grid_point_adapter, only: evaluate_local_wannier_grid_point, evaluate_windowed_kg_grid_point
  use rt_dg_wpw_canonical_faces, only: s_wpw_canonical_face_list, build_wpw_canonical_face_list
  implicit none
  type(s_dg_fragment_rt)::dg
  type(s_dg_fragment_rt)::dgface
  type(s_wpw_canonical_face_list)::faces
  complex(8)::w(1),gw(3,1),p,gp(3)
  integer::info
  real(8),parameter::tol=1d-13

  dg%n_frag=1; dg%nspin=1; dg%ifrag_start=1; dg%ifrag_end=1
  dg%num_fragment=[1,1,1]; dg%lgnum_total=[1,1,1]; dg%hgs=[1d0,1d0,1d0]
  allocate(dg%n_basis(1,1),dg%ixyz_frag(3,1),dg%nxyz_domain(3,1))
  dg%n_basis(1,1)=2; dg%ixyz_frag(:,1)=[1,1,1]; dg%nxyz_domain(:,1)=[1,1,1]
  dg%has_global_wannier_local_basis=.true.; dg%gradient_basis_cache_valid=.true.
  allocate(dg%global_wannier_local_nkeep(1),dg%global_wannier_local_ids(1,1))
  allocate(dg%global_wannier_local_coef(2,1,1,1),dg%phi_frag(1,1,1,2,1))
  allocate(dg%gradient_basis_cache(1,1,1,3,2,1))
  dg%global_wannier_local_nkeep=1; dg%global_wannier_local_ids(1,1)=5
  dg%global_wannier_local_coef(:,1,1,1)=[(1d0,0d0),(0d0,1d0)]
  dg%phi_frag(1,1,1,:,1)=[2d0,3d0]
  dg%gradient_basis_cache=0d0
  dg%gradient_basis_cache(1,1,1,1,1,1)=1d0
  dg%gradient_basis_cache(1,1,1,2,2,1)=2d0
  call evaluate_local_wannier_grid_point(dg,1,1,[5],[1,1,1],w,gw,info)
  if(info/=0 .or. abs(w(1)-(2d0,3d0))>tol .or. abs(gw(1,1)-(1d0,0d0))>tol .or. &
     abs(gw(2,1)-(0d0,2d0))>tol) error stop 1

  dg%has_wpw_window=.true.; dg%n_plane_waves=1
  allocate(dg%k_pw(3,1)); dg%k_pw=0d0
  allocate(dg%wpw_window_box_lo(3,1),dg%wpw_window_box_hi(3,1))
  allocate(dg%wpw_chi(1,1,1,1),dg%wpw_grad_chi(3,1,1,1,1))
  dg%wpw_window_box_lo(:,1)=[1,1,1]; dg%wpw_window_box_hi(:,1)=[1,1,1]
  dg%wpw_chi=1d0; dg%wpw_grad_chi=0d0
  call evaluate_windowed_kg_grid_point(dg,1,1,[1,1,1],4d0,p,gp,info)
  if(info/=0 .or. abs(p-(0.5d0,0d0))>tol .or. maxval(abs(gp))>tol) error stop 2
  dgface%n_frag=2; dgface%ifrag_start=1; dgface%ifrag_end=1
  dgface%num_fragment=[2,1,1]; dgface%lgnum_total=[2,1,1]
  allocate(dgface%ixyz_frag(3,2),dgface%nxyz_domain(3,2))
  dgface%ixyz_frag(:,1)=[1,1,1]; dgface%ixyz_frag(:,2)=[2,1,1]
  dgface%nxyz_domain(:,1)=[1,1,1]; dgface%nxyz_domain(:,2)=[1,1,1]
  call build_wpw_canonical_face_list(dgface,faces,info)
  if(info/=0 .or. faces%nface/=2) error stop 3
  if(any(faces%k_minus/=[1,1]) .or. any(faces%k_plus/=[2,2]) .or. &
     any(faces%axis/=[1,1]) .or. any(faces%side_from_k_minus/=[-1,1])) error stop 4
  write(*,'(a)') 'PASS real-grid W/(K,G) point adapter numerical fixture'
end program test_dg_wpw_grid_point_adapter
