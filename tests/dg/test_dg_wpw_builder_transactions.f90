program test_dg_wpw_builder_transactions
  use mpi
  use dg_wpw_production_context,only:s_dg_wpw_production_context,initialize_dg_wpw_production_context
  use rt_dg_wpw_production_builder,only:build_dg_wpw_rank_local_quadrature,&
    scan_dg_wpw_canonical_faces,build_dg_wpw_production_operator,bind_dg_wpw_hs_callbacks,&
    install_dg_wpw_nonlocal_volume,replace_dg_wpw_potential_volume
  use rt_dg_wpw_trace_halo_provider,only:s_dg_wpw_trace_halo_state,prepare_dg_wpw_trace_halo
  use rt_dg_wpw_face_trace_provider,only:s_wpw_face_trace_provider
  implicit none
  type(s_dg_wpw_production_context)::ctx
  type(s_dg_wpw_trace_halo_state),target::halo
  type(s_wpw_face_trace_provider)::provider
  complex(8)::wh(1),ws(1),ph(1),ps(1),saved(1),x(3),y0(3),y1(3),wp_nl(1),pp_nl(1)
  complex(8)::wm(1,1),wp(1,1),gwm(3,1,1),gwp(3,1,1)
  complex(8)::pm(2,1),pp(2,1),gpm(3,2,1),gpp(3,2,1)
  integer::info,ierr
  call MPI_Init(ierr)
  call initialize_dg_wpw_production_context(ctx,MPI_COMM_SELF,2,1,1,[1],[1],info)
  wh=(1d0,0d0);ws=(0.1d0,0d0);ph=(2d0,0d0);ps=(1d0,0d0)
  call build_dg_wpw_rank_local_quadrature(ctx,[1],[1],wh,ws,[1],[1],ph,ps,info)
  if(info/=0)error stop 1
  wp_nl=(7d0,0d0);pp_nl=(11d0,0d0)
  call install_dg_wpw_nonlocal_volume(ctx,wp_nl,pp_nl,info)
  if(info/=0.or.abs(ctx%wp_h(1)-8d0)>1d-13.or.abs(ctx%pp_h(1)-13d0)>1d-13)error stop 6
  call replace_dg_wpw_potential_volume(ctx,[(3d0,0d0)],[(5d0,0d0)],info)
  if(info/=0.or.abs(ctx%wp_h(1)-10d0)>1d-13.or.abs(ctx%pp_h(1)-16d0)>1d-13)error stop 7
  ctx%halo_epoch=1;ctx%scan_epoch=1;ctx%face_valid=.true.
  call build_dg_wpw_production_operator(ctx,info);call bind_dg_wpw_hs_callbacks(ctx,info)
  x=[(1d0,0d0),(2d0,0d0),(3d0,0d0)];call ctx%apply_h(x,y0,info)
  call build_dg_wpw_rank_local_quadrature(ctx,[1],[1,1],wh,ws,[1],[1],ph,ps,info)
  if(info==0.or..not.ctx%callbacks_bound)error stop 2
  call ctx%apply_h(x,y1,info);if(info/=0.or.any(abs(y1-y0)>1d-13))error stop 3
  wh=(4d0,0d0)
  call build_dg_wpw_rank_local_quadrature(ctx,[1],[1],wh,ws,[1],[1],ph,ps,info)
  if(info/=0.or..not.ctx%quadrature_valid.or.ctx%callbacks_bound.or.abs(ctx%wp_h(1)-wh(1))>1d-13)error stop 4
  saved=ctx%wp_h
  wm=1;wp=2;gwm=0;gwp=0;pm=1;pp=1;gpm=0;gpp=0;gpm(1,:,:)=1;gpp(1,:,:)=1
  call prepare_dg_wpw_trace_halo(ctx,halo,provider,2,1,2,1,1,reshape([1,1,1],[3,1]),&
    wm,wp,gwm,gwp,pm,pp,gpm,gpp,info)
  call scan_dg_wpw_canonical_faces(ctx,2,provider,1,2,1,1,[1,1,1],[1,1,1],[1,1,1],[1,1,1],&
    [1d0,1d0,1d0],[1],[1,2],info)
  if(info==0.or.any(abs(ctx%wp_h-saved)>1d-13))error stop 5
  print '(a)','PASS production builder updates are transactional'
  call MPI_Finalize(ierr)
end program
