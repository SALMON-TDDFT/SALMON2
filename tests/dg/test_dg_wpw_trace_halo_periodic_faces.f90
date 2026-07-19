program test_dg_wpw_trace_halo_periodic_faces
  use mpi
  use dg_wpw_production_context,only:s_dg_wpw_production_context,initialize_dg_wpw_production_context
  use rt_dg_wpw_trace_halo_provider,only:s_dg_wpw_trace_halo_set,prepare_dg_wpw_trace_halo_face
  use rt_dg_wpw_face_trace_provider,only:s_wpw_face_trace_provider
  use rt_dg_wpw_face_trace_scanner,only:assemble_wpw_canonical_face_grid
  implicit none
  type(s_dg_wpw_production_context)::ctx
  type(s_dg_wpw_trace_halo_set),target::set
  type(s_wpw_face_trace_provider)::provider
  complex(8)::wm(1,1),wp(1,1),gwm(3,1,1),gwp(3,1,1),pm(1,1),pp(1,1),gpm(3,1,1),gpp(3,1,1)
  complex(8)::face_h(1,1),minus_h,plus_h
  integer::ierr,info
  call MPI_Init(ierr)
  call initialize_dg_wpw_production_context(ctx,MPI_COMM_SELF,2,1,1,[1],[1],info)
  if(info/=0)error stop 1
  wm=1;wp=2;gwm=0;gwp=0;pm=1;pp=1;gpm=0;gpp=0;gpm(1,1,1)=1;gpp(1,1,1)=1
  call prepare_dg_wpw_trace_halo_face(ctx,set,provider,3,1,2,1,-1,reshape([1,1,1],[3,1]),&
    wm,wp,gwm,gwp,pm,pp,gpm,gpp,info)
  if(info/=0)error stop 2
  call prepare_dg_wpw_trace_halo_face(ctx,set,provider,3,1,2,1,1,reshape([4,1,1],[3,1]),&
    2*wm,3*wp,gwm,gwp,pm,pp,gpm,gpp,info)
  if(info/=0.or.set%nface/=2)error stop 3
  call assemble_wpw_canonical_face_grid(provider,1,2,1,-1,[1,1,1],[1,1,1],[1,1,1],[1,1,1],&
    [1d0,1d0,1d0],[1],[1],face_h,info)
  if(info/=0)error stop 4
  minus_h=face_h(1,1)
  call assemble_wpw_canonical_face_grid(provider,1,2,1,1,[4,1,1],[4,1,1],[4,1,1],[4,1,1],&
    [1d0,1d0,1d0],[1],[1],face_h,info)
  if(info/=0)error stop 5
  plus_h=face_h(1,1)
  if(abs(minus_h-plus_h)<1d-13)error stop 6
  print '(a)','PASS distinct periodic canonical-face trace records'
  call MPI_Finalize(ierr)
end program
