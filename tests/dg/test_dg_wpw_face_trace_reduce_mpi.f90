program test_dg_wpw_face_trace_reduce_mpi
  use mpi
  use dg_wpw_production_context,only:s_dg_wpw_production_context,initialize_dg_wpw_production_context
  use rt_dg_wpw_trace_halo_provider,only:s_dg_wpw_trace_halo_set,reduce_dg_wpw_face_trace_parts
  use rt_dg_wpw_face_trace_provider,only:s_wpw_face_trace_provider,evaluate_wpw_face_traces
  implicit none
  type(s_dg_wpw_production_context)::ctx
  type(s_dg_wpw_trace_halo_set),target::set
  type(s_wpw_face_trace_provider)::provider
  integer::ierr,rank,info,coverage(2),grid(3,2)
  complex(8)::wm(1,2),wp(1,2),gwm(3,1,2),gwp(3,1,2)
  complex(8)::pm(1,2),pp(1,2),gpm(3,1,2),gpp(3,1,2)
  complex(8)::rw(1),rwp(1),rgw(3,1),rgwp(3,1),rp(1),rpp(1),rgp(3,1),rgpp(3,1)

  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  if(rank==0)then
    call initialize_dg_wpw_production_context(ctx,MPI_COMM_SELF,2,1,1,[1],[1],info)
    if(info/=0)error stop 1
  endif
  grid=reshape([1,1,1,1,2,1],[3,2]);coverage=0
  wm=0;wp=0;gwm=0;gwp=0;pm=0;pp=0;gpm=0;gpp=0
  coverage(rank+1)=1
  wm(1,rank+1)=cmplx(rank+1,0d0,8);wp(1,rank+1)=cmplx(rank+3,0d0,8)
  pm(1,rank+1)=1;pp(1,rank+1)=1
  gwm(2,1,rank+1)=cmplx(10+rank,0d0,8);gwp(2,1,rank+1)=-gwm(2,1,rank+1)
  call reduce_dg_wpw_face_trace_parts(MPI_COMM_WORLD,0,ctx,set,provider,3,1,2,1,1,grid,coverage,&
    wm,wp,gwm,gwp,pm,pp,gpm,gpp,info)
  if(info/=0)error stop 2
  if(rank==0)then
    call evaluate_wpw_face_traces(provider,1,2,1,1,grid(:,1),rw,rwp,rgw,rgwp,rp,rpp,rgp,rgpp,info)
    if(info/=0.or.abs(rw(1)-(1d0,0d0))>1d-13)error stop 3
    call evaluate_wpw_face_traces(provider,1,2,1,1,grid(:,2),rw,rwp,rgw,rgwp,rp,rpp,rgp,rgpp,info)
    if(info/=0.or.abs(rw(1)-(2d0,0d0))>1d-13.or.set%nface/=1)error stop 4
  endif
  coverage=0;if(rank==0)coverage(1)=1
  call reduce_dg_wpw_face_trace_parts(MPI_COMM_WORLD,0,ctx,set,provider,4,1,2,1,1,grid,coverage,&
    wm,wp,gwm,gwp,pm,pp,gpm,gpp,info)
  if(info==0)error stop 5
  if(rank==0.and.(set%epoch/=3.or.set%nface/=1))error stop 6
  if(rank==0)print '(a)','PASS fragment-rank face traces publish only after complete coverage'
  call MPI_Finalize(ierr)
end program
