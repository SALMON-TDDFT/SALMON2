program test_dg_wpw_face_trace_binding_mpi
  use mpi
  use dg_wpw_production_context,only:s_dg_wpw_production_context,initialize_dg_wpw_production_context
  use rt_dg_wpw_face_side_halo,only:s_dg_wpw_face_side_send,s_dg_wpw_face_side_state,&
    exchange_dg_wpw_face_side_schedule
  use rt_dg_wpw_face_trace_binding,only:bind_dg_wpw_canonical_face_sides
  use rt_dg_wpw_trace_halo_provider,only:s_dg_wpw_trace_halo_set
  use rt_dg_wpw_face_trace_provider,only:s_wpw_face_trace_provider,evaluate_wpw_face_traces
  implicit none
  type(s_dg_wpw_face_side_send)::send(2)
  type(s_dg_wpw_face_side_state),allocatable::remote(:)
  type(s_dg_wpw_production_context)::ctx
  type(s_dg_wpw_trace_halo_set),target::set
  type(s_wpw_face_trace_provider)::provider
  integer::ierr,rank,info,i,j,side,peer,nid
  complex(8)::wm(2),wp(2),gwm(3,2),gwp(3,2),pm(2),pp(2),gpm(3,2),gpp(3,2)
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);peer=1-rank
  do i=1,2
    side=merge(-1,1,i==1)
    nid=rank+2
    send(i)%peer=peer;send(i)%epoch=7;send(i)%local_fragment=rank+1
    send(i)%neighbor_fragment=peer+1;send(i)%axis=1;send(i)%side_from_local=side
    allocate(send(i)%w_ids(nid),send(i)%p_ids(nid))
    send(i)%w_ids=[(rank+j,j=1,nid)];send(i)%p_ids=send(i)%w_ids
    allocate(send(i)%grid(3,1),send(i)%w(nid,1),send(i)%grad_w(3,nid,1))
    allocate(send(i)%p(nid,1),send(i)%grad_p(3,nid,1))
    send(i)%grid(:,1)=[rank*2+merge(1,2,side<0),1,1]
    send(i)%w(:,1)=cmplx([(100*rank+10*i+j,j=1,nid)],0d0,8);send(i)%grad_w=0
    send(i)%p(:,1)=cmplx([(1000*rank+10*i+j,j=1,nid)],0d0,8);send(i)%grad_p=0
  enddo
  call exchange_dg_wpw_face_side_schedule(MPI_COMM_WORLD,7,send,remote,info)
  if(info/=0)error stop 1
  if(rank==0)then
    call initialize_dg_wpw_production_context(ctx,MPI_COMM_SELF,2,1,1,[1],[1],info)
    if(info/=0)error stop 2
    do i=1,2
      call bind_dg_wpw_canonical_face_sides(ctx,set,provider,send(i),remote(i),info)
      if(info/=0)error stop 3
    enddo
    if(set%nface/=2.or.set%epoch/=7)error stop 4
    do i=1,2
      call evaluate_wpw_face_traces(provider,1,2,1,send(i)%side_from_local,send(i)%grid(:,1),&
        wm,wp,gwm,gwp,pm,pp,gpm,gpp,info)
      if(info/=0.or.any(abs(wm-cmplx([10*i+1,10*i+2],0d0,8))>1d-13).or.&
         abs(wp(1))>1d-13.or.abs(wp(2)-cmplx(100+10*(3-i)+1,0d0,8))>1d-13)error stop 5
    enddo
    print '(a)','PASS canonical owner binds local/remote sides without merging periodic faces'
  endif
  call MPI_Finalize(ierr)
end program
