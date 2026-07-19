module rt_dg_wpw_face_trace_binding
  use dg_wpw_production_context,only:s_dg_wpw_production_context
  use rt_dg_wpw_face_side_halo,only:s_dg_wpw_face_side_send,s_dg_wpw_face_side_state
  use rt_dg_wpw_trace_halo_provider,only:s_dg_wpw_trace_halo_set,prepare_dg_wpw_trace_halo_face
  use rt_dg_wpw_face_trace_provider,only:s_wpw_face_trace_provider
  implicit none
  private
  public::bind_dg_wpw_canonical_face_sides
contains
  subroutine bind_dg_wpw_canonical_face_sides(ctx,set,provider,local_side,remote_side,info)
    type(s_dg_wpw_production_context),intent(inout)::ctx
    type(s_dg_wpw_trace_halo_set),target,intent(inout)::set
    type(s_wpw_face_trace_provider),intent(inout)::provider
    type(s_dg_wpw_face_side_send),intent(in)::local_side
    type(s_dg_wpw_face_side_state),intent(in)::remote_side
    integer,intent(out)::info
    integer::tangent(2),i,match(1)
    complex(8),allocatable::remote_w(:,:),remote_grad_w(:,:,:),remote_p(:,:),remote_grad_p(:,:,:)

    info=1
    if(local_side%local_fragment<=0.or.local_side%local_fragment>=local_side%neighbor_fragment.or.&
       local_side%epoch<=0.or.local_side%epoch/=remote_side%epoch.or..not.remote_side%valid.or.&
       local_side%peer/=remote_side%source_rank.or.&
       local_side%neighbor_fragment/=remote_side%remote_fragment.or.&
       local_side%axis/=remote_side%axis.or.&
       remote_side%side_from_remote/=-local_side%side_from_local.or.&
       .not.allocated(local_side%grid).or..not.allocated(remote_side%grid).or.&
       size(local_side%grid,1)/=3.or.any(shape(remote_side%grid)/=shape(local_side%grid)).or.&
       size(remote_side%w,1)/=size(remote_side%w_ids).or.&
       any(shape(remote_side%grad_w)/=[3,size(remote_side%w_ids),size(remote_side%grid,2)]).or.&
       size(remote_side%p,1)/=size(remote_side%p_ids).or.&
       any(shape(remote_side%grad_p)/=[3,size(remote_side%p_ids),size(remote_side%grid,2)]))return
    select case(local_side%axis)
    case(1);tangent=[2,3]
    case(2);tangent=[1,3]
    case(3);tangent=[1,2]
    case default;return
    end select
    if(any(remote_side%grid(tangent,:)/=local_side%grid(tangent,:)))return
    allocate(remote_w(size(local_side%w_ids),size(local_side%grid,2)),&
      remote_grad_w(3,size(local_side%w_ids),size(local_side%grid,2)))
    allocate(remote_p(size(local_side%p_ids),size(local_side%grid,2)),&
      remote_grad_p(3,size(local_side%p_ids),size(local_side%grid,2)))
    remote_w=0;remote_grad_w=0;remote_p=0;remote_grad_p=0
    do i=1,size(local_side%w_ids)
      match=findloc(remote_side%w_ids,local_side%w_ids(i))
      if(match(1)>0)then
        remote_w(i,:)=remote_side%w(match(1),:)
        remote_grad_w(:,i,:)=remote_side%grad_w(:,match(1),:)
      endif
    enddo
    do i=1,size(local_side%p_ids)
      match=findloc(remote_side%p_ids,local_side%p_ids(i))
      if(match(1)>0)then
        remote_p(i,:)=remote_side%p(match(1),:)
        remote_grad_p(:,i,:)=remote_side%grad_p(:,match(1),:)
      endif
    enddo
    call prepare_dg_wpw_trace_halo_face(ctx,set,provider,local_side%epoch,&
      local_side%local_fragment,local_side%neighbor_fragment,local_side%axis,&
      local_side%side_from_local,local_side%grid,local_side%w,remote_w,&
      local_side%grad_w,remote_grad_w,local_side%p,remote_p,&
      local_side%grad_p,remote_grad_p,info)
  end subroutine
end module
