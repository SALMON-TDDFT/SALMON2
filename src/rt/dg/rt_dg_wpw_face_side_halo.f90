module rt_dg_wpw_face_side_halo
  use mpi
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  type,public::s_dg_wpw_face_side_send
    integer::peer=-1,epoch=-1,local_fragment=0,neighbor_fragment=0,axis=0,side_from_local=0
    integer,allocatable::w_ids(:),p_ids(:)
    integer,allocatable::grid(:,:)
    complex(8),allocatable::w(:,:),grad_w(:,:,:),p(:,:),grad_p(:,:,:)
  end type
  type,public::s_dg_wpw_face_side_state
    integer::epoch=-1,source_rank=-1,remote_fragment=0,axis=0,side_from_remote=0
    logical::valid=.false.
    integer,allocatable::w_ids(:),p_ids(:)
    integer,allocatable::grid(:,:)
    complex(8),allocatable::w(:,:),grad_w(:,:,:),p(:,:),grad_p(:,:,:)
  end type
  public::exchange_dg_wpw_face_side_schedule
  public::reduce_dg_wpw_face_side_parts
contains
  subroutine reduce_dg_wpw_face_side_parts(comm,root,peer,epoch,local_fragment,neighbor_fragment,&
      axis,side_from_local,w_ids,p_ids,grid,coverage,w,grad_w,p,grad_p,side,info)
    integer,intent(in)::comm,root,peer,epoch,local_fragment,neighbor_fragment,axis,side_from_local
    integer,intent(in)::w_ids(:),p_ids(:),grid(:,:),coverage(:)
    complex(8),intent(in)::w(:,:),grad_w(:,:,:),p(:,:),grad_p(:,:,:)
    type(s_dg_wpw_face_side_send),intent(inout)::side
    integer,intent(out)::info
    type(s_dg_wpw_face_side_send)::candidate
    integer::rank,nrank,ierr,local_bad,global_bad,npoint,root_info,metadata(6),reference_metadata(6)
    integer::local_shape(4),reference_shape(4)
    integer,allocatable::reference_w_ids(:),reference_p_ids(:),reference_grid(:,:),total_coverage(:)
    complex(8),allocatable::rw(:,:),rgw(:,:,:),rp(:,:),rgp(:,:,:)

    info=1;local_bad=0;npoint=size(grid,2)
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    metadata=[peer,epoch,local_fragment,neighbor_fragment,axis,side_from_local]
    local_shape=[size(grid,1),npoint,size(w_ids),size(p_ids)];reference_shape=local_shape
    if(root<0.or.root>=nrank.or.peer<0.or.epoch<=0.or.local_fragment<=0.or.&
       neighbor_fragment<=0.or.neighbor_fragment==local_fragment.or.axis<1.or.axis>3.or.&
       abs(side_from_local)/=1.or.npoint<=0.or.size(grid,1)/=3.or.size(coverage)/=npoint.or.&
       any(coverage<0).or.any(coverage>1).or..not.strictly_increasing(w_ids).or.&
       .not.strictly_increasing(p_ids).or.any(shape(w)/=[size(w_ids),npoint]).or.&
       any(shape(grad_w)/=[3,size(w_ids),npoint]).or.any(shape(p)/=[size(p_ids),npoint]).or.&
       any(shape(grad_p)/=[3,size(p_ids),npoint]).or..not.all(finite_complex(w)).or.&
       .not.all(finite_complex(grad_w)).or..not.all(finite_complex(p)).or.&
       .not.all(finite_complex(grad_p)))local_bad=1
    reference_metadata=metadata
    call MPI_Bcast(reference_shape,4,MPI_INTEGER,root,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(reference_shape/=local_shape))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    allocate(reference_w_ids,source=w_ids);allocate(reference_p_ids,source=p_ids)
    allocate(reference_grid,source=grid)
    call MPI_Bcast(reference_metadata,6,MPI_INTEGER,root,comm,ierr)
    call MPI_Bcast(reference_w_ids,size(w_ids),MPI_INTEGER,root,comm,ierr)
    call MPI_Bcast(reference_p_ids,size(p_ids),MPI_INTEGER,root,comm,ierr)
    call MPI_Bcast(reference_grid,size(grid),MPI_INTEGER,root,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(reference_metadata/=metadata).or.any(reference_w_ids/=w_ids).or.&
       any(reference_p_ids/=p_ids).or.any(reference_grid/=grid))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    allocate(total_coverage(npoint));allocate(rw,source=w);allocate(rgw,source=grad_w)
    allocate(rp,source=p);allocate(rgp,source=grad_p)
    call MPI_Reduce(coverage,total_coverage,npoint,MPI_INTEGER,MPI_SUM,root,comm,ierr)
    call MPI_Reduce(w,rw,size(w),MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
    call MPI_Reduce(grad_w,rgw,size(grad_w),MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
    call MPI_Reduce(p,rp,size(p),MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
    call MPI_Reduce(grad_p,rgp,size(grad_p),MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
    root_info=1
    if(rank==root.and.ierr==MPI_SUCCESS.and.all(total_coverage==1))then
      candidate%peer=peer;candidate%epoch=epoch;candidate%local_fragment=local_fragment
      candidate%neighbor_fragment=neighbor_fragment;candidate%axis=axis
      candidate%side_from_local=side_from_local;candidate%w_ids=w_ids;candidate%p_ids=p_ids
      candidate%grid=grid;candidate%w=rw;candidate%grad_w=rgw;candidate%p=rp;candidate%grad_p=rgp
      side=candidate;root_info=0
    endif
    call MPI_Bcast(root_info,1,MPI_INTEGER,root,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.root_info/=0)return
    info=0
  end subroutine reduce_dg_wpw_face_side_parts

  subroutine exchange_dg_wpw_face_side_schedule(comm,epoch,send,remote,info)
    integer,intent(in)::comm,epoch
    type(s_dg_wpw_face_side_send),intent(in)::send(:)
    type(s_dg_wpw_face_side_state),allocatable,intent(out)::remote(:)
    integer,intent(out)::info
    integer::rank,nrank,ierr,nface,i,j,nrequest,local_bad,global_bad
    integer::local_max(3),global_max(3)
    integer,allocatable::send_header(:,:),recv_header(:,:),request(:),status(:,:)

    info=1;local_bad=0;nface=size(send);local_max=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    if(epoch<=0.or.nface<=0)local_bad=1
    do i=1,nface
      if(send(i)%peer<0.or.send(i)%peer>=nrank.or.send(i)%peer==rank.or.send(i)%epoch/=epoch.or.&
         send(i)%local_fragment/=rank+1.or.send(i)%neighbor_fragment/=send(i)%peer+1.or.&
         send(i)%axis<1.or.send(i)%axis>3.or.abs(send(i)%side_from_local)/=1.or.&
         .not.allocated(send(i)%grid).or..not.allocated(send(i)%w).or.&
         .not.allocated(send(i)%grad_w).or..not.allocated(send(i)%p).or.&
         .not.allocated(send(i)%grad_p).or..not.allocated(send(i)%w_ids).or.&
         .not.allocated(send(i)%p_ids))then;local_bad=1;cycle;endif
      if(size(send(i)%grid,1)/=3.or.size(send(i)%grid,2)<=0.or.&
         size(send(i)%w,2)/=size(send(i)%grid,2).or.&
         any(shape(send(i)%grad_w)/=[3,size(send(i)%w,1),size(send(i)%grid,2)]).or.&
         size(send(i)%p,2)/=size(send(i)%grid,2).or.&
         any(shape(send(i)%grad_p)/=[3,size(send(i)%p,1),size(send(i)%grid,2)]).or.&
         size(send(i)%w_ids)/=size(send(i)%w,1).or.size(send(i)%p_ids)/=size(send(i)%p,1).or.&
         .not.strictly_increasing(send(i)%w_ids).or..not.strictly_increasing(send(i)%p_ids).or.&
         .not.all(finite_complex(send(i)%w)).or..not.all(finite_complex(send(i)%grad_w)).or.&
         .not.all(finite_complex(send(i)%p)).or..not.all(finite_complex(send(i)%grad_p)))local_bad=1
      local_max=max(local_max,[size(send(i)%grid,2),size(send(i)%w,1),size(send(i)%p,1)])
      do j=1,i-1
        if(send(j)%peer==send(i)%peer.and.send(j)%axis==send(i)%axis.and.&
           send(j)%side_from_local==send(i)%side_from_local)local_bad=1
      enddo
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;allocate(remote(0));return;endif
    call MPI_Allreduce(local_max,global_max,3,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(global_max<=0))then;allocate(remote(0));return;endif

    allocate(send_header(8,nface),recv_header(8,nface),request(2*nface),status(MPI_STATUS_SIZE,2*nface))
    nrequest=0
    do i=1,nface
      send_header(:,i)=[epoch,send(i)%local_fragment,send(i)%neighbor_fragment,send(i)%axis,&
        send(i)%side_from_local,size(send(i)%grid,2),size(send(i)%w,1),size(send(i)%p,1)]
      nrequest=nrequest+1
      call MPI_Irecv(recv_header(:,i),8,MPI_INTEGER,send(i)%peer,&
        face_tag(1200,send(i)%axis,-send(i)%side_from_local),comm,request(nrequest),ierr)
    enddo
    do i=1,nface
      nrequest=nrequest+1
      call MPI_Isend(send_header(:,i),8,MPI_INTEGER,send(i)%peer,&
        face_tag(1200,send(i)%axis,send(i)%side_from_local),comm,request(nrequest),ierr)
    enddo
    call MPI_Waitall(nrequest,request,status,ierr);local_bad=merge(0,1,ierr==MPI_SUCCESS)
    allocate(remote(nface))
    do i=1,nface
      if(recv_header(1,i)/=epoch.or.recv_header(2,i)/=send(i)%neighbor_fragment.or.&
         recv_header(3,i)/=send(i)%local_fragment.or.recv_header(4,i)/=send(i)%axis.or.&
         recv_header(5,i)/=-send(i)%side_from_local.or.any(recv_header(6:8,i)<=0).or.&
         any(recv_header(6:8,i)>global_max).or.recv_header(6,i)/=size(send(i)%grid,2))then
        local_bad=1;cycle
      endif
      remote(i)%epoch=epoch;remote(i)%source_rank=send(i)%peer
      remote(i)%remote_fragment=recv_header(2,i);remote(i)%axis=send(i)%axis
      remote(i)%side_from_remote=recv_header(5,i)
      allocate(remote(i)%grid(3,recv_header(6,i)),remote(i)%w(recv_header(7,i),recv_header(6,i)),&
        remote(i)%grad_w(3,recv_header(7,i),recv_header(6,i)))
      allocate(remote(i)%w_ids(recv_header(7,i)),remote(i)%p_ids(recv_header(8,i)))
      allocate(remote(i)%p(recv_header(8,i),recv_header(6,i)),&
        remote(i)%grad_p(3,recv_header(8,i),recv_header(6,i)))
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return

    deallocate(request,status);allocate(request(14*nface),status(MPI_STATUS_SIZE,14*nface));nrequest=0
    do i=1,nface
      call post_receives(i)
    enddo
    do i=1,nface
      call post_sends(i)
    enddo
    call MPI_Waitall(nrequest,request,status,ierr);local_bad=merge(0,1,ierr==MPI_SUCCESS)
    do i=1,nface
      if(.not.all(finite_complex(remote(i)%w)).or..not.all(finite_complex(remote(i)%grad_w)).or.&
         .not.all(finite_complex(remote(i)%p)).or..not.all(finite_complex(remote(i)%grad_p)))local_bad=1
      do j=1,3
        if(j/=send(i)%axis.and.any(remote(i)%grid(j,:)/=send(i)%grid(j,:)))local_bad=1
      enddo
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    do i=1,nface;remote(i)%valid=.true.;enddo
    info=0
  contains
    subroutine post_receives(iface)
      integer,intent(in)::iface
      integer::tag_side
      tag_side=-send(iface)%side_from_local
      nrequest=nrequest+1;call MPI_Irecv(remote(iface)%grid,size(remote(iface)%grid),MPI_INTEGER,&
        send(iface)%peer,face_tag(1300,send(iface)%axis,tag_side),comm,request(nrequest),ierr)
      nrequest=nrequest+1;call MPI_Irecv(remote(iface)%w_ids,size(remote(iface)%w_ids),MPI_INTEGER,&
        send(iface)%peer,face_tag(1350,send(iface)%axis,tag_side),comm,request(nrequest),ierr)
      nrequest=nrequest+1;call MPI_Irecv(remote(iface)%w,size(remote(iface)%w),MPI_DOUBLE_COMPLEX,&
        send(iface)%peer,face_tag(1400,send(iface)%axis,tag_side),comm,request(nrequest),ierr)
      nrequest=nrequest+1;call MPI_Irecv(remote(iface)%grad_w,size(remote(iface)%grad_w),MPI_DOUBLE_COMPLEX,&
        send(iface)%peer,face_tag(1500,send(iface)%axis,tag_side),comm,request(nrequest),ierr)
      nrequest=nrequest+1;call MPI_Irecv(remote(iface)%p,size(remote(iface)%p),MPI_DOUBLE_COMPLEX,&
        send(iface)%peer,face_tag(1600,send(iface)%axis,tag_side),comm,request(nrequest),ierr)
      nrequest=nrequest+1;call MPI_Irecv(remote(iface)%p_ids,size(remote(iface)%p_ids),MPI_INTEGER,&
        send(iface)%peer,face_tag(1650,send(iface)%axis,tag_side),comm,request(nrequest),ierr)
      nrequest=nrequest+1;call MPI_Irecv(remote(iface)%grad_p,size(remote(iface)%grad_p),MPI_DOUBLE_COMPLEX,&
        send(iface)%peer,face_tag(1700,send(iface)%axis,tag_side),comm,request(nrequest),ierr)
    end subroutine
    subroutine post_sends(iface)
      integer,intent(in)::iface
      integer::tag_side
      tag_side=send(iface)%side_from_local
      nrequest=nrequest+1;call MPI_Isend(send(iface)%grid,size(send(iface)%grid),MPI_INTEGER,&
        send(iface)%peer,face_tag(1300,send(iface)%axis,tag_side),comm,request(nrequest),ierr)
      nrequest=nrequest+1;call MPI_Isend(send(iface)%w_ids,size(send(iface)%w_ids),MPI_INTEGER,&
        send(iface)%peer,face_tag(1350,send(iface)%axis,tag_side),comm,request(nrequest),ierr)
      nrequest=nrequest+1;call MPI_Isend(send(iface)%w,size(send(iface)%w),MPI_DOUBLE_COMPLEX,&
        send(iface)%peer,face_tag(1400,send(iface)%axis,tag_side),comm,request(nrequest),ierr)
      nrequest=nrequest+1;call MPI_Isend(send(iface)%grad_w,size(send(iface)%grad_w),MPI_DOUBLE_COMPLEX,&
        send(iface)%peer,face_tag(1500,send(iface)%axis,tag_side),comm,request(nrequest),ierr)
      nrequest=nrequest+1;call MPI_Isend(send(iface)%p,size(send(iface)%p),MPI_DOUBLE_COMPLEX,&
        send(iface)%peer,face_tag(1600,send(iface)%axis,tag_side),comm,request(nrequest),ierr)
      nrequest=nrequest+1;call MPI_Isend(send(iface)%p_ids,size(send(iface)%p_ids),MPI_INTEGER,&
        send(iface)%peer,face_tag(1650,send(iface)%axis,tag_side),comm,request(nrequest),ierr)
      nrequest=nrequest+1;call MPI_Isend(send(iface)%grad_p,size(send(iface)%grad_p),MPI_DOUBLE_COMPLEX,&
        send(iface)%peer,face_tag(1700,send(iface)%axis,tag_side),comm,request(nrequest),ierr)
    end subroutine
  end subroutine exchange_dg_wpw_face_side_schedule

  pure integer function face_tag(base,axis,side)result(tag)
    integer,intent(in)::base,axis,side
    tag=base+2*(axis-1)+merge(1,0,side>0)
  end function
  elemental logical function finite_complex(value)result(ok)
    complex(8),intent(in)::value
    ok=ieee_is_finite(real(value,8)).and.ieee_is_finite(aimag(value))
  end function
  logical function strictly_increasing(ids)result(ok)
    integer,intent(in)::ids(:);integer::i
    ok=size(ids)>0
    do i=2,size(ids);if(ids(i)<=ids(i-1))then;ok=.false.;return;endif;enddo
  end function
end module
