module dg_wpw_owner_exchange
  use mpi,only:MPI_Comm_rank,MPI_Comm_size,MPI_Allreduce,MPI_Irecv,MPI_Isend,MPI_Waitall,&
    MPI_Dist_graph_create,MPI_Dist_graph_neighbors_count,MPI_Dist_graph_neighbors,&
    MPI_INTEGER,MPI_MAX,MPI_MIN,MPI_DOUBLE_COMPLEX,MPI_STATUSES_IGNORE,MPI_SUCCESS,&
    MPI_UNWEIGHTED,MPI_INFO_NULL,MPI_COMM_NULL,MPI_Comm_free
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use,intrinsic::iso_fortran_env,only:int64
  implicit none
  private
  type,public::s_dg_wpw_owner_schedule
    integer::comm=MPI_COMM_NULL,rank=-1,nrank=0
    integer,allocatable::peers(:),owned_ids(:),required_ids(:),required_owner(:),required_position(:)
    integer,allocatable::remote_counts(:),remote_offsets(:),remote_ids(:)
    logical::valid=.false.
  end type
  public::initialize_dg_wpw_owner_schedule,fetch_rows_from_owners,reduce_w_partial_to_owners
  public::release_dg_wpw_owner_schedule
  public::peer_sets_equal
  public::exchange_sizes_fit
  public::collective_exchange_sizes_fit
contains
  subroutine initialize_dg_wpw_owner_schedule(schedule,comm,peers,owned_ids,required_ids,info)
    type(s_dg_wpw_owner_schedule),intent(out)::schedule
    integer,intent(in)::comm,peers(:),owned_ids(:),required_ids(:)
    integer,intent(out)::info
    integer::rank,nrank,ierr,i,j,npeer,nrequest,local_bad,global_bad,total,nmatch,before
    integer::graph_comm,indegree,outdegree
    real(8)::stage_clock
    integer::graph_sources(1),graph_degrees(1)
    logical::weighted
    integer,allocatable::request(:),send_counts(:),sources(:),destinations(:),source_weights(:),dest_weights(:)
    info=1;local_bad=0;npeer=size(peers);call cpu_time(stage_clock)
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    if(.not.strictly_increasing(peers).or..not.strictly_increasing(owned_ids).or.&
       .not.strictly_increasing(required_ids).or.any(peers<0).or.any(peers>=nrank).or.any(peers==rank).or.&
       any(owned_ids<=0).or.any(required_ids<=0))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call trace_owner_schedule('graph',rank,stage_clock,size(owned_ids),size(required_ids))
    graph_sources=[rank];graph_degrees=[npeer]
    call MPI_Dist_graph_create(comm,1,graph_sources,graph_degrees,peers,MPI_UNWEIGHTED,&
      MPI_INFO_NULL,.false.,graph_comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    schedule%comm=graph_comm;schedule%rank=rank;schedule%nrank=nrank;schedule%peers=peers
    call MPI_Dist_graph_neighbors_count(graph_comm,indegree,outdegree,weighted,ierr)
    if(ierr/=MPI_SUCCESS.or.indegree/=npeer.or.outdegree/=npeer)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call release_dg_wpw_owner_schedule(schedule);return;endif
    allocate(sources(max(0,indegree)),destinations(max(0,outdegree)),source_weights(max(0,indegree)),&
      dest_weights(max(0,outdegree)))
    call MPI_Dist_graph_neighbors(graph_comm,indegree,sources,source_weights,outdegree,destinations,dest_weights,ierr)
    if(ierr/=MPI_SUCCESS)then
      local_bad=1
    else if(indegree==npeer.and.outdegree==npeer)then
      if(.not.peer_sets_equal(sources,peers).or..not.peer_sets_equal(destinations,peers))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call release_dg_wpw_owner_schedule(schedule);return;endif
    schedule%owned_ids=owned_ids;schedule%required_ids=required_ids
    call trace_owner_schedule('id_exchange',rank,stage_clock,size(owned_ids),size(required_ids))
    allocate(schedule%remote_counts(npeer),send_counts(npeer),request(max(1,2*npeer)))
    send_counts=size(owned_ids);schedule%remote_counts=-1;nrequest=0
    do i=1,npeer
      nrequest=nrequest+1;call MPI_Irecv(schedule%remote_counts(i),1,MPI_INTEGER,peers(i),7301,comm,request(nrequest),ierr)
      if(ierr/=MPI_SUCCESS)local_bad=1
    enddo
    do i=1,npeer
      nrequest=nrequest+1;call MPI_Isend(send_counts(i),1,MPI_INTEGER,peers(i),7301,comm,request(nrequest),ierr)
      if(ierr/=MPI_SUCCESS)local_bad=1
    enddo
    if(nrequest>0)call MPI_Waitall(nrequest,request,MPI_STATUSES_IGNORE,ierr)
    if(ierr/=MPI_SUCCESS.or.any(schedule%remote_counts<0))local_bad=1
    if(.not.exchange_sizes_fit(schedule%remote_counts,send_counts,1))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call release_dg_wpw_owner_schedule(schedule);return;endif
    allocate(schedule%remote_offsets(npeer+1));schedule%remote_offsets(1)=1
    do i=1,npeer;schedule%remote_offsets(i+1)=schedule%remote_offsets(i)+schedule%remote_counts(i);enddo
    total=schedule%remote_offsets(npeer+1)-1;allocate(schedule%remote_ids(max(0,total)))
    deallocate(request);allocate(request(max(1,2*npeer)));nrequest=0
    do i=1,npeer
      nrequest=nrequest+1;call MPI_Irecv(schedule%remote_ids(schedule%remote_offsets(i)),schedule%remote_counts(i),&
        MPI_INTEGER,peers(i),7302,comm,request(nrequest),ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    enddo
    do i=1,npeer
      nrequest=nrequest+1;call MPI_Isend(owned_ids,size(owned_ids),MPI_INTEGER,peers(i),7302,comm,request(nrequest),ierr)
      if(ierr/=MPI_SUCCESS)local_bad=1
    enddo
    if(nrequest>0)call MPI_Waitall(nrequest,request,MPI_STATUSES_IGNORE,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    allocate(schedule%required_owner(size(required_ids)),schedule%required_position(size(required_ids)))
    schedule%required_owner=-2;schedule%required_position=0
    call trace_owner_schedule('owner_resolution',rank,stage_clock,size(owned_ids),size(required_ids))
    do j=1,size(required_ids)
      nmatch=0
      do i=1,size(owned_ids)
        if(owned_ids(i)==required_ids(j))then
          schedule%required_owner(j)=-1;schedule%required_position(j)=i;nmatch=nmatch+1
        endif
      enddo
      do i=1,npeer
        before=nmatch
        call locate(schedule%remote_ids(schedule%remote_offsets(i):schedule%remote_offsets(i+1)-1),&
          required_ids(j),schedule%required_position(j),nmatch)
        if(nmatch>before)schedule%required_owner(j)=i
      enddo
      if(nmatch/=1)local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call release_dg_wpw_owner_schedule(schedule);return;endif
    schedule%valid=.true.;info=0
  end subroutine

  subroutine trace_owner_schedule(stage,rank,stage_clock,nowned,nrequired)
    character(*),intent(in)::stage
    integer,intent(in)::rank,nowned,nrequired
    real(8),intent(inout)::stage_clock
    real(8)::now
    call cpu_time(now)
    write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,f10.3)')'[DG-WPW-OWNER-STAGE] stage=',trim(stage),&
      ' rank=',rank,' owned=',nowned,' required=',nrequired,' cpu_seconds=',now-stage_clock
    flush(6);stage_clock=now
  end subroutine

  subroutine release_dg_wpw_owner_schedule(schedule)
    type(s_dg_wpw_owner_schedule),intent(inout)::schedule
    integer::ierr
    if(schedule%comm/=MPI_COMM_NULL)call MPI_Comm_free(schedule%comm,ierr)
    schedule=s_dg_wpw_owner_schedule()
  end subroutine

  subroutine fetch_rows_from_owners(schedule,owned_x,fetched_x,info)
    type(s_dg_wpw_owner_schedule),intent(in)::schedule
    complex(8),intent(in)::owned_x(:,:)
    complex(8),intent(out)::fetched_x(:,:)
    integer,intent(out)::info
    integer::i,j,nvec,nvec_min,nvec_max,npeer,nrequest,ierr,local_bad,global_bad,total,base
    integer,allocatable::request(:)
    complex(8),allocatable::remote_x(:)
    fetched_x=0;info=1;local_bad=0;nvec=size(owned_x,2);npeer=size(schedule%peers)
    if(.not.schedule%valid.or.size(owned_x,1)/=size(schedule%owned_ids).or.&
       any(shape(fetched_x)/=[size(schedule%required_ids),nvec]).or..not.finite(owned_x))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call MPI_Allreduce(nvec,nvec_min,1,MPI_INTEGER,MPI_MIN,schedule%comm,ierr)
    call MPI_Allreduce(nvec,nvec_max,1,MPI_INTEGER,MPI_MAX,schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.nvec_min/=nvec_max)return
    call collective_exchange_sizes_fit(schedule,nvec,ierr);if(ierr/=0)return
    total=sum(schedule%remote_counts)*nvec;allocate(remote_x(max(1,total)),request(max(1,2*npeer)))
    nrequest=0;base=1
    do i=1,npeer
      nrequest=nrequest+1;call MPI_Irecv(remote_x(base),schedule%remote_counts(i)*nvec,MPI_DOUBLE_COMPLEX,&
        schedule%peers(i),7311,schedule%comm,request(nrequest),ierr);base=base+schedule%remote_counts(i)*nvec
    enddo
    do i=1,npeer
      nrequest=nrequest+1;call MPI_Isend(owned_x,size(owned_x),MPI_DOUBLE_COMPLEX,schedule%peers(i),7311,&
        schedule%comm,request(nrequest),ierr)
    enddo
    if(nrequest>0)call MPI_Waitall(nrequest,request,MPI_STATUSES_IGNORE,ierr)
    if(ierr/=MPI_SUCCESS)return
    do j=1,size(schedule%required_ids)
      if(schedule%required_owner(j)==-1)then
        fetched_x(j,:)=owned_x(schedule%required_position(j),:)
      else
        i=schedule%required_owner(j);base=sum(schedule%remote_counts(1:i-1))*nvec
        do nrequest=1,nvec
          fetched_x(j,nrequest)=remote_x(base+(nrequest-1)*schedule%remote_counts(i)+schedule%required_position(j))
        enddo
      endif
    enddo
    info=0
  end subroutine

  subroutine reduce_w_partial_to_owners(schedule,partial,owned_sum,info)
    type(s_dg_wpw_owner_schedule),intent(in)::schedule
    complex(8),intent(in)::partial(:,:)
    complex(8),intent(out)::owned_sum(:,:)
    integer,intent(out)::info
    integer::i,j,k,nvec,nvec_min,nvec_max,npeer,nrequest,ierr,local_bad,global_bad,send_total,recv_total,base
    integer,allocatable::request(:),send_offsets(:)
    complex(8),allocatable::sendbuf(:),recvbuf(:)
    owned_sum=0;info=1;local_bad=0;nvec=size(partial,2);npeer=size(schedule%peers)
    if(.not.schedule%valid.or.any(shape(partial)/=[size(schedule%required_ids),nvec]).or.&
       any(shape(owned_sum)/=[size(schedule%owned_ids),nvec]).or..not.finite(partial))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call MPI_Allreduce(nvec,nvec_min,1,MPI_INTEGER,MPI_MIN,schedule%comm,ierr)
    call MPI_Allreduce(nvec,nvec_max,1,MPI_INTEGER,MPI_MAX,schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.nvec_min/=nvec_max)return
    call collective_exchange_sizes_fit(schedule,nvec,ierr);if(ierr/=0)return
    allocate(send_offsets(npeer+1));send_offsets(1)=1
    do i=1,npeer;send_offsets(i+1)=send_offsets(i)+schedule%remote_counts(i)*nvec;enddo
    send_total=send_offsets(npeer+1)-1;recv_total=npeer*size(schedule%owned_ids)*nvec
    allocate(sendbuf(max(1,send_total)),recvbuf(max(1,recv_total)),request(max(1,2*npeer)));sendbuf=0
    do j=1,size(schedule%required_ids)
      i=schedule%required_owner(j)
      if(i==-1)then
        owned_sum(schedule%required_position(j),:)=owned_sum(schedule%required_position(j),:)+partial(j,:)
      else
        do k=1,nvec
          base=send_offsets(i)+(k-1)*schedule%remote_counts(i)+schedule%required_position(j)-1
          sendbuf(base)=sendbuf(base)+partial(j,k)
        enddo
      endif
    enddo
    nrequest=0;base=1
    do i=1,npeer
      nrequest=nrequest+1;call MPI_Irecv(recvbuf(base),size(schedule%owned_ids)*nvec,MPI_DOUBLE_COMPLEX,&
        schedule%peers(i),7312,schedule%comm,request(nrequest),ierr);base=base+size(schedule%owned_ids)*nvec
    enddo
    do i=1,npeer
      nrequest=nrequest+1;call MPI_Isend(sendbuf(send_offsets(i)),schedule%remote_counts(i)*nvec,&
        MPI_DOUBLE_COMPLEX,schedule%peers(i),7312,schedule%comm,request(nrequest),ierr)
    enddo
    if(nrequest>0)call MPI_Waitall(nrequest,request,MPI_STATUSES_IGNORE,ierr)
    if(ierr/=MPI_SUCCESS)return
    base=1
    do i=1,npeer
      do k=1,nvec
        owned_sum(:,k)=owned_sum(:,k)+recvbuf(base+(k-1)*size(schedule%owned_ids):base+k*size(schedule%owned_ids)-1)
      enddo
      base=base+size(schedule%owned_ids)*nvec
    enddo
    info=0
  end subroutine

  subroutine locate(ids,target,position,nmatch)
    integer,intent(in)::ids(:),target;integer,intent(inout)::position,nmatch
    integer::i
    do i=1,size(ids);if(ids(i)==target)then;position=i;nmatch=nmatch+1;endif;enddo
  end subroutine
  logical function strictly_increasing(ids)result(ok)
    integer,intent(in)::ids(:);integer::i;ok=.true.
    do i=2,size(ids);if(ids(i)<=ids(i-1))then;ok=.false.;return;endif;enddo
  end function
  logical function peer_sets_equal(left,right)result(equal)
    integer,intent(in)::left(:),right(:)
    integer::i
    equal=size(left)==size(right)
    if(.not.equal)return
    do i=1,size(left)
      if(count(right==left(i))/=1.or.count(left==right(i))/=1)then
        equal=.false.;return
      endif
    enddo
  end function
  logical function exchange_sizes_fit(send_counts,recv_counts,nvec)result(fit)
    integer,intent(in)::send_counts(:),recv_counts(:),nvec
    integer(int64)::limit,total_send,total_recv
    fit=.false.;limit=int(huge(1),int64)
    if(nvec<0.or.any(send_counts<0).or.any(recv_counts<0))return
    total_send=sum(int(send_counts,int64));total_recv=sum(int(recv_counts,int64))
    if(total_send>limit.or.total_recv>limit)return
    if(int(nvec,int64)>0)then
      if(total_send>limit/int(nvec,int64).or.total_recv>limit/int(nvec,int64))return
      if(any(int(send_counts,int64)>limit/int(nvec,int64)))return
      if(any(int(recv_counts,int64)>limit/int(nvec,int64)))return
    endif
    fit=.true.
  end function
  subroutine collective_exchange_sizes_fit(schedule,nvec,info)
    type(s_dg_wpw_owner_schedule),intent(in)::schedule
    integer,intent(in)::nvec
    integer,intent(out)::info
    integer::i,ierr,local_bad,global_bad,npeer
    npeer=size(schedule%peers);local_bad=0
    if(.not.exchange_sizes_fit(schedule%remote_counts,[(size(schedule%owned_ids),i=1,npeer)],nvec))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,schedule%comm,ierr)
    info=merge(0,1,ierr==MPI_SUCCESS.and.global_bad==0)
  end subroutine
  logical function finite(x)result(ok)
    complex(8),intent(in)::x(:,:);ok=all(ieee_is_finite(real(x,8))).and.all(ieee_is_finite(aimag(x)))
  end function
end module
