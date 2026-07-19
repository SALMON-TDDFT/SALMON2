module dg_wpw_nonlocal_projector
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use mpi,only:MPI_Comm_rank,MPI_Comm_size,MPI_Allreduce,MPI_Irecv,MPI_Isend,MPI_Waitall,&
    MPI_INTEGER,MPI_INTEGER8,MPI_MAX,MPI_MIN,MPI_DOUBLE_PRECISION,MPI_DOUBLE_COMPLEX,MPI_STATUSES_IGNORE,&
    MPI_SUCCESS,MPI_Reduce,MPI_SUM,MPI_Gather,MPI_Gatherv,MPI_Bcast
  implicit none
  private
  type,public::s_dg_wpw_projector_overlap
    integer::basis_id=0
    integer::projector_id=0
    complex(8)::value=(0d0,0d0)
  end type
  public::assemble_dg_wpw_nonlocal_blocks
  public::exchange_dg_wpw_projector_overlaps
  public::reduce_dg_wpw_p_projector_partials
  public::build_dg_wpw_projector_overlap_partials
  public::reduce_dg_wpw_projector_dense,validate_dg_wpw_projector_provenance
  public::reduce_dg_wpw_projector_records
  public::validate_dg_wpw_record_counts
contains
  subroutine validate_peer_schedule(comm,peers,info)
    integer,intent(in)::comm,peers(:)
    integer,intent(out)::info
    integer::rank,nrank,ierr,i,j,r,total,local_bad,global_bad,astat,count_info
    integer,allocatable::counts(:),displs(:),all_peers(:)
    logical::reciprocal
    info=1;local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    allocate(counts(nrank),displs(nrank),stat=astat);if(astat/=0)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call MPI_Gather(size(peers),1,MPI_INTEGER,counts,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    total=0;displs=0
    if(rank==0)then
      call validate_dg_wpw_record_counts(counts,count_info);if(count_info/=0)local_bad=1
      if(local_bad==0)then
        do i=1,nrank;displs(i)=total;total=total+counts(i);enddo
      endif
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    allocate(all_peers(max(1,total)),stat=astat);local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call MPI_Gatherv(peers,size(peers),MPI_INTEGER,all_peers,counts,displs,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    global_bad=0
    if(rank==0)then
      do r=0,nrank-1
        do i=displs(r+1)+1,displs(r+1)+counts(r+1)
          reciprocal=.false.
          do j=displs(all_peers(i)+1)+1,displs(all_peers(i)+1)+counts(all_peers(i)+1)
            if(all_peers(j)==r)reciprocal=.true.
          enddo
          if(.not.reciprocal)global_bad=1
        enddo
      enddo
    endif
    call MPI_Bcast(global_bad,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    info=0
  end subroutine validate_peer_schedule

  subroutine validate_dg_wpw_record_counts(counts,info)
    integer,intent(in)::counts(:)
    integer,intent(out)::info
    integer::i
    integer(8)::total
    info=1;total=0_8
    do i=1,size(counts)
      if(counts(i)<0.or.counts(i)>huge(0)/2)return
      total=total+int(counts(i),8)
      if(total>int(huge(0)/2,8))return
    enddo
    info=0
  end subroutine validate_dg_wpw_record_counts

  subroutine validate_dg_wpw_projector_provenance(comm,projector_ids,atom_ids,weights,info,local_info)
    integer,intent(in)::comm,projector_ids(:),atom_ids(:)
    real(8),intent(in)::weights(:)
    integer,intent(out)::info
    integer,intent(in),optional::local_info
    integer::i,ierr,local_bad,global_bad,nlocal,nmin,nmax
    integer(8)::h1,h2,h1min,h1max,h2min,h2max,wbits
    info=1;local_bad=0;nlocal=size(projector_ids)
    if(present(local_info))local_bad=merge(0,1,local_info==0)
    if(size(atom_ids)/=nlocal.or.size(weights)/=nlocal.or.nlocal<=0.or.&
       any(projector_ids<=0).or.any(atom_ids<=0).or..not.all(ieee_is_finite(weights)))local_bad=1
    do i=2,nlocal;if(projector_ids(i)<=projector_ids(i-1))local_bad=1;enddo
    call MPI_Allreduce(nlocal,nmin,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(nlocal,nmax,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    if(nmin/=nmax)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    h1=int(nlocal,8);h2=not(int(nlocal,8))
    do i=1,nlocal
      wbits=transfer(weights(i),wbits)
      h1=ieor(ishftc(h1,7),ieor(int(projector_ids(i),8),ishftc(int(atom_ids(i),8),17)))
      h1=ieor(ishftc(h1,11),wbits)
      h2=ieor(ishftc(h2,13),ieor(int(atom_ids(i),8),ishftc(int(projector_ids(i),8),29)))
      h2=ieor(ishftc(h2,19),ishftc(wbits,31))
    enddo
    call MPI_Allreduce(h1,h1min,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(h1,h1max,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(h2,h2min,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(h2,h2max,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    local_bad=merge(0,1,h1min==h1max.and.h2min==h2max)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    info=0
  end subroutine validate_dg_wpw_projector_provenance

  subroutine reduce_dg_wpw_projector_dense(comm,root,local_values,reduced_values,local_info,info)
    integer,intent(in)::comm,root,local_info
    complex(8),intent(in)::local_values(:,:)
    complex(8),intent(out)::reduced_values(:,:)
    integer,intent(out)::info
    integer::rank,nrank,ierr,local_bad,global_bad,dims(2),dmin(2),dmax(2)
    reduced_values=0;info=1;local_bad=merge(0,1,local_info==0)
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    if(root<0.or.root>=nrank.or.any(shape(reduced_values)/=shape(local_values)).or.&
       .not.all_finite_2d(local_values))local_bad=1
    dims=shape(local_values)
    call MPI_Allreduce(dims,dmin,2,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(dims,dmax,2,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    if(any(dmin/=dmax))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call MPI_Reduce(local_values,reduced_values,size(local_values),MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    info=0
  end subroutine reduce_dg_wpw_projector_dense

  subroutine reduce_dg_wpw_projector_records(comm,root,local_records,reduced_records,local_info,info)
    integer,intent(in)::comm,root,local_info
    type(s_dg_wpw_projector_overlap),intent(in)::local_records(:)
    type(s_dg_wpw_projector_overlap),allocatable,intent(out)::reduced_records(:)
    integer,intent(out)::info
    integer::rank,nrank,ierr,local_bad,global_bad,nlocal,total,i,k,astat
    integer(8)::total8
    integer,allocatable::counts(:),displs(:),counts2(:),displs2(:),send_i(:),recv_i(:)
    complex(8),allocatable::send_z(:),recv_z(:)
    type(s_dg_wpw_projector_overlap),allocatable::work(:)
    info=1;local_bad=merge(0,1,local_info==0);nlocal=size(local_records)
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    if(root<0.or.root>=nrank)local_bad=1
    do i=1,nlocal
      if(local_records(i)%basis_id<=0.or.local_records(i)%projector_id<=0.or.&
         .not.ieee_is_finite(real(local_records(i)%value,8)).or.&
         .not.ieee_is_finite(aimag(local_records(i)%value)))local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;allocate(reduced_records(0));return;endif
    allocate(counts(nrank),displs(nrank),counts2(nrank),displs2(nrank),stat=astat)
    local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;allocate(reduced_records(0));return;endif
    call MPI_Gather(nlocal,1,MPI_INTEGER,counts,1,MPI_INTEGER,root,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;allocate(reduced_records(0));return;endif
    total=0;total8=0_8;counts2=0;displs=0;displs2=0
    if(nlocal>huge(0)/2)local_bad=1
    if(rank==root)then
      do i=1,nrank
        total8=total8+int(counts(i),8)
        if(counts(i)>huge(0)/2.or.total8>int(huge(0)/2,8))local_bad=1
      enddo
      if(local_bad==0)then
        total=0
        do i=1,nrank
          displs(i)=total;displs2(i)=2*total;counts2(i)=2*counts(i);total=total+counts(i)
        enddo
      endif
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;allocate(reduced_records(0));return;endif
    allocate(send_i(max(1,2*nlocal)),send_z(max(1,nlocal)),stat=astat)
    local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;allocate(reduced_records(0));return;endif
    do i=1,nlocal
      send_i(2*i-1)=local_records(i)%basis_id;send_i(2*i)=local_records(i)%projector_id
      send_z(i)=local_records(i)%value
    enddo
    allocate(recv_i(max(1,2*total)),recv_z(max(1,total)),stat=astat)
    local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;allocate(reduced_records(0));return;endif
    call MPI_Gatherv(send_i,2*nlocal,MPI_INTEGER,recv_i,counts2,displs2,MPI_INTEGER,root,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;allocate(reduced_records(0));return;endif
    call MPI_Gatherv(send_z,nlocal,MPI_DOUBLE_COMPLEX,recv_z,counts,displs,MPI_DOUBLE_COMPLEX,root,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;allocate(reduced_records(0));return;endif
    if(rank==root)then
      allocate(work(total),reduced_records(total),stat=astat)
    else
      allocate(work(0),reduced_records(0),stat=astat)
    endif
    local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    if(rank/=root)then;info=0;return;endif
    do i=1,total
      work(i)=s_dg_wpw_projector_overlap(recv_i(2*i-1),recv_i(2*i),recv_z(i))
    enddo
    call sort_overlap_records(work);k=0
    do i=1,total
      if(k>0)then
        if(reduced_records(k)%basis_id==work(i)%basis_id.and.&
           reduced_records(k)%projector_id==work(i)%projector_id)then
          reduced_records(k)%value=reduced_records(k)%value+work(i)%value;cycle
        endif
      endif
      k=k+1;reduced_records(k)=work(i)
    enddo
    if(k<total)reduced_records=reduced_records(1:k)
    info=0
  end subroutine reduce_dg_wpw_projector_records

  subroutine build_dg_wpw_projector_overlap_partials(grid_ids,point_weights,w_ids,w_values,p_ids,p_values,&
      sample_grid_ids,sample_projector_ids,sample_values,n_w_namespace,records,info)
    integer,intent(in)::grid_ids(:),w_ids(:),p_ids(:),sample_grid_ids(:),sample_projector_ids(:),n_w_namespace
    real(8),intent(in)::point_weights(:)
    complex(8),intent(in)::w_values(:,:),p_values(:,:),sample_values(:)
    type(s_dg_wpw_projector_overlap),allocatable,intent(out)::records(:)
    integer,intent(out)::info
    type(s_dg_wpw_projector_overlap),allocatable::work(:)
    complex(8)::overlap
    integer::ibasis,iprojector,isample,point,nprojector,nrecord,first_sample,last_sample,basis_id,nbasis
    integer(8)::nbasis8
    integer,allocatable::grid_order(:)

    info=1
    if(size(grid_ids)<=0.or.size(point_weights)/=size(grid_ids).or.&
       any(shape(w_values)/=[size(w_ids),size(grid_ids)]).or.&
       any(shape(p_values)/=[size(p_ids),size(grid_ids)]).or.&
       size(sample_grid_ids)/=size(sample_projector_ids).or.&
       size(sample_values)/=size(sample_grid_ids).or.size(sample_grid_ids)<=0.or.&
       n_w_namespace<=0.or.any(w_ids<=0).or.any(w_ids>n_w_namespace).or.any(p_ids<=0).or.&
       .not.strictly_increasing(w_ids).or..not.strictly_increasing(p_ids).or.&
       .not.all(ieee_is_finite(point_weights)).or..not.all_finite_2d(w_values).or.&
       .not.all_finite_2d(p_values).or..not.all_finite(sample_values))return
    do isample=1,size(sample_projector_ids)
      if(sample_grid_ids(isample)<=0.or.sample_projector_ids(isample)<=0)return
      if(isample>1)then
        if(sample_projector_ids(isample)<sample_projector_ids(isample-1))return
      endif
    enddo
    allocate(grid_order(size(grid_ids)));grid_order=[(point,point=1,size(grid_ids))]
    call sort_grid_order(grid_ids,grid_order)
    do point=1,size(grid_ids)
      if(grid_ids(grid_order(point))<=0)return
      if(point>1)then
        if(grid_ids(grid_order(point))==grid_ids(grid_order(point-1)))return
      endif
    enddo
    nprojector=1
    do isample=2,size(sample_projector_ids)
      if(sample_projector_ids(isample)/=sample_projector_ids(isample-1))nprojector=nprojector+1
    enddo
    nbasis8=int(size(w_ids),8)+int(size(p_ids),8)
    if(nbasis8>int(huge(0),8).or.nbasis8*int(nprojector,8)>int(huge(0),8))return
    nbasis=int(nbasis8);allocate(work(nbasis*nprojector));nrecord=0;first_sample=1
    do iprojector=1,nprojector
      last_sample=first_sample
      do while(last_sample<size(sample_projector_ids))
        if(sample_projector_ids(last_sample+1)/=sample_projector_ids(first_sample))exit
        last_sample=last_sample+1
      enddo
      do ibasis=1,nbasis
        overlap=(0d0,0d0)
        do isample=first_sample,last_sample
          point=find_grid_point(grid_ids,grid_order,sample_grid_ids(isample));if(point==0)cycle
          if(ibasis<=size(w_ids))then
            overlap=overlap+conjg(sample_values(isample))*w_values(ibasis,point)*point_weights(point)
          else
            overlap=overlap+conjg(sample_values(isample))*p_values(ibasis-size(w_ids),point)*point_weights(point)
          endif
        enddo
        if(overlap==(0d0,0d0))cycle
        nrecord=nrecord+1
        if(ibasis<=size(w_ids))then;basis_id=w_ids(ibasis)
        else;basis_id=n_w_namespace+p_ids(ibasis-size(w_ids));endif
        work(nrecord)=s_dg_wpw_projector_overlap(basis_id,sample_projector_ids(first_sample),overlap)
      enddo
      first_sample=last_sample+1
    enddo
    if(.not.all_finite_records(work(1:nrecord)))return
    call sort_overlap_records(work(1:nrecord))
    allocate(records(nrecord));records=work(1:nrecord);info=0
  end subroutine build_dg_wpw_projector_overlap_partials

  integer function find_grid_point(grid_ids,order,target)result(point)
    integer,intent(in)::grid_ids(:),order(:),target
    integer::lo,hi,mid
    point=0;lo=1;hi=size(order)
    do while(lo<=hi)
      mid=(lo+hi)/2
      if(grid_ids(order(mid))<target)then;lo=mid+1
      else if(grid_ids(order(mid))>target)then;hi=mid-1
      else;point=order(mid);return;endif
    enddo
  end function find_grid_point

  subroutine sort_grid_order(grid_ids,order)
    integer,intent(in)::grid_ids(:)
    integer,intent(inout)::order(:)
    integer,allocatable::work(:)
    integer::width,left,mid,right,i,j,k,n
    n=size(order);if(n<2)return;allocate(work(n));width=1
    do while(width<n)
      left=1
      do while(left<=n)
        mid=min(left+width,n+1);right=min(left+2*width,n+1);i=left;j=mid;k=left
        do while(i<mid.and.j<right)
          if(grid_ids(order(i))<=grid_ids(order(j)))then;work(k)=order(i);i=i+1
          else;work(k)=order(j);j=j+1;endif
          k=k+1
        enddo
        do while(i<mid);work(k)=order(i);i=i+1;k=k+1;enddo
        do while(j<right);work(k)=order(j);j=j+1;k=k+1;enddo
        left=right
      enddo
      order=work;width=2*width
    enddo
  end subroutine sort_grid_order

  subroutine reduce_dg_wpw_p_projector_partials(comm,n_g,peers,partial,owned,info,basis_offset)
    integer,intent(in)::comm,n_g,peers(:)
    type(s_dg_wpw_projector_overlap),intent(in)::partial(:)
    type(s_dg_wpw_projector_overlap),allocatable,intent(out)::owned(:)
    integer,intent(out)::info
    integer,intent(in),optional::basis_offset
    integer::rank,nrank,ierr,local_bad,global_bad,npeer,nrequest,i,j,k,total,nlocal,owner,&
      max_send,max_recv,offset,count_info,astat
    integer(8)::delta,owner8
    integer,allocatable::send_counts(:),recv_counts(:),requests(:),send_i(:,:,:),recv_i(:,:,:)
    complex(8),allocatable::send_z(:,:),recv_z(:,:)
    type(s_dg_wpw_projector_overlap),allocatable::work(:),reduced(:)

    info=1;local_bad=0;npeer=size(peers);offset=0;if(present(basis_offset))offset=basis_offset
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    if(n_g<=0.or..not.strictly_increasing(peers).or.any(peers<0).or.any(peers>=nrank).or.&
       any(peers==rank).or..not.valid_sorted_records(partial))local_bad=1
    if(npeer>huge(0)/4)local_bad=1
    do i=1,size(partial)
      if(partial(i)%basis_id<=offset)local_bad=1
      delta=int(partial(i)%basis_id,8)-int(offset,8)-1_8
      owner8=delta/int(max(1,n_g),8)
      if(delta<0_8.or.owner8<0_8.or.owner8>=int(nrank,8))local_bad=1
      owner=int(max(0_8,min(owner8,int(huge(0),8))))
      if(owner/=rank.and..not.any(peers==owner))local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call validate_peer_schedule(comm,peers,count_info);if(count_info/=0)return

    allocate(send_counts(npeer),recv_counts(npeer),stat=astat)
    local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    send_counts=0;recv_counts=-1
    do i=1,size(partial)
      owner=int((int(partial(i)%basis_id,8)-int(offset,8)-1_8)/int(n_g,8))
      do j=1,npeer;if(peers(j)==owner)send_counts(j)=send_counts(j)+1;enddo
    enddo
    allocate(requests(max(1,2*npeer)),stat=astat);local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    nrequest=0
    do i=1,npeer
      nrequest=nrequest+1;call MPI_Irecv(recv_counts(i),1,MPI_INTEGER,peers(i),9611,comm,requests(nrequest),ierr)
      nrequest=nrequest+1;call MPI_Isend(send_counts(i),1,MPI_INTEGER,peers(i),9611,comm,requests(nrequest),ierr)
    enddo
    if(nrequest>0)call MPI_Waitall(nrequest,requests,MPI_STATUSES_IGNORE,ierr)
    if(ierr/=MPI_SUCCESS.or.any(recv_counts<0))local_bad=1
    call validate_dg_wpw_record_counts(send_counts,count_info);if(count_info/=0)local_bad=1
    call validate_dg_wpw_record_counts(recv_counts,count_info);if(count_info/=0)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    max_send=1;max_recv=1
    if(npeer>0)then
      max_send=max(1,maxval(send_counts));max_recv=max(1,maxval(recv_counts))
    endif
    if(2_8*int(max_send,8)*int(max(1,npeer),8)>int(huge(0),8).or.&
       2_8*int(max_recv,8)*int(max(1,npeer),8)>int(huge(0),8))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    allocate(send_i(2,max_send,max(1,npeer)),send_z(max_send,max(1,npeer)),&
      recv_i(2,max_recv,max(1,npeer)),recv_z(max_recv,max(1,npeer)),stat=astat)
    local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    send_i=0;send_z=0;recv_i=0;recv_z=0;send_counts=0
    do i=1,size(partial)
      owner=int((int(partial(i)%basis_id,8)-int(offset,8)-1_8)/int(n_g,8))
      do j=1,npeer
        if(peers(j)==owner)then
          send_counts(j)=send_counts(j)+1;k=send_counts(j)
          send_i(:,k,j)=[partial(i)%basis_id,partial(i)%projector_id];send_z(k,j)=partial(i)%value
        endif
      enddo
    enddo
    deallocate(requests);allocate(requests(max(1,4*npeer)),stat=astat)
    local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    nrequest=0
    do i=1,npeer
      nrequest=nrequest+1;call MPI_Irecv(recv_i(1,1,i),2*recv_counts(i),MPI_INTEGER,peers(i),9612,comm,requests(nrequest),ierr)
      nrequest=nrequest+1;call MPI_Irecv(recv_z(1,i),recv_counts(i),MPI_DOUBLE_COMPLEX,peers(i),9613,comm,requests(nrequest),ierr)
      nrequest=nrequest+1;call MPI_Isend(send_i(1,1,i),2*send_counts(i),MPI_INTEGER,peers(i),9612,comm,requests(nrequest),ierr)
      nrequest=nrequest+1;call MPI_Isend(send_z(1,i),send_counts(i),MPI_DOUBLE_COMPLEX,peers(i),9613,comm,requests(nrequest),ierr)
    enddo
    if(nrequest>0)call MPI_Waitall(nrequest,requests,MPI_STATUSES_IGNORE,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1

    nlocal=0
    do i=1,size(partial)
      if((int(partial(i)%basis_id,8)-int(offset,8)-1_8)/int(n_g,8)==int(rank,8))nlocal=nlocal+1
    enddo
    call validate_dg_wpw_record_counts([nlocal,recv_counts],count_info)
    if(count_info/=0)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    total=nlocal+sum(recv_counts);allocate(work(total),reduced(total),stat=astat)
    local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    k=0
    do i=1,size(partial)
      if((int(partial(i)%basis_id,8)-int(offset,8)-1_8)/int(n_g,8)==int(rank,8))then
        k=k+1;work(k)=partial(i)
      endif
    enddo
    do i=1,npeer;do j=1,recv_counts(i)
      k=k+1;work(k)%basis_id=recv_i(1,j,i);work(k)%projector_id=recv_i(2,j,i);work(k)%value=recv_z(j,i)
      if((int(work(k)%basis_id,8)-int(offset,8)-1_8)/int(n_g,8)/=int(rank,8))local_bad=1
    enddo;enddo
    call sort_overlap_records(work)
    k=0
    do i=1,total
      if(k>0)then
        if(reduced(k)%basis_id==work(i)%basis_id.and.&
           reduced(k)%projector_id==work(i)%projector_id)then
          reduced(k)%value=reduced(k)%value+work(i)%value
          cycle
        endif
      endif
      k=k+1;reduced(k)=work(i)
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    allocate(owned(k),stat=astat);local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    owned=reduced(1:k);info=0
  end subroutine reduce_dg_wpw_p_projector_partials

  subroutine exchange_dg_wpw_projector_overlaps(comm,peers,owned,support,info)
    integer,intent(in)::comm,peers(:)
    type(s_dg_wpw_projector_overlap),intent(in)::owned(:)
    type(s_dg_wpw_projector_overlap),allocatable,intent(out)::support(:)
    integer,intent(out)::info
    integer::rank,nrank,ierr,local_bad,global_bad,npeer,nrequest,i,j,k,total,send_count,buffer_offset,count_info,astat
    integer,allocatable::recv_counts(:),requests(:),send_i(:),recv_i(:,:),offsets(:)
    complex(8),allocatable::send_z(:),recv_z(:,:)
    type(s_dg_wpw_projector_overlap),allocatable::work(:)
    info=1;local_bad=0;npeer=size(peers)
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    if(.not.strictly_increasing(peers).or.any(peers<0).or.any(peers>=nrank).or.any(peers==rank).or.&
       .not.valid_sorted_records(owned))local_bad=1
    if(npeer>huge(0)/4)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call validate_peer_schedule(comm,peers,count_info);if(count_info/=0)return
    send_count=size(owned)
    allocate(recv_counts(npeer),requests(max(1,2*npeer)),stat=astat)
    local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    recv_counts=-1;nrequest=0
    do i=1,npeer
      nrequest=nrequest+1;call MPI_Irecv(recv_counts(i),1,MPI_INTEGER,peers(i),9601,comm,requests(nrequest),ierr)
      if(ierr/=MPI_SUCCESS)local_bad=1
    enddo
    do i=1,npeer
      nrequest=nrequest+1;call MPI_Isend(send_count,1,MPI_INTEGER,peers(i),9601,comm,requests(nrequest),ierr)
      if(ierr/=MPI_SUCCESS)local_bad=1
    enddo
    if(nrequest>0)call MPI_Waitall(nrequest,requests,MPI_STATUSES_IGNORE,ierr)
    if(ierr/=MPI_SUCCESS.or.any(recv_counts<0))local_bad=1
    call validate_dg_wpw_record_counts([send_count,recv_counts],count_info)
    if(count_info/=0)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    total=sum(recv_counts);allocate(offsets(npeer+1),send_i(2*size(owned)),send_z(size(owned)),&
      recv_i(2,max(1,total)),recv_z(1,max(1,total)),stat=astat)
    local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    offsets(1)=1
    do i=1,npeer;offsets(i+1)=offsets(i)+recv_counts(i);enddo
    do i=1,size(owned)
      send_i(2*i-1:2*i)=[owned(i)%basis_id,owned(i)%projector_id];send_z(i)=owned(i)%value
    enddo
    recv_i=0;recv_z=0
    deallocate(requests);allocate(requests(max(1,4*npeer)),stat=astat)
    local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    nrequest=0
    do i=1,npeer
      buffer_offset=max(1,min(offsets(i),max(1,total)))
      nrequest=nrequest+1;call MPI_Irecv(recv_i(1,buffer_offset),2*recv_counts(i),MPI_INTEGER,&
        peers(i),9602,comm,requests(nrequest),ierr)
      if(ierr/=MPI_SUCCESS)local_bad=1
      nrequest=nrequest+1;call MPI_Irecv(recv_z(1,buffer_offset),recv_counts(i),MPI_DOUBLE_COMPLEX,&
        peers(i),9603,comm,requests(nrequest),ierr)
      if(ierr/=MPI_SUCCESS)local_bad=1
    enddo
    do i=1,npeer
      nrequest=nrequest+1;call MPI_Isend(send_i,size(send_i),MPI_INTEGER,peers(i),9602,comm,requests(nrequest),ierr)
      if(ierr/=MPI_SUCCESS)local_bad=1
      nrequest=nrequest+1;call MPI_Isend(send_z,size(send_z),MPI_DOUBLE_COMPLEX,peers(i),9603,comm,requests(nrequest),ierr)
      if(ierr/=MPI_SUCCESS)local_bad=1
    enddo
    if(nrequest>0)call MPI_Waitall(nrequest,requests,MPI_STATUSES_IGNORE,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    allocate(work(size(owned)+total),stat=astat);local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    work(1:size(owned))=owned;k=size(owned)
    do i=1,npeer
      do j=offsets(i),offsets(i+1)-1
        k=k+1;work(k)%basis_id=recv_i(1,j);work(k)%projector_id=recv_i(2,j);work(k)%value=recv_z(1,j)
      enddo
    enddo
    call sort_overlap_records(work)
    if(.not.valid_sorted_records(work))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call move_alloc(work,support);info=0
  end subroutine

  subroutine assemble_dg_wpw_nonlocal_blocks(records,projector_ids,weights,rows,columns,values,info)
    type(s_dg_wpw_projector_overlap),intent(in)::records(:)
    integer,intent(in)::projector_ids(:),rows(:),columns(:)
    real(8),intent(in)::weights(:)
    complex(8),intent(out)::values(:)
    integer,intent(out)::info
    integer::i,j,row_first,row_last,column_first,column_last,row_record,column_record
    values=(0d0,0d0);info=1
    if(size(projector_ids)/=size(weights).or.size(rows)/=size(columns).or.size(rows)/=size(values).or.&
       any(projector_ids<=0).or.any(rows<=0).or.any(columns<=0).or..not.all(ieee_is_finite(weights)))return
    do i=2,size(projector_ids)
      if(projector_ids(i)<=projector_ids(i-1))return
    enddo
    do i=1,size(records)
      if(records(i)%basis_id<=0.or.records(i)%projector_id<=0.or..not.finite_complex(records(i)%value))return
      if(find_integer_sorted(projector_ids,records(i)%projector_id)==0)return
      if(i>1)then
        if(records(i)%basis_id<records(i-1)%basis_id.or.&
           (records(i)%basis_id==records(i-1)%basis_id.and.&
            records(i)%projector_id<=records(i-1)%projector_id))return
      endif
    enddo
    do i=1,size(values)
      call find_basis_range(records,rows(i),row_first,row_last)
      call find_basis_range(records,columns(i),column_first,column_last)
      row_record=row_first;column_record=column_first
      do while(row_record<=row_last.and.column_record<=column_last)
        if(records(row_record)%projector_id<records(column_record)%projector_id)then
          row_record=row_record+1
        else if(records(row_record)%projector_id>records(column_record)%projector_id)then
          column_record=column_record+1
        else
          j=find_integer_sorted(projector_ids,records(row_record)%projector_id)
          if(j<=0)then;values=0;return;endif
          values(i)=values(i)+conjg(records(row_record)%value)*weights(j)*records(column_record)%value
          row_record=row_record+1;column_record=column_record+1
        endif
      enddo
    enddo
    if(.not.all_finite(values))then;values=0;return;endif
    info=0
  end subroutine

  subroutine find_basis_range(records,basis_id,first,last)
    type(s_dg_wpw_projector_overlap),intent(in)::records(:)
    integer,intent(in)::basis_id
    integer,intent(out)::first,last
    integer::lo,hi,mid
    first=1;last=0;lo=1;hi=size(records)
    do while(lo<=hi)
      mid=(lo+hi)/2
      if(records(mid)%basis_id<basis_id)then;lo=mid+1
      else;hi=mid-1;endif
    enddo
    if(lo>size(records))return
    if(records(lo)%basis_id/=basis_id)return
    first=lo;hi=size(records)
    do while(lo<=hi)
      mid=(lo+hi)/2
      if(records(mid)%basis_id<=basis_id)then;lo=mid+1
      else;hi=mid-1;endif
    enddo
    last=hi
  end subroutine find_basis_range

  integer function find_integer_sorted(values,target)result(position)
    integer,intent(in)::values(:),target
    integer::lo,hi,mid
    position=0;lo=1;hi=size(values)
    do while(lo<=hi)
      mid=(lo+hi)/2
      if(values(mid)<target)then;lo=mid+1
      else if(values(mid)>target)then;hi=mid-1
      else;position=mid;return;endif
    enddo
  end function find_integer_sorted

  integer function find_record(records,basis_id,projector_id)result(position)
    type(s_dg_wpw_projector_overlap),intent(in)::records(:)
    integer,intent(in)::basis_id,projector_id
    integer::i
    position=0
    do i=1,size(records)
      if(records(i)%basis_id>basis_id)return
      if(records(i)%basis_id==basis_id.and.records(i)%projector_id==projector_id)then
        position=i;return
      endif
    enddo
  end function

  integer function find_integer(values,target)result(position)
    integer,intent(in)::values(:),target
    integer::i
    position=0
    do i=1,size(values)
      if(values(i)==target)then;position=i;return;endif
    enddo
  end function

  logical function finite_complex(value)result(ok)
    complex(8),intent(in)::value
    ok=ieee_is_finite(real(value,8)).and.ieee_is_finite(aimag(value))
  end function

  logical function all_finite(values)result(ok)
    complex(8),intent(in)::values(:)
    integer::i
    ok=.true.
    do i=1,size(values)
      if(.not.finite_complex(values(i)))then;ok=.false.;return;endif
    enddo
  end function

  logical function all_finite_2d(values)result(ok)
    complex(8),intent(in)::values(:,:)
    ok=all(ieee_is_finite(real(values,8))).and.all(ieee_is_finite(aimag(values)))
  end function

  logical function all_finite_records(records)result(ok)
    type(s_dg_wpw_projector_overlap),intent(in)::records(:)
    integer::i
    ok=.true.
    do i=1,size(records);if(.not.finite_complex(records(i)%value))then;ok=.false.;return;endif;enddo
  end function

  logical function valid_sorted_records(records)result(ok)
    type(s_dg_wpw_projector_overlap),intent(in)::records(:)
    integer::i
    ok=.true.
    do i=1,size(records)
      if(records(i)%basis_id<=0.or.records(i)%projector_id<=0.or..not.finite_complex(records(i)%value))then
        ok=.false.;return
      endif
      if(i>1)then
        if(records(i)%basis_id<records(i-1)%basis_id.or.&
           (records(i)%basis_id==records(i-1)%basis_id.and.&
            records(i)%projector_id<=records(i-1)%projector_id))then;ok=.false.;return;endif
      endif
    enddo
  end function

  logical function strictly_increasing(values)result(ok)
    integer,intent(in)::values(:)
    integer::i
    ok=.true.
    do i=2,size(values);if(values(i)<=values(i-1))then;ok=.false.;return;endif;enddo
  end function

  subroutine sort_overlap_records(records)
    type(s_dg_wpw_projector_overlap),intent(inout)::records(:)
    type(s_dg_wpw_projector_overlap),allocatable::work(:)
    integer::width,left,mid,right,i,j,k,n
    n=size(records);if(n<2)return;allocate(work(n));width=1
    do while(width<n)
      left=1
      do while(left<=n)
        mid=min(left+width,n+1);right=min(left+2*width,n+1);i=left;j=mid;k=left
        do while(i<mid.and.j<right)
          if(record_less_equal(records(i),records(j)))then;work(k)=records(i);i=i+1
          else;work(k)=records(j);j=j+1;endif
          k=k+1
        enddo
        do while(i<mid);work(k)=records(i);i=i+1;k=k+1;enddo
        do while(j<right);work(k)=records(j);j=j+1;k=k+1;enddo
        left=right
      enddo
      records=work;width=2*width
    enddo
  end subroutine

  logical function record_less_equal(a,b)result(le)
    type(s_dg_wpw_projector_overlap),intent(in)::a,b
    le=a%basis_id<b%basis_id.or.(a%basis_id==b%basis_id.and.a%projector_id<=b%projector_id)
  end function
end module
