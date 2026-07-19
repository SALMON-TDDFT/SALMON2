module rt_dg_wpw_volume_halo_provider
  use mpi
  use,intrinsic::iso_fortran_env,only:int64
  implicit none
  private

  type,public::s_dg_wpw_volume_halo_state
    integer::epoch=-1,source_rank=-1
    integer::image_shift(3)=0
    integer::box_lo(3)=0,box_hi(3)=-1
    logical::valid=.false.
    integer(int64)::checksum=0_int64
    integer,allocatable::w_ids(:)
    complex(8),allocatable::values(:,:),gradients(:,:,:)
  end type
  type,public::s_dg_wpw_volume_halo_send
    integer::peer=-1
    integer::image_shift(3)=0
    integer::box_lo(3)=0,box_hi(3)=-1
    integer,allocatable::w_ids(:)
    complex(8),allocatable::values(:,:),gradients(:,:,:)
  end type
  public::exchange_dg_wpw_volume_halo,exchange_dg_wpw_volume_halo_schedule
  public::read_dg_wpw_volume_halo
  public::dg_wpw_volume_payload_checksum
  public::map_dg_wpw_volume_halo_box_to_canonical
  public::pack_dg_wpw_volume_halo_send
  public::broadcast_dg_wpw_volume_halos
  public::dg_wpw_volume_header_is_bounded

contains

  pure logical function dg_wpw_volume_header_is_bounded(header,max_w,max_extent,max_points)result(ok)
    integer,intent(in)::header(8),max_w,max_extent(3)
    integer(int64),intent(in)::max_points
    integer::extent(3),i
    integer(int64)::points
    ok=.false.
    if(max_w<=0.or.any(max_extent<=0).or.max_points<=0_int64.or.header(2)<=0.or.header(2)>max_w.or.&
       any(header(6:8)<header(3:5)))return
    extent=header(6:8)-header(3:5)+1
    if(any(extent<=0).or.any(extent>max_extent))return
    points=1_int64
    do i=1,3
      if(points>max_points/int(extent(i),int64))return
      points=points*int(extent(i),int64)
    enddo
    ok=points<=max_points.and.points<=int(huge(1),int64)
  end function dg_wpw_volume_header_is_bounded

  subroutine broadcast_dg_wpw_volume_halos(comm,root,halos,info)
    integer,intent(in)::comm,root
    type(s_dg_wpw_volume_halo_state),allocatable,intent(inout)::halos(:)
    integer,intent(out)::info
    type(s_dg_wpw_volume_halo_state),allocatable::candidate(:)
    integer::rank,nrank,ierr,nhalo,i,header(12),npoint,local_bad,global_bad,extent(3)
    integer::max_w,max_extent(3),header8(8)
    integer(int64)::checksum,max_points,points64

    info=1;local_bad=0;max_w=0;max_extent=0;max_points=0_int64
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    if(root<0.or.root>=nrank)local_bad=1
    nhalo=0
    if(rank==root)then
      if(.not.allocated(halos))then;local_bad=1
      else
        nhalo=size(halos);if(nhalo<0.or.nhalo>26)local_bad=1
        do i=1,nhalo
          if(.not.halos(i)%valid.or.halos(i)%epoch<=0.or.halos(i)%source_rank<0.or.&
             any(abs(halos(i)%image_shift)>1).or.any(halos(i)%box_hi<halos(i)%box_lo).or.&
             .not.allocated(halos(i)%w_ids).or..not.allocated(halos(i)%values).or.&
             .not.allocated(halos(i)%gradients))then;local_bad=1;cycle;endif
          extent=halos(i)%box_hi-halos(i)%box_lo+1;points64=product(int(extent,int64))
          if(points64<=0.or.points64>int(huge(npoint),int64))then;local_bad=1;cycle;endif
          npoint=int(points64);max_w=max(max_w,size(halos(i)%w_ids));max_extent=max(max_extent,extent)
          max_points=max(max_points,points64)
          if(.not.strictly_increasing(halos(i)%w_ids).or.&
             any(shape(halos(i)%values)/=[size(halos(i)%w_ids),npoint]).or.&
             any(shape(halos(i)%gradients)/=[3,size(halos(i)%w_ids),npoint]))local_bad=1
        enddo
      endif
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call MPI_Bcast(nhalo,1,MPI_INTEGER,root,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.nhalo<0.or.nhalo>26)return
    call MPI_Bcast(max_w,1,MPI_INTEGER,root,comm,ierr)
    call MPI_Bcast(max_extent,3,MPI_INTEGER,root,comm,ierr)
    call MPI_Bcast(max_points,1,MPI_INTEGER8,root,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.(nhalo>0.and.(max_w<=0.or.any(max_extent<=0).or.max_points<=0)))return
    if(rank==root)then;allocate(candidate,source=halos)
    else;allocate(candidate(nhalo));endif
    do i=1,nhalo
      if(rank==root)header=[halos(i)%epoch,halos(i)%source_rank,halos(i)%image_shift,&
        halos(i)%box_lo,halos(i)%box_hi,size(halos(i)%w_ids)]
      call MPI_Bcast(header,12,MPI_INTEGER,root,comm,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
      header8=[header(1),header(12),header(6:8),header(9:11)]
      if(.not.dg_wpw_volume_header_is_bounded(header8,max_w,max_extent,max_points))local_bad=1
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
        if(allocated(candidate))then
          do npoint=1,size(candidate);candidate(npoint)%valid=.false.;enddo
          call move_alloc(candidate,halos)
        endif
        return
      endif
      if(rank/=root)then
        candidate(i)%epoch=header(1);candidate(i)%source_rank=header(2);candidate(i)%image_shift=header(3:5)
        candidate(i)%box_lo=header(6:8);candidate(i)%box_hi=header(9:11)
        npoint=product(candidate(i)%box_hi-candidate(i)%box_lo+1)
        allocate(candidate(i)%w_ids(header(12)),candidate(i)%values(header(12),npoint),&
          candidate(i)%gradients(3,header(12),npoint))
      endif
      checksum=candidate(i)%checksum
      call MPI_Bcast(checksum,1,MPI_INTEGER8,root,comm,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
      candidate(i)%checksum=checksum
      call MPI_Bcast(candidate(i)%w_ids,size(candidate(i)%w_ids),MPI_INTEGER,root,comm,ierr)
      if(ierr/=MPI_SUCCESS)local_bad=1
      call MPI_Bcast(candidate(i)%values,size(candidate(i)%values),MPI_DOUBLE_COMPLEX,root,comm,ierr)
      if(ierr/=MPI_SUCCESS)local_bad=1
      call MPI_Bcast(candidate(i)%gradients,size(candidate(i)%gradients),MPI_DOUBLE_COMPLEX,root,comm,ierr)
      if(ierr/=MPI_SUCCESS)local_bad=1
      checksum=dg_wpw_volume_payload_checksum(candidate(i)%w_ids,candidate(i)%box_lo,candidate(i)%box_hi,&
        candidate(i)%values,candidate(i)%gradients)
      if(checksum/=candidate(i)%checksum.or..not.all(finite_complex(candidate(i)%values)).or.&
         .not.all(finite_complex(candidate(i)%gradients)))local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      do i=1,nhalo;candidate(i)%valid=.false.;enddo
      call move_alloc(candidate,halos);return
    endif
    do i=1,nhalo;candidate(i)%valid=.true.;enddo
    call move_alloc(candidate,halos);info=0
  end subroutine broadcast_dg_wpw_volume_halos

  pure subroutine map_dg_wpw_volume_halo_box_to_canonical(unwrapped_lo,unwrapped_hi,cell_extent,&
      image_shift,canonical_lo,canonical_hi,info)
    integer,intent(in)::unwrapped_lo(3),unwrapped_hi(3),cell_extent(3),image_shift(3)
    integer,intent(out)::canonical_lo(3),canonical_hi(3),info
    canonical_lo=0;canonical_hi=-1;info=1
    if(any(unwrapped_hi<unwrapped_lo).or.any(cell_extent<=0).or.any(abs(image_shift)>1))return
    canonical_lo=unwrapped_lo-image_shift*cell_extent
    canonical_hi=unwrapped_hi-image_shift*cell_extent
    info=0
  end subroutine map_dg_wpw_volume_halo_box_to_canonical

  subroutine pack_dg_wpw_volume_halo_send(peer,route_shift,periodic_shift,unwrapped_lo,unwrapped_hi,&
      cell_extent,buffer_lo,buffer_hi,w_ids,buffer_values,buffer_gradients,send,info)
    integer,intent(in)::peer,route_shift(3),periodic_shift(3),unwrapped_lo(3),unwrapped_hi(3),cell_extent(3)
    integer,intent(in)::buffer_lo(3),buffer_hi(3),w_ids(:)
    complex(8),intent(in)::buffer_values(:,:),buffer_gradients(:,:,:)
    type(s_dg_wpw_volume_halo_send),intent(out)::send
    integer,intent(out)::info
    complex(8),allocatable::packed_values(:,:),packed_gradients(:,:,:)
    integer::buffer_extent(3),box_extent(3),grid(3),source_point,destination_point,map_info,ix,iy,iz

    info=1
    if(peer<0.or.any(abs(route_shift)>1).or.all(route_shift==0).or.any(abs(periodic_shift)>1).or.&
       any(unwrapped_hi<unwrapped_lo).or.any(buffer_hi<buffer_lo).or.&
       any(unwrapped_lo<buffer_lo).or.any(unwrapped_hi>buffer_hi).or..not.strictly_increasing(w_ids))return
    buffer_extent=buffer_hi-buffer_lo+1;box_extent=unwrapped_hi-unwrapped_lo+1
    if(any(shape(buffer_values)/=[size(w_ids),product(buffer_extent)]).or.&
       any(shape(buffer_gradients)/=[3,size(w_ids),product(buffer_extent)]))return
    allocate(packed_values(size(w_ids),product(box_extent)),&
      packed_gradients(3,size(w_ids),product(box_extent)))
    destination_point=0
    do iz=unwrapped_lo(3),unwrapped_hi(3)
      do iy=unwrapped_lo(2),unwrapped_hi(2)
        do ix=unwrapped_lo(1),unwrapped_hi(1)
          grid=[ix,iy,iz]
          destination_point=destination_point+1
          source_point=grid(1)-buffer_lo(1)+1+(grid(2)-buffer_lo(2))*buffer_extent(1)+&
            (grid(3)-buffer_lo(3))*buffer_extent(1)*buffer_extent(2)
          packed_values(:,destination_point)=buffer_values(:,source_point)
          packed_gradients(:,:,destination_point)=buffer_gradients(:,:,source_point)
        enddo
      enddo
    enddo
    if(.not.all(finite_complex(packed_values)).or..not.all(finite_complex(packed_gradients)))return
    call map_dg_wpw_volume_halo_box_to_canonical(unwrapped_lo,unwrapped_hi,cell_extent,periodic_shift,&
      send%box_lo,send%box_hi,map_info)
    if(map_info/=0)return
    send%peer=peer;send%image_shift=route_shift
    allocate(send%w_ids,source=w_ids);call move_alloc(packed_values,send%values)
    call move_alloc(packed_gradients,send%gradients);info=0
  end subroutine pack_dg_wpw_volume_halo_send

  subroutine exchange_dg_wpw_volume_halo_schedule(comm,epoch,send,halos,info)
    integer,intent(in)::comm,epoch
    type(s_dg_wpw_volume_halo_send),intent(in)::send(:)
    type(s_dg_wpw_volume_halo_state),allocatable,intent(out)::halos(:)
    integer,intent(out)::info
    integer::rank,nrank,ierr,i,j,nneighbor,npoint,recv_npoint,extent(3)
    integer::local_bad,global_bad,nrequest
    integer::local_max_w,global_max_w,local_max_extent(3),global_max_extent(3)
    integer(int64)::local_max_points,global_max_points,points64
    integer,allocatable::send_header(:,:),recv_header(:,:),request(:),statuses(:,:)
    integer(int64),allocatable::send_checksum(:),recv_checksum(:)

    info=0;nneighbor=size(send);local_bad=0;local_max_w=0;local_max_extent=0;local_max_points=0_int64
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    if(epoch<=0.or.nneighbor<=0)local_bad=1
    do i=1,nneighbor
      if(send(i)%peer<0.or.send(i)%peer>=nrank.or.send(i)%peer==rank)local_bad=1
      if(any(abs(send(i)%image_shift)>1))local_bad=1
      do j=1,i-1
        if(send(j)%peer==send(i)%peer.and.all(send(j)%image_shift==send(i)%image_shift))local_bad=1
      enddo
      if(.not.allocated(send(i)%w_ids).or..not.allocated(send(i)%values).or.&
         .not.allocated(send(i)%gradients).or.any(send(i)%box_hi<send(i)%box_lo))then
        local_bad=1;cycle
      endif
      extent=send(i)%box_hi-send(i)%box_lo+1;points64=product(int(extent,int64))
      if(points64<=0_int64.or.points64>int(huge(npoint),int64))then;local_bad=1;cycle;endif
      npoint=int(points64);local_max_w=max(local_max_w,size(send(i)%w_ids))
      local_max_extent=max(local_max_extent,extent);local_max_points=max(local_max_points,points64)
      if(.not.strictly_increasing(send(i)%w_ids).or.&
         any(shape(send(i)%values)/=[size(send(i)%w_ids),npoint]).or.&
         any(shape(send(i)%gradients)/=[3,size(send(i)%w_ids),npoint]).or.&
         .not.all(finite_complex(send(i)%values)).or..not.all(finite_complex(send(i)%gradients)))local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;info=1;allocate(halos(0));return;endif
    call MPI_Allreduce(local_max_w,global_max_w,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_max_extent,global_max_extent,3,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_max_points,global_max_points,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_max_w<=0.or.any(global_max_extent<=0).or.global_max_points<=0)then
      info=1;allocate(halos(0));return
    endif

    allocate(send_header(8,nneighbor),recv_header(8,nneighbor),send_checksum(nneighbor),&
      recv_checksum(nneighbor),request(4*nneighbor),statuses(MPI_STATUS_SIZE,4*nneighbor))
    nrequest=0
    do i=1,nneighbor
      send_header(:,i)=[epoch,size(send(i)%w_ids),send(i)%box_lo,send(i)%box_hi]
      send_checksum(i)=dg_wpw_volume_payload_checksum(send(i)%w_ids,send(i)%box_lo,send(i)%box_hi,&
        send(i)%values,send(i)%gradients)
      nrequest=nrequest+1
      call MPI_Irecv(recv_header(:,i),8,MPI_INTEGER,send(i)%peer,&
        600+encode_image_shift(-send(i)%image_shift),comm,request(nrequest),ierr)
      nrequest=nrequest+1
      call MPI_Irecv(recv_checksum(i),1,MPI_INTEGER8,send(i)%peer,&
        700+encode_image_shift(-send(i)%image_shift),comm,request(nrequest),ierr)
    enddo
    do i=1,nneighbor
      nrequest=nrequest+1
      call MPI_Isend(send_header(:,i),8,MPI_INTEGER,send(i)%peer,&
        600+encode_image_shift(send(i)%image_shift),comm,request(nrequest),ierr)
      nrequest=nrequest+1
      call MPI_Isend(send_checksum(i),1,MPI_INTEGER8,send(i)%peer,&
        700+encode_image_shift(send(i)%image_shift),comm,request(nrequest),ierr)
    enddo
    call MPI_Waitall(nrequest,request,statuses,ierr)
    local_bad=merge(0,1,ierr==MPI_SUCCESS)
    allocate(halos(nneighbor))
    do i=1,nneighbor
      if(recv_header(1,i)/=epoch.or..not.dg_wpw_volume_header_is_bounded(recv_header(:,i),&
         global_max_w,global_max_extent,global_max_points))then
        local_bad=1;cycle
      endif
      recv_npoint=int(product(int(recv_header(6:8,i)-recv_header(3:5,i)+1,int64)))
      allocate(halos(i)%w_ids(recv_header(2,i)),halos(i)%values(recv_header(2,i),recv_npoint),&
        halos(i)%gradients(3,recv_header(2,i),recv_npoint))
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;info=2;return;endif

    deallocate(request,statuses);allocate(request(6*nneighbor),statuses(MPI_STATUS_SIZE,6*nneighbor));nrequest=0
    do i=1,nneighbor
      nrequest=nrequest+1;call MPI_Irecv(halos(i)%w_ids,size(halos(i)%w_ids),MPI_INTEGER,&
        send(i)%peer,800+encode_image_shift(-send(i)%image_shift),comm,request(nrequest),ierr)
      nrequest=nrequest+1;call MPI_Irecv(halos(i)%values,size(halos(i)%values),MPI_DOUBLE_COMPLEX,&
        send(i)%peer,900+encode_image_shift(-send(i)%image_shift),comm,request(nrequest),ierr)
      nrequest=nrequest+1;call MPI_Irecv(halos(i)%gradients,size(halos(i)%gradients),MPI_DOUBLE_COMPLEX,&
        send(i)%peer,1000+encode_image_shift(-send(i)%image_shift),comm,request(nrequest),ierr)
    enddo
    do i=1,nneighbor
      nrequest=nrequest+1;call MPI_Isend(send(i)%w_ids,size(send(i)%w_ids),MPI_INTEGER,&
        send(i)%peer,800+encode_image_shift(send(i)%image_shift),comm,request(nrequest),ierr)
      nrequest=nrequest+1;call MPI_Isend(send(i)%values,size(send(i)%values),MPI_DOUBLE_COMPLEX,&
        send(i)%peer,900+encode_image_shift(send(i)%image_shift),comm,request(nrequest),ierr)
      nrequest=nrequest+1;call MPI_Isend(send(i)%gradients,size(send(i)%gradients),MPI_DOUBLE_COMPLEX,&
        send(i)%peer,1000+encode_image_shift(send(i)%image_shift),comm,request(nrequest),ierr)
    enddo
    call MPI_Waitall(nrequest,request,statuses,ierr)
    local_bad=merge(0,1,ierr==MPI_SUCCESS)
    do i=1,nneighbor
      if(.not.strictly_increasing(halos(i)%w_ids).or..not.all(finite_complex(halos(i)%values)).or.&
         .not.all(finite_complex(halos(i)%gradients)))local_bad=1
      halos(i)%checksum=dg_wpw_volume_payload_checksum(halos(i)%w_ids,recv_header(3:5,i),&
        recv_header(6:8,i),halos(i)%values,halos(i)%gradients)
      if(halos(i)%checksum/=recv_checksum(i))local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;info=3;return;endif
    do i=1,nneighbor
      halos(i)%epoch=epoch;halos(i)%source_rank=send(i)%peer
      halos(i)%image_shift=send(i)%image_shift
      halos(i)%box_lo=recv_header(3:5,i);halos(i)%box_hi=recv_header(6:8,i);halos(i)%valid=.true.
    enddo
  end subroutine exchange_dg_wpw_volume_halo_schedule

  subroutine exchange_dg_wpw_volume_halo(comm,epoch,peer,send_w_ids,send_box_lo,send_box_hi, &
      send_values,send_gradients,halo,info)
    integer,intent(in)::comm,epoch,peer,send_w_ids(:),send_box_lo(3),send_box_hi(3)
    complex(8),intent(in)::send_values(:,:),send_gradients(:,:,:)
    type(s_dg_wpw_volume_halo_state),intent(out)::halo
    integer,intent(out)::info
    integer::rank,nrank,ierr,local_bad,global_bad,npoint,recv_npoint
    integer::send_header(8),recv_header(8),status(MPI_STATUS_SIZE),i
    integer(int64)::send_checksum,recv_checksum

    halo%valid=.false.;info=0;local_bad=0
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Comm_size(comm,nrank,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    npoint=0
    if(any(send_box_hi<send_box_lo))then
      local_bad=1
    else
      npoint=product(send_box_hi-send_box_lo+1)
    endif
    if(epoch<=0.or.peer<0.or.peer>=nrank.or.peer==rank.or.size(send_w_ids)<=0.or.&
       .not.strictly_increasing(send_w_ids).or.any(shape(send_values)/=[size(send_w_ids),npoint]).or.&
       any(shape(send_gradients)/=[3,size(send_w_ids),npoint]))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;info=1;return;endif

    send_header=[epoch,size(send_w_ids),send_box_lo,send_box_hi]
    send_checksum=dg_wpw_volume_payload_checksum(send_w_ids,send_box_lo,send_box_hi,send_values,send_gradients)
    call MPI_Sendrecv(send_header,8,MPI_INTEGER,peer,510,recv_header,8,MPI_INTEGER,peer,510, &
      comm,status,ierr)
    if(ierr/=MPI_SUCCESS)then;info=2;return;endif
    call MPI_Sendrecv(send_checksum,1,MPI_INTEGER8,peer,514,recv_checksum,1,MPI_INTEGER8,peer,514,&
      comm,status,ierr)
    if(ierr/=MPI_SUCCESS)then;info=2;return;endif
    recv_npoint=0;local_bad=0
    if(recv_header(1)/=epoch.or.recv_header(2)<=0.or.any(recv_header(6:8)<recv_header(3:5)))then
      local_bad=1
    else
      recv_npoint=product(recv_header(6:8)-recv_header(3:5)+1)
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;info=2;return;endif

    allocate(halo%w_ids(recv_header(2)),halo%values(recv_header(2),recv_npoint), &
      halo%gradients(3,recv_header(2),recv_npoint))
    call MPI_Sendrecv(send_w_ids,size(send_w_ids),MPI_INTEGER,peer,511,halo%w_ids,size(halo%w_ids), &
      MPI_INTEGER,peer,511,comm,status,ierr)
    if(ierr/=MPI_SUCCESS)then;info=3;return;endif
    call MPI_Sendrecv(send_values,size(send_values),MPI_DOUBLE_COMPLEX,peer,512,halo%values, &
      size(halo%values),MPI_DOUBLE_COMPLEX,peer,512,comm,status,ierr)
    if(ierr/=MPI_SUCCESS)then;info=3;return;endif
    call MPI_Sendrecv(send_gradients,size(send_gradients),MPI_DOUBLE_COMPLEX,peer,513,halo%gradients, &
      size(halo%gradients),MPI_DOUBLE_COMPLEX,peer,513,comm,status,ierr)
    if(ierr/=MPI_SUCCESS)then;info=3;return;endif
    local_bad=merge(0,1,strictly_increasing(halo%w_ids).and.all(finite_complex(halo%values)).and.&
      all(finite_complex(halo%gradients)))
    halo%checksum=dg_wpw_volume_payload_checksum(halo%w_ids,recv_header(3:5),recv_header(6:8),&
      halo%values,halo%gradients)
    if(halo%checksum/=recv_checksum)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;info=4;return;endif
    halo%epoch=epoch;halo%source_rank=peer
    halo%box_lo=recv_header(3:5);halo%box_hi=recv_header(6:8);halo%valid=.true.
  end subroutine

  subroutine read_dg_wpw_volume_halo(halo,w_id,grid,epoch,value,gradient,info)
    use,intrinsic::ieee_arithmetic,only:ieee_is_finite
    type(s_dg_wpw_volume_halo_state),intent(in)::halo
    integer,intent(in)::w_id,grid(3),epoch
    complex(8),intent(out)::value,gradient(3)
    integer,intent(out)::info
    integer::location(1),extent(3),point
    value=(0d0,0d0);gradient=(0d0,0d0);info=1
    if(.not.halo%valid.or.epoch/=halo%epoch.or.any(grid<halo%box_lo).or.any(grid>halo%box_hi))return
    location=findloc(halo%w_ids,w_id)
    if(location(1)==0)return
    extent=halo%box_hi-halo%box_lo+1
    point=grid(1)-halo%box_lo(1)+1+(grid(2)-halo%box_lo(2))*extent(1)+&
      (grid(3)-halo%box_lo(3))*extent(1)*extent(2)
    value=halo%values(location(1),point);gradient=halo%gradients(:,location(1),point)
    if(.not.ieee_is_finite(real(value,8)).or..not.ieee_is_finite(aimag(value)))return
    info=0
  end subroutine

  logical function strictly_increasing(ids)result(ok)
    integer,intent(in)::ids(:);integer::i
    ok=size(ids)>0
    do i=2,size(ids);if(ids(i)<=ids(i-1))then;ok=.false.;return;endif;enddo
  end function

  elemental logical function finite_complex(value)result(ok)
    use,intrinsic::ieee_arithmetic,only:ieee_is_finite
    complex(8),intent(in)::value
    ok=ieee_is_finite(real(value,8)).and.ieee_is_finite(aimag(value))
  end function

  pure integer function encode_image_shift(shift)result(tag_offset)
    integer,intent(in)::shift(3)
    tag_offset=(shift(1)+1)*9+(shift(2)+1)*3+shift(3)+1
  end function encode_image_shift

  integer(int64) function dg_wpw_volume_payload_checksum(w_ids,box_lo,box_hi,values,gradients)result(checksum)
    integer,intent(in)::w_ids(:),box_lo(3),box_hi(3)
    complex(8),intent(in)::values(:,:),gradients(:,:,:)
    integer::i
    integer(int64)::bits
    complex(8),allocatable::flat_values(:),flat_gradients(:)
    checksum=int(z'6A09E667F3BCC909',int64)
    do i=1,size(w_ids);checksum=ieor(ishftc(checksum,7),int(w_ids(i),int64));enddo
    do i=1,3
      checksum=ieor(ishftc(checksum,7),int(box_lo(i),int64))
      checksum=ieor(ishftc(checksum,7),int(box_hi(i),int64))
    enddo
    flat_values=reshape(values,[size(values)])
    flat_gradients=reshape(gradients,[size(gradients)])
    do i=1,size(flat_values)
      bits=transfer(real(flat_values(i),8),bits);checksum=ieor(ishftc(checksum,7),bits)
      bits=transfer(aimag(flat_values(i)),bits);checksum=ieor(ishftc(checksum,7),bits)
    enddo
    do i=1,size(flat_gradients)
      bits=transfer(real(flat_gradients(i),8),bits);checksum=ieor(ishftc(checksum,7),bits)
      bits=transfer(aimag(flat_gradients(i)),bits);checksum=ieor(ishftc(checksum,7),bits)
    enddo
  end function dg_wpw_volume_payload_checksum
end module rt_dg_wpw_volume_halo_provider
