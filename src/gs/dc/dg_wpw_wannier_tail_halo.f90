module dg_wpw_wannier_tail_halo
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mpi
  implicit none
  private
  public::validate_sawf_wannier_tail_schedule
  public::exchange_sawf_wannier_tail_values
  public::exchange_sawf_discovered_wannier_tails
  public::locate_sawf_wannier_tail_core
  public::locate_sawf_wannier_tail_rank
  public::qualify_sawf_wannier_buffer_tail
  public::is_sawf_outer_buffer_shell
contains
  logical function is_sawf_outer_buffer_shell(point,box_shape,total_shape)result(is_outer)
    integer,intent(in)::point(3),box_shape(3),total_shape(3)
    integer::axis
    is_outer=.false.
    if(any(box_shape<=0).or.any(total_shape<=0).or.any(box_shape>total_shape).or.&
      any(point<1).or.any(point>box_shape))return
    do axis=1,3
      if(box_shape(axis)<total_shape(axis).and.&
        (point(axis)==1.or.point(axis)==box_shape(axis)))then
        is_outer=.true.;return
      endif
    enddo
  end function is_sawf_outer_buffer_shell

  subroutine qualify_sawf_wannier_buffer_tail(total_norm,omitted_norm,outer_shell_norm,&
      omitted_charge_bound,tolerance,outer_shell_ratio,omitted_tail_ratio,info)
    real(8),intent(in)::total_norm,omitted_norm,outer_shell_norm,omitted_charge_bound,tolerance
    real(8),intent(out)::outer_shell_ratio,omitted_tail_ratio
    integer,intent(out)::info

    info=1;outer_shell_ratio=huge(1d0);omitted_tail_ratio=huge(1d0)
    if(.not.all(ieee_is_finite([total_norm,omitted_norm,outer_shell_norm,&
      omitted_charge_bound,tolerance])).or.total_norm<=0d0.or.omitted_norm<0d0.or.&
      outer_shell_norm<0d0.or.omitted_charge_bound<0d0.or.tolerance<=0d0)return
    outer_shell_ratio=outer_shell_norm/total_norm
    omitted_tail_ratio=omitted_norm/total_norm
    if(max(outer_shell_ratio,max(omitted_tail_ratio,omitted_charge_bound))>tolerance)return
    info=0
  end subroutine qualify_sawf_wannier_buffer_tail

  subroutine locate_sawf_wannier_tail_rank(destination_fragment,destination_local,rank_fragment,&
      rank_orbital_lane,rank_grid_lo,rank_grid_hi,destination_rank,info)
    integer,intent(in)::destination_fragment,destination_local(3),rank_fragment(:),&
      rank_orbital_lane(:),rank_grid_lo(:,:),rank_grid_hi(:,:)
    integer,intent(out)::destination_rank,info
    integer::candidate,count_owner

    info=1;destination_rank=-1;count_owner=0
    if(destination_fragment<=0.or.any(destination_local<=0).or.size(rank_fragment)<=0.or.&
      size(rank_orbital_lane)/=size(rank_fragment).or.any(shape(rank_grid_lo)/=[3,size(rank_fragment)]).or.&
      any(shape(rank_grid_hi)/=[3,size(rank_fragment)]).or.any(rank_grid_lo<=0).or.&
      any(rank_grid_hi<rank_grid_lo))return
    do candidate=1,size(rank_fragment)
      if(rank_fragment(candidate)==destination_fragment.and.rank_orbital_lane(candidate)==0.and.&
        all(destination_local>=rank_grid_lo(:,candidate)).and.&
        all(destination_local<=rank_grid_hi(:,candidate)))then
        count_owner=count_owner+1;destination_rank=candidate-1
      endif
    enddo
    if(count_owner/=1)then;destination_rank=-1;return;endif
    info=0
  end subroutine locate_sawf_wannier_tail_rank

  subroutine locate_sawf_wannier_tail_core(global_point,total_shape,fragment_start,fragment_shape,&
      destination_fragment,destination_local,info)
    integer,intent(in)::global_point(3),total_shape(3),fragment_start(:,:),fragment_shape(:,:)
    integer,intent(out)::destination_fragment,destination_local(3),info
    integer::wrapped(3),local(3),fragment,owner_count

    info=1;destination_fragment=0;destination_local=0;owner_count=0
    if(any(total_shape<=0).or.size(fragment_start,1)/=3.or.size(fragment_shape,1)/=3.or.&
      size(fragment_start,2)<=0.or.size(fragment_shape,2)/=size(fragment_start,2).or.&
      any(fragment_start<1).or.any(fragment_start>spread(total_shape,2,size(fragment_start,2))).or.&
      any(fragment_shape<=0).or.any(fragment_shape>spread(total_shape,2,size(fragment_shape,2))))return
    wrapped=modulo(global_point-1,total_shape)+1
    do fragment=1,size(fragment_start,2)
      local=modulo(wrapped-fragment_start(:,fragment),total_shape)+1
      if(all(local<=fragment_shape(:,fragment)))then
        owner_count=owner_count+1;destination_fragment=fragment;destination_local=local
      endif
    enddo
    if(owner_count/=1)then;destination_fragment=0;destination_local=0;return;endif
    info=0
  end subroutine locate_sawf_wannier_tail_core

  subroutine exchange_sawf_discovered_wannier_tails(comm,send_source,send_destination,send_point,&
      send_value,received_source,received_point,received_value,info)
    integer,intent(in)::comm,send_source(:,:),send_destination(:),send_point(:)
    complex(8),intent(in)::send_value(:)
    integer,allocatable,intent(out)::received_source(:,:),received_point(:)
    complex(8),allocatable,intent(out)::received_value(:)
    integer,intent(out)::info
    integer::rank,nrank,i,j,p,position,total_send,total_recv,ierr,local_bad,global_bad
    integer,allocatable::send_count(:),recv_count(:),send_displ(:),recv_displ(:),cursor(:),&
      packed_source(:),transport_source(:),packed_point(:)
    complex(8),allocatable::packed_value(:)

    info=1;local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)return
    if(size(send_source,1)/=5.or.size(send_source,2)/=size(send_destination).or.&
      size(send_point)/=size(send_destination).or.size(send_value)/=size(send_destination).or.&
      any(send_destination<0).or.any(send_destination>=nrank).or.any(send_point<=0).or.&
      .not.all(ieee_is_finite(real(send_value,8))).or.&
      .not.all(ieee_is_finite(aimag(send_value))))local_bad=1
    if(local_bad/=0)write(*,'(1x,a,i0,5(a,i0),3(a,l1))')&
      '[DG-WPW-TAIL-EXCHANGE-FAIL] validation rank=',rank,&
      ' nsource=',size(send_source,2),' ndestination=',size(send_destination),&
      ' npoint=',size(send_point),' nvalue=',size(send_value),' source_rows=',size(send_source,1),&
      ' bad_destination=',any(send_destination<0).or.any(send_destination>=nrank),&
      ' bad_point=',any(send_point<=0),' bad_value=',&
      .not.all(ieee_is_finite(real(send_value,8))).or..not.all(ieee_is_finite(aimag(send_value)))
    if(local_bad==0)then
      do i=1,size(send_destination);do j=i+1,size(send_destination)
        if(send_destination(i)==send_destination(j).and.send_point(i)==send_point(j).and.&
          all(send_source(:,i)==send_source(:,j)))then
          write(*,'(1x,a,i0,a,2(i0,1x),a,i0,a,i0)')&
            '[DG-WPW-TAIL-EXCHANGE-FAIL] duplicate_send rank=',rank,' records=',i,j,&
            ' destination=',send_destination(i),' point=',send_point(i)
          local_bad=1
        endif
      enddo;enddo
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    allocate(send_count(nrank),recv_count(nrank),send_displ(nrank),recv_displ(nrank),cursor(nrank))
    send_count=0
    do i=1,size(send_destination)
      send_count(send_destination(i)+1)=send_count(send_destination(i)+1)+1
    enddo
    call MPI_Alltoall(send_count,1,MPI_INTEGER,recv_count,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      write(*,'(1x,a,i0,a,i0)')'[DG-WPW-TAIL-EXCHANGE-FAIL] count_exchange rank=',rank,' mpi=',ierr
      return
    endif
    send_displ=0;recv_displ=0
    do p=2,nrank
      send_displ(p)=send_displ(p-1)+send_count(p-1)
      recv_displ(p)=recv_displ(p-1)+recv_count(p-1)
    enddo
    total_send=sum(send_count);total_recv=sum(recv_count);cursor=send_displ
    allocate(packed_source(5*total_send),transport_source(5*total_recv),packed_point(total_send),&
      packed_value(total_send),received_source(5,total_recv),received_point(total_recv),&
      received_value(total_recv))
    do i=1,size(send_destination)
      p=send_destination(i)+1;position=cursor(p)+1;cursor(p)=cursor(p)+1
      packed_source(5*(position-1)+1:5*position)=send_source(:,i)
      packed_point(position)=send_point(i);packed_value(position)=send_value(i)
    enddo
    call MPI_Alltoallv(packed_source,5*send_count,5*send_displ,MPI_INTEGER,&
      transport_source,5*recv_count,5*recv_displ,MPI_INTEGER,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Alltoallv(packed_point,send_count,send_displ,MPI_INTEGER,&
      received_point,recv_count,recv_displ,MPI_INTEGER,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Alltoallv(packed_value,send_count,send_displ,MPI_DOUBLE_COMPLEX,&
      received_value,recv_count,recv_displ,MPI_DOUBLE_COMPLEX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      write(*,'(1x,a,i0,a,i0)')'[DG-WPW-TAIL-EXCHANGE-FAIL] payload_exchange rank=',rank,' mpi=',ierr
      return
    endif
    if(total_recv>0)received_source=reshape(transport_source,[5,total_recv])
    local_bad=0
    do i=1,total_recv;do j=i+1,total_recv
      if(received_point(i)==received_point(j).and.&
        all(received_source(:,i)==received_source(:,j)))then
        write(*,'(1x,a,i0,a,2(i0,1x),a,i0)')&
          '[DG-WPW-TAIL-EXCHANGE-FAIL] duplicate_receive rank=',rank,' records=',i,j,&
          ' point=',received_point(i)
        local_bad=1
      endif
    enddo;enddo
    if(.not.all(ieee_is_finite(real(received_value,8))).or.&
      .not.all(ieee_is_finite(aimag(received_value))))then
      write(*,'(1x,a,i0)')'[DG-WPW-TAIL-EXCHANGE-FAIL] nonfinite_receive rank=',rank
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS.and.global_bad==0)then
      info=0
    else
      write(*,'(1x,a,i0,a,i0,a,i0)')'[DG-WPW-TAIL-EXCHANGE-FAIL] final_reduce rank=',rank,&
        ' global_bad=',global_bad,' mpi=',ierr
    endif
  end subroutine exchange_sawf_discovered_wannier_tails

  subroutine validate_sawf_wannier_tail_schedule(comm,send_source,send_destination,send_point,&
      expected_source,expected_point,info)
    integer,intent(in)::comm
    integer,intent(in)::send_source(:,:),send_destination(:),send_point(:)
    integer,intent(in)::expected_source(:,:),expected_point(:)
    integer,intent(out)::info
    integer::rank,nrank,i,j,ierr,local_bad,global_bad,total_send,total_recv,position
    integer,allocatable::send_count(:),recv_count(:),send_displ(:),recv_displ(:),cursor(:)
    integer,allocatable::packed_send(:),packed_recv(:)
    logical,allocatable::matched(:)

    info=1;local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)return
    if(size(send_source,1)/=5.or.size(expected_source,1)/=5.or.&
        size(send_source,2)/=size(send_destination).or.size(send_point)/=size(send_destination).or.&
        size(expected_source,2)/=size(expected_point).or.any(send_destination<0).or.&
        any(send_destination>=nrank).or.any(send_point<=0).or.any(expected_point<=0))local_bad=1
    if(local_bad==0)then
      do i=1,size(send_destination);do j=i+1,size(send_destination)
        if(send_destination(i)==send_destination(j).and.send_point(i)==send_point(j).and.&
          all(send_source(:,i)==send_source(:,j)))local_bad=1
      enddo;enddo
      do i=1,size(expected_point);do j=i+1,size(expected_point)
        if(expected_point(i)==expected_point(j).and.&
          all(expected_source(:,i)==expected_source(:,j)))local_bad=1
      enddo;enddo
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return

    allocate(send_count(nrank),recv_count(nrank),send_displ(nrank),recv_displ(nrank),cursor(nrank))
    send_count=0
    do i=1,size(send_destination);send_count(send_destination(i)+1)=send_count(send_destination(i)+1)+1;enddo
    call MPI_Alltoall(send_count,1,MPI_INTEGER,recv_count,1,MPI_INTEGER,comm,ierr);if(ierr/=MPI_SUCCESS)return
    send_displ=0;recv_displ=0
    do i=2,nrank
      send_displ(i)=send_displ(i-1)+send_count(i-1)
      recv_displ(i)=recv_displ(i-1)+recv_count(i-1)
    enddo
    total_send=sum(send_count);total_recv=sum(recv_count)
    allocate(packed_send(6*total_send),packed_recv(6*total_recv));cursor=send_displ
    do i=1,size(send_destination)
      j=send_destination(i)+1;position=cursor(j)+1;cursor(j)=cursor(j)+1
      packed_send(6*(position-1)+1:6*position-1)=send_source(:,i)
      packed_send(6*position)=send_point(i)
    enddo
    call MPI_Alltoallv(packed_send,6*send_count,6*send_displ,MPI_INTEGER,&
      packed_recv,6*recv_count,6*recv_displ,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    local_bad=merge(0,1,total_recv==size(expected_point));allocate(matched(size(expected_point)));matched=.false.
    if(local_bad==0)then
      do i=1,total_recv
        position=0
        do j=1,size(expected_point)
          if(.not.matched(j).and.expected_point(j)==packed_recv(6*i).and.&
            all(expected_source(:,j)==packed_recv(6*(i-1)+1:6*i-1)))then;position=j;exit;endif
        enddo
        if(position==0)then;local_bad=1;exit;endif
        matched(position)=.true.
      enddo
      if(.not.all(matched))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS.and.global_bad==0)info=0
  end subroutine validate_sawf_wannier_tail_schedule

  subroutine exchange_sawf_wannier_tail_values(comm,send_source,send_destination,send_point,send_value,&
      expected_source,expected_point,received_value,info)
    integer,intent(in)::comm
    integer,intent(in)::send_source(:,:),send_destination(:),send_point(:)
    complex(8),intent(in)::send_value(:)
    integer,intent(in)::expected_source(:,:),expected_point(:)
    complex(8),intent(out)::received_value(:)
    integer,intent(out)::info
    integer::nrank,i,j,ierr,total_send,total_recv,position,local_bad,global_bad
    integer,allocatable::send_count(:),recv_count(:),send_displ(:),recv_displ(:),cursor(:)
    integer,allocatable::packed_source(:),received_source(:)
    complex(8),allocatable::packed_value(:),transport_value(:)

    info=1;received_value=(0d0,0d0)
    if(size(send_value)/=size(send_destination).or.size(received_value)/=size(expected_point))return
    call validate_sawf_wannier_tail_schedule(comm,send_source,send_destination,send_point,&
      expected_source,expected_point,info)
    if(info/=0)return
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)then;info=1;return;endif
    allocate(send_count(nrank),recv_count(nrank),send_displ(nrank),recv_displ(nrank),cursor(nrank))
    send_count=0
    do i=1,size(send_destination);send_count(send_destination(i)+1)=send_count(send_destination(i)+1)+1;enddo
    call MPI_Alltoall(send_count,1,MPI_INTEGER,recv_count,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;info=1;return;endif
    send_displ=0;recv_displ=0
    do i=2,nrank
      send_displ(i)=send_displ(i-1)+send_count(i-1)
      recv_displ(i)=recv_displ(i-1)+recv_count(i-1)
    enddo
    total_send=sum(send_count);total_recv=sum(recv_count);cursor=send_displ
    allocate(packed_source(6*total_send),received_source(6*total_recv),&
      packed_value(total_send),transport_value(total_recv))
    do i=1,size(send_destination)
      j=send_destination(i)+1;position=cursor(j)+1;cursor(j)=cursor(j)+1
      packed_source(6*(position-1)+1:6*position-1)=send_source(:,i)
      packed_source(6*position)=send_point(i);packed_value(position)=send_value(i)
    enddo
    call MPI_Alltoallv(packed_source,6*send_count,6*send_displ,MPI_INTEGER,&
      received_source,6*recv_count,6*recv_displ,MPI_INTEGER,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Alltoallv(packed_value,send_count,send_displ,MPI_DOUBLE_COMPLEX,&
      transport_value,recv_count,recv_displ,MPI_DOUBLE_COMPLEX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;info=1;return;endif
    local_bad=0
    do i=1,total_recv
      position=0
      do j=1,size(expected_point)
        if(expected_point(j)==received_source(6*i).and.&
          all(expected_source(:,j)==received_source(6*(i-1)+1:6*i-1)))then;position=j;exit;endif
      enddo
      if(position==0)then;local_bad=1;exit;endif
      received_value(position)=transport_value(i)
    enddo
    if(.not.all(ieee_is_finite(real(received_value,8))).or.&
      .not.all(ieee_is_finite(aimag(received_value))))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS.and.global_bad==0)then;info=0;else;info=1;endif
  end subroutine exchange_sawf_wannier_tail_values
end module dg_wpw_wannier_tail_halo
