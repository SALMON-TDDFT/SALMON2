program test_lcfo_halo_direction_tag_mpi
  use mpi_f08
#ifndef LEGACY_TAG
  use dc_fragment_geometry, only: fragment_direction_tag
#endif
  implicit none
  integer :: ierr, rank, nfrag, destination, source, fragment, source_fragment, h, tag, expected
  integer :: direction(3,2), send_value(2), recv_value(2)
  type(MPI_Request) :: send_request(2), recv_request(2)
#ifndef LEGACY_TAG
  integer :: sender, dx, dy, dz, max_tag
  integer(kind=MPI_ADDRESS_KIND) :: mpi_tag_upper
  logical :: flag, valid
  logical, allocatable :: seen(:)
#endif

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nfrag,ierr)
  if(nfrag/=2.and.nfrag/=3)error stop 'fixture requires two or three ranks'
  fragment=rank+1
  direction(:,1)=[-1,0,0];direction(:,2)=[1,0,0]

  do h=1,2
    source=modulo(rank-direction(1,h),nfrag)
    source_fragment=source+1
    send_value(h)=100*fragment+h
    call make_tag(source_fragment,direction(:,h),tag)
    call MPI_Irecv(recv_value(h),1,MPI_INTEGER,source,tag,MPI_COMM_WORLD,recv_request(h),ierr)
  end do
  do h=1,2
    destination=modulo(rank+direction(1,h),nfrag)
    call make_tag(fragment,direction(:,h),tag)
    call MPI_Isend(send_value(h),1,MPI_INTEGER,destination,tag,MPI_COMM_WORLD,send_request(h),ierr)
  end do
  call MPI_Waitall(2,recv_request,MPI_STATUSES_IGNORE,ierr)
  call MPI_Waitall(2,send_request,MPI_STATUSES_IGNORE,ierr)
  do h=1,2
    source=modulo(rank-direction(1,h),nfrag)
    expected=100*(source+1)+h
    if(recv_value(h)/=expected)then
      write(*,'(a,i0,a,i0,a,i0,a,i0)')'rank ',rank,' direction slot ',h,&
        ' received ',recv_value(h),' expected same-dvec ',expected
      error stop 'LCFO halo direction payload did not retain its dvec'
    end if
  end do

#ifndef LEGACY_TAG
  call MPI_Comm_get_attr(MPI_COMM_WORLD,MPI_TAG_UB,mpi_tag_upper,flag,ierr)
  if(.not.flag)error stop 'MPI_TAG_UB is unavailable'
  allocate(seen(0:27000))
  seen=.false.;max_tag=-1
  do sender=1,1000
    do dx=-1,1;do dy=-1,1;do dz=-1,1
      if(dx==0.and.dy==0.and.dz==0)cycle
      call fragment_direction_tag(sender,1000,[dx,dy,dz],tag,valid,int(mpi_tag_upper))
      if(.not.valid.or.tag<0.or.int(tag,MPI_ADDRESS_KIND)>mpi_tag_upper) &
        error stop 'realistic fragment tag exceeds MPI_TAG_UB'
      if(tag>ubound(seen,1))error stop 'unexpected realistic tag range'
      if(seen(tag))error stop 'fragment/direction tag collision'
      seen(tag)=.true.;max_tag=max(max_tag,tag)
    end do;end do;end do
  end do
  deallocate(seen)
  if(rank==0)write(*,'(a,i0,a,i0)')'MPI_TAG_UB=',mpi_tag_upper,' maximum tested tag=',max_tag
#endif

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if(rank==0)write(*,'(a,i0,a)')'PASS LCFO same-direction tags on ',nfrag,' MPI ranks'
  call MPI_Finalize(ierr)

contains

  subroutine make_tag(sender,dvec,result_tag)
    integer,intent(in)::sender,dvec(3)
    integer,intent(out)::result_tag
#ifndef LEGACY_TAG
    logical :: tag_valid
    call fragment_direction_tag(sender,nfrag,dvec,result_tag,tag_valid)
    if(.not.tag_valid)error stop 'production tag construction failed'
#else
    result_tag=sender
#endif
  end subroutine make_tag

end program test_lcfo_halo_direction_tag_mpi
