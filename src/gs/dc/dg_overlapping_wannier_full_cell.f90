#include "config.h"
module dg_overlapping_wannier_full_cell
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::project_dg_full_cell_hamiltonian_tiles
  public::s_dg_full_cell_redistribution_schedule,initialize_dg_full_cell_redistribution,&
    apply_dg_full_cell_redistribution_forward,apply_dg_full_cell_redistribution_reverse,&
    clear_dg_full_cell_redistribution

  type s_dg_full_cell_redistribution_schedule
    logical::initialized=.false.
    integer::comm=0,nproc=0,global_count=0,source_count=0,destination_count=0
    integer,allocatable::send_counts(:),send_displs(:),recv_counts(:),recv_displs(:)
    integer,allocatable::source_pack_positions(:),destination_unpack_positions(:)
    integer(int64)::source_fingerprint=0_int64,destination_fingerprint=0_int64
    integer(int64)::workspace_bytes=0_int64
  end type

  abstract interface
    subroutine dg_full_cell_tile_operator(tile_in,tile_out,ok)
      import real64
      complex(real64),intent(in)::tile_in(:,:)
      complex(real64),intent(out)::tile_out(:,:)
      logical,intent(out)::ok
    end subroutine
  end interface
contains
  subroutine clear_dg_full_cell_redistribution(schedule)
    type(s_dg_full_cell_redistribution_schedule),intent(inout)::schedule
    if(allocated(schedule%send_counts))deallocate(schedule%send_counts)
    if(allocated(schedule%send_displs))deallocate(schedule%send_displs)
    if(allocated(schedule%recv_counts))deallocate(schedule%recv_counts)
    if(allocated(schedule%recv_displs))deallocate(schedule%recv_displs)
    if(allocated(schedule%source_pack_positions))deallocate(schedule%source_pack_positions)
    if(allocated(schedule%destination_unpack_positions))deallocate(schedule%destination_unpack_positions)
    schedule%initialized=.false.;schedule%comm=0;schedule%nproc=0;schedule%global_count=0
    schedule%source_count=0;schedule%destination_count=0
    schedule%source_fingerprint=0_int64;schedule%destination_fingerprint=0_int64
    schedule%workspace_bytes=0_int64
  end subroutine

  subroutine initialize_dg_full_cell_redistribution(comm,global_count,source_ids,destination_ids,&
      schedule,workspace_bytes,ok,message)
    integer,intent(in)::comm,global_count
    integer(int64),intent(in)::source_ids(:),destination_ids(:)
    type(s_dg_full_cell_redistribution_schedule),intent(inout)::schedule
    integer(int64),intent(out)::workspace_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,local_bad,global_bad,allocation_status,p,r,id_index,broker,bucket_count
    integer::source_total,destination_total,record_count,response_count,target,slot
    integer,allocatable::source_route_counts(:),source_route_displs(:),destination_route_counts(:),&
      destination_route_displs(:),source_record_counts(:),source_record_displs(:),&
      destination_record_counts(:),destination_record_displs(:),cursor(:),response_counts(:),&
      response_displs(:),response_record_counts(:),response_record_displs(:),target_owner(:),&
      target_position(:),seen(:)
    integer(int64),allocatable::source_send(:),source_recv(:),destination_send(:),destination_recv(:),&
      response_send(:),response_recv(:),send_destination_positions(:),recv_destination_positions(:)
    integer,allocatable::bucket_source_owner(:),bucket_source_position(:),bucket_destination_owner(:),&
      bucket_destination_position(:)
    integer(int64)::id,integer_elements
    call clear_dg_full_cell_redistribution(schedule)
    ok=.false.;message='';workspace_bytes=0_int64;local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)then;message='redistribution rank query failed';return;endif
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)then;message='redistribution size query failed';return;endif
    if(global_count<=0)local_bad=1
    if(any(source_ids<1_int64).or.any(source_ids>int(global_count,int64)))local_bad=1
    if(any(destination_ids<1_int64).or.any(destination_ids>int(global_count,int64)))local_bad=1
    call MPI_Allreduce(size(source_ids),source_total,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(size(destination_ids),destination_total,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    if(source_total/=global_count.or.destination_total/=global_count)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid full-cell redistribution ID extent';return;endif

    allocate(source_route_counts(0:nproc-1),source_route_displs(0:nproc-1),&
      destination_route_counts(0:nproc-1),destination_route_displs(0:nproc-1),&
      source_record_counts(0:nproc-1),source_record_displs(0:nproc-1),&
      destination_record_counts(0:nproc-1),destination_record_displs(0:nproc-1),cursor(0:nproc-1),&
      response_counts(0:nproc-1),response_displs(0:nproc-1),response_record_counts(0:nproc-1),&
      response_record_displs(0:nproc-1),schedule%send_counts(0:nproc-1),schedule%send_displs(0:nproc-1),&
      schedule%recv_counts(0:nproc-1),schedule%recv_displs(0:nproc-1),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call clear_dg_full_cell_redistribution(schedule);message='cannot allocate redistribution count workspace';return
    endif
    source_route_counts=0;destination_route_counts=0
    do p=1,size(source_ids);broker=int(modulo(source_ids(p)-1_int64,int(nproc,int64)));source_route_counts(broker)=source_route_counts(broker)+1;enddo
    do p=1,size(destination_ids);broker=int(modulo(destination_ids(p)-1_int64,int(nproc,int64)));destination_route_counts(broker)=destination_route_counts(broker)+1;enddo
    call prefix_counts(source_route_counts,source_route_displs)
    call prefix_counts(destination_route_counts,destination_route_displs)
    source_record_counts=3*source_route_counts;source_record_displs=3*source_route_displs
    destination_record_counts=3*destination_route_counts;destination_record_displs=3*destination_route_displs
    allocate(source_send(3*size(source_ids)),destination_send(3*size(destination_ids)),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call clear_dg_full_cell_redistribution(schedule);message='cannot allocate redistribution route workspace';return
    endif
    cursor=source_route_displs
    do p=1,size(source_ids)
      broker=int(modulo(source_ids(p)-1_int64,int(nproc,int64)));slot=cursor(broker);cursor(broker)=slot+1
      source_send(3*slot+1:3*slot+3)=[source_ids(p),int(rank,int64),int(p,int64)]
    enddo
    cursor=destination_route_displs
    do p=1,size(destination_ids)
      broker=int(modulo(destination_ids(p)-1_int64,int(nproc,int64)));slot=cursor(broker);cursor(broker)=slot+1
      destination_send(3*slot+1:3*slot+3)=[destination_ids(p),int(rank,int64),int(p,int64)]
    enddo
    call MPI_Alltoall(source_route_counts,1,MPI_INTEGER,response_counts,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    source_record_counts=3*source_route_counts;source_record_displs=3*source_route_displs
    response_record_counts=3*response_counts;call prefix_counts(response_counts,response_displs)
    response_record_displs=3*response_displs;record_count=sum(response_counts)
    allocate(source_recv(3*record_count),stat=allocation_status);if(allocation_status/=0)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='cannot allocate routed source records';return;endif
    call MPI_Alltoallv(source_send,source_record_counts,source_record_displs,MPI_INTEGER8,&
      source_recv,response_record_counts,response_record_displs,MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1

    call MPI_Alltoall(destination_route_counts,1,MPI_INTEGER,response_counts,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    response_record_counts=3*response_counts;call prefix_counts(response_counts,response_displs)
    response_record_displs=3*response_displs;response_count=sum(response_counts)
    allocate(destination_recv(3*response_count),stat=allocation_status);if(allocation_status/=0)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='cannot allocate routed destination records';return;endif
    call MPI_Alltoallv(destination_send,destination_record_counts,destination_record_displs,MPI_INTEGER8,&
      destination_recv,response_record_counts,response_record_displs,MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1

    bucket_count=(global_count+nproc-1-rank)/nproc
    allocate(bucket_source_owner(bucket_count),bucket_source_position(bucket_count),&
      bucket_destination_owner(bucket_count),bucket_destination_position(bucket_count),stat=allocation_status)
    if(allocation_status/=0)local_bad=1
    if(local_bad==0)then
      bucket_source_owner=-1;bucket_source_position=0;bucket_destination_owner=-1;bucket_destination_position=0
      do p=1,record_count
        id=source_recv(3*p-2);id_index=int((id-1_int64)/int(nproc,int64))+1
        if(id_index<1.or.id_index>bucket_count.or.bucket_source_owner(id_index)/=-1)then
          local_bad=1
        else
          bucket_source_owner(id_index)=int(source_recv(3*p-1));bucket_source_position(id_index)=int(source_recv(3*p))
        endif
      enddo
      do p=1,response_count
        id=destination_recv(3*p-2);id_index=int((id-1_int64)/int(nproc,int64))+1
        if(id_index<1.or.id_index>bucket_count.or.bucket_destination_owner(id_index)/=-1)then
          local_bad=1
        else
          bucket_destination_owner(id_index)=int(destination_recv(3*p-1));&
          bucket_destination_position(id_index)=int(destination_recv(3*p))
        endif
      enddo
      if(any(bucket_source_owner<0).or.any(bucket_destination_owner<0))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='duplicate, missing, or mismatched redistribution physical IDs';return;endif

    response_counts=0
    do id_index=1,bucket_count
      response_counts(bucket_source_owner(id_index))=response_counts(bucket_source_owner(id_index))+1
    enddo
    call prefix_counts(response_counts,response_displs);response_record_counts=3*response_counts
    response_record_displs=3*response_displs
    allocate(response_send(3*bucket_count),stat=allocation_status);if(allocation_status/=0)local_bad=1
    cursor=response_displs
    if(local_bad==0)then
      do id_index=1,bucket_count
        target=bucket_source_owner(id_index);slot=cursor(target);cursor(target)=slot+1
        response_send(3*slot+1:3*slot+3)=[int(bucket_source_position(id_index),int64),&
          int(bucket_destination_owner(id_index),int64),int(bucket_destination_position(id_index),int64)]
      enddo
    endif
    call MPI_Alltoall(response_counts,1,MPI_INTEGER,source_route_counts,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call prefix_counts(source_route_counts,source_route_displs);source_record_counts=3*source_route_counts
    source_record_displs=3*source_route_displs
    allocate(response_recv(3*sum(source_route_counts)),target_owner(size(source_ids)),&
      target_position(size(source_ids)),seen(size(source_ids)),stat=allocation_status)
    if(allocation_status/=0)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='cannot allocate redistribution match response';return;endif
    call MPI_Alltoallv(response_send,response_record_counts,response_record_displs,MPI_INTEGER8,&
      response_recv,source_record_counts,source_record_displs,MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    seen=0;target_owner=-1;target_position=0
    do p=1,sum(source_route_counts)
      slot=int(response_recv(3*p-2))
      if(slot<1.or.slot>size(source_ids).or.seen(slot)/=0)then;local_bad=1;cycle;endif
      seen(slot)=1;target_owner(slot)=int(response_recv(3*p-1));target_position(slot)=int(response_recv(3*p))
    enddo
    if(any(seen/=1).or.any(target_owner<0).or.any(target_owner>=nproc))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid redistribution owner response';return;endif

    schedule%send_counts=0
    do p=1,size(source_ids);schedule%send_counts(target_owner(p))=schedule%send_counts(target_owner(p))+1;enddo
    call prefix_counts(schedule%send_counts,schedule%send_displs)
    call MPI_Alltoall(schedule%send_counts,1,MPI_INTEGER,schedule%recv_counts,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='redistribution value-count exchange failed';return;endif
    call prefix_counts(schedule%recv_counts,schedule%recv_displs)
    deallocate(seen)
    allocate(schedule%source_pack_positions(size(source_ids)),send_destination_positions(size(source_ids)),&
      recv_destination_positions(size(destination_ids)),schedule%destination_unpack_positions(size(destination_ids)),&
      seen(size(destination_ids)),stat=allocation_status)
    if(allocation_status/=0)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='cannot allocate redistribution permutation';return;endif
    cursor=schedule%send_displs
    do p=1,size(source_ids)
      target=target_owner(p);slot=cursor(target)+1;cursor(target)=cursor(target)+1
      schedule%source_pack_positions(slot)=p;send_destination_positions(slot)=int(target_position(p),int64)
    enddo
    call MPI_Alltoallv(send_destination_positions,schedule%send_counts,schedule%send_displs,MPI_INTEGER8,&
      recv_destination_positions,schedule%recv_counts,schedule%recv_displs,MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    seen=0
    do p=1,size(destination_ids)
      slot=int(recv_destination_positions(p))
      if(slot<1.or.slot>size(destination_ids).or.seen(slot)/=0)then;local_bad=1;cycle;endif
      seen(slot)=1;schedule%destination_unpack_positions(p)=slot
    enddo
    if(any(seen/=1))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid redistribution destination permutation';return;endif

    schedule%initialized=.true.;schedule%comm=comm;schedule%nproc=nproc;schedule%global_count=global_count
    schedule%source_count=size(source_ids);schedule%destination_count=size(destination_ids)
    schedule%source_fingerprint=id_fingerprint(source_ids,rank)
    schedule%destination_fingerprint=id_fingerprint(destination_ids,rank)
    call MPI_Allreduce(MPI_IN_PLACE,schedule%source_fingerprint,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,schedule%destination_fingerprint,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
    integer_elements=int(6*nproc+2*size(source_ids)+size(destination_ids),int64)
    schedule%workspace_bytes=4_int64*integer_elements
    workspace_bytes=schedule%workspace_bytes;ok=.true.;message=''
#else
    call clear_dg_full_cell_redistribution(schedule)
    workspace_bytes=0_int64;ok=.false.;message='full-cell redistribution requires MPI'
#endif
  end subroutine

  subroutine apply_dg_full_cell_redistribution_forward(schedule,source_values,destination_values,ok,message)
    type(s_dg_full_cell_redistribution_schedule),intent(in)::schedule
    complex(real64),intent(in)::source_values(:,:)
    complex(real64),intent(out)::destination_values(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::width,p,i,ierr,local_bad,global_bad,allocation_status
    integer,allocatable::send_counts(:),send_displs(:),recv_counts(:),recv_displs(:)
    complex(real64),allocatable::send_values(:),recv_values(:)
    ok=.false.;message='';local_bad=0;width=size(source_values,1)
    if(.not.schedule%initialized.or.width<=0.or.size(source_values,2)/=schedule%source_count.or.&
      size(destination_values,1)/=width.or.size(destination_values,2)/=schedule%destination_count)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid forward redistribution value shape';return;endif
    allocate(send_counts(0:schedule%nproc-1),send_displs(0:schedule%nproc-1),&
      recv_counts(0:schedule%nproc-1),recv_displs(0:schedule%nproc-1),&
      send_values(width*schedule%source_count),recv_values(width*schedule%destination_count),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='cannot allocate forward redistribution values';return;endif
    send_counts=width*schedule%send_counts;send_displs=width*schedule%send_displs
    recv_counts=width*schedule%recv_counts;recv_displs=width*schedule%recv_displs
    do p=1,schedule%source_count;do i=1,width
      send_values(width*(p-1)+i)=source_values(i,schedule%source_pack_positions(p))
    enddo;enddo
    call MPI_Alltoallv(send_values,send_counts,send_displs,MPI_DOUBLE_COMPLEX,&
      recv_values,recv_counts,recv_displs,MPI_DOUBLE_COMPLEX,schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='forward redistribution MPI Alltoallv failed';return;endif
    do p=1,schedule%destination_count;do i=1,width
      destination_values(i,schedule%destination_unpack_positions(p))=recv_values(width*(p-1)+i)
    enddo;enddo
    ok=finite_complex(destination_values);if(ok)then;message='';else;message='nonfinite forward redistribution values';endif
#else
    ok=.false.;message='full-cell redistribution requires MPI'
#endif
  end subroutine

  subroutine apply_dg_full_cell_redistribution_reverse(schedule,destination_values,source_values,ok,message)
    type(s_dg_full_cell_redistribution_schedule),intent(in)::schedule
    complex(real64),intent(in)::destination_values(:,:)
    complex(real64),intent(out)::source_values(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::width,p,i,ierr,local_bad,global_bad,allocation_status
    integer,allocatable::send_counts(:),send_displs(:),recv_counts(:),recv_displs(:)
    complex(real64),allocatable::send_values(:),recv_values(:)
    ok=.false.;message='';local_bad=0;width=size(destination_values,1)
    if(.not.schedule%initialized.or.width<=0.or.size(destination_values,2)/=schedule%destination_count.or.&
      size(source_values,1)/=width.or.size(source_values,2)/=schedule%source_count)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid reverse redistribution value shape';return;endif
    allocate(send_counts(0:schedule%nproc-1),send_displs(0:schedule%nproc-1),&
      recv_counts(0:schedule%nproc-1),recv_displs(0:schedule%nproc-1),&
      send_values(width*schedule%destination_count),recv_values(width*schedule%source_count),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='cannot allocate reverse redistribution values';return;endif
    send_counts=width*schedule%recv_counts;send_displs=width*schedule%recv_displs
    recv_counts=width*schedule%send_counts;recv_displs=width*schedule%send_displs
    do p=1,schedule%destination_count;do i=1,width
      send_values(width*(p-1)+i)=destination_values(i,schedule%destination_unpack_positions(p))
    enddo;enddo
    call MPI_Alltoallv(send_values,send_counts,send_displs,MPI_DOUBLE_COMPLEX,&
      recv_values,recv_counts,recv_displs,MPI_DOUBLE_COMPLEX,schedule%comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='reverse redistribution MPI Alltoallv failed';return;endif
    do p=1,schedule%source_count;do i=1,width
      source_values(i,schedule%source_pack_positions(p))=recv_values(width*(p-1)+i)
    enddo;enddo
    ok=finite_complex(source_values);if(ok)then;message='';else;message='nonfinite reverse redistribution values';endif
#else
    ok=.false.;message='full-cell redistribution requires MPI'
#endif
  end subroutine

#ifdef USE_MPI
  subroutine prefix_counts(counts,displs)
    integer,intent(in)::counts(0:)
    integer,intent(out)::displs(0:)
    integer::r
    displs(0)=0
    do r=1,ubound(counts,1);displs(r)=displs(r-1)+counts(r-1);enddo
  end subroutine

  integer(int64) function id_fingerprint(ids,rank)
    integer(int64),intent(in)::ids(:)
    integer,intent(in)::rank
    integer::p
    id_fingerprint=int(rank+1,int64)*104729_int64
    do p=1,size(ids)
      id_fingerprint=id_fingerprint+ids(p)*int(2*p+1,int64)
    enddo
  end function
#endif

  subroutine project_dg_full_cell_hamiltonian_tiles(comm,global_spatial_count,spatial_ids,weights,&
      basis_values,row_ids,tile_width,apply_tile,matrix_rows,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm,global_spatial_count,tile_width
    integer(int64),intent(in)::spatial_ids(:),row_ids(:)
    real(real64),intent(in)::weights(:)
    complex(real64),intent(in)::basis_values(:,:)
    procedure(dg_full_cell_tile_operator)::apply_tile
    complex(real64),allocatable,intent(out)::matrix_rows(:,:)
    integer(int64),intent(out)::workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::i,j,j0,j1,width,nstate,nlocal,rank,nproc,ierr,local_bad,global_bad,allocation_status
    integer::minimum_value,maximum_value,root,local_position
    integer,allocatable::spatial_count(:),row_count(:),row_owner(:),row_position(:)
    integer(int64)::complex_elements,integer_elements
    complex(real64),allocatable::tile_in(:,:),tile_out(:,:),local_row(:),reduced_row(:)
    logical::callback_ok
    ok=.false.;message='';workspace_peak_bytes=0_int64;local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    nstate=size(basis_values,1);nlocal=size(spatial_ids)
    call agree_integer(global_spatial_count,minimum_value,maximum_value,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_value/=maximum_value)then
      message='inconsistent full-cell spatial extent';return
    endif
    call agree_integer(nstate,minimum_value,maximum_value,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_value/=maximum_value)then
      message='inconsistent full-cell orbital extent';return
    endif
    call agree_integer(tile_width,minimum_value,maximum_value,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_value/=maximum_value)then
      message='inconsistent full-cell tile width';return
    endif
    if(global_spatial_count<=0.or.nstate<=0.or.tile_width<=0)local_bad=1
    if(size(weights)/=nlocal.or.size(basis_values,2)/=nlocal)local_bad=1
    if(any(spatial_ids<1_int64).or.any(spatial_ids>int(global_spatial_count,int64)))local_bad=1
    if(any(row_ids<1_int64).or.any(row_ids>int(nstate,int64)))local_bad=1
    if(.not.all(ieee_is_finite(weights)).or..not.finite_complex(basis_values))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid full-cell tiled projection contract';return
    endif
    allocate(spatial_count(global_spatial_count),row_count(nstate),row_owner(nstate),row_position(nstate),&
      stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(spatial_count))deallocate(spatial_count)
      if(allocated(row_count))deallocate(row_count)
      if(allocated(row_owner))deallocate(row_owner)
      if(allocated(row_position))deallocate(row_position)
      message='cannot allocate full-cell ownership workspace';return
    endif
    spatial_count=0;row_count=0;row_owner=-1;row_position=0
    do i=1,nlocal;spatial_count(int(spatial_ids(i)))=spatial_count(int(spatial_ids(i)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,spatial_count,global_spatial_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    do i=1,size(row_ids)
      row_count(int(row_ids(i)))=row_count(int(row_ids(i)))+1
      row_owner(int(row_ids(i)))=rank;row_position(int(row_ids(i)))=i
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,row_count,nstate,MPI_INTEGER,MPI_SUM,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,row_owner,nstate,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,row_position,nstate,MPI_INTEGER,MPI_MAX,comm,ierr)
    local_bad=merge(0,1,all(spatial_count==1).and.all(row_count==1))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='duplicate or missing full-cell spatial/orbital owner';return
    endif
    ! Peak owned workspace: output rows, two bounded complex tiles, two row buffers,
    ! and four integer ownership arrays.  Accumulate in int64 before allocating.
    complex_elements=int(size(row_ids),int64)*int(nstate,int64)
    if(2_int64*int(tile_width,int64)>huge(complex_elements)/max(1_int64,int(nlocal,int64)))local_bad=1
    if(local_bad==0)complex_elements=complex_elements+2_int64*int(tile_width,int64)*int(nlocal,int64)+&
      2_int64*int(tile_width,int64)
    integer_elements=2_int64*int(global_spatial_count,int64)+2_int64*int(nstate,int64)
    if(complex_elements>huge(workspace_peak_bytes)/16_int64)local_bad=1
    if(integer_elements>huge(workspace_peak_bytes)/4_int64)local_bad=1
    if(local_bad==0.and.16_int64*complex_elements>huge(workspace_peak_bytes)-4_int64*integer_elements)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='full-cell tiled projection workspace overflow';return
    endif
    workspace_peak_bytes=16_int64*complex_elements+4_int64*integer_elements
    allocate(matrix_rows(size(row_ids),nstate),tile_in(tile_width,nlocal),tile_out(tile_width,nlocal),&
      local_row(tile_width),reduced_row(tile_width),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup_outputs();message='cannot allocate full-cell tiled projection arrays';return
    endif
    matrix_rows=(0d0,0d0)
    do j0=1,nstate,tile_width
      j1=min(nstate,j0+tile_width-1);width=j1-j0+1
      tile_in(1:width,:)=basis_values(j0:j1,:)
      call apply_tile(tile_in(1:width,:),tile_out(1:width,:),callback_ok)
      local_bad=merge(0,1,callback_ok.and.finite_complex(tile_out(1:width,:)))
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
        call cleanup_outputs();message='full-cell tile Hamiltonian callback failed';return
      endif
      do i=1,nstate
        do j=1,width
          local_row(j)=sum(weights*conjg(basis_values(i,:))*tile_out(j,:))
        enddo
        root=row_owner(i);reduced_row(1:width)=(0d0,0d0)
        call MPI_Reduce(local_row,reduced_row,width,MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
        if(ierr/=MPI_SUCCESS)then
          call cleanup_outputs();message='full-cell projected-row reduction failed';return
        endif
        if(rank==root)then
          local_position=row_position(i)
          matrix_rows(local_position,j0:j1)=reduced_row(1:width)
        endif
      enddo
    enddo
    ok=.true.
#else
    ok=.false.;message='full-cell tiled projection requires MPI';workspace_peak_bytes=0_int64
#endif
  contains
#ifdef USE_MPI
    subroutine cleanup_outputs()
      if(allocated(matrix_rows))deallocate(matrix_rows)
      if(allocated(tile_in))deallocate(tile_in)
      if(allocated(tile_out))deallocate(tile_out)
      if(allocated(local_row))deallocate(local_row)
      if(allocated(reduced_row))deallocate(reduced_row)
    end subroutine
#endif
  end subroutine

  logical function finite_complex(values)
    complex(real64),intent(in)::values(:,:)
    finite_complex=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function

#ifdef USE_MPI
  subroutine agree_integer(value,minimum_value,maximum_value,comm,ierr)
    integer,intent(in)::value,comm
    integer,intent(out)::minimum_value,maximum_value,ierr
    call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER,MPI_MAX,comm,ierr)
  end subroutine
#endif
end module dg_overlapping_wannier_full_cell
