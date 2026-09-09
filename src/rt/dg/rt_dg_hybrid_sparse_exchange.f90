#include "config.h"
module rt_dg_hybrid_sparse_exchange
  use,intrinsic::iso_fortran_env,only:int64,real64
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  type,public::s_rt_dg_sparse_exchange
    logical::valid=.false.
    integer::comm=0,nproc=0,local_count=0
    integer(int64)::workspace_peak_bytes=0_int64,catalog_fingerprint=0_int64
    integer,allocatable::send_counts(:),send_displacements(:),receive_counts(:),receive_displacements(:)
    integer,allocatable::send_positions(:),value_slots(:)
    complex(real64),allocatable::send_values(:),receive_values(:)
  end type s_rt_dg_sparse_exchange
  public::build_rt_dg_sparse_exchange,exchange_rt_dg_sparse_values,&
    apply_rt_dg_sparse_rows_tiled,accumulate_rt_dg_sparse_density_tiled,clear_rt_dg_sparse_exchange,&
    checked_rt_dg_sparse_extent_product
contains
  subroutine build_rt_dg_sparse_exchange(comm,global_count,catalog_fingerprint,owned_row_ids,needed_ids,plan,ok,message)
    integer,intent(in)::comm,global_count,needed_ids(:)
    integer(int64),intent(in)::catalog_fingerprint
    integer(int64),intent(in)::owned_row_ids(:)
    type(s_rt_dg_sparse_exchange),intent(out)::plan
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,i,j,target,total_send,total_receive,allocation_status,local_bad,global_bad,&
      minimum_count,maximum_count
    integer(int64)::integer_elements,complex_elements,minimum_fingerprint,maximum_fingerprint
    integer,allocatable::counts(:),displacements(:),received_counts(:),received_displacements(:),cursor(:)
    integer,allocatable::send_ids(:),send_positions(:),received_ids(:),received_positions(:)
    integer,allocatable::directory_owner(:),directory_count(:),sorted_ids(:),sorted_positions(:)
    integer,allocatable::query_counts(:),query_displacements(:),query_received_counts(:),query_received_displacements(:)
    integer,allocatable::query_ids(:),query_edges(:),received_queries(:),query_answers(:),received_answers(:),edge_owner(:)
    integer,allocatable::request_ids(:),request_edges(:),received_requests(:)
    integer,allocatable::remote_ids(:),remote_edges(:),edge_unique(:),unique_slots(:)
    integer::remote_count,unique_count
    ok=.false.;message='';call clear_rt_dg_sparse_exchange(plan)
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    if(nproc>huge(0)/2)then;message='sparse exchange communicator extent overflow';return;endif
    call MPI_Allreduce(global_count,minimum_count,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(global_count,maximum_count,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(minimum_count/=maximum_count)then;message='rank-disagreeing sparse global extent';return;endif
    call MPI_Allreduce(catalog_fingerprint,minimum_fingerprint,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(catalog_fingerprint,maximum_fingerprint,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(minimum_fingerprint/=maximum_fingerprint.or.catalog_fingerprint==0_int64)then
      message='rank-disagreeing sparse exchange provenance';return
    endif
    local_bad=0
    if(global_count<1.or.any(owned_row_ids<1_int64).or.any(owned_row_ids>int(global_count,int64)).or.&
      any(needed_ids<1).or.any(needed_ids>global_count))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid sparse exchange contract';return;endif
    allocate(counts(nproc),displacements(nproc),received_counts(nproc),received_displacements(nproc),cursor(nproc),&
      sorted_ids(size(owned_row_ids)),sorted_positions(size(owned_row_ids)),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0);call consensus(local_bad,global_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='cannot allocate sparse ownership counts';return;endif
    counts=0
    do i=1,size(owned_row_ids);sorted_ids(i)=int(owned_row_ids(i));sorted_positions(i)=i;enddo
    if(size(sorted_ids)>1)call sort_pairs(sorted_ids,sorted_positions,1,size(sorted_ids))
    do i=1,size(owned_row_ids);target=mod(int(owned_row_ids(i))-1,nproc)+1;counts(target)=counts(target)+1;enddo
    call offsets(counts,displacements,total_send,local_bad)
    call consensus(local_bad,global_bad,ierr);if(ierr/=MPI_SUCCESS.or.global_bad/=0)goto 900
    call MPI_Alltoall(counts,1,MPI_INTEGER,received_counts,1,MPI_INTEGER,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call offsets(received_counts,received_displacements,total_receive,local_bad)
    call consensus(local_bad,global_bad,ierr);if(ierr/=MPI_SUCCESS.or.global_bad/=0)goto 900
    allocate(send_ids(total_send),send_positions(total_send),received_ids(total_receive),received_positions(total_receive),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0);call consensus(local_bad,global_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='cannot allocate sparse ownership directory';return;endif
    cursor=displacements+1
    do i=1,size(owned_row_ids)
      target=mod(int(owned_row_ids(i))-1,nproc)+1;j=cursor(target);cursor(target)=j+1
      send_ids(j)=int(owned_row_ids(i));send_positions(j)=i
    enddo
    call MPI_Alltoallv(send_ids,counts,displacements,MPI_INTEGER,received_ids,received_counts,received_displacements,&
      MPI_INTEGER,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    ! The source rank of each registration is reconstructed from receive blocks.
    do target=1,nproc
      do j=received_displacements(target)+1,received_displacements(target)+received_counts(target)
        received_positions(j)=target-1
      enddo
    enddo
    if(rank>=global_count)then;i=0;else;i=(global_count-1-rank)/nproc+1;endif
    allocate(directory_owner(i),directory_count(i),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0);call consensus(local_bad,global_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='cannot allocate sparse ownership shard';return;endif
    directory_owner=-1;directory_count=0
    do j=1,total_receive
      if(mod(received_ids(j)-1,nproc)/=rank)then;local_bad=1;cycle;endif
      i=(received_ids(j)-1)/nproc+1
      if(i<1.or.i>size(directory_count))then;local_bad=1;cycle;endif
      directory_count(i)=directory_count(i)+1;directory_owner(i)=received_positions(j)
    enddo
    if(any(directory_count/=1))local_bad=1
    call consensus(local_bad,global_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='sparse row ownership is not exactly once';return;endif

    allocate(plan%value_slots(size(needed_ids)),edge_owner(size(needed_ids)),edge_unique(size(needed_ids)),&
      remote_ids(size(needed_ids)),remote_edges(size(needed_ids)),query_counts(nproc),query_displacements(nproc),&
      query_received_counts(nproc),query_received_displacements(nproc),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0);call consensus(local_bad,global_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='cannot allocate sparse query metadata';return;endif
    plan%value_slots=0;edge_owner=-1;edge_unique=0;query_counts=0;remote_count=0
    do i=1,size(needed_ids)
      j=find_local(needed_ids(i))
      if(j>0)then
        plan%value_slots(i)=j
      else
        remote_count=remote_count+1;remote_ids(remote_count)=needed_ids(i);remote_edges(remote_count)=i
      endif
    enddo
    if(remote_count>1)call sort_pairs(remote_ids,remote_edges,1,remote_count)
    unique_count=0
    do i=1,remote_count
      if(i==1)then
        unique_count=unique_count+1;remote_ids(unique_count)=remote_ids(i)
      else if(remote_ids(i)/=remote_ids(i-1))then
        unique_count=unique_count+1;remote_ids(unique_count)=remote_ids(i)
      endif
      edge_unique(remote_edges(i))=unique_count
    enddo
    do i=1,unique_count;target=mod(remote_ids(i)-1,nproc)+1;query_counts(target)=query_counts(target)+1;enddo
    call offsets(query_counts,query_displacements,total_send,local_bad)
    call consensus(local_bad,global_bad,ierr);if(ierr/=MPI_SUCCESS.or.global_bad/=0)goto 900
    call MPI_Alltoall(query_counts,1,MPI_INTEGER,query_received_counts,1,MPI_INTEGER,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call offsets(query_received_counts,query_received_displacements,total_receive,local_bad)
    call consensus(local_bad,global_bad,ierr);if(ierr/=MPI_SUCCESS.or.global_bad/=0)goto 900
    allocate(query_ids(total_send),query_edges(total_send),received_queries(total_receive),query_answers(total_receive),&
      received_answers(total_send),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0);call consensus(local_bad,global_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='cannot allocate sparse directory queries';return;endif
    cursor=query_displacements+1
    do i=1,unique_count
      target=mod(remote_ids(i)-1,nproc)+1;j=cursor(target);cursor(target)=j+1
      query_ids(j)=remote_ids(i);query_edges(j)=i
    enddo
    call MPI_Alltoallv(query_ids,query_counts,query_displacements,MPI_INTEGER,received_queries,query_received_counts,&
      query_received_displacements,MPI_INTEGER,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    do i=1,total_receive
      j=(received_queries(i)-1)/nproc+1
      if(j<1.or.j>size(directory_owner))then;local_bad=1;query_answers(i)=-1;else;query_answers(i)=directory_owner(j);endif
    enddo
    call consensus(local_bad,global_bad,ierr);if(ierr/=MPI_SUCCESS.or.global_bad/=0)goto 900
    call MPI_Alltoallv(query_answers,query_received_counts,query_received_displacements,MPI_INTEGER,received_answers,&
      query_counts,query_displacements,MPI_INTEGER,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    do i=1,total_send;edge_owner(query_edges(i))=received_answers(i);enddo

    allocate(unique_slots(unique_count),plan%send_counts(nproc),plan%send_displacements(nproc),plan%receive_counts(nproc),&
      plan%receive_displacements(nproc),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0);call consensus(local_bad,global_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='cannot allocate sparse exchange counts';return;endif
    plan%receive_counts=0
    do i=1,unique_count;plan%receive_counts(edge_owner(i)+1)=plan%receive_counts(edge_owner(i)+1)+1;enddo
    call offsets(plan%receive_counts,plan%receive_displacements,total_receive,local_bad)
    call consensus(local_bad,global_bad,ierr);if(ierr/=MPI_SUCCESS.or.global_bad/=0)goto 900
    call MPI_Alltoall(plan%receive_counts,1,MPI_INTEGER,plan%send_counts,1,MPI_INTEGER,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call offsets(plan%send_counts,plan%send_displacements,total_send,local_bad)
    call consensus(local_bad,global_bad,ierr);if(ierr/=MPI_SUCCESS.or.global_bad/=0)goto 900
    allocate(request_ids(total_receive),request_edges(total_receive),received_requests(total_send),&
      plan%send_positions(total_send),plan%send_values(total_send),plan%receive_values(total_receive),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0);call consensus(local_bad,global_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='cannot allocate sparse halo payload';return;endif
    cursor=plan%receive_displacements+1
    do i=1,unique_count
      target=edge_owner(i)+1;j=cursor(target);cursor(target)=j+1
      request_ids(j)=remote_ids(i);request_edges(j)=i
    enddo
    call MPI_Alltoallv(request_ids,plan%receive_counts,plan%receive_displacements,MPI_INTEGER,received_requests,&
      plan%send_counts,plan%send_displacements,MPI_INTEGER,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    do i=1,total_send
      plan%send_positions(i)=find_local(received_requests(i))
      if(plan%send_positions(i)==0)local_bad=1
    enddo
    do i=1,total_receive;unique_slots(request_edges(i))=-i;enddo
    do i=1,size(needed_ids)
      if(plan%value_slots(i)==0)plan%value_slots(i)=unique_slots(edge_unique(i))
    enddo
    call consensus(local_bad,global_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='sparse halo owner lookup failed';return;endif
    integer_elements=0_int64;complex_elements=0_int64
    call add_size(integer_elements,counts,local_bad);call add_size(integer_elements,displacements,local_bad)
    call add_size(integer_elements,received_counts,local_bad);call add_size(integer_elements,received_displacements,local_bad)
    call add_size(integer_elements,cursor,local_bad);call add_size(integer_elements,send_ids,local_bad)
    call add_size(integer_elements,send_positions,local_bad);call add_size(integer_elements,received_ids,local_bad)
    call add_size(integer_elements,received_positions,local_bad);call add_size(integer_elements,directory_owner,local_bad)
    call add_size(integer_elements,directory_count,local_bad);call add_size(integer_elements,sorted_ids,local_bad)
    call add_size(integer_elements,sorted_positions,local_bad)
    call add_size(integer_elements,query_counts,local_bad);call add_size(integer_elements,query_displacements,local_bad)
    call add_size(integer_elements,query_received_counts,local_bad);call add_size(integer_elements,query_received_displacements,local_bad)
    call add_size(integer_elements,query_ids,local_bad);call add_size(integer_elements,query_edges,local_bad)
    call add_size(integer_elements,received_queries,local_bad);call add_size(integer_elements,query_answers,local_bad)
    call add_size(integer_elements,received_answers,local_bad);call add_size(integer_elements,edge_owner,local_bad)
    call add_size(integer_elements,request_ids,local_bad);call add_size(integer_elements,request_edges,local_bad)
    call add_size(integer_elements,received_requests,local_bad);call add_size(integer_elements,plan%send_counts,local_bad)
    call add_size(integer_elements,remote_ids,local_bad);call add_size(integer_elements,remote_edges,local_bad)
    call add_size(integer_elements,edge_unique,local_bad);call add_size(integer_elements,unique_slots,local_bad)
    call add_size(integer_elements,plan%send_displacements,local_bad);call add_size(integer_elements,plan%receive_counts,local_bad)
    call add_size(integer_elements,plan%receive_displacements,local_bad);call add_size(integer_elements,plan%send_positions,local_bad)
    call add_size(integer_elements,plan%value_slots,local_bad)
    complex_elements=int(size(plan%send_values),int64)+int(size(plan%receive_values),int64)
    if(integer_elements>huge(plan%workspace_peak_bytes)/4_int64.or.&
      complex_elements>huge(plan%workspace_peak_bytes)/16_int64)local_bad=1
    if(local_bad==0)then
      if(4_int64*integer_elements>huge(plan%workspace_peak_bytes)-16_int64*complex_elements)local_bad=1
    endif
    call consensus(local_bad,global_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();call clear_rt_dg_sparse_exchange(plan);message='sparse exchange receipt overflow';return;endif
    plan%workspace_peak_bytes=4_int64*integer_elements+16_int64*complex_elements
    plan%comm=comm;plan%nproc=nproc;plan%local_count=size(owned_row_ids)
    plan%catalog_fingerprint=catalog_fingerprint;plan%valid=.true.;ok=.true.;call cleanup();return
900 call cleanup();call clear_rt_dg_sparse_exchange(plan);message='sparse exchange MPI setup failed';return
#else
    ok=.false.;message='sparse exchange requires MPI'
#endif
  contains
#ifdef USE_MPI
    subroutine consensus(local_value,global_value,status)
      integer,intent(in)::local_value;integer,intent(out)::global_value,status
      call MPI_Allreduce(local_value,global_value,1,MPI_INTEGER,MPI_MAX,comm,status)
    end subroutine consensus
    subroutine offsets(input_counts,output_displacements,total,bad)
      integer,intent(in)::input_counts(:);integer,intent(out)::output_displacements(:),total
      integer,intent(inout)::bad;integer::q
      total=0
      do q=1,size(input_counts)
        output_displacements(q)=total
        if(input_counts(q)<0.or.input_counts(q)>huge(0)-total)then;bad=1;return;endif
        total=total+input_counts(q)
      enddo
    end subroutine offsets
    integer function find_local(id)
      integer,intent(in)::id;integer::left,right,middle
      find_local=0;left=1;right=size(sorted_ids)
      do while(left<=right)
        middle=left+(right-left)/2
        if(sorted_ids(middle)==id)then;find_local=sorted_positions(middle);return
        else if(sorted_ids(middle)<id)then;left=middle+1
        else;right=middle-1;endif
      enddo
    end function find_local
    recursive subroutine sort_pairs(keys,positions,left,right)
      integer,intent(inout)::keys(:),positions(:);integer,intent(in)::left,right
      integer::i1,j1,pivot,temp
      if(left>=right)return;i1=left;j1=right;pivot=keys(left+(right-left)/2)
      do
        do while(keys(i1)<pivot);i1=i1+1;enddo
        do while(keys(j1)>pivot);j1=j1-1;enddo
        if(i1<=j1)then
          temp=keys(i1);keys(i1)=keys(j1);keys(j1)=temp
          temp=positions(i1);positions(i1)=positions(j1);positions(j1)=temp;i1=i1+1;j1=j1-1
        endif
        if(i1>j1)exit
      enddo
      if(left<j1)call sort_pairs(keys,positions,left,j1)
      if(i1<right)call sort_pairs(keys,positions,i1,right)
    end subroutine sort_pairs
    subroutine add_size(total,array,bad)
      integer(int64),intent(inout)::total;integer,intent(in)::array(:);integer,intent(inout)::bad
      integer(int64)::amount
      if(bad/=0)return;amount=int(size(array),int64)
      if(amount>huge(total)-total)then;bad=1;else;total=total+amount;endif
    end subroutine add_size
    subroutine cleanup()
      if(allocated(counts))deallocate(counts);if(allocated(displacements))deallocate(displacements)
      if(allocated(received_counts))deallocate(received_counts);if(allocated(received_displacements))deallocate(received_displacements)
      if(allocated(cursor))deallocate(cursor);if(allocated(send_ids))deallocate(send_ids)
      if(allocated(send_positions))deallocate(send_positions);if(allocated(received_ids))deallocate(received_ids)
      if(allocated(received_positions))deallocate(received_positions);if(allocated(directory_owner))deallocate(directory_owner)
      if(allocated(directory_count))deallocate(directory_count)
      if(allocated(sorted_ids))deallocate(sorted_ids);if(allocated(sorted_positions))deallocate(sorted_positions)
      if(allocated(query_counts))deallocate(query_counts);if(allocated(query_displacements))deallocate(query_displacements)
      if(allocated(query_received_counts))deallocate(query_received_counts)
      if(allocated(query_received_displacements))deallocate(query_received_displacements)
      if(allocated(query_ids))deallocate(query_ids);if(allocated(query_edges))deallocate(query_edges)
      if(allocated(received_queries))deallocate(received_queries);if(allocated(query_answers))deallocate(query_answers)
      if(allocated(received_answers))deallocate(received_answers)
      if(allocated(edge_owner))deallocate(edge_owner)
      if(allocated(remote_ids))deallocate(remote_ids);if(allocated(remote_edges))deallocate(remote_edges)
      if(allocated(edge_unique))deallocate(edge_unique);if(allocated(unique_slots))deallocate(unique_slots)
      if(allocated(request_ids))deallocate(request_ids);if(allocated(request_edges))deallocate(request_edges)
      if(allocated(received_requests))deallocate(received_requests)
    end subroutine cleanup
#endif
  end subroutine build_rt_dg_sparse_exchange

  subroutine exchange_rt_dg_sparse_values(comm,plan,local_values,needed_values,ierr)
    integer,intent(in)::comm
    type(s_rt_dg_sparse_exchange),intent(inout)::plan
    complex(real64),intent(in)::local_values(:)
    complex(real64),intent(out)::needed_values(:)
    integer,intent(out)::ierr
#ifdef USE_MPI
    integer::i,local_bad,global_bad,comparison,actual_nproc
    local_bad=0
    if(.not.plan%valid.or.plan%nproc<1.or.plan%local_count/=size(local_values))local_bad=1
    if(.not.allocated(plan%value_slots).or..not.allocated(plan%send_positions).or.&
      .not.allocated(plan%send_values).or..not.allocated(plan%receive_values))local_bad=1
    if(local_bad==0)then
      call MPI_Comm_compare(comm,plan%comm,comparison,ierr)
      if(ierr/=MPI_SUCCESS.or.(comparison/=MPI_IDENT.and.comparison/=MPI_CONGRUENT))local_bad=1
      call MPI_Comm_size(comm,actual_nproc,ierr)
      if(ierr/=MPI_SUCCESS.or.actual_nproc/=plan%nproc)local_bad=1
      if(.not.allocated(plan%send_counts).or..not.allocated(plan%send_displacements).or.&
        .not.allocated(plan%receive_counts).or..not.allocated(plan%receive_displacements))local_bad=1
    endif
    if(local_bad==0)then
      if(size(plan%send_counts)/=plan%nproc.or.size(plan%send_displacements)/=plan%nproc.or.&
        size(plan%receive_counts)/=plan%nproc.or.size(plan%receive_displacements)/=plan%nproc)local_bad=1
      if(size(needed_values)/=size(plan%value_slots).or.size(plan%send_positions)/=size(plan%send_values).or.&
        size(plan%receive_values)/=sum(plan%receive_counts))local_bad=1
      if(.not.valid_layout(plan%send_counts,plan%send_displacements,size(plan%send_values)).or.&
        .not.valid_layout(plan%receive_counts,plan%receive_displacements,size(plan%receive_values)))local_bad=1
      if(any(plan%send_positions<1).or.any(plan%send_positions>size(local_values)))local_bad=1
      if(any(plan%value_slots==0).or.any(plan%value_slots>size(local_values)).or.&
        any(plan%value_slots < -size(plan%receive_values)))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;ierr=MPI_ERR_OTHER;return;endif
    do i=1,size(plan%send_positions);plan%send_values(i)=local_values(plan%send_positions(i));enddo
    call MPI_Alltoallv(plan%send_values,plan%send_counts,plan%send_displacements,MPI_DOUBLE_COMPLEX,&
      plan%receive_values,plan%receive_counts,plan%receive_displacements,MPI_DOUBLE_COMPLEX,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    do i=1,size(needed_values)
      if(plan%value_slots(i)>0)then;needed_values(i)=local_values(plan%value_slots(i));else;needed_values(i)=plan%receive_values(-plan%value_slots(i));endif
    enddo
#else
    ierr=1
#endif
  contains
#ifdef USE_MPI
    logical function valid_layout(counts,displacements,total)
      integer,intent(in)::counts(:),displacements(:),total;integer::q,expected
      valid_layout=.false.;expected=0
      do q=1,size(counts)
        if(counts(q)<0.or.displacements(q)/=expected.or.counts(q)>huge(0)-expected)return
        expected=expected+counts(q)
      enddo
      valid_layout=expected==total
    end function valid_layout
#endif
  end subroutine exchange_rt_dg_sparse_values

  subroutine apply_rt_dg_sparse_rows_tiled(comm,plan,row_offsets,matrix_values,column_slots,&
      local_values,result,tile_width,workspace_peak_bytes,payload_collective_count,ok,message)
    integer,intent(in)::comm,row_offsets(:),column_slots(:),tile_width
    type(s_rt_dg_sparse_exchange),intent(in)::plan
    complex(real64),intent(in)::matrix_values(:),local_values(:,:)
    complex(real64),intent(out)::result(:,:)
    integer(int64),intent(out)::workspace_peak_bytes
    integer,intent(out)::payload_collective_count
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::first,width,i,edge,j,ierr,local_bad,global_bad
    integer(int64)::tile_bytes
    complex(real64),allocatable::received(:)
    ok=.false.;message='';workspace_peak_bytes=0_int64;payload_collective_count=0;local_bad=0
    call validate_tiled_plan_collective(comm,plan,size(local_values,1),size(local_values,2),tile_width,&
      ok,message)
    if(.not.ok)return
    if(tile_width<1.or.size(row_offsets)/=size(result,1)+1.or.size(matrix_values)/=size(column_slots).or.&
      size(local_values,1)/=plan%local_count.or.size(result,2)/=size(local_values,2))local_bad=1
    if(local_bad==0)then
      if(row_offsets(1)/=1.or.row_offsets(size(row_offsets))/=size(matrix_values)+1.or.&
        any(row_offsets(2:)<row_offsets(:size(row_offsets)-1)).or.any(column_slots<1).or.&
        any(column_slots>size(plan%value_slots)))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='tiled sparse action validation reduction failed';return;endif
    if(global_bad/=0)then;message='invalid tiled sparse multi-RHS action contract';return;endif
    result=(0d0,0d0)
    do first=1,size(local_values,2),tile_width
      width=min(tile_width,size(local_values,2)-first+1)
      call exchange_unique_rows_tile(comm,plan,local_values,first,width,received,tile_bytes,ok,message)
      if(.not.ok)return
      workspace_peak_bytes=max(workspace_peak_bytes,tile_bytes)
      payload_collective_count=payload_collective_count+1
      do i=1,size(result,1)
        do edge=row_offsets(i),row_offsets(i+1)-1
          if(plan%value_slots(column_slots(edge))>0)then
            do j=1,width
              result(i,first+j-1)=result(i,first+j-1)+matrix_values(edge)*&
                local_values(plan%value_slots(column_slots(edge)),first+j-1)
            enddo
          else
            do j=1,width
              result(i,first+j-1)=result(i,first+j-1)+matrix_values(edge)*&
                received((-plan%value_slots(column_slots(edge))-1)*width+j)
            enddo
          endif
        enddo
      enddo
      deallocate(received)
    enddo
    ok=.true.;message=''
#else
    ok=.false.;message='tiled sparse multi-RHS action requires MPI'
    workspace_peak_bytes=0_int64;payload_collective_count=0
#endif
  end subroutine apply_rt_dg_sparse_rows_tiled

  subroutine accumulate_rt_dg_sparse_density_tiled(comm,plan,point_offsets,support_slots,support_values,&
      local_values,occupations,density,tile_width,workspace_peak_bytes,payload_collective_count,ok,message)
    integer,intent(in)::comm,point_offsets(:),support_slots(:),tile_width
    type(s_rt_dg_sparse_exchange),intent(in)::plan
    complex(real64),intent(in)::support_values(:),local_values(:,:)
    real(real64),intent(in)::occupations(:)
    real(real64),intent(out)::density(:)
    integer(int64),intent(out)::workspace_peak_bytes
    integer,intent(out)::payload_collective_count
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::first,width,p,edge,j,ierr,local_bad,global_bad
    integer(int64)::tile_bytes
    complex(real64),allocatable::received(:),orbitals(:)
    ok=.false.;message='';workspace_peak_bytes=0_int64;payload_collective_count=0;local_bad=0
    call validate_tiled_plan_collective(comm,plan,size(local_values,1),size(local_values,2),tile_width,&
      ok,message)
    if(.not.ok)return
    if(tile_width<1.or.size(point_offsets)/=size(density)+1.or.size(support_slots)/=size(support_values).or.&
      size(occupations)/=size(local_values,2).or.size(local_values,1)/=plan%local_count)local_bad=1
    if(local_bad==0)then
      if(point_offsets(1)/=1.or.point_offsets(size(point_offsets))/=size(support_values)+1.or.&
        any(point_offsets(2:)<point_offsets(:size(point_offsets)-1)).or.any(support_slots<1).or.&
        any(support_slots>size(plan%value_slots)))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='tiled point-CSR validation reduction failed';return;endif
    if(global_bad/=0)then;message='invalid tiled point-CSR density contract';return;endif
    density=0d0
    do first=1,size(local_values,2),tile_width
      width=min(tile_width,size(local_values,2)-first+1)
      call exchange_unique_rows_tile(comm,plan,local_values,first,width,received,tile_bytes,ok,message)
      if(.not.ok)return
      allocate(orbitals(width),stat=local_bad)
      local_bad=merge(0,1,local_bad==0)
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='tiled point-CSR allocation reduction failed';return;endif
      if(global_bad/=0)then;message='cannot allocate tiled point-CSR orbital workspace';return;endif
      do p=1,size(density)
        orbitals=(0d0,0d0)
        do edge=point_offsets(p),point_offsets(p+1)-1
          if(plan%value_slots(support_slots(edge))>0)then
            orbitals=orbitals+support_values(edge)*&
              local_values(plan%value_slots(support_slots(edge)),first:first+width-1)
          else
            do j=1,width
              orbitals(j)=orbitals(j)+support_values(edge)*&
                received((-plan%value_slots(support_slots(edge))-1)*width+j)
            enddo
          endif
        enddo
        do j=1,width;density(p)=density(p)+occupations(first+j-1)*abs(orbitals(j))**2;enddo
      enddo
      workspace_peak_bytes=max(workspace_peak_bytes,tile_bytes+16_int64*int(size(orbitals),int64))
      payload_collective_count=payload_collective_count+1
      deallocate(received,orbitals)
    enddo
    ok=.true.;message=''
#else
    ok=.false.;message='tiled point-CSR density requires MPI'
    workspace_peak_bytes=0_int64;payload_collective_count=0
#endif
  end subroutine accumulate_rt_dg_sparse_density_tiled

#ifdef USE_MPI
  subroutine validate_tiled_plan_collective(comm,plan,nlocal,nrhs,tile_width,ok,message)
    integer,intent(in)::comm,nlocal,nrhs,tile_width
    type(s_rt_dg_sparse_exchange),intent(in)::plan
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::nproc,ierr,local_bad,global_bad,minimum_value,maximum_value,comm_relation
    integer(int64)::minimum_fingerprint,maximum_fingerprint
    ok=.false.;message='';local_bad=0
    call MPI_Comm_size(comm,nproc,ierr)
    if(ierr/=MPI_SUCCESS)then;message='tiled sparse plan communicator query failed';return;endif
    comm_relation=MPI_UNEQUAL
    if(plan%valid)call MPI_Comm_compare(plan%comm,comm,comm_relation,ierr)
    if(ierr/=MPI_SUCCESS)then;message='tiled sparse plan communicator comparison failed';return;endif
    if(.not.plan%valid.or.(comm_relation/=MPI_IDENT.and.comm_relation/=MPI_CONGRUENT).or.&
      plan%nproc/=nproc.or.plan%local_count/=nlocal.or.&
      nlocal<0.or.nrhs<1.or.tile_width<1.or.plan%catalog_fingerprint==0_int64)local_bad=1
    if(.not.allocated(plan%send_counts).or..not.allocated(plan%send_displacements).or.&
      .not.allocated(plan%receive_counts).or..not.allocated(plan%receive_displacements).or.&
      .not.allocated(plan%send_positions).or..not.allocated(plan%value_slots).or.&
      .not.allocated(plan%send_values).or..not.allocated(plan%receive_values))then
      local_bad=1
    else
      if(size(plan%send_counts)/=nproc.or.size(plan%send_displacements)/=nproc.or.&
        size(plan%receive_counts)/=nproc.or.size(plan%receive_displacements)/=nproc)then
        local_bad=1
      else
        if(.not.valid_counts_layout(plan%send_counts,plan%send_displacements,size(plan%send_positions)).or.&
          .not.valid_counts_layout(plan%receive_counts,plan%receive_displacements,size(plan%receive_values)))local_bad=1
      endif
      if(size(plan%send_values)/=size(plan%send_positions))local_bad=1
      if(any(plan%send_positions<1).or.any(plan%send_positions>nlocal))local_bad=1
      if(any(plan%value_slots==0).or.any(plan%value_slots>nlocal).or.&
        any(plan%value_slots < -size(plan%receive_values)))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='tiled sparse plan validation reduction failed';return;endif
    if(global_bad/=0)then;message='invalid tiled sparse exchange plan';return;endif
    call MPI_Allreduce(tile_width,minimum_value,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='tiled sparse tile-width agreement failed';return;endif
    call MPI_Allreduce(tile_width,maximum_value,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='tiled sparse tile-width agreement failed';return;endif
    if(minimum_value/=maximum_value)then;message='rank-disagreeing tiled sparse tile width';return;endif
    call MPI_Allreduce(nrhs,minimum_value,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='tiled sparse RHS agreement failed';return;endif
    call MPI_Allreduce(nrhs,maximum_value,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='tiled sparse RHS agreement failed';return;endif
    if(minimum_value/=maximum_value)then;message='rank-disagreeing tiled sparse RHS count';return;endif
    call MPI_Allreduce(plan%catalog_fingerprint,minimum_fingerprint,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='tiled sparse provenance agreement failed';return;endif
    call MPI_Allreduce(plan%catalog_fingerprint,maximum_fingerprint,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='tiled sparse provenance agreement failed';return;endif
    if(minimum_fingerprint/=maximum_fingerprint)then
      message='rank-disagreeing tiled sparse exchange provenance';return
    endif
    ok=.true.
  contains
    logical function valid_counts_layout(counts,displacements,total)
      integer,intent(in)::counts(:),displacements(:),total
      integer::q,expected
      valid_counts_layout=.false.;expected=0
      do q=1,size(counts)
        if(counts(q)<0.or.displacements(q)/=expected.or.counts(q)>huge(0)-expected)return
        expected=expected+counts(q)
      enddo
      valid_counts_layout=expected==total
    end function valid_counts_layout
  end subroutine validate_tiled_plan_collective

  subroutine exchange_unique_rows_tile(comm,plan,local_values,first,width,received,workspace_bytes,ok,message)
    integer,intent(in)::comm,first,width
    type(s_rt_dg_sparse_exchange),intent(in)::plan
    complex(real64),intent(in)::local_values(:,:)
    complex(real64),allocatable,intent(out)::received(:)
    integer(int64),intent(out)::workspace_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::i,j,ierr,local_bad,global_bad,allocation_status,total_send,total_receive
    integer,allocatable::send_counts(:),send_displacements(:),receive_counts(:),receive_displacements(:)
    complex(real64),allocatable::send_values(:)
    ok=.false.;message='';workspace_bytes=0_int64;local_bad=0
    if(.not.plan%valid.or.width<1.or.first<1.or.first>size(local_values,2)-width+1) local_bad=1
    call checked_default_product(size(plan%send_positions),width,total_send,local_bad)
    call checked_default_product(size(plan%receive_values),width,total_receive,local_bad)
    if(local_bad==0)then
      if(int(total_send,int64)>huge(0_int64)/16_int64.or.int(total_receive,int64)>huge(0_int64)/16_int64)local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='tiled sparse extent reduction failed';return;endif
    if(global_bad/=0)then;message='tiled sparse MPI extent overflow';return;endif
    allocate(send_counts(plan%nproc),send_displacements(plan%nproc),receive_counts(plan%nproc),&
      receive_displacements(plan%nproc),send_values(total_send),received(total_receive),&
      stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='unique-row tiled allocation reduction failed';return;endif
    if(global_bad/=0)then;message='cannot allocate unique-row tiled halo';return;endif
    do i=1,plan%nproc
      call checked_default_product(plan%send_counts(i),width,send_counts(i),local_bad)
      call checked_default_product(plan%send_displacements(i),width,send_displacements(i),local_bad)
      call checked_default_product(plan%receive_counts(i),width,receive_counts(i),local_bad)
      call checked_default_product(plan%receive_displacements(i),width,receive_displacements(i),local_bad)
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='tiled sparse count reduction failed';return;endif
    if(global_bad/=0)then;message='tiled sparse count/displacement overflow';return;endif
    do i=1,size(plan%send_positions);do j=1,width
      send_values((i-1)*width+j)=local_values(plan%send_positions(i),first+j-1)
    enddo;enddo
    call MPI_Alltoallv(send_values,send_counts,send_displacements,MPI_DOUBLE_COMPLEX,received,&
      receive_counts,receive_displacements,MPI_DOUBLE_COMPLEX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='unique-row tiled halo payload exchange failed';return;endif
    workspace_bytes=16_int64*(int(total_send,int64)+int(total_receive,int64))+&
      16_int64*int(plan%nproc,int64)
    ok=.true.;message=''
  end subroutine exchange_unique_rows_tile

  subroutine checked_default_product(left,right,value,bad)
    integer,intent(in)::left,right
    integer,intent(out)::value
    integer,intent(inout)::bad
    integer(int64)::product64
    value=0
    if(bad/=0)return
    if(left<0.or.right<0)then;bad=1;return;endif
    if(left/=0.and.int(right,int64)>int(huge(0),int64)/int(left,int64))then;bad=1;return;endif
    product64=int(left,int64)*int(right,int64)
    if(product64>int(huge(0),int64))then;bad=1;return;endif
    value=int(product64)
  end subroutine checked_default_product
#endif

  subroutine checked_rt_dg_sparse_extent_product(left,right,value,ok)
    integer(int64),intent(in)::left,right
    integer(int64),intent(out)::value
    logical,intent(out)::ok
    value=0_int64;ok=left>=0_int64.and.right>=0_int64
    if(.not.ok)return
    if(left/=0_int64.and.right>int(huge(0),int64)/left)then;ok=.false.;return;endif
    value=left*right
  end subroutine checked_rt_dg_sparse_extent_product

  subroutine clear_rt_dg_sparse_exchange(plan)
    type(s_rt_dg_sparse_exchange),intent(inout)::plan
    if(allocated(plan%send_counts))deallocate(plan%send_counts)
    if(allocated(plan%send_displacements))deallocate(plan%send_displacements)
    if(allocated(plan%receive_counts))deallocate(plan%receive_counts)
    if(allocated(plan%receive_displacements))deallocate(plan%receive_displacements)
    if(allocated(plan%send_positions))deallocate(plan%send_positions)
    if(allocated(plan%value_slots))deallocate(plan%value_slots)
    if(allocated(plan%send_values))deallocate(plan%send_values)
    if(allocated(plan%receive_values))deallocate(plan%receive_values)
    plan%valid=.false.
    plan%workspace_peak_bytes=0_int64;plan%catalog_fingerprint=0_int64;plan%local_count=0
  end subroutine clear_rt_dg_sparse_exchange
end module rt_dg_hybrid_sparse_exchange
