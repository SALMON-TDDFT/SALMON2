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
    integer,allocatable::requests(:)
  end type s_rt_dg_sparse_exchange
  public::build_rt_dg_sparse_exchange,exchange_rt_dg_sparse_values,exchange_rt_dg_sparse_matrix,&
    clear_rt_dg_sparse_exchange
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
      plan%send_positions(total_send),plan%send_values(total_send),plan%receive_values(total_receive),&
      plan%requests(2*nproc),stat=allocation_status)
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
    call add_size(integer_elements,plan%requests,local_bad)
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
    integer::i,rank,local_bad,global_bad,nrequest,comparison,actual_nproc,post_status
    local_bad=0
    if(.not.plan%valid.or.plan%nproc<1.or.plan%local_count/=size(local_values))local_bad=1
    if(.not.allocated(plan%value_slots).or..not.allocated(plan%send_positions).or.&
      .not.allocated(plan%send_values).or..not.allocated(plan%receive_values).or..not.allocated(plan%requests))local_bad=1
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
        size(plan%requests)<2*plan%nproc) local_bad=1
      if(.not.valid_layout(plan%send_counts,plan%send_displacements,size(plan%send_values)).or.&
        .not.valid_layout(plan%receive_counts,plan%receive_displacements,size(plan%receive_values)))local_bad=1
      if(any(plan%send_positions<1).or.any(plan%send_positions>size(local_values)))local_bad=1
      if(any(plan%value_slots==0).or.any(plan%value_slots>size(local_values)).or.&
        any(plan%value_slots < -size(plan%receive_values)))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;ierr=MPI_ERR_OTHER;return;endif
    do i=1,size(plan%send_positions);plan%send_values(i)=local_values(plan%send_positions(i));enddo
    nrequest=0
    do rank=0,plan%nproc-1
      if(plan%receive_counts(rank+1)>0)then
        nrequest=nrequest+1
        call MPI_Irecv(plan%receive_values(plan%receive_displacements(rank+1)+1),plan%receive_counts(rank+1),&
          MPI_DOUBLE_COMPLEX,rank,3817,comm,plan%requests(nrequest),post_status)
        if(post_status/=MPI_SUCCESS)then;call cancel_requests(nrequest-1);ierr=post_status;return;endif
      endif
    enddo
    do rank=0,plan%nproc-1
      if(plan%send_counts(rank+1)>0)then
        nrequest=nrequest+1
        call MPI_Isend(plan%send_values(plan%send_displacements(rank+1)+1),plan%send_counts(rank+1),&
          MPI_DOUBLE_COMPLEX,rank,3817,comm,plan%requests(nrequest),post_status)
        if(post_status/=MPI_SUCCESS)then;call cancel_requests(nrequest-1);ierr=post_status;return;endif
      endif
    enddo
    if(nrequest>0)call MPI_Waitall(nrequest,plan%requests,MPI_STATUSES_IGNORE,ierr)
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
    subroutine cancel_requests(count)
      integer,intent(in)::count;integer::q,cancel_status,wait_status
      do q=1,count;call MPI_Cancel(plan%requests(q),cancel_status);enddo
      if(count>0)call MPI_Waitall(count,plan%requests,MPI_STATUSES_IGNORE,wait_status)
    end subroutine cancel_requests
#endif
  end subroutine exchange_rt_dg_sparse_values

  subroutine exchange_rt_dg_sparse_matrix(comm,plan,local_values,needed_values,workspace_peak_bytes,&
      payload_collective_count,ok,message)
    integer,intent(in)::comm
    type(s_rt_dg_sparse_exchange),intent(in)::plan
    complex(real64),intent(in)::local_values(:,:)
    complex(real64),intent(out)::needed_values(:,:)
    integer(int64),intent(out)::workspace_peak_bytes
    integer,intent(out)::payload_collective_count
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nrhs,i,j,ierr,local_bad,global_bad,comparison,actual_nproc,allocation_status
    integer,allocatable::send_counts(:),send_displacements(:),receive_counts(:),receive_displacements(:)
    complex(real64),allocatable::send_values(:),receive_values(:)
    ok=.false.;message='';workspace_peak_bytes=0_int64;payload_collective_count=0
    nrhs=size(local_values,2);local_bad=0
    if(.not.plan%valid.or.plan%nproc<1.or.plan%local_count/=size(local_values,1).or.nrhs<1.or.&
      size(needed_values,1)/=size(plan%value_slots).or.size(needed_values,2)/=nrhs)local_bad=1
    if(.not.allocated(plan%value_slots).or..not.allocated(plan%send_positions).or.&
      .not.allocated(plan%send_counts).or..not.allocated(plan%send_displacements).or.&
      .not.allocated(plan%receive_counts).or..not.allocated(plan%receive_displacements))local_bad=1
    if(local_bad==0)then
      call MPI_Comm_compare(comm,plan%comm,comparison,ierr)
      if(ierr/=MPI_SUCCESS.or.(comparison/=MPI_IDENT.and.comparison/=MPI_CONGRUENT))local_bad=1
      call MPI_Comm_size(comm,actual_nproc,ierr)
      if(ierr/=MPI_SUCCESS.or.actual_nproc/=plan%nproc)local_bad=1
      if(size(plan%send_counts)/=plan%nproc.or.size(plan%send_displacements)/=plan%nproc.or.&
        size(plan%receive_counts)/=plan%nproc.or.size(plan%receive_displacements)/=plan%nproc)local_bad=1
    endif
    if(local_bad==0)then
      if(any(plan%send_counts>huge(0)/nrhs).or.any(plan%send_displacements>huge(0)/nrhs).or.&
        any(plan%receive_counts>huge(0)/nrhs).or.any(plan%receive_displacements>huge(0)/nrhs))local_bad=1
      if(any(plan%send_positions<1).or.any(plan%send_positions>size(local_values,1)))local_bad=1
      if(any(plan%value_slots==0).or.any(plan%value_slots>size(local_values,1)).or.&
        any(plan%value_slots < -sum(plan%receive_counts)))local_bad=1
      if(int(size(plan%send_positions),int64)>huge(0_int64)/int(nrhs,int64).or.&
        int(sum(plan%receive_counts),int64)>huge(0_int64)/int(nrhs,int64))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid sparse coefficient-matrix exchange contract';return;endif
    allocate(send_counts(plan%nproc),send_displacements(plan%nproc),receive_counts(plan%nproc),&
      receive_displacements(plan%nproc),send_values(size(plan%send_positions)*nrhs),&
      receive_values(sum(plan%receive_counts)*nrhs),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='cannot allocate sparse coefficient-matrix exchange';return;endif
    send_counts=plan%send_counts*nrhs;send_displacements=plan%send_displacements*nrhs
    receive_counts=plan%receive_counts*nrhs;receive_displacements=plan%receive_displacements*nrhs
    do i=1,size(plan%send_positions);do j=1,nrhs
      send_values((i-1)*nrhs+j)=local_values(plan%send_positions(i),j)
    enddo;enddo
    call MPI_Alltoallv(send_values,send_counts,send_displacements,MPI_DOUBLE_COMPLEX,&
      receive_values,receive_counts,receive_displacements,MPI_DOUBLE_COMPLEX,comm,ierr)
    payload_collective_count=1
    if(ierr/=MPI_SUCCESS)then;message='sparse coefficient-matrix payload exchange failed';return;endif
    do i=1,size(plan%value_slots);do j=1,nrhs
      if(plan%value_slots(i)>0)then
        needed_values(i,j)=local_values(plan%value_slots(i),j)
      else
        needed_values(i,j)=receive_values((-plan%value_slots(i)-1)*nrhs+j)
      endif
    enddo;enddo
    workspace_peak_bytes=16_int64*int(nrhs,int64)*&
      int(size(plan%send_positions)+sum(plan%receive_counts),int64)
    ok=.true.;message=''
#else
    ok=.false.;message='sparse coefficient-matrix exchange requires MPI'
    workspace_peak_bytes=0_int64;payload_collective_count=0
#endif
  end subroutine exchange_rt_dg_sparse_matrix

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
    if(allocated(plan%requests))deallocate(plan%requests)
    plan%valid=.false.
    plan%workspace_peak_bytes=0_int64;plan%catalog_fingerprint=0_int64;plan%local_count=0
  end subroutine clear_rt_dg_sparse_exchange
end module rt_dg_hybrid_sparse_exchange
