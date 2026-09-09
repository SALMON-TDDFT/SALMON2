#include "config.h"
module rt_dg_hybrid_structural_graph
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::build_rt_dg_hybrid_structural_graph,checked_rt_dg_hybrid_structural_capacity,&
    collective_rt_dg_hybrid_structural_capacity_status
#ifdef USE_MPI
  type::key_set
    integer(int64),allocatable::slot(:)
    integer(int64)::count=0_int64
    integer(int64)::peak_capacity=0_int64
    logical::failed=.false.
  end type key_set
#endif
contains
  subroutine build_rt_dg_hybrid_structural_graph(comm,global_count,row_ids,basis_values,metric_rows,kinetic_rows,&
      nonlocal_rows,local_rows,sipg_rows,hamiltonian_rows,position_rows,metric_offsets,metric_columns,&
      operator_offsets,operator_columns,ok,message,local_unique_candidates,peak_workspace_keys)
    integer,intent(in)::comm,global_count
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::basis_values(:,:),metric_rows(:,:),kinetic_rows(:,:),nonlocal_rows(:,:),&
      local_rows(:,:),sipg_rows(:,:),hamiltonian_rows(:,:),position_rows(:,:,:)
    integer,allocatable,intent(out)::metric_offsets(:),metric_columns(:),operator_offsets(:),operator_columns(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer(int64),intent(out),optional::local_unique_candidates,peak_workspace_keys
#ifdef USE_MPI
    integer::i,j,p,row,nowned,npoint,nactive,ierr,local_bad,global_bad,rank,nproc,metric_requested,operator_requested
    integer,allocatable::active(:),owners(:),owner_marks(:)
    integer(int64)::workspace_peak
    type(key_set)::point_support,metric_support,operator_support,closure_support
    logical::local_capacity_ok,operator_capacity_ok
    ok=.false.;message='';nowned=size(row_ids);npoint=size(basis_values,2);local_bad=0;workspace_peak=0_int64
    if(present(local_unique_candidates))local_unique_candidates=0_int64
    if(present(peak_workspace_keys))peak_workspace_keys=0_int64
    if(global_count<1.or.size(basis_values,1)/=global_count.or.&
      any(shape(metric_rows)/=[nowned,global_count]).or.any(shape(kinetic_rows)/=[nowned,global_count]).or.&
      any(shape(nonlocal_rows)/=[nowned,global_count]).or.any(shape(local_rows)/=[nowned,global_count]).or.&
      any(shape(sipg_rows)/=[nowned,global_count]).or.any(shape(hamiltonian_rows)/=[nowned,global_count]).or.&
      any(shape(position_rows)/=[3,nowned,global_count]))local_bad=1
    if(any(row_ids<1_int64).or.any(row_ids>int(global_count,int64)))local_bad=1
    if(.not.finite_matrix(basis_values).or..not.finite_matrix(metric_rows).or..not.finite_matrix(kinetic_rows).or.&
      .not.finite_matrix(nonlocal_rows).or..not.finite_matrix(local_rows).or..not.finite_matrix(sipg_rows).or.&
      .not.finite_matrix(hamiltonian_rows).or..not.finite_rank3(position_rows))local_bad=1
    if(int(global_count,int64)>huge(0_int64)/int(global_count,int64))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid Hybrid structural graph inputs';return;endif
    call checked_rt_dg_hybrid_structural_capacity(2,nowned,metric_requested,local_capacity_ok)
    call checked_rt_dg_hybrid_structural_capacity(4,nowned,operator_requested,operator_capacity_ok)
    local_capacity_ok=local_capacity_ok.and.operator_capacity_ok
    call collective_rt_dg_hybrid_structural_capacity_status(comm,local_capacity_ok,ok,message)
    if(.not.ok)return
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    allocate(owners(global_count),owner_marks(global_count));owners=0;owner_marks=0
    do i=1,nowned
      owners(int(row_ids(i)))=rank+1;owner_marks(int(row_ids(i)))=owner_marks(int(row_ids(i)))+1
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,owners,global_count,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(MPI_IN_PLACE,owner_marks,global_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(owner_marks/=1))then;message='Hybrid structural rows do not have unique owners';return;endif

    allocate(active(global_count));call initialize_set(metric_support,metric_requested)
    call initialize_set(operator_support,operator_requested)
    ! Deduplicate the support induced by this rank's local grid once, then
    ! route each key directly to its row owner in one packed collective.
    ! Communication count is independent of MPI size and payload remains
    ! proportional to the locally generated sparse support.
    call reset_set(point_support,64)
    do p=1,npoint
      nactive=0
      do i=1,global_count
        if(basis_values(i,p)/=(0d0,0d0))then;nactive=nactive+1;active(nactive)=i;endif
      enddo
      do i=1,nactive;do j=1,nactive
        call insert_key(point_support,directed_key(active(i),active(j),global_count))
      enddo;enddo
    enddo
    call collective_rt_dg_hybrid_structural_capacity_status(comm,.not.point_support%failed,ok,message)
    if(.not.ok)return
    if(present(local_unique_candidates))local_unique_candidates=point_support%count
    workspace_peak=max(workspace_peak,point_support%peak_capacity)
    call route_set_to_row_owners(comm,global_count,owners,point_support,operator_support,.false.,ierr,workspace_peak)
    if(ierr/=MPI_SUCCESS)then;message='Hybrid point-support owner exchange failed';return;endif
    do i=1,nowned
      row=int(row_ids(i))
      do j=1,global_count
        if(metric_rows(i,j)/=(0d0,0d0))call insert_key(metric_support,directed_key(row,j,global_count))
        if(metric_rows(i,j)/=(0d0,0d0).or.kinetic_rows(i,j)/=(0d0,0d0).or.&
          nonlocal_rows(i,j)/=(0d0,0d0).or.local_rows(i,j)/=(0d0,0d0).or.sipg_rows(i,j)/=(0d0,0d0).or.&
          hamiltonian_rows(i,j)/=(0d0,0d0).or.any(position_rows(:,i,j)/=(0d0,0d0)))&
          call insert_key(operator_support,directed_key(row,j,global_count))
      enddo
    enddo
    local_capacity_ok=.not.metric_support%failed.and..not.operator_support%failed
    call collective_rt_dg_hybrid_structural_capacity_status(comm,local_capacity_ok,ok,message)
    if(.not.ok)return
    call clone_set(metric_support,closure_support)
    call route_set_to_row_owners(comm,global_count,owners,closure_support,metric_support,.true.,ierr,workspace_peak)
    if(ierr==MPI_SUCCESS)call merge_set(metric_support,operator_support)
    local_capacity_ok=ierr==MPI_SUCCESS.and..not.operator_support%failed
    call collective_rt_dg_hybrid_structural_capacity_status(comm,local_capacity_ok,ok,message)
    if(.not.ok)return
    if(ierr==MPI_SUCCESS)call clone_set(operator_support,closure_support)
    if(ierr==MPI_SUCCESS)call route_set_to_row_owners(comm,global_count,owners,closure_support,operator_support,&
      .true.,ierr,workspace_peak)
    if(ierr/=MPI_SUCCESS)then;message='Hybrid Hermitian support closure exchange failed';return;endif
    local_bad=merge(1,0,metric_support%count>int(huge(0),int64).or.&
      operator_support%count>int(huge(0),int64))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='Hybrid CSR count exceeds MPI integer extent';return;endif
    call build_owned_csr(global_count,row_ids,metric_support,metric_offsets,metric_columns)
    call build_owned_csr(global_count,row_ids,operator_support,operator_offsets,operator_columns)
    if(present(peak_workspace_keys))peak_workspace_keys=max(workspace_peak,point_support%peak_capacity,&
      metric_support%peak_capacity,operator_support%peak_capacity)
    ok=.true.;message=''
#else
    ok=.false.;message='Hybrid structural graph requires MPI'
#endif
  end subroutine build_rt_dg_hybrid_structural_graph

  pure subroutine checked_rt_dg_hybrid_structural_capacity(multiplier,local_count,requested,ok)
    integer,intent(in)::multiplier,local_count
    integer,intent(out)::requested
    logical,intent(out)::ok
    if(multiplier<=0.or.local_count<0)then;requested=0;ok=.false.;return;endif
    if(local_count>huge(0)/multiplier)then;requested=0;ok=.false.;return;endif
    requested=max(64,multiplier*local_count);ok=.true.
  end subroutine checked_rt_dg_hybrid_structural_capacity

  subroutine collective_rt_dg_hybrid_structural_capacity_status(comm,local_ok,ok,message)
    integer,intent(in)::comm
    logical,intent(in)::local_ok
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::local_bad,global_bad,ierr
    local_bad=merge(0,1,local_ok)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='Hybrid structural capacity status reduction failed';return;endif
    ok=global_bad==0
    if(ok)then;message='';else;message='Hybrid structural hash capacity exceeds integer extent';endif
#else
    ok=.false.;message='Hybrid structural capacity status requires MPI'
#endif
  end subroutine collective_rt_dg_hybrid_structural_capacity_status

#ifdef USE_MPI
  pure integer(int64) function directed_key(row,column,n) result(key)
    integer,intent(in)::row,column,n
    key=int(row-1,int64)*int(n,int64)+int(column,int64)
  end function directed_key
  subroutine initialize_set(set,requested)
    type(key_set),intent(inout)::set
    integer,intent(in)::requested
    integer::capacity
    capacity=64;do while(capacity<requested.and.capacity<=huge(capacity)/2);capacity=capacity*2;enddo
    allocate(set%slot(capacity));set%slot=0_int64;set%count=0_int64;set%peak_capacity=int(capacity,int64)
    set%failed=.false.
  end subroutine initialize_set
  subroutine reset_set(set,requested)
    type(key_set),intent(inout)::set
    integer,intent(in)::requested
    if(allocated(set%slot))deallocate(set%slot)
    call initialize_set(set,requested)
  end subroutine reset_set
  subroutine clone_set(source,target)
    type(key_set),intent(in)::source
    type(key_set),intent(inout)::target
    if(allocated(target%slot))deallocate(target%slot)
    allocate(target%slot,source=source%slot);target%count=source%count;target%peak_capacity=source%peak_capacity
    target%failed=source%failed
  end subroutine clone_set
  subroutine insert_key(set,key)
    type(key_set),intent(inout)::set
    integer(int64),intent(in)::key
    integer::position,capacity,next_capacity
    logical::growth_ok
    if(set%failed)return
    if(key<=0_int64)then;set%failed=.true.;return;endif
    if(.not.allocated(set%slot))call initialize_set(set,64)
    capacity=size(set%slot)
    if(set%count*10_int64>=int(capacity,int64)*7_int64)then
      call checked_rt_dg_hybrid_structural_capacity(2,capacity,next_capacity,growth_ok)
      if(.not.growth_ok)then;set%failed=.true.;return;endif
      call rehash_set(set,next_capacity);capacity=size(set%slot)
    endif
    position=int(modulo(key-1_int64,int(capacity,int64)))+1
    do
      if(set%slot(position)==0_int64)then;set%slot(position)=key;set%count=set%count+1_int64;return
      elseif(set%slot(position)==key)then;return;endif
      position=position+1;if(position>capacity)position=1
    enddo
  end subroutine insert_key
  subroutine rehash_set(set,new_capacity)
    type(key_set),intent(inout)::set
    integer,intent(in)::new_capacity
    integer(int64),allocatable::old(:)
    integer::i
    call move_alloc(set%slot,old);allocate(set%slot(new_capacity));set%slot=0_int64;set%count=0_int64
    set%peak_capacity=max(set%peak_capacity,int(new_capacity,int64))
    do i=1,size(old);if(old(i)>0_int64)call insert_key(set,old(i));enddo
  end subroutine rehash_set
  subroutine extract_sorted(set,keys)
    type(key_set),intent(in)::set
    integer(int64),allocatable,intent(out)::keys(:)
    integer::i,q
    if(set%failed.or.set%count>int(huge(0),int64))then;allocate(keys(0));return;endif
    allocate(keys(int(set%count)));q=0
    do i=1,size(set%slot);if(set%slot(i)>0_int64)then;q=q+1;keys(q)=set%slot(i);endif;enddo
    if(size(keys)>1)call sort_int64(keys,1,size(keys))
  end subroutine extract_sorted
  subroutine route_set_to_row_owners(comm,n,owners,source,target,reverse,ierr,workspace_peak)
    integer,intent(in)::comm,n,owners(:)
    type(key_set),intent(in)::source
    type(key_set),intent(inout)::target
    logical,intent(in)::reverse
    integer,intent(out)::ierr
    integer(int64),intent(inout)::workspace_peak
    integer::nproc,p,q,row,column,destination,total_send,total_recv,local_bad,global_bad,send_bad,recv_bad
    integer(int64)::total_send64,total_recv64
    integer,allocatable::send_counts(:),recv_counts(:),send_displacements(:),recv_displacements(:),cursor(:)
    integer(int64),allocatable::keys(:),send_keys(:),recv_keys(:)
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    allocate(send_counts(nproc),recv_counts(nproc),send_displacements(nproc),recv_displacements(nproc),cursor(nproc))
    local_bad=merge(1,0,source%failed.or.target%failed.or.source%count>int(huge(0),int64))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;ierr=-2;return;endif
    send_counts=0;call extract_sorted(source,keys)
    do q=1,size(keys)
      row=int((keys(q)-1_int64)/int(n,int64))+1;column=int(modulo(keys(q)-1_int64,int(n,int64)))+1
      if(reverse)then;p=row;row=column;column=p;endif
      destination=owners(row);if(destination<1.or.destination>nproc)then;ierr=1;return;endif
      if(send_counts(destination)==huge(0))then;local_bad=1;else;send_counts(destination)=send_counts(destination)+1;endif
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;ierr=-2;return;endif
    call MPI_Alltoall(send_counts,1,MPI_INTEGER,recv_counts,1,MPI_INTEGER,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call make_checked_displacements(send_counts,send_displacements,total_send,total_send64,send_bad)
    call make_checked_displacements(recv_counts,recv_displacements,total_recv,total_recv64,recv_bad)
    local_bad=max(send_bad,recv_bad)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;ierr=-2;return;endif
    allocate(send_keys(total_send),recv_keys(total_recv));cursor=send_displacements
    do q=1,size(keys)
      row=int((keys(q)-1_int64)/int(n,int64))+1;column=int(modulo(keys(q)-1_int64,int(n,int64)))+1
      if(reverse)then;p=row;row=column;column=p;endif
      destination=owners(row);cursor(destination)=cursor(destination)+1;send_keys(cursor(destination))=directed_key(row,column,n)
    enddo
    call MPI_Alltoallv(send_keys,send_counts,send_displacements,MPI_INTEGER8,recv_keys,recv_counts,&
      recv_displacements,MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    do q=1,total_recv;call insert_key(target,recv_keys(q));enddo
    local_bad=merge(0,1,.not.target%failed)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    if(global_bad/=0)then;ierr=-3;return;endif
    workspace_peak=max(workspace_peak,int(total_send,int64),int(total_recv,int64))
  end subroutine route_set_to_row_owners
  subroutine make_checked_displacements(counts,displacements,total,total64,bad)
    integer,intent(in)::counts(:)
    integer,intent(out)::displacements(:),total,bad
    integer(int64),intent(out)::total64
    integer::p
    integer(int64)::running
    bad=0;running=0_int64
    do p=1,size(counts)
      if(running>int(huge(0),int64))then;bad=1;displacements(p)=0;else;displacements(p)=int(running);endif
      running=running+int(counts(p),int64)
    enddo
    total64=running
    if(running>int(huge(0),int64))then;bad=1;total=0;else;total=int(running);endif
  end subroutine make_checked_displacements
  subroutine merge_set(source,target)
    type(key_set),intent(in)::source
    type(key_set),intent(inout)::target
    integer::i
    do i=1,size(source%slot);if(source%slot(i)>0_int64)call insert_key(target,source%slot(i));enddo
  end subroutine merge_set
  subroutine build_owned_csr(n,row_ids,set,offsets,columns)
    integer,intent(in)::n
    integer(int64),intent(in)::row_ids(:)
    type(key_set),intent(in)::set
    integer,allocatable,intent(out)::offsets(:),columns(:)
    integer(int64),allocatable::keys(:)
    integer::i,q,row,key_row,count
    call extract_sorted(set,keys);allocate(offsets(size(row_ids)+1));offsets(1)=1;count=0
    do i=1,size(row_ids)
      row=int(row_ids(i));do q=1,size(keys);key_row=int((keys(q)-1_int64)/int(n,int64))+1;if(key_row==row)count=count+1;enddo
      offsets(i+1)=count+1
    enddo
    allocate(columns(count));count=0
    do i=1,size(row_ids)
      row=int(row_ids(i))
      do q=1,size(keys);key_row=int((keys(q)-1_int64)/int(n,int64))+1
        if(key_row==row)then;count=count+1;columns(count)=int(modulo(keys(q)-1_int64,int(n,int64)))+1;endif
      enddo
    enddo
  end subroutine build_owned_csr
  recursive subroutine sort_int64(values,left,right)
    integer(int64),intent(inout)::values(:)
    integer,intent(in)::left,right
    integer::i,j
    integer(int64)::pivot,temp
    if(left>=right)return
    i=left;j=right;pivot=values((left+right)/2)
    do
      do while(values(i)<pivot);i=i+1;enddo;do while(values(j)>pivot);j=j-1;enddo
      if(i<=j)then;temp=values(i);values(i)=values(j);values(j)=temp;i=i+1;j=j-1;endif
      if(i>j)exit
    enddo
    if(left<j)call sort_int64(values,left,j);if(i<right)call sort_int64(values,i,right)
  end subroutine sort_int64
  logical function finite_matrix(values)
    complex(real64),intent(in)::values(:,:)
    finite_matrix=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_matrix
  logical function finite_rank3(values)
    complex(real64),intent(in)::values(:,:,:)
    finite_rank3=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_rank3
#endif
end module rt_dg_hybrid_structural_graph
