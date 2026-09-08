#include "config.h"
module rt_dg_hybrid_sparse_projection
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::project_rt_dg_hybrid_sparse_edges,validate_rt_dg_hybrid_sparse_hermiticity
#ifdef USE_MPI
  type::key_value_set
    integer(int64),allocatable::keys(:)
    complex(real64),allocatable::values(:)
    integer(int64)::count=0_int64
  end type key_value_set
#endif
contains
  subroutine project_rt_dg_hybrid_sparse_edges(comm,global_count,row_ids,row_offsets,column_ids,grid_ids,&
      grid_weights,basis_values,potential_values,local_values,ok,message)
    integer,intent(in)::comm,global_count,row_offsets(:),column_ids(:)
    integer(int64),intent(in)::row_ids(:),grid_ids(:)
    real(real64),intent(in)::grid_weights(:),potential_values(:)
    complex(real64),intent(in)::basis_values(:,:)
    complex(real64),intent(out)::local_values(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::p,i,j,row,nactive,ierr,local_bad,global_bad,destination,nproc
    integer,allocatable::active(:),owners(:)
    type(key_value_set)::contributions
    complex(real64)::factor
    ok=.false.;message='';local_bad=0
    call validate_contract(global_count,row_ids,row_offsets,column_ids,grid_ids,grid_weights,basis_values,&
      potential_values,local_values,local_bad)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid sparse Hybrid projection contract';return;endif
    call build_row_owners(comm,global_count,row_ids,owners,ok,ierr)
    if(ierr/=MPI_SUCCESS.or..not.ok)then;message='sparse projection rows do not have unique owners';return;endif
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    allocate(active(global_count));local_values=(0d0,0d0)
    ! Bound the temporary projection map by one owner's rows.  Re-scanning the
    ! local grid avoids ever materializing a global edge catalog on a rank.
    do destination=1,nproc
      call reset_map(contributions,64)
      do p=1,size(grid_ids)
        nactive=0
        do i=1,global_count
          if(basis_values(i,p)/=(0d0,0d0))then;nactive=nactive+1;active(nactive)=i;endif
        enddo
        factor=cmplx(grid_weights(p)*potential_values(p),0d0,real64)
        do i=1,nactive
          if(owners(active(i))/=destination)cycle
          do j=1,nactive
            call add_value(contributions,directed_key(active(i),active(j),global_count),&
              factor*conjg(basis_values(active(i),p))*basis_values(active(j),p))
          enddo
        enddo
      enddo
      call route_contributions_to_single_owner(comm,global_count,destination,contributions,row_ids,row_offsets,&
        column_ids,local_values,ierr)
      if(ierr/=MPI_SUCCESS)exit
    enddo
    ok=ierr==MPI_SUCCESS
    if(ok)then;message='';else;message='sparse local-potential owner exchange failed';endif
#else
    ok=.false.;message='sparse Hybrid projection requires MPI'
#endif
  end subroutine project_rt_dg_hybrid_sparse_edges

  subroutine validate_rt_dg_hybrid_sparse_hermiticity(comm,global_count,row_ids,row_offsets,column_ids,values,&
      tolerance,ok,message)
    integer,intent(in)::comm,global_count,row_offsets(:),column_ids(:)
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::values(:)
    real(real64),intent(in)::tolerance
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::i,edge,row,nproc,p,destination,ierr,local_bad,global_bad,total_send,total_recv,location
    integer,allocatable::owners(:),send_counts(:),recv_counts(:),send_displacements(:),recv_displacements(:),cursor(:)
    integer(int64),allocatable::send_keys(:),recv_keys(:)
    complex(real64),allocatable::send_values(:),recv_values(:)
    real(real64)::local_defect,global_defect,local_scale,global_scale
    ok=.false.;message='';local_bad=0
    if(size(row_offsets)/=size(row_ids)+1.or.size(column_ids)/=size(values))local_bad=1
    if(size(row_offsets)>0)then
      if(row_offsets(1)/=1.or.row_offsets(size(row_offsets))/=size(values)+1)local_bad=1
    endif
    if(.not.ieee_is_finite(tolerance).or.tolerance<0d0.or.&
      .not.all(ieee_is_finite(real(values))).or..not.all(ieee_is_finite(aimag(values))))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid sparse Hermiticity contract';return;endif
    call build_row_owners(comm,global_count,row_ids,owners,ok,ierr)
    if(ierr/=MPI_SUCCESS.or..not.ok)then;message='sparse Hermitian rows do not have unique owners';return;endif
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    allocate(send_counts(nproc),recv_counts(nproc),send_displacements(nproc),recv_displacements(nproc),cursor(nproc))
    send_counts=0
    do i=1,size(row_ids)
      do edge=row_offsets(i),row_offsets(i+1)-1
        destination=owners(column_ids(edge));send_counts(destination)=send_counts(destination)+1
      enddo
    enddo
    call MPI_Alltoall(send_counts,1,MPI_INTEGER,recv_counts,1,MPI_INTEGER,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call make_displacements(send_counts,send_displacements,total_send)
    call make_displacements(recv_counts,recv_displacements,total_recv)
    allocate(send_keys(total_send),recv_keys(total_recv),send_values(total_send),recv_values(total_recv));cursor=send_displacements
    do i=1,size(row_ids)
      row=int(row_ids(i))
      do edge=row_offsets(i),row_offsets(i+1)-1
        destination=owners(column_ids(edge));cursor(destination)=cursor(destination)+1;p=cursor(destination)
        send_keys(p)=directed_key(column_ids(edge),row,global_count);send_values(p)=conjg(values(edge))
      enddo
    enddo
    call MPI_Alltoallv(send_keys,send_counts,send_displacements,MPI_INTEGER8,recv_keys,recv_counts,&
      recv_displacements,MPI_INTEGER8,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Alltoallv(send_values,send_counts,send_displacements,MPI_DOUBLE_COMPLEX,&
      recv_values,recv_counts,recv_displacements,MPI_DOUBLE_COMPLEX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='sparse Hermitian partner owner exchange failed';return;endif
    local_defect=0d0;local_scale=1d0
    do p=1,total_recv
      row=int((recv_keys(p)-1_int64)/int(global_count,int64))+1
      i=find_row_position(row_ids,row)
      if(i==0)then;local_bad=1;cycle;endif
      location=find_column(column_ids,row_offsets(i),row_offsets(i+1)-1,&
        int(modulo(recv_keys(p)-1_int64,int(global_count,int64)))+1)
      if(location==0)then;local_bad=1;cycle;endif
      local_scale=max(local_scale,abs(values(location)),abs(recv_values(p)))
      local_defect=max(local_defect,abs(values(location)-recv_values(p)))
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_defect,global_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_scale,global_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    ok=ierr==MPI_SUCCESS.and.global_bad==0.and.global_defect<=tolerance*global_scale
    if(ok)then;message='';else;message='sparse operator is non-Hermitian or missing conjugate partner';endif
#else
    ok=.false.;message='sparse Hermiticity validation requires MPI'
#endif
  end subroutine validate_rt_dg_hybrid_sparse_hermiticity

#ifdef USE_MPI
  subroutine validate_contract(n,row_ids,offsets,columns,grid_ids,weights,basis,potential,values,bad)
    integer,intent(in)::n,offsets(:),columns(:)
    integer(int64),intent(in)::row_ids(:),grid_ids(:)
    real(real64),intent(in)::weights(:),potential(:)
    complex(real64),intent(in)::basis(:,:),values(:)
    integer,intent(out)::bad
    integer::p
    bad=0
    if(n<1.or.size(offsets)/=size(row_ids)+1.or.size(values)/=size(columns).or.&
      size(weights)/=size(grid_ids).or.size(potential)/=size(grid_ids).or.any(shape(basis)/=[n,size(grid_ids)]))bad=1
    if(size(offsets)>0)then;if(offsets(1)/=1.or.offsets(size(offsets))/=size(columns)+1)bad=1;endif
    if(any(row_ids<1_int64).or.any(row_ids>int(n,int64)).or.any(columns<1).or.any(columns>n))bad=1
    do p=1,size(row_ids)
      if(offsets(p)<1.or.offsets(p+1)<offsets(p).or.offsets(p+1)>size(columns)+1)bad=1
      if(offsets(p+1)-offsets(p)>1)then
        if(any(columns(offsets(p)+1:offsets(p+1)-1)<=columns(offsets(p):offsets(p+1)-2)))bad=1
      endif
    enddo
    if(.not.all(ieee_is_finite(weights)).or..not.all(ieee_is_finite(potential)).or.&
      .not.all(ieee_is_finite(real(basis))).or..not.all(ieee_is_finite(aimag(basis))))bad=1
  end subroutine validate_contract
  pure integer(int64) function directed_key(row,column,n) result(key)
    integer,intent(in)::row,column,n
    key=int(row-1,int64)*int(n,int64)+int(column,int64)
  end function directed_key
  subroutine build_row_owners(comm,n,row_ids,owners,ok,ierr)
    integer,intent(in)::comm,n
    integer(int64),intent(in)::row_ids(:)
    integer,allocatable,intent(out)::owners(:)
    logical,intent(out)::ok
    integer,intent(out)::ierr
    integer::rank,i
    integer,allocatable::marks(:)
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    allocate(owners(n),marks(n));owners=0;marks=0
    do i=1,size(row_ids);owners(int(row_ids(i)))=rank+1;marks(int(row_ids(i)))=1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,owners,n,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(MPI_IN_PLACE,marks,n,MPI_INTEGER,MPI_SUM,comm,ierr)
    ok=ierr==MPI_SUCCESS.and.all(marks==1)
  end subroutine build_row_owners
  subroutine initialize_map(map,requested)
    type(key_value_set),intent(inout)::map
    integer,intent(in)::requested
    integer::capacity
    capacity=64;do while(capacity<requested.and.capacity<=huge(capacity)/2);capacity=capacity*2;enddo
    allocate(map%keys(capacity),map%values(capacity));map%keys=0_int64;map%values=(0d0,0d0);map%count=0_int64
  end subroutine initialize_map
  subroutine reset_map(map,requested)
    type(key_value_set),intent(inout)::map
    integer,intent(in)::requested
    if(allocated(map%keys))deallocate(map%keys,map%values)
    call initialize_map(map,requested)
  end subroutine reset_map
  subroutine add_value(map,key,value)
    type(key_value_set),intent(inout)::map
    integer(int64),intent(in)::key
    complex(real64),intent(in)::value
    integer::position,capacity
    capacity=size(map%keys)
    if(map%count*10_int64>=int(capacity,int64)*7_int64)then;call rehash_map(map,2*capacity);capacity=size(map%keys);endif
    position=int(modulo(key-1_int64,int(capacity,int64)))+1
    do
      if(map%keys(position)==0_int64)then
        map%keys(position)=key;map%values(position)=value;map%count=map%count+1_int64;return
      elseif(map%keys(position)==key)then;map%values(position)=map%values(position)+value;return;endif
      position=position+1;if(position>capacity)position=1
    enddo
  end subroutine add_value
  subroutine rehash_map(map,new_capacity)
    type(key_value_set),intent(inout)::map
    integer,intent(in)::new_capacity
    integer(int64),allocatable::old_keys(:)
    complex(real64),allocatable::old_values(:)
    integer::i
    call move_alloc(map%keys,old_keys);call move_alloc(map%values,old_values)
    allocate(map%keys(new_capacity),map%values(new_capacity));map%keys=0_int64;map%values=(0d0,0d0);map%count=0_int64
    do i=1,size(old_keys);if(old_keys(i)>0_int64)call add_value(map,old_keys(i),old_values(i));enddo
  end subroutine rehash_map
  subroutine route_contributions_to_single_owner(comm,n,destination,map,row_ids,offsets,columns,local_values,ierr)
    integer,intent(in)::comm,n,destination,offsets(:),columns(:)
    type(key_value_set),intent(in)::map
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(inout)::local_values(:)
    integer,intent(out)::ierr
    integer::rank,nproc,sender,count,i,q,status(MPI_STATUS_SIZE)
    integer(int64),allocatable::keys(:),received_keys(:)
    complex(real64),allocatable::values(:),received_values(:)
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    if(map%count>int(huge(0),int64))then;ierr=1;return;endif
    allocate(keys(int(map%count)),values(int(map%count)));q=0
    do i=1,size(map%keys);if(map%keys(i)>0_int64)then;q=q+1;keys(q)=map%keys(i);values(q)=map%values(i);endif;enddo
    do sender=0,nproc-1
      count=merge(size(keys),0,rank==sender)
      call MPI_Bcast(count,1,MPI_INTEGER,sender,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(sender==destination-1)then
        if(rank==destination-1)call accumulate_owner_values(n,keys,values,row_ids,offsets,columns,local_values,ierr)
      elseif(rank==sender)then
        call MPI_Send(keys,count,MPI_INTEGER8,destination-1,1800+sender,comm,ierr)
        if(ierr==MPI_SUCCESS)call MPI_Send(values,count,MPI_DOUBLE_COMPLEX,destination-1,1900+sender,comm,ierr)
      elseif(rank==destination-1)then
        allocate(received_keys(count),received_values(count))
        call MPI_Recv(received_keys,count,MPI_INTEGER8,sender,1800+sender,comm,status,ierr)
        if(ierr==MPI_SUCCESS)call MPI_Recv(received_values,count,MPI_DOUBLE_COMPLEX,sender,1900+sender,comm,status,ierr)
        if(ierr==MPI_SUCCESS)call accumulate_owner_values(n,received_keys,received_values,row_ids,offsets,columns,&
          local_values,ierr)
        deallocate(received_keys,received_values)
      endif
      if(ierr/=MPI_SUCCESS)return
    enddo
  end subroutine route_contributions_to_single_owner
  subroutine accumulate_owner_values(n,keys,values,row_ids,offsets,columns,local_values,ierr)
    integer,intent(in)::n,offsets(:),columns(:)
    integer(int64),intent(in)::keys(:),row_ids(:)
    complex(real64),intent(in)::values(:)
    complex(real64),intent(inout)::local_values(:)
    integer,intent(out)::ierr
    integer::p,row,i,location
    ierr=MPI_SUCCESS
    do p=1,size(keys)
      row=int((keys(p)-1_int64)/int(n,int64))+1;i=find_row_position(row_ids,row)
      if(i==0)then;ierr=1;return;endif
      location=find_column(columns,offsets(i),offsets(i+1)-1,int(modulo(keys(p)-1_int64,int(n,int64)))+1)
      if(location==0)then;ierr=1;return;endif
      local_values(location)=local_values(location)+values(p)
    enddo
  end subroutine accumulate_owner_values
  subroutine make_displacements(counts,displacements,total)
    integer,intent(in)::counts(:)
    integer,intent(out)::displacements(:),total
    integer::p
    displacements(1)=0;do p=2,size(counts);displacements(p)=displacements(p-1)+counts(p-1);enddo;total=sum(counts)
  end subroutine make_displacements
  pure integer function find_row_position(rows,target) result(location)
    integer(int64),intent(in)::rows(:)
    integer,intent(in)::target
    integer::i
    location=0;do i=1,size(rows);if(rows(i)==int(target,int64))then;location=i;return;endif;enddo
  end function find_row_position
  pure integer function find_column(columns,left,right,target) result(location)
    integer,intent(in)::columns(:),left,right,target
    integer::lo,hi,middle
    lo=left;hi=right;location=0
    do while(lo<=hi)
      middle=lo+(hi-lo)/2
      if(columns(middle)==target)then;location=middle;return
      elseif(columns(middle)<target)then;lo=middle+1;else;hi=middle-1;endif
    enddo
  end function find_column
#endif
end module rt_dg_hybrid_sparse_projection
