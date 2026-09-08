#include "config.h"
module rt_dg_hybrid_structural_graph
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::build_rt_dg_hybrid_structural_graph
contains
  subroutine build_rt_dg_hybrid_structural_graph(comm,global_count,row_ids,basis_values,metric_rows,kinetic_rows,&
      nonlocal_rows,local_rows,sipg_rows,hamiltonian_rows,position_rows,metric_offsets,metric_columns,&
      operator_offsets,operator_columns,ok,message)
    integer,intent(in)::comm,global_count
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::basis_values(:,:),metric_rows(:,:),kinetic_rows(:,:),nonlocal_rows(:,:),&
      local_rows(:,:),sipg_rows(:,:),hamiltonian_rows(:,:),position_rows(:,:,:)
    integer,allocatable,intent(out)::metric_offsets(:),metric_columns(:),operator_offsets(:),operator_columns(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::i,j,p,row,nowned,npoint,nactive,metric_raw_count,operator_raw_count,metric_unique_count,&
      operator_unique_count,ierr,local_bad,global_bad
    integer,allocatable::active(:)
    integer(int64),allocatable::metric_keys(:),operator_keys(:)
    ok=.false.;message='';nowned=size(row_ids);npoint=size(basis_values,2);local_bad=0
    if(global_count<1.or.size(basis_values,1)/=global_count.or.&
      any(shape(metric_rows)/=[nowned,global_count]).or.any(shape(kinetic_rows)/=[nowned,global_count]).or.&
      any(shape(nonlocal_rows)/=[nowned,global_count]).or.any(shape(local_rows)/=[nowned,global_count]).or.&
      any(shape(sipg_rows)/=[nowned,global_count]).or.any(shape(hamiltonian_rows)/=[nowned,global_count]).or.&
      any(shape(position_rows)/=[3,nowned,global_count]))local_bad=1
    if(any(row_ids<1_int64).or.any(row_ids>int(global_count,int64)))local_bad=1
    if(.not.finite_matrix(basis_values).or..not.finite_matrix(metric_rows).or..not.finite_matrix(kinetic_rows).or.&
      .not.finite_matrix(nonlocal_rows).or..not.finite_matrix(local_rows).or..not.finite_matrix(sipg_rows).or.&
      .not.finite_matrix(hamiltonian_rows).or..not.finite_rank3(position_rows))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid Hybrid structural graph inputs';return;endif
    if(int(global_count,int64)*int(global_count,int64)>huge(0_int64))then
      message='Hybrid structural graph key extent overflow';return
    endif
    allocate(active(global_count));metric_raw_count=0;operator_raw_count=0
    do p=1,npoint
      nactive=0
      do i=1,global_count
        if(basis_values(i,p)/=(0d0,0d0))then;nactive=nactive+1;active(nactive)=i;endif
      enddo
      operator_raw_count=operator_raw_count+nactive*(nactive+1)/2
    enddo
    do i=1,nowned;do j=1,global_count
      if(metric_rows(i,j)/=(0d0,0d0))metric_raw_count=metric_raw_count+1
      if(metric_rows(i,j)/=(0d0,0d0).or.kinetic_rows(i,j)/=(0d0,0d0).or.&
        nonlocal_rows(i,j)/=(0d0,0d0).or.local_rows(i,j)/=(0d0,0d0).or.sipg_rows(i,j)/=(0d0,0d0).or.&
        hamiltonian_rows(i,j)/=(0d0,0d0).or.any(position_rows(:,i,j)/=(0d0,0d0)))&
        operator_raw_count=operator_raw_count+1
    enddo;enddo
    allocate(metric_keys(metric_raw_count),operator_keys(operator_raw_count));metric_raw_count=0;operator_raw_count=0
    do p=1,npoint
      nactive=0
      do i=1,global_count
        if(basis_values(i,p)/=(0d0,0d0))then;nactive=nactive+1;active(nactive)=i;endif
      enddo
      do i=1,nactive;do j=i,nactive
        operator_raw_count=operator_raw_count+1
        operator_keys(operator_raw_count)=pair_key(active(i),active(j),global_count)
      enddo;enddo
    enddo
    do i=1,nowned
      row=int(row_ids(i))
      do j=1,global_count
        if(metric_rows(i,j)/=(0d0,0d0))then
          metric_raw_count=metric_raw_count+1;metric_keys(metric_raw_count)=pair_key(row,j,global_count)
        endif
        if(metric_rows(i,j)/=(0d0,0d0).or.kinetic_rows(i,j)/=(0d0,0d0).or.&
          nonlocal_rows(i,j)/=(0d0,0d0).or.local_rows(i,j)/=(0d0,0d0).or.sipg_rows(i,j)/=(0d0,0d0).or.&
          hamiltonian_rows(i,j)/=(0d0,0d0).or.any(position_rows(:,i,j)/=(0d0,0d0)))then
          operator_raw_count=operator_raw_count+1;operator_keys(operator_raw_count)=pair_key(row,j,global_count)
        endif
      enddo
    enddo
    call global_unique_keys(comm,metric_keys,metric_unique_count,ierr)
    if(ierr==MPI_SUCCESS)call global_unique_keys(comm,operator_keys,operator_unique_count,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Hybrid structural support exchange failed';return;endif
    call merge_key_sets(metric_keys,metric_unique_count,operator_keys,operator_unique_count)
    call build_owned_csr(global_count,row_ids,metric_keys,metric_unique_count,metric_offsets,metric_columns)
    call build_owned_csr(global_count,row_ids,operator_keys,operator_unique_count,operator_offsets,operator_columns)
    ok=.true.;message=''
#else
    ok=.false.;message='Hybrid structural graph requires MPI'
#endif
  end subroutine build_rt_dg_hybrid_structural_graph

#ifdef USE_MPI
  pure integer(int64) function pair_key(row,column,n) result(key)
    integer,intent(in)::row,column,n
    integer::lower,upper
    lower=min(row,column);upper=max(row,column);key=int(lower-1,int64)*int(n,int64)+int(upper,int64)
  end function pair_key

  subroutine global_unique_keys(comm,keys,unique_count,ierr)
    integer,intent(in)::comm
    integer(int64),allocatable,intent(inout)::keys(:)
    integer,intent(out)::unique_count,ierr
    integer::rank_count,p,local_count,global_count
    integer,allocatable::counts(:),displacements(:)
    integer(int64),allocatable::global_keys(:)
    call sort_int64(keys,1,size(keys));call unique_prefix(keys,local_count)
    call MPI_Comm_size(comm,rank_count,ierr);if(ierr/=MPI_SUCCESS)return
    allocate(counts(rank_count),displacements(rank_count))
    call MPI_Allgather(local_count,1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,ierr);if(ierr/=MPI_SUCCESS)return
    displacements(1)=0
    do p=2,rank_count;displacements(p)=displacements(p-1)+counts(p-1);enddo
    global_count=sum(counts);allocate(global_keys(global_count))
    call MPI_Allgatherv(keys,local_count,MPI_INTEGER8,global_keys,counts,displacements,MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call sort_int64(global_keys,1,size(global_keys));call unique_prefix(global_keys,unique_count)
    deallocate(keys);allocate(keys(unique_count));if(unique_count>0)keys=global_keys(:unique_count)
  end subroutine global_unique_keys

  subroutine merge_key_sets(metric_keys,metric_count,operator_keys,operator_count)
    integer(int64),intent(in)::metric_keys(:)
    integer,intent(in)::metric_count
    integer(int64),allocatable,intent(inout)::operator_keys(:)
    integer,intent(inout)::operator_count
    integer(int64),allocatable::combined(:)
    integer::combined_count
    allocate(combined(metric_count+operator_count))
    if(metric_count>0)combined(:metric_count)=metric_keys(:metric_count)
    if(operator_count>0)combined(metric_count+1:)=operator_keys(:operator_count)
    call sort_int64(combined,1,size(combined));call unique_prefix(combined,combined_count)
    deallocate(operator_keys);allocate(operator_keys(combined_count));if(combined_count>0)operator_keys=combined(:combined_count)
    operator_count=combined_count
  end subroutine merge_key_sets

  subroutine build_owned_csr(n,row_ids,keys,key_count,offsets,columns)
    integer,intent(in)::n,key_count
    integer(int64),intent(in)::row_ids(:),keys(:)
    integer,allocatable,intent(out)::offsets(:),columns(:)
    integer::i,q,row,lower,upper,count,edge
    allocate(offsets(size(row_ids)+1));offsets(1)=1;count=0
    do i=1,size(row_ids)
      row=int(row_ids(i))
      do q=1,key_count
        lower=int((keys(q)-1_int64)/int(n,int64))+1;upper=int(mod(keys(q)-1_int64,int(n,int64)))+1
        if(lower==row)count=count+1
        if(upper==row.and.upper/=lower)count=count+1
      enddo
      offsets(i+1)=count+1
    enddo
    allocate(columns(count));edge=0
    do i=1,size(row_ids)
      row=int(row_ids(i))
      do q=1,key_count
        lower=int((keys(q)-1_int64)/int(n,int64))+1;upper=int(mod(keys(q)-1_int64,int(n,int64)))+1
        if(upper==row.and.upper/=lower)then;edge=edge+1;columns(edge)=lower;endif
      enddo
      do q=1,key_count
        lower=int((keys(q)-1_int64)/int(n,int64))+1;upper=int(mod(keys(q)-1_int64,int(n,int64)))+1
        if(lower==row)then;edge=edge+1;columns(edge)=upper;endif
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
      do while(values(i)<pivot);i=i+1;enddo
      do while(values(j)>pivot);j=j-1;enddo
      if(i<=j)then;temp=values(i);values(i)=values(j);values(j)=temp;i=i+1;j=j-1;endif
      if(i>j)exit
    enddo
    if(left<j)call sort_int64(values,left,j)
    if(i<right)call sort_int64(values,i,right)
  end subroutine sort_int64

  subroutine unique_prefix(values,count)
    integer(int64),intent(inout)::values(:)
    integer,intent(out)::count
    integer::i
    count=0
    do i=1,size(values)
      if(i==1)then
        count=count+1;values(count)=values(i)
      elseif(values(i)/=values(i-1))then
        count=count+1;values(count)=values(i)
      endif
    enddo
  end subroutine unique_prefix

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
