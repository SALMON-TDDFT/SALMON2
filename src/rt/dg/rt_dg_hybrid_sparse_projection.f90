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
    integer::p,edge,row,local_edge_count,global_edge_count,nproc,ierr,local_bad,global_bad
    integer,allocatable::edge_counts(:),edge_displacements(:),global_edge_columns(:)
    integer(int64),allocatable::local_edge_rows(:),global_edge_rows(:)
    complex(real64),allocatable::partial_values(:)
    ok=.false.;message='';local_edge_count=size(column_ids);local_bad=0
    if(global_count<1.or.size(row_offsets)/=size(row_ids)+1.or.size(local_values)/=local_edge_count.or.&
      size(grid_weights)/=size(grid_ids).or.size(potential_values)/=size(grid_ids).or.&
      any(shape(basis_values)/=[global_count,size(grid_ids)]))local_bad=1
    if(size(row_offsets)>0)then
      if(row_offsets(1)/=1.or.row_offsets(size(row_offsets))/=local_edge_count+1)local_bad=1
    endif
    if(any(row_ids<1_int64).or.any(row_ids>int(global_count,int64)).or.&
      any(column_ids<1).or.any(column_ids>global_count))local_bad=1
    do p=1,size(row_ids)
      if(row_offsets(p)<1.or.row_offsets(p+1)<row_offsets(p).or.row_offsets(p+1)>local_edge_count+1)local_bad=1
      if(row_offsets(p+1)-row_offsets(p)>1)then
        if(any(column_ids(row_offsets(p)+1:row_offsets(p+1)-1)<=&
          column_ids(row_offsets(p):row_offsets(p+1)-2)))local_bad=1
      endif
    enddo
    if(.not.all(ieee_is_finite(grid_weights)).or..not.all(ieee_is_finite(potential_values)).or.&
      .not.all(ieee_is_finite(real(basis_values))).or..not.all(ieee_is_finite(aimag(basis_values))))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid sparse Hybrid projection contract';return;endif
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    allocate(edge_counts(nproc),edge_displacements(nproc))
    call MPI_Allgather(local_edge_count,1,MPI_INTEGER,edge_counts,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='sparse projection edge-count exchange failed';return;endif
    edge_displacements(1)=0
    do p=2,nproc;edge_displacements(p)=edge_displacements(p-1)+edge_counts(p-1);enddo
    global_edge_count=sum(edge_counts)
    allocate(local_edge_rows(local_edge_count),global_edge_rows(global_edge_count),&
      global_edge_columns(global_edge_count),partial_values(global_edge_count))
    do p=1,size(row_ids)
      local_edge_rows(row_offsets(p):row_offsets(p+1)-1)=row_ids(p)
    enddo
    call MPI_Allgatherv(local_edge_rows,local_edge_count,MPI_INTEGER8,global_edge_rows,edge_counts,&
      edge_displacements,MPI_INTEGER8,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allgatherv(column_ids,local_edge_count,MPI_INTEGER,global_edge_columns,&
      edge_counts,edge_displacements,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='sparse projection graph exchange failed';return;endif
    partial_values=(0d0,0d0)
    do edge=1,global_edge_count
      row=int(global_edge_rows(edge))
      do p=1,size(grid_ids)
        partial_values(edge)=partial_values(edge)+grid_weights(p)*conjg(basis_values(row,p))*&
          basis_values(global_edge_columns(edge),p)*potential_values(p)
      enddo
    enddo
    call MPI_Reduce_scatter(partial_values,local_values,edge_counts,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    ok=ierr==MPI_SUCCESS
    if(ok)then;message='';else;message='sparse local-potential projection failed';endif
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
    integer::nproc,p,i,edge,q,row,total,ierr,local_bad,global_bad
    integer,allocatable::counts(:),displacements(:),all_columns(:)
    integer(int64),allocatable::local_rows(:),all_rows(:),directed_keys(:)
    complex(real64),allocatable::all_values(:)
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
    call MPI_Comm_size(comm,nproc,ierr);allocate(counts(nproc),displacements(nproc))
    call MPI_Allgather(size(values),1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,ierr);if(ierr/=MPI_SUCCESS)return
    displacements(1)=0;do p=2,nproc;displacements(p)=displacements(p-1)+counts(p-1);enddo;total=sum(counts)
    allocate(local_rows(size(values)),all_rows(total),all_columns(total),all_values(total))
    do i=1,size(row_ids);local_rows(row_offsets(i):row_offsets(i+1)-1)=row_ids(i);enddo
    call MPI_Allgatherv(local_rows,size(values),MPI_INTEGER8,all_rows,counts,displacements,MPI_INTEGER8,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allgatherv(column_ids,size(values),MPI_INTEGER,all_columns,counts,displacements,&
      MPI_INTEGER,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allgatherv(values,size(values),MPI_DOUBLE_COMPLEX,all_values,counts,displacements,&
      MPI_DOUBLE_COMPLEX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='sparse Hermitian partner exchange failed';return;endif
    allocate(directed_keys(total))
    do edge=1,total
      directed_keys(edge)=(all_rows(edge)-1_int64)*int(global_count,int64)+int(all_columns(edge),int64)
    enddo
    call sort_directed_edges(directed_keys,all_values,1,total)
    do edge=2,total
      if(directed_keys(edge)==directed_keys(edge-1))local_bad=1
    enddo
    local_defect=0d0;local_scale=1d0
    do edge=1,total
      row=int((directed_keys(edge)-1_int64)/int(global_count,int64))+1
      i=int(mod(directed_keys(edge)-1_int64,int(global_count,int64)))+1
      q=find_directed_edge(directed_keys,(int(i-1,int64)*int(global_count,int64)+int(row,int64)))
      if(q==0)then;local_bad=1;cycle;endif
      local_scale=max(local_scale,abs(all_values(edge)),abs(all_values(q)))
      local_defect=max(local_defect,abs(all_values(edge)-conjg(all_values(q))))
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
  recursive subroutine sort_directed_edges(keys,values,left,right)
    integer(int64),intent(inout)::keys(:)
    complex(real64),intent(inout)::values(:)
    integer,intent(in)::left,right
    integer::i,j
    integer(int64)::pivot,key_temp
    complex(real64)::value_temp
    if(left>=right)return
    i=left;j=right;pivot=keys((left+right)/2)
    do
      do while(keys(i)<pivot);i=i+1;enddo
      do while(keys(j)>pivot);j=j-1;enddo
      if(i<=j)then
        key_temp=keys(i);keys(i)=keys(j);keys(j)=key_temp
        value_temp=values(i);values(i)=values(j);values(j)=value_temp
        i=i+1;j=j-1
      endif
      if(i>j)exit
    enddo
    if(left<j)call sort_directed_edges(keys,values,left,j)
    if(i<right)call sort_directed_edges(keys,values,i,right)
  end subroutine sort_directed_edges

  pure integer function find_directed_edge(keys,target) result(location)
    integer(int64),intent(in)::keys(:),target
    integer::left,right,middle
    left=1;right=size(keys);location=0
    do while(left<=right)
      middle=left+(right-left)/2
      if(keys(middle)==target)then
        location=middle;return
      elseif(keys(middle)<target)then
        left=middle+1
      else
        right=middle-1
      endif
    enddo
  end function find_directed_edge
#endif
end module rt_dg_hybrid_sparse_projection
