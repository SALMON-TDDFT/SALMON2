#include "config.h"
module dg_hybrid_sparse_operators
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi, only: MPI_Allreduce, MPI_BXOR, MPI_INTEGER, MPI_INTEGER8, MPI_MAX, MPI_SUCCESS, MPI_SUM
#endif
  implicit none
  private
  type,public::s_dg_hybrid_sparse_operators
    logical::valid=.false.;integer::global_count=0
    integer(int64)::selection_fingerprint=0_int64,window_fingerprint=0_int64,packet_fingerprint=0_int64,&
      complement_fingerprint=0_int64,metric_fingerprint=0_int64,position_convention_fingerprint=0_int64,&
      fingerprint=0_int64,persistent_bytes=0_int64,transient_peak_bytes=0_int64
    integer(int64),allocatable::owned_row_ids(:)
    integer,allocatable::row_offsets(:),column_ids(:)
    complex(real64),allocatable::metric_values(:),hamiltonian_values(:),position_values(:,:)
  end type s_dg_hybrid_sparse_operators
  type,public::s_dg_hybrid_coupling_envelope
    logical::valid=.false.;integer::global_count=0
    integer(int64)::structure_fingerprint=0_int64
    integer(int64),allocatable::owned_row_ids(:)
    integer,allocatable::row_offsets(:),column_ids(:)
    complex(real64),allocatable::component_values(:,:)
  contains
    procedure::find_edge=>find_envelope_edge
  end type s_dg_hybrid_coupling_envelope
  public::build_dg_hybrid_operator_envelope,set_dg_hybrid_operator_hermitian_edge
contains
  subroutine build_dg_hybrid_operator_envelope(icomm,global_count,owned_row_ids,basis_edges,nonlocal_edges,face_edges,&
      component_count,envelope,ok,message)
    integer,intent(in)::icomm,global_count,basis_edges(:,:),nonlocal_edges(:,:),face_edges(:,:),component_count
    integer(int64),intent(in)::owned_row_ids(:)
    type(s_dg_hybrid_coupling_envelope),intent(out)::envelope
    logical,intent(out)::ok;character(*),intent(out)::message
#ifdef USE_MPI
    integer::i,j,row,edge,nraw,nunique,local_bad,global_bad,ierr
    integer,allocatable::raw_offsets(:),raw_columns(:),cursor(:),unique_counts(:)
    integer(int64)::local_xor,global_xor,local_sum,global_sum,pair_hash
    ok=.false.;message='';local_bad=0
    if(global_count<1.or.component_count<1.or.size(basis_edges,1)/=2.or.size(nonlocal_edges,1)/=2.or.&
        size(face_edges,1)/=2.or.any(owned_row_ids<1_int64).or.any(owned_row_ids>int(global_count,int64)))local_bad=1
    do i=1,size(owned_row_ids);if(count(owned_row_ids==owned_row_ids(i))/=1)local_bad=1;enddo
    call validate_edges(basis_edges);call validate_edges(nonlocal_edges);call validate_edges(face_edges)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid row-owned operator coupling envelope';return;endif
    allocate(raw_offsets(size(owned_row_ids)+1),cursor(size(owned_row_ids)),unique_counts(size(owned_row_ids)))
    raw_offsets(1)=1
    do i=1,size(owned_row_ids)
      raw_offsets(i+1)=raw_offsets(i)+1+count(basis_edges(1,:)==int(owned_row_ids(i)))+&
        count(nonlocal_edges(1,:)==int(owned_row_ids(i)))+count(face_edges(1,:)==int(owned_row_ids(i)))
    enddo
    nraw=raw_offsets(size(raw_offsets))-1;allocate(raw_columns(nraw));cursor=raw_offsets(:size(owned_row_ids))
    do i=1,size(owned_row_ids);raw_columns(cursor(i))=int(owned_row_ids(i));cursor(i)=cursor(i)+1;enddo
    call insert_edges(basis_edges);call insert_edges(nonlocal_edges);call insert_edges(face_edges)
    do i=1,size(owned_row_ids)
      call sort_segment(raw_columns,raw_offsets(i),raw_offsets(i+1)-1)
      unique_counts(i)=1
      do j=raw_offsets(i)+1,raw_offsets(i+1)-1
        if(raw_columns(j)/=raw_columns(j-1))unique_counts(i)=unique_counts(i)+1
      enddo
    enddo
    allocate(envelope%owned_row_ids(size(owned_row_ids)),envelope%row_offsets(size(owned_row_ids)+1))
    envelope%owned_row_ids=owned_row_ids;envelope%row_offsets(1)=1
    do i=1,size(owned_row_ids);envelope%row_offsets(i+1)=envelope%row_offsets(i)+unique_counts(i);enddo
    nunique=envelope%row_offsets(size(envelope%row_offsets))-1;allocate(envelope%column_ids(nunique));edge=0
    do i=1,size(owned_row_ids)
      do j=raw_offsets(i),raw_offsets(i+1)-1
        if(j==raw_offsets(i))then;edge=edge+1;envelope%column_ids(edge)=raw_columns(j)
        else if(raw_columns(j)/=raw_columns(j-1))then;edge=edge+1;envelope%column_ids(edge)=raw_columns(j)
        endif
      enddo
    enddo
    allocate(envelope%component_values(component_count,nunique));envelope%component_values=(0d0,0d0)
    local_xor=0_int64;local_sum=0_int64
    do i=1,size(owned_row_ids);row=int(owned_row_ids(i))
      do edge=envelope%row_offsets(i),envelope%row_offsets(i+1)-1
        pair_hash=ieor(ishftc(int(row,int64),29),int(envelope%column_ids(edge),int64))
        pair_hash=ieor(pair_hash,ishftc(pair_hash,17));local_xor=ieor(local_xor,pair_hash);local_sum=local_sum+pair_hash
      enddo
    enddo
    call MPI_Allreduce(local_xor,global_xor,1,MPI_INTEGER8,MPI_BXOR,icomm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(local_sum,global_sum,1,MPI_INTEGER8,MPI_SUM,icomm,ierr);if(ierr/=MPI_SUCCESS)return
    envelope%structure_fingerprint=ieor(ishftc(global_xor,13),global_sum)
    envelope%structure_fingerprint=ieor(envelope%structure_fingerprint,int(global_count,int64))
    if(envelope%structure_fingerprint==0_int64)envelope%structure_fingerprint=1_int64
    envelope%global_count=global_count;envelope%valid=.true.;ok=.true.
  contains
    subroutine validate_edges(edges)
      integer,intent(in)::edges(:,:);integer::q
      do q=1,size(edges,2)
        if(edges(1,q)<1.or.edges(1,q)>global_count.or.edges(2,q)<1.or.edges(2,q)>global_count.or.&
            find_owned(edges(1,q),owned_row_ids)==0)local_bad=1
      enddo
    end subroutine validate_edges
    subroutine insert_edges(edges)
      integer,intent(in)::edges(:,:);integer::q,p
      do q=1,size(edges,2);p=find_owned(edges(1,q),owned_row_ids);raw_columns(cursor(p))=edges(2,q);cursor(p)=cursor(p)+1;enddo
    end subroutine insert_edges
#else
    ok=.false.;message='row-owned operator envelope requires MPI'
#endif
  end subroutine build_dg_hybrid_operator_envelope

  subroutine set_dg_hybrid_operator_hermitian_edge(icomm,envelope,component,row,column,value,ok,message)
    integer,intent(in)::icomm,component,row,column
    type(s_dg_hybrid_coupling_envelope),intent(inout)::envelope
    complex(real64),intent(in)::value
    logical,intent(out)::ok;character(*),intent(out)::message
#ifdef USE_MPI
    integer::position,local_updates,global_updates,ierr
    ok=.false.;message='';local_updates=0
    if(.not.envelope%valid.or.component<1.or.component>size(envelope%component_values,1).or.&
        .not.ieee_is_finite(real(value,real64)).or..not.ieee_is_finite(aimag(value)))then
      message='invalid Hermitian operator component';return
    endif
    position=envelope%find_edge(row,column)
    if(position>0)then;envelope%component_values(component,position)=value;local_updates=local_updates+1;endif
    if(row/=column)then;position=envelope%find_edge(column,row)
      if(position>0)then;envelope%component_values(component,position)=conjg(value);local_updates=local_updates+1;endif
    endif
    call MPI_Allreduce(local_updates,global_updates,1,MPI_INTEGER,MPI_SUM,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_updates/=merge(1,2,row==column))then
      message='Hermitian operator edge or reverse is outside distributed envelope';return
    endif
    ok=.true.
#else
    ok=.false.;message='Hermitian operator update requires MPI'
#endif
  end subroutine set_dg_hybrid_operator_hermitian_edge

  integer function find_envelope_edge(envelope,row,column) result(position)
    class(s_dg_hybrid_coupling_envelope),intent(in)::envelope;integer,intent(in)::row,column
    integer::local_row,edge
    position=0;if(.not.envelope%valid)return;local_row=find_owned(row,envelope%owned_row_ids);if(local_row==0)return
    do edge=envelope%row_offsets(local_row),envelope%row_offsets(local_row+1)-1
      if(envelope%column_ids(edge)==column)then;position=edge;return;endif
    enddo
  end function find_envelope_edge
  integer function find_owned(row,owned_ids) result(position)
    integer,intent(in)::row;integer(int64),intent(in)::owned_ids(:);integer::i
    position=0;do i=1,size(owned_ids);if(owned_ids(i)==int(row,int64))then;position=i;return;endif;enddo
  end function find_owned
  subroutine sort_segment(values,first,last)
    integer,intent(inout)::values(:);integer,intent(in)::first,last;integer::i,j,key
    do i=first+1,last;key=values(i);j=i-1
      do while(j>=first);if(values(j)<=key)exit;values(j+1)=values(j);j=j-1;enddo;values(j+1)=key
    enddo
  end subroutine sort_segment
end module dg_hybrid_sparse_operators
