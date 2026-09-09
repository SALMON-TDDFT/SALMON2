#include "config.h"
program test_rt_dg_hybrid_sparse_projection_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use rt_dg_hybrid_structural_graph,only:build_rt_dg_hybrid_structural_graph,&
    checked_rt_dg_hybrid_structural_capacity,collective_rt_dg_hybrid_structural_capacity_status,&
    checked_rt_dg_hybrid_global_square_extent
  use rt_dg_hybrid_sparse_projection,only:project_rt_dg_hybrid_sparse_edges,&
    validate_rt_dg_hybrid_sparse_hermiticity,checked_rt_dg_hybrid_projection_capacity
  implicit none
  integer::comm,rank,nproc,ierr
  logical::ok
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  call exercise_capacity_boundary
  call exercise_structural_capacity_boundary
  call exercise_zero_global_count
  call exercise_structural_graph
  call exercise_repeated_support_scaling
  call exercise_sparse_projection
  call exercise_invalid_projection_contracts
  call exercise_missing_edge_collective_failure
  if(rank==0)write(*,'(a,i0,a)')'PASS structural Hybrid sparse projection on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine exercise_capacity_boundary
    integer::next_capacity
    logical::capacity_ok
    call checked_rt_dg_hybrid_projection_capacity(64,next_capacity,capacity_ok)
    call require(capacity_ok.and.next_capacity==128,'valid sparse projection hash growth was rejected')
    call checked_rt_dg_hybrid_projection_capacity(huge(0)/2+1,next_capacity,capacity_ok)
    call require(.not.capacity_ok,'sparse projection hash capacity overflow was not rejected')
  end subroutine exercise_capacity_boundary
  subroutine exercise_structural_capacity_boundary
    integer::requested
    logical::capacity_ok
    character(256)::capacity_message
    call checked_rt_dg_hybrid_structural_capacity(4,50,requested,capacity_ok)
    call require(capacity_ok.and.requested==200,'valid structural hash initial capacity was rejected')
    call checked_rt_dg_hybrid_structural_capacity(4,huge(0)/4+1,requested,capacity_ok)
    call require(.not.capacity_ok,'structural hash initial capacity overflow was not rejected')
    call checked_rt_dg_hybrid_structural_capacity(0,50,requested,capacity_ok)
    call require(.not.capacity_ok,'zero structural capacity multiplier was not rejected safely')
    call checked_rt_dg_hybrid_structural_capacity(-1,50,requested,capacity_ok)
    call require(.not.capacity_ok,'negative structural capacity multiplier was not rejected safely')
    call checked_rt_dg_hybrid_global_square_extent(0,capacity_ok)
    call require(.not.capacity_ok,'zero structural global extent was not rejected before division')
    call collective_rt_dg_hybrid_structural_capacity_status(comm,rank/=0,capacity_ok,capacity_message)
    call require(.not.capacity_ok.and.index(capacity_message,'capacity')>0,&
      'one-rank structural capacity failure was not collectively rejected')
  end subroutine exercise_structural_capacity_boundary
  subroutine exercise_zero_global_count
    integer(int64)::empty_rows(0)
    integer,allocatable::metric_offsets(:),metric_columns(:),operator_offsets(:),operator_columns(:)
    complex(real64)::empty_basis(0,0),empty_matrix(0,0),empty_position(3,0,0)
    logical::graph_ok
    character(256)::graph_message
    call build_rt_dg_hybrid_structural_graph(comm,0,empty_rows,empty_basis,empty_matrix,empty_matrix,empty_matrix,&
      empty_matrix,empty_matrix,empty_matrix,empty_position,metric_offsets,metric_columns,operator_offsets,&
      operator_columns,graph_ok,graph_message)
    call require(.not.graph_ok.and.index(graph_message,'invalid Hybrid structural graph inputs')>0,&
      'zero global_count was not safely rejected before division')
  end subroutine exercise_zero_global_count
  subroutine exercise_invalid_projection_contracts
    integer::n,owned,row
    integer(int64),allocatable::rows(:),empty_grid(:)
    integer,allocatable::offsets(:),columns(:)
    real(real64),allocatable::empty_real(:)
    complex(real64),allocatable::empty_basis(:,:),values(:)
    logical::projection_ok
    character(256)::projection_message
    n=nproc+1;owned=merge(2,1,rank==0);row=rank+2
    allocate(rows(owned),offsets(owned+1),columns(owned),values(owned),empty_grid(0),empty_real(0),empty_basis(n,0))
    if(rank==0)then;rows=[1_int64,2_int64];offsets=[1,2,3];columns=[1,2]
    else;rows=[int(row,int64)];offsets=[1,2];columns=[row];endif
    values=(0d0,0d0)
    if(rank==0)offsets(2)=0
    call validate_rt_dg_hybrid_sparse_hermiticity(comm,n,rows,offsets,columns,values,1d-12,&
      projection_ok,projection_message)
    call require(.not.projection_ok.and.index(projection_message,'contract')>0,&
      'Hermiticity entry accepted one-rank corrupt CSR offset')
    call project_rt_dg_hybrid_sparse_edges(comm,n,rows,offsets,columns,empty_grid,empty_real,empty_basis,&
      empty_real,values,projection_ok,projection_message)
    call require(.not.projection_ok.and.index(projection_message,'contract')>0,&
      'one-rank corrupt CSR offset was not collectively rejected')
    if(rank==0)then;offsets=[1,2,3];columns=[1,n+1]
    else;offsets=[1,2];columns=[row];endif
    call validate_rt_dg_hybrid_sparse_hermiticity(comm,n,rows,offsets,columns,values,1d-12,&
      projection_ok,projection_message)
    call require(.not.projection_ok.and.index(projection_message,'contract')>0,&
      'Hermiticity entry accepted one-rank out-of-range CSR column')
    call project_rt_dg_hybrid_sparse_edges(comm,n,rows,offsets,columns,empty_grid,empty_real,empty_basis,&
      empty_real,values,projection_ok,projection_message)
    call require(.not.projection_ok.and.index(projection_message,'contract')>0,&
      'one-rank out-of-range CSR column was not collectively rejected')
    deallocate(rows,offsets,columns,values)
    owned=merge(3,1,rank==0);allocate(rows(owned),offsets(owned+1),columns(owned),values(owned))
    if(rank==0)then;rows=[1_int64,1_int64,2_int64];offsets=[1,2,3,4];columns=[1,1,2]
    else;rows=[int(rank+2,int64)];offsets=[1,2];columns=[rank+2];endif
    values=(0d0,0d0)
    call validate_rt_dg_hybrid_sparse_hermiticity(comm,n,rows,offsets,columns,values,1d-12,&
      projection_ok,projection_message)
    call require(.not.projection_ok.and.index(projection_message,'unique owners')>0,&
      'Hermiticity entry accepted duplicate owned CSR row')
    call project_rt_dg_hybrid_sparse_edges(comm,n,rows,offsets,columns,empty_grid,empty_real,empty_basis,&
      empty_real,values,projection_ok,projection_message)
    call require(.not.projection_ok.and.index(projection_message,'unique owners')>0,&
      'duplicate owned CSR row was not collectively rejected')
    deallocate(rows,offsets,columns,values)
    owned=1;allocate(rows(owned),offsets(owned+1),columns(owned),values(owned))
    offsets=1;values=(0d0,0d0)
    if(owned>0)then
      offsets=[1,2]
      rows(1)=merge(int(n,int64),int(rank+2,int64),rank==0)
      columns(1)=int(rows(1))
    endif
    call project_rt_dg_hybrid_sparse_edges(comm,n,rows,offsets,columns,empty_grid,empty_real,empty_basis,&
      empty_real,values,projection_ok,projection_message)
    call require(.not.projection_ok.and.index(projection_message,'unique owners')>0,&
      'missing/mismatched CSR row owner was not collectively rejected')
  end subroutine exercise_invalid_projection_contracts
  subroutine exercise_repeated_support_scaling
    integer,parameter::n=400,np_global=16000
    integer::nowned,npoint,row,point,i,j,k,slot,fragment,first,nactive,owner,local_nnz,min_nnz,max_nnz,edge
    integer(int64)::local_unique,peak_workspace,global_unique,global_peak
    integer(int64),allocatable::rows(:),grid_ids(:)
    integer,allocatable::mo(:),mc(:),oo(:),oc(:)
    real(real64),allocatable::weights(:),potential(:)
    complex(real64),allocatable::basis(:,:),metric(:,:),zero2(:,:),position(:,:,:),actual(:),dense_local(:,:),dense(:,:)
    logical::scale_ok,projection_ok,large_values_ok,zero_row_ok
    nowned=0
    do row=1,n
      fragment=(row-1)/50;owner=min(nproc-1,fragment*nproc/8)
      if(owner==rank)nowned=nowned+1
    enddo
    npoint=0
    do point=1,np_global
      fragment=mod(point-1,8);owner=min(nproc-1,fragment*nproc/8)
      if(owner==rank)npoint=npoint+1
    enddo
    allocate(rows(nowned),grid_ids(npoint),weights(npoint),potential(npoint),basis(n,npoint),metric(nowned,n),&
      zero2(nowned,n),position(3,nowned,n))
    i=0
    do row=1,n
      fragment=(row-1)/50;owner=min(nproc-1,fragment*nproc/8)
      if(owner==rank)then;i=i+1;rows(i)=row;endif
    enddo
    if(nproc==8)then
      call require(nowned==50,'eight-rank scaling does not own one 50-basis fragment per rank')
      call require(all(rows==[(int(rank*50+i,int64),i=1,50)]),'eight-rank fragment rows are not contiguous')
    endif
    basis=(0d0,0d0);metric=(0d0,0d0);zero2=(0d0,0d0);position=(0d0,0d0)
    do i=1,nowned;metric(i,int(rows(i)))=(1d0,0d0);enddo
    slot=0
    do point=1,np_global
      fragment=mod(point-1,8);owner=min(nproc-1,fragment*nproc/8)
      if(owner/=rank)cycle
      slot=slot+1;grid_ids(slot)=point;weights(slot)=0.25d0+1d-5*point
      potential(slot)=(-1d0)**point*(0.4d0+2d-5*point)
      first=fragment*50+mod((point-1)/8,46)+1
      nactive=merge(1,4,fragment==7)
      do j=first,first+nactive-1
        basis(j,slot)=cmplx(1d0+0.01d0*j+1d-6*point,0.02d0*j-2d-6*point,real64)
      enddo
    enddo
    call build_rt_dg_hybrid_structural_graph(comm,n,rows,basis,metric,zero2,zero2,zero2,zero2,zero2,position,&
      mo,mc,oo,oc,scale_ok,message,local_unique_candidates=local_unique,peak_workspace_keys=peak_workspace)
    call require(scale_ok,'repeated-support structural graph failed: '//trim(message))
    local_nnz=size(oc)
    call MPI_Allreduce(local_unique,global_unique,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
    call MPI_Allreduce(peak_workspace,global_peak,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    call require(global_unique<20000_int64,'localized support graph grew with repeated grid points')
    call require(global_peak<50000_int64,'structural graph workspace is not bounded by unique sparse support')
    call require(int(local_nnz,int64)<=global_unique+int(n,int64),&
      'owner-local CSR exceeds point support plus metric diagonal')
    call MPI_Allreduce(local_nnz,min_nnz,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(local_nnz,max_nnz,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(nproc==8)call require(min_nnz<max_nnz,'unequal fragment sparsity was not exercised')

    allocate(actual(size(oc)),dense_local(n,n),dense(n,n));dense_local=(0d0,0d0)
    do slot=1,npoint
      do j=1,n
        if(basis(j,slot)==(0d0,0d0))cycle
        do k=1,n
          if(basis(k,slot)==(0d0,0d0))cycle
          dense_local(j,k)=dense_local(j,k)+weights(slot)*potential(slot)*conjg(basis(j,slot))*basis(k,slot)
        enddo
      enddo
    enddo
    call MPI_Allreduce(dense_local,dense,n*n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call project_rt_dg_hybrid_sparse_edges(comm,n,rows,oo,oc,grid_ids,weights,basis,potential,actual,&
      projection_ok,message)
    call require(projection_ok,'large owner-local sparse projection failed: '//trim(message))
    large_values_ok=.true.
    do i=1,nowned
      row=int(rows(i))
      do edge=oo(i),oo(i+1)-1
        large_values_ok=large_values_ok.and.abs(actual(edge)-dense(row,oc(edge)))<1d-11
      enddo
    enddo
    call require(large_values_ok,'large sparse projection differs from dense oracle')
    zero_row_ok=.true.
    if(any(rows==400_int64))then
      i=0;do j=1,nowned;if(rows(j)==400_int64)i=j;enddo
      zero_row_ok=all(abs(actual(oo(i):oo(i+1)-1))<1d-15)
    endif
    call require(zero_row_ok,'zero-support basis row projected a nonzero value')
  end subroutine exercise_repeated_support_scaling

  subroutine exercise_structural_graph
    integer,parameter::r=3,g=4
    integer::i,j,row,nowned,npoint,point,metric_nnz,operator_nnz,edge
    integer(int64),allocatable::row_ids(:),grid_ids(:)
    integer,allocatable::metric_offsets(:),metric_columns(:),operator_offsets(:),operator_columns(:)
    complex(real64),allocatable::basis(:,:),metric(:,:),kinetic(:,:),nonlocal(:,:),local(:,:),sipg(:,:),h(:,:),&
      position(:,:,:),values(:),dense_local(:,:),dense(:,:),x(:),sparse_action(:),metric_action(:)
    real(real64),allocatable::weights(:)
    logical::found_12,found_21,found_23,found_32,found_13,found_31,sorted_ok,operator_action_ok,metric_action_ok
    nowned=count([(mod(row-1,nproc)==rank,row=1,r)])
    npoint=count([(mod(point-1,nproc)==rank,point=1,g)])
    allocate(row_ids(nowned),grid_ids(npoint),weights(npoint),basis(r,npoint),metric(nowned,r),kinetic(nowned,r),&
      nonlocal(nowned,r),local(nowned,r),sipg(nowned,r),h(nowned,r),position(3,nowned,r))
    i=0
    do row=1,r;if(mod(row-1,nproc)==rank)then;i=i+1;row_ids(i)=row;endif;enddo
    basis=(0d0,0d0);weights=1d0;i=0
    do point=1,g
      if(mod(point-1,nproc)/=rank)cycle
      i=i+1;grid_ids(i)=point
      select case(point)
      case(1);basis(:,i)=[(1d0,0d0),(0d0,2d0),(0d0,0d0)]
      case(2);basis(:,i)=[(1d0,0d0),(0d0,-2d0),(0d0,0d0)]
      case(3);basis(:,i)=[(0d0,0d0),(1d-30,0d0),(3d0,0d0)]
      case(4);basis(:,i)=[(0d0,0d0),(0d0,0d0),(1d0,0d0)]
      end select
    enddo
    metric=(0d0,0d0);kinetic=(0d0,0d0);nonlocal=(0d0,0d0);local=(0d0,0d0)
    sipg=(0d0,0d0);h=(0d0,0d0);position=(0d0,0d0)
    do i=1,nowned
      row=int(row_ids(i));metric(i,row)=(1d0,0d0)
      if(row==1)kinetic(i,3)=(1d-30,0d0)
    enddo
    call build_rt_dg_hybrid_structural_graph(comm,r,row_ids,basis,metric,kinetic,nonlocal,local,sipg,h,position,&
      metric_offsets,metric_columns,operator_offsets,operator_columns,ok,message)
    call require(ok,'structural graph construction failed: '//trim(message))
    metric_nnz=size(metric_columns);operator_nnz=size(operator_columns)
    call require(metric_nnz==nowned,'identity metric CSR changed extent')
    found_12=.false.;found_21=.false.;found_23=.false.;found_32=.false.;found_13=.false.;found_31=.false.;sorted_ok=.true.
    do i=1,nowned
      row=int(row_ids(i))
      do edge=operator_offsets(i),operator_offsets(i+1)-1
        j=operator_columns(edge)
        if(row==1.and.j==2)found_12=.true.;if(row==2.and.j==1)found_21=.true.
        if(row==2.and.j==3)found_23=.true.;if(row==3.and.j==2)found_32=.true.
        if(row==1.and.j==3)found_13=.true.;if(row==3.and.j==1)found_31=.true.
        if(edge>operator_offsets(i))sorted_ok=sorted_ok.and.j>operator_columns(edge-1)
      enddo
    enddo
    call require(sorted_ok,'operator CSR is not sorted')
    call require(any_rank(found_12).and.any_rank(found_21),'cancelled pointwise basis edge lost Hermitian closure')
    call require(any_rank(found_23).and.any_rank(found_32),'tiny basis support below old cutoff was dropped')
    call require(any_rank(found_13).and.any_rank(found_31),'asymmetric tiny fixed edge lost reverse partner')
    allocate(values(operator_nnz),dense_local(r,r),dense(r,r),x(r),sparse_action(nowned),metric_action(nowned))
    values=(0d0,0d0);dense_local=(0d0,0d0)
    x=[(0.25d0,-0.5d0),(-0.75d0,0.125d0),(1.5d0,0.375d0)]
    do i=1,nowned
      row=int(row_ids(i))
      do edge=operator_offsets(i),operator_offsets(i+1)-1
        j=operator_columns(edge);values(edge)=cmplx(real(row+j,real64),real(row-j,real64),real64)
        dense_local(row,j)=values(edge)
      enddo
    enddo
    call validate_rt_dg_hybrid_sparse_hermiticity(comm,r,row_ids,operator_offsets,operator_columns,values,&
      1d-14,ok,message)
    call require(ok,'Hermitian structural values were rejected: '//trim(message))
    call MPI_Allreduce(dense_local,dense,r*r,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    sparse_action=(0d0,0d0);metric_action=(0d0,0d0);operator_action_ok=.true.;metric_action_ok=.true.
    do i=1,nowned
      row=int(row_ids(i))
      do edge=operator_offsets(i),operator_offsets(i+1)-1
        sparse_action(i)=sparse_action(i)+values(edge)*x(operator_columns(edge))
      enddo
      do edge=metric_offsets(i),metric_offsets(i+1)-1
        metric_action(i)=metric_action(i)+metric(i,metric_columns(edge))*x(metric_columns(edge))
      enddo
      operator_action_ok=operator_action_ok.and.&
        abs(sparse_action(i)-dot_product(conjg(dense(row,:)),x))<1d-13
      metric_action_ok=metric_action_ok.and.abs(metric_action(i)-dot_product(conjg(metric(i,:)),x))<1d-13
    enddo
    call require(operator_action_ok,'sparse structural operator action differs from dense oracle')
    call require(metric_action_ok,'sparse structural metric action differs from dense oracle')
    do i=1,nowned
      if(row_ids(i)==1_int64)then
        edge=find_edge(operator_offsets,operator_columns,i,3);values(edge)=values(edge)+(1d-12,0d0)
      endif
    enddo
    call validate_rt_dg_hybrid_sparse_hermiticity(comm,r,row_ids,operator_offsets,operator_columns,values,&
      1d-14,ok,message)
    call require(.not.ok,'asymmetric near-threshold operator pair was accepted')
  end subroutine exercise_structural_graph

  subroutine exercise_sparse_projection
    integer,parameter::r=3,g=5
    integer::row,i,edge,nowned,npoint,point,nnz
    integer(int64),allocatable::row_ids(:),grid_ids(:)
    integer,allocatable::offsets(:),columns(:)
    real(real64),allocatable::weights(:),potential(:)
    complex(real64),allocatable::basis(:,:),actual(:),dense_local(:,:),dense(:,:)
    logical::values_ok
    nowned=count([(mod(row-1,nproc)==rank,row=1,r)])
    npoint=count([(mod(point-1,nproc)==rank,point=1,g)])
    allocate(row_ids(nowned),offsets(nowned+1));i=0;nnz=0;offsets(1)=1
    do row=1,r
      if(mod(row-1,nproc)/=rank)cycle
      i=i+1;row_ids(i)=row;if(row<3)nnz=nnz+2;offsets(i+1)=nnz+1
    enddo
    allocate(columns(nnz));edge=0
    do i=1,nowned
      if(row_ids(i)<3)then;columns(edge+1:edge+2)=[1,2];edge=edge+2;endif
    enddo
    allocate(grid_ids(npoint),weights(npoint),potential(npoint),basis(r,npoint),actual(nnz),&
      dense_local(r,r),dense(r,r))
    i=0
    do point=1,g
      if(mod(point-1,nproc)/=rank)cycle
      i=i+1;grid_ids(i)=point;weights(i)=0.5d0+0.1d0*point;potential(i)=(-1d0)**point*(0.2d0+0.07d0*point)
      basis(1,i)=cmplx(0.3d0*point,0.11d0*point,real64)
      basis(2,i)=cmplx(-0.17d0*point,0.05d0*(g-point),real64)
      basis(3,i)=(0d0,0d0)
    enddo
    dense_local=(0d0,0d0)
    do point=1,npoint;do row=1,r;do i=1,r
      dense_local(row,i)=dense_local(row,i)+weights(point)*conjg(basis(row,point))*basis(i,point)*potential(point)
    enddo;enddo;enddo
    call MPI_Allreduce(dense_local,dense,r*r,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call project_rt_dg_hybrid_sparse_edges(comm,r,row_ids,offsets,columns,grid_ids,weights,basis,potential,&
      actual,ok,message)
    call require(ok,'production sparse edge projection failed: '//trim(message))
    values_ok=.true.
    do i=1,nowned
      row=int(row_ids(i))
      do edge=offsets(i),offsets(i+1)-1
        values_ok=values_ok.and.abs(actual(edge)-dense(row,columns(edge)))<1d-13
      enddo
    enddo
    call require(values_ok,'sparse projection differs from dense oracle')
    if(nproc==4)call require(any_rank(nowned==0),'four-rank projection did not exercise a zero-owned rank')
    call require(any_rank(any(row_ids==3_int64)),'zero-degree operator row was not owned')
  end subroutine exercise_sparse_projection

  subroutine exercise_missing_edge_collective_failure
    integer,parameter::r=2,g=1
    integer::row,i,nowned,npoint
    integer(int64),allocatable::row_ids(:),grid_ids(:)
    integer,allocatable::offsets(:),columns(:)
    real(real64),allocatable::weights(:),potential(:)
    complex(real64),allocatable::basis(:,:),actual(:)
    logical::missing_ok
    nowned=count([(mod(row-1,nproc)==rank,row=1,r)])
    npoint=merge(1,0,rank==0)
    allocate(row_ids(nowned),offsets(nowned+1),columns(nowned),grid_ids(npoint),weights(npoint),&
      potential(npoint),basis(r,npoint),actual(nowned))
    i=0;offsets(1)=1
    do row=1,r
      if(mod(row-1,nproc)/=rank)cycle
      i=i+1;row_ids(i)=row;columns(i)=row;offsets(i+1)=i+1
    enddo
    if(npoint==1)then
      grid_ids(1)=1_int64;weights(1)=1d0;potential(1)=2d0
      basis(:,1)=[(1d0,0d0),(0.5d0,0.25d0)]
    endif
    call project_rt_dg_hybrid_sparse_edges(comm,r,row_ids,offsets,columns,grid_ids,weights,basis,potential,&
      actual,missing_ok,message)
    call require(.not.missing_ok,'projection accepted contributions missing from structural CSR')
    call require(index(message,'missing')>0,'missing-edge collective failure has no named diagnostic')
  end subroutine exercise_missing_edge_collective_failure

  integer function find_edge(offsets,columns,row_position,column) result(position)
    integer,intent(in)::offsets(:),columns(:),row_position,column
    integer::q
    position=0
    do q=offsets(row_position),offsets(row_position+1)-1;if(columns(q)==column)then;position=q;return;endif;enddo
  end function find_edge
  logical function any_rank(local_value)
    logical,intent(in)::local_value
    integer::a,b
    a=merge(1,0,local_value);call MPI_Allreduce(a,b,1,MPI_INTEGER,MPI_MAX,comm,ierr);any_rank=b==1
  end function any_rank
  subroutine require(condition,text)
    logical,intent(in)::condition;character(*),intent(in)::text
    integer::bad,global_bad
    bad=merge(0,1,condition);call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;if(rank==0)write(0,'(a)')trim(text);call MPI_Abort(comm,1,ierr);endif
  end subroutine require
end program test_rt_dg_hybrid_sparse_projection_mpi
