#include "config.h"
program test_rt_dg_hybrid_distributed_v4_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use rt_dg_hybrid_sparse_exchange,only:s_rt_dg_sparse_exchange,build_rt_dg_sparse_exchange,&
    exchange_rt_dg_sparse_matrix
  use rt_dg_hybrid_point_density,only:reconstruct_rt_dg_point_csr_density
  use rt_dg_hybrid_sparse_projection,only:project_rt_dg_hybrid_point_csr_edges
  implicit none
  type(s_rt_dg_sparse_exchange)::plan
  integer::comm,rank,nproc,ierr,n,nlocal,nrhs,i,j,k,collective_count
  integer(int64)::fingerprint,workspace_peak
  integer(int64),allocatable::owned(:)
  integer,allocatable::needed(:)
  complex(real64),allocatable::local_values(:,:),needed_values(:,:)
  logical::ok
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  nlocal=50;n=nlocal*nproc;nrhs=7
  allocate(owned(nlocal),needed(4*nlocal),local_values(nlocal,nrhs),needed_values(4*nlocal,nrhs))
  do i=1,nlocal
    owned(i)=int(rank*nlocal+i,int64)
    do j=1,nrhs;local_values(i,j)=canonical_value(int(owned(i)),j);enddo
    needed(4*i-3)=int(owned(i))
    needed(4*i-2)=mod(int(owned(i)),n)+1
    needed(4*i-1)=mod(int(owned(i))+nlocal-1,n)+1
    needed(4*i)=needed(4*i-2)
  enddo
  fingerprint=7717_int64
  call build_rt_dg_sparse_exchange(comm,n,fingerprint,owned,needed,plan,ok,message)
  call require(ok,'distributed coefficient halo setup failed: '//trim(message))
  call exchange_rt_dg_sparse_matrix(comm,plan,local_values,needed_values,workspace_peak,collective_count,ok,message)
  call require(ok,'distributed coefficient matrix exchange failed: '//trim(message))
  call require(collective_count==1,'coefficient halo did not use exactly one payload collective')
  call require(workspace_peak<=16_int64*int(nrhs,int64)*&
    int(size(plan%send_positions)+size(plan%receive_values),int64),&
    'coefficient halo workspace exceeds local sparse payload')
  do i=1,size(needed);do j=1,nrhs
    call require(abs(needed_values(i,j)-canonical_value(needed(i),j))<1d-14,&
      'packed distributed coefficient halo disagrees with dense oracle')
  enddo;enddo
  call exercise_point_density
  call exercise_point_projection
  call exercise_zero_owned
  if(rank==0)write(*,'(a,i0,a,i0,a,i0)')'PASS distributed-v4 coefficient halo ranks=',nproc,&
    ' local_basis=',nlocal,' repeated_edges=',size(needed)
  call MPI_Finalize(ierr)
contains
  pure complex(real64) function canonical_value(row,column)
    integer,intent(in)::row,column
    canonical_value=cmplx(0.125d0*row+0.03125d0*column,-0.0625d0*row+0.015625d0*column,real64)
  end function canonical_value
  subroutine exercise_point_density
    type(s_rt_dg_sparse_exchange)::density_plan
    integer,parameter::npoint=4000,degree=3
    integer::nhalo,remote_rank,point,slot,status,payload_count
    integer,allocatable::halo_ids(:),point_offsets(:),support_slots(:)
    complex(real64),allocatable::support_values(:)
    real(real64),allocatable::occupations(:),density(:),oracle(:)
    integer(int64)::density_workspace
    logical::density_ok
    character(256)::density_message
    remote_rank=mod(rank+1,nproc);nhalo=merge(2*nlocal,nlocal,nproc>1)
    allocate(halo_ids(nhalo),point_offsets(npoint+1),support_slots(degree*npoint),&
      support_values(degree*npoint),occupations(nrhs),density(npoint),oracle(npoint))
    do i=1,nlocal;halo_ids(i)=int(rank*nlocal+i);enddo
    if(nproc>1)then;do i=1,nlocal;halo_ids(nlocal+i)=remote_rank*nlocal+i;enddo;endif
    do j=1,nrhs;occupations(j)=0.25d0+0.125d0*j;enddo
    point_offsets(1)=1
    do point=1,npoint
      support_slots(degree*(point-1)+1)=mod(point-1,nlocal)+1
      support_slots(degree*(point-1)+2)=mod(3*point-1,nhalo)+1
      support_slots(degree*(point-1)+3)=mod(7*point-1,nhalo)+1
      support_values(degree*(point-1)+1)=cmplx(0.5d0,0.125d0,real64)
      support_values(degree*(point-1)+2)=cmplx(-0.25d0,0.375d0,real64)
      support_values(degree*(point-1)+3)=cmplx(0.0625d0,-0.125d0,real64)
      point_offsets(point+1)=degree*point+1
      oracle(point)=0d0
      do j=1,nrhs
        oracle(point)=oracle(point)+occupations(j)*abs(sum([&
          (support_values(degree*(point-1)+slot)*canonical_value(&
            halo_ids(support_slots(degree*(point-1)+slot)),j),slot=1,degree)]))**2
      enddo
    enddo
    call build_rt_dg_sparse_exchange(comm,n,fingerprint,owned,halo_ids,density_plan,density_ok,density_message)
    call require(density_ok,'point-CSR coefficient halo setup failed: '//trim(density_message))
    call reconstruct_rt_dg_point_csr_density(comm,density_plan,point_offsets,support_slots,support_values,&
      local_values,occupations,density,density_workspace,payload_count,density_ok,density_message)
    call require(density_ok,'point-CSR local density failed: '//trim(density_message))
    call require(maxval(abs(density-oracle))<1d-12,'point-CSR density disagrees with dense oracle')
    call require(payload_count==1,'point-CSR density payload collective count depends on rank count')
    call require(density_workspace<=16_int64*int(nrhs,int64)*int(3*nhalo+1,int64),&
      'point-CSR density workspace scales with repeated grid support')
  end subroutine exercise_point_density
  subroutine exercise_point_projection
    integer::small_n,row,previous,next_row,q,e,p
    integer(int64),allocatable::projection_rows(:),projection_grid_ids(:)
    integer,allocatable::projection_offsets(:),projection_columns(:),point_offsets(:),support_ids(:)
    complex(real64),allocatable::support_values(:),projected(:),oracle_values(:)
    real(real64),allocatable::weights(:),potential(:)
    logical::projection_ok
    character(256)::projection_message
    small_n=3*nproc
    allocate(projection_rows(3),projection_offsets(4),projection_columns(9),projection_grid_ids(3),&
      point_offsets(4),support_ids(6),support_values(6),weights(3),potential(3),projected(9),oracle_values(9))
    projection_offsets=[1,4,7,10];point_offsets=[1,3,5,7];oracle_values=(0d0,0d0)
    do i=1,3
      row=3*rank+i;projection_rows(i)=int(row,int64);projection_grid_ids(i)=int(row,int64)
      previous=mod(row-2+small_n,small_n)+1;next_row=mod(row,small_n)+1
      projection_columns(3*(i-1)+1:3*i)=[previous,row,next_row]
      call sort_three(projection_columns(3*(i-1)+1:3*i))
      support_ids(2*(i-1)+1:2*i)=[row,next_row]
      support_values(2*(i-1)+1)=cmplx(0.5d0+0.01d0*row,0.125d0,real64)
      support_values(2*i)=cmplx(-0.25d0,0.0625d0*row,real64)
      weights(i)=0.5d0+0.01d0*row;potential(i)=-0.75d0+0.02d0*row
    enddo
    call project_rt_dg_hybrid_point_csr_edges(comm,small_n,projection_rows,projection_offsets,&
      projection_columns,projection_grid_ids,weights,point_offsets,support_ids,support_values,potential,&
      projected,projection_ok,projection_message)
    call require(projection_ok,'point-CSR sparse projection failed: '//trim(projection_message))
    do i=1,3
      row=int(projection_rows(i))
      do e=projection_offsets(i),projection_offsets(i+1)-1
        do q=1,small_n
          next_row=mod(q,small_n)+1
          if((row==q.or.row==next_row).and.&
             (projection_columns(e)==q.or.projection_columns(e)==next_row))then
            oracle_values(e)=oracle_values(e)+(0.5d0+0.01d0*q)*(-0.75d0+0.02d0*q)*&
              conjg(projection_point_value(q,row,small_n))*&
              projection_point_value(q,projection_columns(e),small_n)
          endif
        enddo
      enddo
    enddo
    call require(maxval(abs(projected-oracle_values))<1d-12,&
      'production point-CSR projection disagrees with dense oracle')
  end subroutine exercise_point_projection
  subroutine exercise_zero_owned
    type(s_rt_dg_sparse_exchange)::zero_plan
    integer::zero_n,zero_local,zero_payload_count
    integer(int64),allocatable::zero_rows(:)
    integer,allocatable::zero_needed(:)
    complex(real64),allocatable::zero_values(:,:),zero_received(:,:)
    integer(int64)::zero_workspace
    logical::zero_ok,received_ok
    character(256)::zero_message
    if(nproc<2)return
    zero_n=50*(nproc-1);zero_local=merge(0,50,rank==nproc-1)
    allocate(zero_rows(zero_local),zero_values(zero_local,nrhs))
    do i=1,zero_local
      zero_rows(i)=int(rank*50+i,int64)
      do j=1,nrhs;zero_values(i,j)=canonical_value(int(zero_rows(i)),j);enddo
    enddo
    if(rank==nproc-1)then
      allocate(zero_needed(1),zero_received(1,nrhs));zero_needed=[1]
    else
      allocate(zero_needed(0),zero_received(0,nrhs))
    endif
    call build_rt_dg_sparse_exchange(comm,zero_n,8181_int64,zero_rows,zero_needed,zero_plan,zero_ok,zero_message)
    call require(zero_ok,'zero-owned halo setup failed: '//trim(zero_message))
    call exchange_rt_dg_sparse_matrix(comm,zero_plan,zero_values,zero_received,zero_workspace,&
      zero_payload_count,zero_ok,zero_message)
    call require(zero_ok,'zero-owned halo exchange failed: '//trim(zero_message))
    received_ok=.true.
    if(rank==nproc-1)then
      do j=1,nrhs;received_ok=received_ok.and.&
        abs(zero_received(1,j)-canonical_value(1,j))<1d-14;enddo
    endif
    call require(received_ok,'zero-owned rank did not receive requested coefficient')
  end subroutine exercise_zero_owned
  pure complex(real64) function projection_point_value(point,row_id,nsize)
    integer,intent(in)::point,row_id,nsize
    if(row_id==point)then
      projection_point_value=cmplx(0.5d0+0.01d0*point,0.125d0,real64)
    else if(row_id==mod(point,nsize)+1)then
      projection_point_value=cmplx(-0.25d0,0.0625d0*point,real64)
    else
      projection_point_value=(0d0,0d0)
    endif
  end function projection_point_value
  subroutine sort_three(values)
    integer,intent(inout)::values(3);integer::a,b,t
    do a=1,2;do b=a+1,3;if(values(b)<values(a))then;t=values(a);values(a)=values(b);values(b)=t;endif;enddo;enddo
  end subroutine sort_three
  subroutine require(condition,text)
    logical,intent(in)::condition;character(*),intent(in)::text
    integer::bad,total_bad,status
    bad=merge(0,1,condition);call MPI_Allreduce(bad,total_bad,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.total_bad/=0)then
      if(.not.condition)write(0,'(a)')trim(text)
      error stop trim(text)
    endif
  end subroutine require
end program test_rt_dg_hybrid_distributed_v4_mpi
