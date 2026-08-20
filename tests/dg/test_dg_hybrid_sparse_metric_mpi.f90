#include "config.h"
program test_dg_hybrid_sparse_metric_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_sparse_metric,only:s_dg_hybrid_sparse_metric,build_dg_hybrid_sparse_metric,&
    apply_dg_hybrid_sparse_metric
  implicit none
  integer,parameter::n=4
  integer::comm,rank,nproc,ierr,nowned,i,j,pos,nnz
  integer(int64),allocatable::row_ids(:)
  integer,allocatable::offsets(:),columns(:),packets(:),active_packets(:),rejected_packets(:)
  complex(real64),allocatable::values(:),x(:),y(:)
  complex(real64)::dense(n,n),reference(n)
  type(s_dg_hybrid_sparse_metric)::metric
  integer(int64)::fingerprint,reference_fingerprint,workspace
  real(real64)::condition,apply_defect,tolerance
  logical::ok
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  dense=reshape([(1d0,0d0),(0.2d0,-0.1d0),(0.05d0,0d0),(0d0,0d0),&
    (0.2d0,0.1d0),(1.1d0,0d0),(0d0,0d0),(0.04d0,0.01d0),&
    (0.05d0,0d0),(0d0,0d0),(0.9d0,0d0),(0.1d0,-0.02d0),&
    (0d0,0d0),(0.04d0,-0.01d0),(0.1d0,0.02d0),(1.05d0,0d0)],[n,n])
  packets=[1,1,2,2];tolerance=1d-10;call distribute_dense(dense,row_ids,offsets,columns,values)
  call build_dg_hybrid_sparse_metric(comm,n,row_ids,offsets,columns,values,packets,8,&
    101_int64,202_int64,303_int64,tolerance,metric,active_packets,rejected_packets,&
    condition,workspace,fingerprint,ok,message)
  call require(ok,trim(message));reference_fingerprint=fingerprint
  call require(all(active_packets==[1,2]).and.size(rejected_packets)==0,'valid packets were rejected')
  call require(metric%numerical_rank==n.and.condition>1d0,'metric rank/condition receipt is invalid')
  allocate(x(n));x=[(0.3d0,0.2d0),(-0.4d0,0.1d0),(0.7d0,-0.3d0),(0.2d0,0.5d0)]
  call apply_dg_hybrid_sparse_metric(metric,x,y,ok,message);call require(ok,trim(message))
  reference=matmul(dense,x)
  apply_defect=0d0
  do i=1,size(row_ids)
    apply_defect=max(apply_defect,abs(y(i)-reference(int(row_ids(i)))))
  enddo
  call require(apply_defect<2d-13,'sparse metric apply differs from dense oracle')

  dense(3,3)=1d0;dense(4,4)=1d0;dense(3,4)=1d0;dense(4,3)=1d0
  call distribute_dense(dense,row_ids,offsets,columns,values)
  call build_dg_hybrid_sparse_metric(comm,n,row_ids,offsets,columns,values,packets,8,&
    101_int64,202_int64,303_int64,tolerance,metric,active_packets,rejected_packets,&
    condition,workspace,fingerprint,ok,message)
  call require(ok,'rank-deficient packet selection failed collectively')
  call require(all(active_packets==[1]).and.all(rejected_packets==[2]),'rank-deficient packet was split or retained')
  call require(metric%numerical_rank==2,'packet-atomic numerical rank is incorrect')

  dense(3,3)=1d0;dense(4,4)=1d-13;dense(3,4)=0d0;dense(4,3)=0d0
  call distribute_dense(dense,row_ids,offsets,columns,values)
  call build_dg_hybrid_sparse_metric(comm,n,row_ids,offsets,columns,values,packets,8,&
    101_int64,202_int64,303_int64,tolerance,metric,active_packets,rejected_packets,&
    condition,workspace,fingerprint,ok,message)
  call require(ok.and.all(rejected_packets==[2]),'near-threshold packet was partially retained')

  dense(4,4)=-0.1d0
  call distribute_dense(dense,row_ids,offsets,columns,values)
  call build_dg_hybrid_sparse_metric(comm,n,row_ids,offsets,columns,values,packets,8,&
    101_int64,202_int64,303_int64,tolerance,metric,active_packets,rejected_packets,&
    condition,workspace,fingerprint,ok,message)
  call require(.not.ok,'indefinite metric packet was accepted')

  dense(3,3)=1d0;dense(4,4)=1d0;dense(3,4)=(0.2d0,0d0);dense(4,3)=(0d0,0d0)
  call distribute_dense(dense,row_ids,offsets,columns,values)
  call build_dg_hybrid_sparse_metric(comm,n,row_ids,offsets,columns,values,packets,8,&
    101_int64,202_int64,303_int64,tolerance,metric,active_packets,rejected_packets,&
    condition,workspace,fingerprint,ok,message)
  call require(.not.ok,'missing reverse metric edge was accepted')

  dense(4,3)=conjg(dense(3,4));call distribute_dense(dense,row_ids,offsets,columns,values)
  if(nproc>1)then;if(rank==0)tolerance=2d-10;else;tolerance=-1d0;endif
  call build_dg_hybrid_sparse_metric(comm,n,row_ids,offsets,columns,values,packets,8,&
    101_int64,202_int64,303_int64,tolerance,metric,active_packets,rejected_packets,&
    condition,workspace,fingerprint,ok,message)
  call require(.not.ok,'rank-disagreeing metric tolerance was accepted')

  if(rank==0)then
    write(*,'(a,i0,a,i0)')'HYBRID_METRIC ranks=',nproc,' fingerprint=',reference_fingerprint
    write(*,'(a,i0,a)')'PASS hybrid sparse metric on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine distribute_dense(matrix,ids,row_offsets,column_ids,csr_values)
    complex(real64),intent(in)::matrix(:,:)
    integer(int64),allocatable,intent(out)::ids(:)
    integer,allocatable,intent(out)::row_offsets(:),column_ids(:)
    complex(real64),allocatable,intent(out)::csr_values(:)
    integer::row,column,k
    nowned=count([(mod(row-1,nproc)==rank,row=1,n)])
    allocate(ids(nowned),row_offsets(nowned+1));nnz=0;pos=0;row_offsets(1)=1
    do row=n,1,-1
      if(mod(row-1,nproc)/=rank)cycle
      pos=pos+1;ids(pos)=row
      do column=1,n;if(abs(matrix(row,column))>0d0)nnz=nnz+1;enddo
      row_offsets(pos+1)=nnz+1
    enddo
    allocate(column_ids(nnz),csr_values(nnz));k=0
    do pos=1,nowned
      row=int(ids(pos))
      do column=1,n
        if(abs(matrix(row,column))==0d0)cycle
        k=k+1;column_ids(k)=column;csr_values(k)=matrix(row,column)
      enddo
    enddo
  end subroutine distribute_dense
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label
    integer::lf,gf
    lf=merge(0,1,condition);call MPI_Allreduce(lf,gf,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gf/=0)error stop label
  end subroutine require
end program test_dg_hybrid_sparse_metric_mpi
