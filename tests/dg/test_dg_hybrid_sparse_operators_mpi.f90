#include "config.h"
program test_dg_hybrid_sparse_operators_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_sparse_operators,only:s_dg_hybrid_sparse_operators
  use dg_hybrid_full_cell_operator_adapter,only:project_dg_hybrid_full_cell_sparse_operators
  implicit none
  integer,parameter::ngrid=9,nbasis=4
  integer::comm,rank,nproc,ierr,nlocal,nowned,i,j,p,pos,k,max_materialized_width
  integer(int64),allocatable::spatial_ids(:),row_ids(:)
  integer,allocatable::offsets(:),columns(:)
  real(real64),allocatable::weights(:),coordinates(:,:)
  complex(real64),allocatable::basis(:,:),global_basis(:,:),global_hbasis(:,:)
  complex(real64),allocatable::expected_metric(:)
  complex(real64)::reference_s(nbasis,nbasis),reference_h(nbasis,nbasis),reference_z(3,nbasis,nbasis)
  complex(real64)::phase(nbasis),expected
  type(s_dg_hybrid_sparse_operators)::operators,baseline
  integer(int64)::persistent_bytes,transient_bytes,fingerprint,reference_fingerprint
  real(real64)::x,pi,defect
  logical::ok
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  max_materialized_width=0
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr);pi=acos(-1d0)
  nlocal=count([(mod(p-1,nproc)==rank,p=1,ngrid)])
  nowned=count([(mod(i-1,nproc)==rank,i=1,nbasis)])
  allocate(spatial_ids(nlocal),weights(nlocal),coordinates(3,nlocal),basis(nbasis,nlocal),&
    row_ids(nowned),offsets(nowned+1),columns(nowned*nbasis),expected_metric(nowned*nbasis),&
    global_basis(nbasis,ngrid),global_hbasis(nbasis,ngrid))
  pos=0
  do p=1,ngrid
    x=2d0*pi*real(p-1,real64)/real(ngrid,real64)
    do i=1,nbasis;global_basis(i,p)=exp(cmplx(0d0,real(i-1,real64)*x,real64))/sqrt(real(ngrid,real64));enddo
    if(mod(p-1,nproc)/=rank)cycle
    pos=pos+1;spatial_ids(pos)=p;weights(pos)=1d0;coordinates(:,pos)=[x,0.2d0*sin(x),0.1d0*cos(x)]
    basis(:,pos)=global_basis(:,p)
  enddo
  pos=0;k=0;offsets(1)=1
  do i=nbasis,1,-1
    if(mod(i-1,nproc)/=rank)cycle
    pos=pos+1;row_ids(pos)=i
    do j=1,nbasis;k=k+1;columns(k)=j;enddo
    offsets(pos+1)=k+1
  enddo
  call apply_reference(global_basis,global_hbasis)
  do i=1,nbasis;do j=1,nbasis
    reference_s(i,j)=sum(conjg(global_basis(i,:))*global_basis(j,:))
    reference_h(i,j)=sum(conjg(global_basis(i,:))*global_hbasis(j,:))
    do k=1,3
      reference_z(k,i,j)=sum(conjg(global_basis(i,:))*coordinate_component(k)*global_basis(j,:))
    enddo
  enddo;enddo
  do i=1,nowned;do k=offsets(i),offsets(i+1)-1
    expected_metric(k)=reference_s(int(row_ids(i)),columns(k))
  enddo;enddo
  call project_dg_hybrid_full_cell_sparse_operators(comm,ngrid,nbasis,spatial_ids,weights,coordinates,materialize_basis_tile,&
    row_ids,offsets,columns,expected_metric,2,apply_tile,11_int64,12_int64,13_int64,14_int64,15_int64,1d-11,&
    operators,persistent_bytes,transient_bytes,fingerprint,ok,message)
  call require(ok,trim(message));reference_fingerprint=fingerprint;baseline=operators
  defect=0d0
  do i=1,nowned
    do k=offsets(i),offsets(i+1)-1
      j=columns(k);defect=max(defect,abs(operators%metric_values(k)-reference_s(int(row_ids(i)),j)))
      defect=max(defect,abs(operators%hamiltonian_values(k)-reference_h(int(row_ids(i)),j)))
      defect=max(defect,maxval(abs(operators%position_values(:,k)-reference_z(:,int(row_ids(i)),j))))
    enddo
  enddo
  call require(defect<3d-12,'sparse S/H/Z differs from direct full-cell oracle')
  call require(persistent_bytes>0_int64.and.transient_bytes>0_int64,'operator memory receipts are zero')
  call require(max_materialized_width<=2,'operator adapter requested an unbounded basis tile')
  call require(abs(reference_s(1,3))<1d-12.and.abs(reference_z(1,1,3))>1d-3,&
    'W-P metric/position cross-block contract is not exercised')
  call require(abs(reference_h(2,2)-sum(conjg(global_basis(2,:))*&
    ((0.4d0+0.1d0*cos([(2d0*pi*real(p-1,real64)/ngrid,p=1,ngrid)]))*global_basis(2,:))))>1d-4,&
    'Hamiltonian callback did not include the nonlocal rank-one term')

  phase=[(1d0,0d0),(0d0,1d0),exp(cmplx(0d0,0.37d0,real64)),exp(cmplx(0d0,-0.22d0,real64))]
  do i=1,nbasis;basis(i,:)=phase(i)*basis(i,:);enddo
  do i=1,nowned;do k=offsets(i),offsets(i+1)-1
    expected_metric(k)=conjg(phase(int(row_ids(i))))*phase(columns(k))*reference_s(int(row_ids(i)),columns(k))
  enddo;enddo
  call project_dg_hybrid_full_cell_sparse_operators(comm,ngrid,nbasis,spatial_ids,weights,coordinates,materialize_basis_tile,&
    row_ids,offsets,columns,expected_metric,2,apply_tile,11_int64,12_int64,13_int64,14_int64,15_int64,1d-11,&
    operators,persistent_bytes,transient_bytes,fingerprint,ok,message)
  call require(ok,trim(message));defect=0d0
  do i=1,nowned;do k=offsets(i),offsets(i+1)-1
    j=columns(k);expected=conjg(phase(int(row_ids(i))))*phase(j)
    defect=max(defect,abs(operators%metric_values(k)-expected*baseline%metric_values(k)))
    defect=max(defect,abs(operators%hamiltonian_values(k)-expected*baseline%hamiltonian_values(k)))
    defect=max(defect,maxval(abs(operators%position_values(:,k)-expected*baseline%position_values(:,k))))
  enddo;enddo
  call require(defect<3d-12,'hybrid sparse operators are not gauge covariant')
  if(rank==0)then
    write(*,'(a,i0,a,i0)')'HYBRID_OPERATORS ranks=',nproc,' fingerprint=',reference_fingerprint
    write(*,'(a,i0,a)')'PASS hybrid sparse operators on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine materialize_basis_tile(first_column,column_count,tile_values,tile_ok)
    integer,intent(in)::first_column,column_count
    complex(real64),intent(out)::tile_values(:,:)
    logical,intent(out)::tile_ok
    tile_ok=first_column>=1.and.column_count>=1.and.first_column+column_count-1<=nbasis
    max_materialized_width=max(max_materialized_width,column_count)
    tile_ok=tile_ok.and.size(tile_values,1)==column_count.and.size(tile_values,2)==nlocal
    if(tile_ok)tile_values=basis(first_column:first_column+column_count-1,:)
  end subroutine materialize_basis_tile
  function coordinate_component(component) result(values)
    integer,intent(in)::component
    real(real64)::values(ngrid),angle
    integer::q
    do q=1,ngrid
      angle=2d0*pi*real(q-1,real64)/real(ngrid,real64)
      select case(component);case(1);values(q)=angle;case(2);values(q)=0.2d0*sin(angle);case default;values(q)=0.1d0*cos(angle);end select
    enddo
  end function coordinate_component
  subroutine apply_tile(tile_in,tile_out,callback_ok)
    complex(real64),intent(in)::tile_in(:,:);complex(real64),intent(out)::tile_out(:,:)
    logical,intent(out)::callback_ok
    complex(real64),allocatable::local_full(:,:),full(:,:),projection(:)
    real(real64)::angle,u(ngrid)
    integer::q,t,loc
    allocate(local_full(size(tile_in,1),ngrid),full(size(tile_in,1),ngrid),projection(size(tile_in,1)))
    local_full=(0d0,0d0)
    do loc=1,nlocal;local_full(:,int(spatial_ids(loc)))=tile_in(:,loc);enddo
    call MPI_Allreduce(local_full,full,size(full),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    do q=1,ngrid;angle=2d0*pi*real(q-1,real64)/real(ngrid,real64);u(q)=sin(angle)/sqrt(real(ngrid,real64));enddo
    do t=1,size(tile_in,1);projection(t)=sum(u*full(t,:));enddo
    do loc=1,nlocal
      q=int(spatial_ids(loc));angle=2d0*pi*real(q-1,real64)/real(ngrid,real64)
      do t=1,size(tile_in,1)
        tile_out(t,loc)=(0.4d0+0.1d0*cos(angle))*tile_in(t,loc)+0.17d0*u(q)*projection(t)
      enddo
    enddo
    callback_ok=ierr==MPI_SUCCESS
  end subroutine apply_tile
  subroutine apply_reference(input,output)
    complex(real64),intent(in)::input(:,:);complex(real64),intent(out)::output(:,:)
    real(real64)::angle,u(ngrid);complex(real64)::projection
    integer::q,t
    do q=1,ngrid;angle=2d0*pi*real(q-1,real64)/real(ngrid,real64);u(q)=sin(angle)/sqrt(real(ngrid,real64));enddo
    do t=1,size(input,1)
      projection=sum(u*input(t,:))
      do q=1,ngrid;angle=2d0*pi*real(q-1,real64)/real(ngrid,real64)
        output(t,q)=(0.4d0+0.1d0*cos(angle))*input(t,q)+0.17d0*u(q)*projection
      enddo
    enddo
  end subroutine apply_reference
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label;integer::lf,gf
    lf=merge(0,1,condition);call MPI_Allreduce(lf,gf,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gf/=0)error stop label
  end subroutine require
end program test_dg_hybrid_sparse_operators_mpi
