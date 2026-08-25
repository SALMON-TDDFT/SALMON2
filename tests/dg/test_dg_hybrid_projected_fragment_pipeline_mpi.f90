#include "config.h"
program test_dg_hybrid_projected_fragment_pipeline_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_windowed_pw_types,only:s_dg_hybrid_basis_catalog
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
  use dg_hybrid_projected_fragment_pipeline,only:build_dg_hybrid_projected_fragment_basis
  implicit none
  integer::comm,rank,nproc,ierr,fragment_id,nlocal,i,j
  integer(int64),allocatable::core_ids(:),buffer_ids(:)
  real(real64),allocatable::weights(:),core_coordinates(:,:),core_windows(:,:),buffer_coordinates(:,:),buffer_windows(:,:)
  real(real64),allocatable::g_vectors(:,:)
  complex(real64),allocatable::core_wannier(:,:),buffer_wannier(:,:),local_full(:,:),global_full(:,:)
  integer::wannier_owner(2)
  type(s_dg_hybrid_basis_catalog)::catalog
  type(s_dg_hybrid_fragment_basis)::basis
  integer(int64)::workspace,fingerprint
  logical::ok,values_ok
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  fragment_id=0;if(rank<2)fragment_id=rank+1
  nlocal=0;if(rank<2)nlocal=2
  allocate(core_ids(nlocal),weights(nlocal),core_coordinates(3,nlocal),core_windows(2,nlocal),core_wannier(2,nlocal))
  if(rank==0)core_ids=[1_int64,2_int64]
  if(rank==1)core_ids=[3_int64,4_int64]
  weights=1d0;core_coordinates=0d0;core_windows=0d0;core_wannier=(0d0,0d0)
  do i=1,nlocal
    core_coordinates(1,i)=real(core_ids(i)-1_int64,real64)
    if(core_ids(i)<=2)core_windows(1,i)=1d0
    if(core_ids(i)>=3)core_windows(2,i)=1d0
    if(core_ids(i)==1)core_wannier(1,i)=(1d0,0d0)
    if(core_ids(i)==3)core_wannier(2,i)=(1d0,0d0)
  enddo
  if(fragment_id>0)then
    allocate(buffer_ids(4),source=[1_int64,2_int64,3_int64,4_int64])
    allocate(buffer_coordinates(3,4),buffer_windows(2,4),buffer_wannier(2,4))
    buffer_coordinates=0d0;buffer_windows=0d0;buffer_wannier=(0d0,0d0)
    do i=1,4;buffer_coordinates(1,i)=real(i-1,real64);enddo
    buffer_windows(1,1:2)=1d0;buffer_windows(2,3:4)=1d0
    buffer_wannier(1,1)=(1d0,0d0);buffer_wannier(2,3)=(1d0,0d0)
  else
    allocate(buffer_ids(0),buffer_coordinates(3,0),buffer_windows(2,0),buffer_wannier(2,0))
  endif
  allocate(g_vectors(3,1));g_vectors=0d0;allocate(catalog%packets(2));catalog%valid=.true.
  do i=1,2
    catalog%packets(i)%fragment_id=i;catalog%packets(i)%star_id=1;catalog%packets(i)%owner_rank=i-1
    allocate(catalog%packets(i)%g_indices(1),source=[1])
  enddo
  catalog%packet_fingerprint=31_int64;catalog%catalog_fingerprint=37_int64;wannier_owner=[1,2]
  call build_dg_hybrid_projected_fragment_basis(comm,4,2,fragment_id,core_ids,weights,core_wannier,&
    core_coordinates,core_windows,buffer_ids,buffer_wannier,buffer_coordinates,buffer_windows,catalog,&
    g_vectors,wannier_owner,1,1d-12,41_int64,basis,workspace,fingerprint,ok,message)
  call require(ok,'projected fragment pipeline failed: '//trim(message))
  call require(size(basis%global_ids)==merge(2,0,fragment_id>0),'projected fragment owner count mismatch')
  allocate(local_full(4,4),global_full(4,4));local_full=(0d0,0d0)
  do j=1,size(basis%global_ids);do i=1,size(buffer_ids)
    local_full(int(buffer_ids(i)),int(basis%global_ids(j)))=basis%buffer_values(i,j)
  enddo;enddo
  call MPI_Allreduce(local_full,global_full,16,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  values_ok=abs(sum(conjg(global_full(:,1))*global_full(:,3)))<1d-12.and.&
    abs(sum(conjg(global_full(:,2))*global_full(:,4)))<1d-12.and.&
    abs(global_full(2,3)-1d0)<1d-12.and.abs(global_full(4,4)-1d0)<1d-12
  call require(values_ok,'projected fragment PW is not orthogonal to Wannier space')
  call require(workspace<=512_int64,'projected fragment pipeline workspace is not tile bounded')
  call require(fingerprint/=0_int64,'projected fragment pipeline fingerprint is empty')
  if(rank==0)write(*,'(a,i0,a,i0)')'PROJECTED_FRAGMENT ranks=',nproc,' fingerprint=',fingerprint
  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid projected fragment pipeline on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,text)
    logical,intent(in)::condition;character(*),intent(in)::text;logical::global_condition
    call MPI_Allreduce(condition,global_condition,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(.not.global_condition)then;if(rank==0)write(0,'(a)')trim(text);call MPI_Abort(comm,1,ierr);endif
  end subroutine require
end program test_dg_hybrid_projected_fragment_pipeline_mpi
