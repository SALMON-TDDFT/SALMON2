#include "config.h"
program test_dg_hybrid_production_pw_basis_mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_windowed_pw_types,only:s_dg_hybrid_basis_catalog
  use dg_hybrid_windowed_pw_basis,only:materialize_dg_hybrid_windowed_pw_columns
  use dg_hybrid_production_pw_basis,only:build_dg_hybrid_production_pw_basis
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  integer::comm,rank,nproc,ierr,nowned,i,p
  integer,allocatable::fragment_ids(:),core_fragment_ids(:),row_action(:,:)
  integer(int64),allocatable::box_ids(:),core_ids(:)
  real(real64),allocatable::box_windows(:,:),coordinates(:,:),windows(:,:),g_vectors(:,:)
  real(real64)::reciprocal_lattice(3,3),reciprocal_rotation(3,3,2)
  complex(real64),allocatable::tile(:,:)
  type(s_dg_hybrid_basis_catalog)::catalog
  integer(int64)::fingerprint,workspace
  logical::ok,values_ok
  character(256)::message
#ifdef USE_MPI
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
#else
  comm=0;rank=0;nproc=1
#endif
  nowned=count([(mod(i-1,nproc)==rank,i=1,2)])
  allocate(fragment_ids(nowned),box_ids(4),box_windows(nowned,4),row_action(4,2))
  p=0
  do i=1,2
    if(mod(i-1,nproc)/=rank)cycle
    p=p+1;fragment_ids(p)=i
  enddo
  box_ids=[1_int64,2_int64,3_int64,4_int64]
  do p=1,nowned
    if(fragment_ids(p)==1)box_windows(p,:)=[4d0,3d0,2d0,1d0]
    if(fragment_ids(p)==2)box_windows(p,:)=[1d0,2d0,3d0,4d0]
  enddo
  row_action(:,1)=[1,2,3,4];row_action(:,2)=[4,3,2,1]
  if(rank==0)then
    allocate(core_ids(2),source=[1_int64,2_int64]);allocate(core_fragment_ids(2),source=[1,1])
  elseif(rank==1)then
    allocate(core_ids(2),source=[3_int64,4_int64]);allocate(core_fragment_ids(2),source=[2,2])
  else
    allocate(core_ids(0),core_fragment_ids(0))
  endif
  if(nproc==1)then
    deallocate(core_ids,core_fragment_ids);allocate(core_ids(4),source=box_ids)
    allocate(core_fragment_ids(4),source=[1,1,2,2])
  endif
  allocate(coordinates(3,size(core_ids)));coordinates=0d0
  do p=1,size(core_ids);coordinates(1,p)=real(core_ids(p)-1_int64,real64);enddo
  reciprocal_lattice=0d0;reciprocal_rotation=0d0
  do i=1,3
    reciprocal_lattice(i,i)=1d0;reciprocal_rotation(i,i,1)=1d0;reciprocal_rotation(i,i,2)=-1d0
  enddo
  call build_dg_hybrid_production_pw_basis(comm,4,2,fragment_ids,box_ids,box_windows,core_ids,&
    core_fragment_ids,coordinates,row_action,reciprocal_lattice,reciprocal_rotation,0d0,2,1d-12,&
    windows,g_vectors,catalog,workspace,fingerprint,ok,message)
  call require(ok,'production PW basis rejected: '//trim(message))
  call require(catalog%valid.and.size(catalog%packets)==2,'production packet catalog mismatch')
  call require(size(g_vectors,2)==1,'zero-cutoff production catalog must contain only G=0')
  allocate(tile(1,size(core_ids)))
  call materialize_dg_hybrid_windowed_pw_columns(catalog,g_vectors,coordinates,windows,1,1,tile,ok,message)
  call require(ok,'production PW packet materialization failed: '//trim(message))
  values_ok=.true.
  do p=1,size(core_ids)
    values_ok=values_ok.and.abs(tile(1,p)-cmplx(windows(1,p),0d0,real64))<1d-12
  enddo
  call require(values_ok,'materialized production PW values mismatch')
  call require(fingerprint/=0_int64.and.workspace>0_int64,'production PW receipts are missing')
  if(rank==0)write(*,'(a,i0,a,i0)')'PRODUCTION_PW ranks=',nproc,' fingerprint=',fingerprint
  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid production PW basis on ',nproc,' ranks'
#ifdef USE_MPI
  call MPI_Finalize(ierr)
#endif
contains
  subroutine require(condition,text)
    logical,intent(in)::condition
    character(*),intent(in)::text
    integer::bad,global_bad
    bad=merge(0,1,condition)
#ifdef USE_MPI
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
#else
    global_bad=bad
#endif
    if(global_bad/=0)then
      if(rank==0)write(0,'(a)')trim(text)
#ifdef USE_MPI
      call MPI_Abort(comm,1,ierr)
#else
      error stop 1
#endif
    endif
  end subroutine require
end program test_dg_hybrid_production_pw_basis_mpi
