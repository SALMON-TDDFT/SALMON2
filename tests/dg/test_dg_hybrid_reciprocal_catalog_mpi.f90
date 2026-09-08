#include "config.h"
program test_dg_hybrid_reciprocal_catalog_mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_reciprocal_catalog,only:build_dg_hybrid_reciprocal_catalog
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  integer::comm,rank,nproc,ierr,i
  integer,allocatable::g_integer(:,:),g_action(:,:),g_star(:),g_conjugate(:)
  real(real64)::reciprocal_lattice(3,3),reciprocal_rotation(3,3,2),&
    missing_identity_rotation(3,3,1),nonclosed_rotation(3,3,3),cutoff,effective_cutoff
  real(real64),allocatable::g_vectors(:,:)
  integer(int64)::fingerprint,reference_fingerprint
  integer::shell_added,orbit_added
  logical::ok
  character(256)::message
#ifdef USE_MPI
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
#else
  comm=0;rank=0;nproc=1
#endif
  reciprocal_lattice=0d0
  reciprocal_lattice(1,1)=1d0
  reciprocal_lattice(2,2)=1d0+5d-14
  reciprocal_lattice(3,3)=1d0+2d-14
  reciprocal_rotation=0d0
  do i=1,3
    reciprocal_rotation(i,i,1)=1d0
  enddo
  reciprocal_rotation(1,2,2)=1d0;reciprocal_rotation(2,1,2)=1d0
  reciprocal_rotation(3,3,2)=1d0
  cutoff=0.5d0
  call build_dg_hybrid_reciprocal_catalog(comm,reciprocal_lattice,reciprocal_rotation,cutoff,1d-12,&
    g_integer,g_vectors,g_action,g_star,g_conjugate,fingerprint,effective_cutoff,shell_added,orbit_added,&
    ok,message)
  call require(ok,'valid reciprocal catalog rejected: '//trim(message))
  call require(size(g_integer,2)==7,'cutoff catalog must contain G=0 and six axial modes')
  call require(effective_cutoff>cutoff.and.shell_added==2.and.orbit_added==2,&
    'boundary-shell and authoritative-orbit completion receipts are wrong')
  call require(all(g_integer(:,1)==[-1,0,0]),'catalog ordering is not deterministic lexicographic order')
  call require(all(g_integer(:,4)==[0,0,0]),'zero mode is not in the deterministic center position')
  call require(all(g_integer(:,7)==[1,0,0]),'positive final axial mode is missing')
  do i=1,size(g_integer,2)
    call require(all(g_integer(:,g_conjugate(i))==-g_integer(:,i)),'conjugate mode mismatch')
    call require(g_conjugate(g_conjugate(i))==i,'conjugate catalog is not involutive')
    call require(all(g_integer(:,g_action(i,2))==[g_integer(2,i),g_integer(1,i),g_integer(3,i)]),&
      'nonidentity reciprocal action mismatch')
    call require(g_star(g_conjugate(i))==g_star(i),'conjugates must share one G-star')
  enddo
  reference_fingerprint=fingerprint

  missing_identity_rotation=reciprocal_rotation(:,:,2:2)
  call build_dg_hybrid_reciprocal_catalog(comm,reciprocal_lattice,missing_identity_rotation,cutoff,1d-12,&
    g_integer,g_vectors,g_action,g_star,g_conjugate,fingerprint,effective_cutoff,shell_added,orbit_added,&
    ok,message)
  call require(.not.ok.and.index(message,'identity')>0,&
    'reciprocal catalog accepted an operation list without identity')

  nonclosed_rotation=0d0
  do i=1,3;nonclosed_rotation(i,i,1)=1d0;enddo
  nonclosed_rotation(1,2,2)=-1d0;nonclosed_rotation(2,1,2)=1d0;nonclosed_rotation(3,3,2)=1d0
  nonclosed_rotation(1,1,3)=-1d0;nonclosed_rotation(2,2,3)=-1d0;nonclosed_rotation(3,3,3)=1d0
  call build_dg_hybrid_reciprocal_catalog(comm,reciprocal_lattice,nonclosed_rotation,cutoff,1d-12,&
    g_integer,g_vectors,g_action,g_star,g_conjugate,fingerprint,effective_cutoff,shell_added,orbit_added,&
    ok,message)
  call require(.not.ok.and.index(message,'closed')>0,&
    'reciprocal catalog accepted a nonclosed operation list')

  cutoff=0.5d0
  if(rank==0)cutoff=-1d0
  call build_dg_hybrid_reciprocal_catalog(comm,reciprocal_lattice,reciprocal_rotation,cutoff,1d-12,&
    g_integer,g_vectors,g_action,g_star,g_conjugate,fingerprint,effective_cutoff,shell_added,orbit_added,&
    ok,message)
  call require(.not.ok,'rank-inconsistent cutoff must fail collectively')
  if(rank==0)write(*,'(a,i0,a,i0)')'RECIPROCAL_CATALOG ranks=',nproc,' fingerprint=',reference_fingerprint
  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid reciprocal catalog on ',nproc,' ranks'
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
end program test_dg_hybrid_reciprocal_catalog_mpi
