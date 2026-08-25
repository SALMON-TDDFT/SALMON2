#include "config.h"
program test_dg_hybrid_lcfo_support_redistribution_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
  use dg_hybrid_lcfo_support_redistribution,only:redistribute_dg_hybrid_lcfo_support_tile
  implicit none
  integer,parameter::nglobal_point=6,nglobal_basis=4
  integer::comm,rank,nproc,ierr,nlocal,p,f,nb_local,j
  integer(int64),allocatable::spatial_ids(:)
  type(s_dg_hybrid_fragment_basis),allocatable::bases(:)
  complex(real64),allocatable::tile(:,:)
  complex(real64)::reference(nglobal_basis,nglobal_point)
  integer(int64)::peak_elements,fingerprint
  logical::ok,tile_matches
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  reference=(0d0,0d0)
  reference(1,1:4)=[(1d0,0d0),(0.5d0,0d0),(0.2d0,0d0),(0.1d0,0d0)]
  reference(2,2:4)=[(0.3d0,0d0),(0.7d0,0d0),(0.4d0,0d0)]
  reference(3,3:6)=[(0.6d0,0d0),(0.8d0,0d0),(0.2d0,0d0),(0.1d0,0d0)]
  reference(4,3:5)=[(0.2d0,0d0),(0.9d0,0d0),(0.5d0,0d0)]
  nb_local=count([(mod(f-1,nproc)==rank,f=1,2)]);allocate(bases(nb_local));j=0
  do f=1,2
    if(mod(f-1,nproc)/=rank)cycle
    j=j+1;bases(j)%fragment_id=f;bases(j)%generation=1
    allocate(bases(j)%global_ids(2),bases(j)%sector(2),bases(j)%buffer_point_ids(4))
    bases(j)%global_ids=[int(2*f-1,int64),int(2*f,int64)];bases(j)%sector=[1,2]
    if(f==1)then
      bases(j)%buffer_point_ids=[1_int64,2_int64,3_int64,4_int64]
    else
      bases(j)%buffer_point_ids=[3_int64,4_int64,5_int64,6_int64]
    endif
    allocate(bases(j)%buffer_values(4,2))
    do p=1,4
      bases(j)%buffer_values(p,1)=reference(2*f-1,int(bases(j)%buffer_point_ids(p)))
      bases(j)%buffer_values(p,2)=reference(2*f,int(bases(j)%buffer_point_ids(p)))
    enddo
    bases(j)%provenance_fingerprint=1000_int64+f
  enddo
  nlocal=count([(mod(p-1,nproc)==rank,p=1,nglobal_point)]);allocate(spatial_ids(nlocal));j=0
  do p=1,nglobal_point;if(mod(p-1,nproc)==rank)then;j=j+1;spatial_ids(j)=p;endif;enddo
  do f=1,nglobal_basis,2
    call redistribute_dg_hybrid_lcfo_support_tile(comm,bases,nglobal_point,nglobal_basis,spatial_ids,f,&
      min(2,nglobal_basis-f+1),tile,peak_elements,fingerprint,ok,message)
    call require(ok,trim(message));call require(all(shape(tile)==[min(2,nglobal_basis-f+1),nlocal]),&
      'LCFO support tile shape mismatch')
    tile_matches=.true.
    do j=1,nlocal
      tile_matches=tile_matches.and.&
        maxval(abs(tile(:,j)-reference(f:f+size(tile,1)-1,int(spatial_ids(j)))))<1d-14
    enddo
    call require(tile_matches,'LCFO support tile value mismatch')
    call require(peak_elements<=int(size(tile)+2*size(tile,1),int64),&
      'LCFO support redistribution retained a global-point tile')
    call require(fingerprint/=0_int64,'LCFO support tile fingerprint is empty')
  enddo
  if(rank==0)write(*,'(a,i0,a,i0)')'HYBRID_LCFO_SUPPORT ranks=',nproc,' fingerprint=',fingerprint
  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid LCFO support redistribution on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,text)
    logical,intent(in)::condition;character(*),intent(in)::text;logical::global_condition
    call MPI_Allreduce(condition,global_condition,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(.not.global_condition)then;if(rank==0)write(0,'(a)')trim(text);call MPI_Abort(comm,1,ierr);endif
  end subroutine require
end program test_dg_hybrid_lcfo_support_redistribution_mpi
