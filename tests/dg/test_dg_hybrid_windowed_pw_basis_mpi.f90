#include "config.h"
program test_dg_hybrid_windowed_pw_basis_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_windowed_pw_types,only:s_dg_hybrid_basis_catalog
  use dg_hybrid_windowed_pw_basis,only:build_dg_hybrid_windowed_pw_basis,&
    materialize_dg_hybrid_windowed_pw_columns
  implicit none
  integer,parameter::nrow=8,nfragment=4,ng=3,noperation=4
  integer::comm,rank,nproc,ierr,i,f,g,op,nlocal,position,column,first,width
  integer(int64),allocatable::row_ids(:)
  integer::fragment_action(nfragment,noperation),row_action(nrow,noperation)
  integer::g_action(ng,noperation),g_star(ng),g_conjugate(ng)
  real(real64),allocatable::coordinates(:,:),raw_windows(:,:),windows(:,:)
  real(real64)::g_vectors(3,ng),reciprocal_rotation(3,3,noperation),pi,partition_defect
  real(real64),allocatable::unsafe_coordinates(:,:)
  complex(real64),allocatable::tile(:,:)
  type(s_dg_hybrid_basis_catalog)::catalog
  integer(int64)::workspace,fingerprint,reference_fingerprint
  logical::ok
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  pi=acos(-1d0);nlocal=count([(mod(i-1,nproc)==rank,i=1,nrow)])
  allocate(row_ids(nlocal),coordinates(3,nlocal),raw_windows(nfragment,nlocal))
  position=0
  do i=1,nrow
    if(mod(i-1,nproc)/=rank)cycle
    position=position+1;row_ids(position)=i
    coordinates(:,position)=[2d0*pi*real(i-1,real64)/real(nrow,real64),0d0,0d0]
    raw_windows(:,position)=0d0;raw_windows(1+mod(i-1,nfragment),position)=1d0
  enddo
  do op=1,noperation
    do f=1,nfragment;fragment_action(f,op)=1+ieor(f-1,op-1);enddo
    do i=1,nrow;row_action(i,op)=1+ieor(mod(i-1,nfragment),op-1)+nfragment*((i-1)/nfragment);enddo
    reciprocal_rotation(:,:,op)=0d0
    if(op==2.or.op==4)then
      do i=1,3;reciprocal_rotation(i,i,op)=-1d0;enddo
      g_action(:,op)=[3,2,1]
    else
      do i=1,3;reciprocal_rotation(i,i,op)=1d0;enddo
      g_action(:,op)=[1,2,3]
    endif
  enddo
  g_vectors=0d0;g_vectors(1,:)=[-1d0,0d0,1d0]
  g_star=[1,2,1];g_conjugate=[3,2,1]
  call build_dg_hybrid_windowed_pw_basis(comm,nrow,row_ids,coordinates,raw_windows,&
    fragment_action,row_action,g_vectors,reciprocal_rotation,g_action,g_star,g_conjugate,2,1d-12,&
    windows,catalog,workspace,fingerprint,ok,message)
  call require(ok,trim(message));reference_fingerprint=fingerprint
  partition_defect=maxval(abs(sum(windows**2,dim=1)-1d0))
  call require(partition_defect<1d-14,'partition windows are not normalized')
  call require(size(catalog%packets)==nfragment*2,'incorrect fragment/star packet count')
  call require(workspace>0_int64,'basis workspace receipt is zero')
  column=0
  do first=1,nfragment*ng,2
    width=min(2,nfragment*ng-first+1)
    call materialize_dg_hybrid_windowed_pw_columns(catalog,g_vectors,coordinates,windows,&
      first,width,tile,ok,message)
    call require(ok,trim(message));call require(size(tile,1)==width,'bounded tile width is incorrect')
    call require(all(abs(tile)<=1d0+1d-14),'windowed PW magnitude exceeds its window')
    column=column+width
  enddo
  call require(column==nfragment*ng,'bounded materialization omitted a PW column')
  allocate(unsafe_coordinates,source=coordinates)
  if(size(unsafe_coordinates,2)>0)unsafe_coordinates(1,1)=huge(1d0)/2d0
  call materialize_dg_hybrid_windowed_pw_columns(catalog,g_vectors,unsafe_coordinates,windows,&
    1,1,tile,ok,message)
  call require(.not.ok,'unsafe finite PW phase magnitude was accepted')
  deallocate(unsafe_coordinates)

  if(nlocal>0)raw_windows(:,1)=1d300*raw_windows(:,1)
  call build_dg_hybrid_windowed_pw_basis(comm,nrow,row_ids,coordinates,raw_windows,&
    fragment_action,row_action,g_vectors,reciprocal_rotation,g_action,g_star,g_conjugate,2,1d-12,&
    windows,catalog,workspace,fingerprint,ok,message)
  call require(ok,'scaled finite windows must normalize without overflow')
  call require(fingerprint==reference_fingerprint,'window normalization depends on finite input scale')
  if(nlocal>0)raw_windows(:,1)=raw_windows(:,1)/1d300

  g_star=[1,2,3]
  call build_dg_hybrid_windowed_pw_basis(comm,nrow,row_ids,coordinates,raw_windows,&
    fragment_action,row_action,g_vectors,reciprocal_rotation,g_action,g_star,g_conjugate,2,1d-12,&
    windows,catalog,workspace,fingerprint,ok,message)
  call require(.not.ok,'incomplete conjugate G star was accepted')
  g_star=[1,2,1]

  if(nlocal>0.and.rank==0)raw_windows(1+mod(int(row_ids(1)),nfragment),1)=0.5d0
  call build_dg_hybrid_windowed_pw_basis(comm,nrow,row_ids,coordinates,raw_windows,&
    fragment_action,row_action,g_vectors,reciprocal_rotation,g_action,g_star,g_conjugate,2,1d-12,&
    windows,catalog,workspace,fingerprint,ok,message)
  call require(.not.ok,'noncovariant partition windows were accepted')
  if(nlocal>0.and.rank==0)then
    raw_windows(:,1)=0d0;raw_windows(1+mod(int(row_ids(1)-1_int64),nfragment),1)=1d0
  endif

  if(nproc>1.and.rank==0)fragment_action(1,1)=2
  if(nproc==1)fragment_action(1,1)=0
  call build_dg_hybrid_windowed_pw_basis(comm,nrow,row_ids,coordinates,raw_windows,&
    fragment_action,row_action,g_vectors,reciprocal_rotation,g_action,g_star,g_conjugate,2,1d-12,&
    windows,catalog,workspace,fingerprint,ok,message)
  call require(.not.ok,'rank-disagreeing fragment action was accepted')

  if(rank==0)then
    write(*,'(a,i0,a,i0)')'HYBRID_BASIS ranks=',nproc,' fingerprint=',reference_fingerprint
    write(*,'(a,i0,a)')'PASS hybrid windowed PW basis on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label
    integer::local_failure,global_failure
    local_failure=merge(0,1,condition)
    call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_failure/=0)error stop label
  end subroutine require
end program test_dg_hybrid_windowed_pw_basis_mpi
