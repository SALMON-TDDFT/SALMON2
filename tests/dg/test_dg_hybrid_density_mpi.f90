#include "config.h"
program test_dg_hybrid_density_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use dg_hybrid_density,only:reconstruct_dg_hybrid_density
  implicit none
  integer,parameter::ngrid=8,nbasis=4,nocc=2
  integer::comm,rank,nproc,ierr,nlocal,nowned,p,row,i,position
  integer(int64),allocatable::spatial_ids(:),row_ids(:)
  real(real64),allocatable::weights(:),density(:),rotated_density(:)
  real(real64)::occupations(nocc),electron_count,rotated_electron_count,defect
  complex(real64),allocatable::coefficients(:,:),rotated_coefficients(:,:)
  complex(real64)::basis(nbasis,ngrid),full_coefficients(nbasis,nocc),rotation(nocc,nocc),phase
  integer(int64)::workspace,fingerprint,rotated_fingerprint,reference_fingerprint
  logical::ok,provider_finite
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  do row=1,nbasis;do p=1,ngrid
    phase=exp(cmplx(0d0,2d0*acos(-1d0)*real((row-1)*(p-1),real64)/real(ngrid,real64),real64))
    basis(row,p)=phase/sqrt(real(ngrid,real64))
  enddo;enddo
  full_coefficients(:,1)=[(0.7d0,0.1d0),(-0.2d0,0.3d0),(0.4d0,-0.1d0),(0.1d0,0.2d0)]
  full_coefficients(:,2)=[(-0.1d0,0.2d0),(0.5d0,0.1d0),(0.2d0,0.3d0),(-0.3d0,0.2d0)]
  occupations=[1d0,1d0];rotation=reshape([cmplx(1d0,0d0,real64),cmplx(0d0,1d0,real64),&
    cmplx(0d0,1d0,real64),cmplx(1d0,0d0,real64)],[nocc,nocc])/sqrt(2d0)
  nlocal=count([(mod(p-1,nproc)==rank,p=1,ngrid)]);nowned=count([(mod(row-1,nproc)==rank,row=1,nbasis)])
  allocate(spatial_ids(nlocal),weights(nlocal),row_ids(nowned),coefficients(nowned,nocc),rotated_coefficients(nowned,nocc))
  position=0
  do p=ngrid,1,-1;if(mod(p-1,nproc)==rank)then;position=position+1;spatial_ids(position)=p;weights(position)=1d0;endif;enddo
  position=0
  do row=nbasis,1,-1
    if(mod(row-1,nproc)/=rank)cycle
    position=position+1;row_ids(position)=row;coefficients(position,:)=full_coefficients(row,:)
    rotated_coefficients(position,:)=matmul(full_coefficients(row,:),rotation)
  enddo
  provider_finite=.true.
  call reconstruct_dg_hybrid_density(comm,ngrid,spatial_ids,weights,nbasis,row_ids,coefficients,occupations,&
    materialize_basis,2,1,7117_int64,1d-12,density,electron_count,workspace,fingerprint,ok,message)
  call require(ok,trim(message));reference_fingerprint=fingerprint
  call reconstruct_dg_hybrid_density(comm,ngrid,spatial_ids,weights,nbasis,row_ids,rotated_coefficients,occupations,&
    materialize_basis,2,2,7117_int64,1d-12,rotated_density,rotated_electron_count,workspace,rotated_fingerprint,ok,message)
  call require(ok,trim(message));defect=maxval(abs(density-rotated_density))
  call MPI_Allreduce(MPI_IN_PLACE,defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  call require(defect<2d-12.and.abs(electron_count-rotated_electron_count)<2d-12,&
    'occupied unitary rotation changed reconstructed density')
  call require(fingerprint==rotated_fingerprint,'density fingerprint changed under occupied rotation')
  provider_finite=.false.
  call reconstruct_dg_hybrid_density(comm,ngrid,spatial_ids,weights,nbasis,row_ids,coefficients,occupations,&
    materialize_basis,2,1,7117_int64,1d-12,density,electron_count,workspace,fingerprint,ok,message)
  call require(.not.ok.and..not.allocated(density),'nonfinite basis callback was accepted')
  if(rank==0)then
    write(*,'(a,i0,a,i0)')'HYBRID_DENSITY ranks=',nproc,' fingerprint=',reference_fingerprint
    write(*,'(a,i0,a)')'PASS hybrid density on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine materialize_basis(first_column,column_count,tile_values,callback_ok)
    integer,intent(in)::first_column,column_count
    complex(real64),intent(out)::tile_values(:,:)
    logical,intent(out)::callback_ok
    integer::column,point
    callback_ok=.true.
    do column=1,column_count;do point=1,nlocal
      tile_values(column,point)=basis(first_column+column-1,int(spatial_ids(point)))
    enddo;enddo
    if(.not.provider_finite.and.rank==0.and.column_count>0.and.nlocal>0)&
      tile_values(1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,real64)
  end subroutine materialize_basis
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label;integer::local_bad,global_bad
    local_bad=merge(0,1,condition);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)error stop label
  end subroutine require
end program test_dg_hybrid_density_mpi
