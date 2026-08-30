#include "config.h"
module dg_hybrid_total_energy
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi,only:MPI_Allreduce,MPI_IN_PLACE,MPI_DOUBLE_COMPLEX,MPI_DOUBLE_PRECISION,MPI_INTEGER,&
    MPI_MAX,MPI_SUM,MPI_SUCCESS
#endif
  implicit none
  private
  public::evaluate_dg_hybrid_fixed_energy
contains
  subroutine evaluate_dg_hybrid_fixed_energy(comm,row_ids,coefficients,occupations,&
      broken_kinetic_rows,sipg_rows,nonlocal_rows,kinetic_energy,nonlocal_energy,ok,message)
    integer,intent(in)::comm
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::coefficients(:,:),broken_kinetic_rows(:,:),sipg_rows(:,:),nonlocal_rows(:,:)
    real(real64),intent(in)::occupations(:)
    real(real64),intent(out)::kinetic_energy,nonlocal_energy
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    complex(real64),allocatable::local_coefficients(:,:),global_coefficients(:,:)
    real(real64)::local_energy(2),global_energy(2)
    integer,allocatable::ownership(:)
    integer::n,nocc,i,state,ierr,local_bad,global_bad
    n=size(broken_kinetic_rows,2);nocc=size(coefficients,2);local_bad=0
    if(n<1.or.nocc<1.or.size(row_ids)/=size(coefficients,1).or.&
      any(shape(broken_kinetic_rows)/=[size(row_ids),n]).or.&
      any(shape(sipg_rows)/=shape(broken_kinetic_rows)).or.&
      any(shape(nonlocal_rows)/=shape(broken_kinetic_rows)).or.size(occupations)/=nocc.or.&
      any(row_ids<1_int64).or.any(row_ids>int(n,int64)).or..not.finite_matrix(coefficients).or.&
      .not.finite_matrix(broken_kinetic_rows).or..not.finite_matrix(sipg_rows).or.&
      .not.finite_matrix(nonlocal_rows).or..not.all(ieee_is_finite(occupations)).or.any(occupations<0d0))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      kinetic_energy=huge(1d0);nonlocal_energy=huge(1d0);ok=.false.
      message='invalid distributed DG energy payload';return
    endif
    allocate(ownership(n));ownership=0
    do i=1,size(row_ids);ownership(int(row_ids(i)))=ownership(int(row_ids(i)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,n,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership/=1))then
      kinetic_energy=huge(1d0);nonlocal_energy=huge(1d0);ok=.false.
      message='DG energy rows are not owned exactly once';return
    endif
    allocate(local_coefficients(n,nocc),global_coefficients(n,nocc));local_coefficients=(0d0,0d0)
    do i=1,size(row_ids);local_coefficients(int(row_ids(i)),:)=coefficients(i,:);enddo
    call MPI_Allreduce(local_coefficients,global_coefficients,n*nocc,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      kinetic_energy=huge(1d0);nonlocal_energy=huge(1d0);ok=.false.
      message='DG energy coefficient redistribution failed';return
    endif
    local_energy=0d0
    do i=1,size(row_ids);do state=1,nocc
      local_energy(1)=local_energy(1)+occupations(state)*real(conjg(coefficients(i,state))*&
        sum((broken_kinetic_rows(i,:)+sipg_rows(i,:))*global_coefficients(:,state)),real64)
      local_energy(2)=local_energy(2)+occupations(state)*real(conjg(coefficients(i,state))*&
        sum(nonlocal_rows(i,:)*global_coefficients(:,state)),real64)
    enddo;enddo
    call MPI_Allreduce(local_energy,global_energy,2,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    kinetic_energy=global_energy(1);nonlocal_energy=global_energy(2)
    ok=ierr==MPI_SUCCESS.and.all(ieee_is_finite(global_energy));message=''
    if(.not.ok)message='DG fixed energy reduction failed'
#else
    kinetic_energy=huge(1d0);nonlocal_energy=huge(1d0);ok=.false.;message='DG fixed energy requires MPI'
#endif
  end subroutine
  logical function finite_matrix(values) result(ok)
    complex(real64),intent(in)::values(:,:)
    ok=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function
end module dg_hybrid_total_energy
