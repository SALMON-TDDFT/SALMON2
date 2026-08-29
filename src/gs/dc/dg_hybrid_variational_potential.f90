#include "config.h"
module dg_hybrid_variational_potential
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::assemble_dg_hybrid_total_density,combine_dg_hybrid_fragment_local_potential,&
    project_dg_hybrid_fragment_local_rows
contains
  subroutine assemble_dg_hybrid_total_density(comm,global_count,core_ids,core_density,total_density,ok,message)
    integer,intent(in)::comm,global_count
    integer(int64),intent(in)::core_ids(:)
    real(real64),intent(in)::core_density(:)
    real(real64),allocatable,intent(out)::total_density(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::i,ierr,local_bad,global_bad
    integer,allocatable::ownership(:)
    real(real64),allocatable::local_density(:)
    ok=.false.;message='';local_bad=0
    if(global_count<1.or.size(core_ids)/=size(core_density).or.any(core_ids<1_int64).or.&
        any(core_ids>int(max(0,global_count),int64)).or..not.all(ieee_is_finite(core_density)).or.&
        any(core_density<0d0))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid divided density contract';return;endif
    allocate(ownership(global_count),local_density(global_count),total_density(global_count))
    ownership=0;local_density=0d0
    do i=1,size(core_ids)
      ownership(int(core_ids(i)))=ownership(int(core_ids(i)))+1
      local_density(int(core_ids(i)))=core_density(i)
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,global_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_density,total_density,global_count,&
      MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership/=1))then
      message='divided density core points are not owned exactly once';return
    endif
    ok=all(ieee_is_finite(total_density)).and.all(total_density>=0d0)
    if(ok)then;message='';else;message='assembled total density is invalid';endif
#else
    ok=.false.;message='MPI is required for divided density assembly'
#endif
  end subroutine assemble_dg_hybrid_total_density

  subroutine combine_dg_hybrid_fragment_local_potential(comm,global_count,fragment_point_ids,&
      hartree,exchange_correlation,ionic_local,semilocal_halo_complete,local_potential,ok,message)
    integer,intent(in)::comm,global_count
    integer(int64),intent(in)::fragment_point_ids(:)
    real(real64),intent(in)::hartree(:),exchange_correlation(:),ionic_local(:)
    logical,intent(in)::semilocal_halo_complete
    real(real64),intent(out)::local_potential(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::ierr,local_bad,global_bad
    local_bad=0;ok=.false.;message='';local_potential=0d0
    if(global_count<1.or.size(hartree)/=size(fragment_point_ids).or.&
        size(exchange_correlation)/=size(fragment_point_ids).or.size(ionic_local)/=size(fragment_point_ids).or.&
        size(local_potential)/=size(fragment_point_ids).or.any(fragment_point_ids<1_int64).or.&
        any(fragment_point_ids>int(max(0,global_count),int64)).or..not.semilocal_halo_complete.or.&
        .not.all(ieee_is_finite(hartree)).or..not.all(ieee_is_finite(exchange_correlation)).or.&
        .not.all(ieee_is_finite(ionic_local)))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(.not.semilocal_halo_complete)then
        message='semilocal XC fragment halo is incomplete'
      else
        message='invalid fragment local-potential contract'
      endif
      return
    endif
    local_potential=hartree+exchange_correlation+ionic_local
    ok=all(ieee_is_finite(local_potential))
    if(ok)then;message='';else;message='combined fragment local potential is invalid';endif
#else
    ok=.false.;message='MPI is required for fragment local-potential assembly'
#endif
  end subroutine combine_dg_hybrid_fragment_local_potential

  subroutine project_dg_hybrid_fragment_local_rows(comm,weights,basis_values,local_potential,local_rows,ok,message)
    integer,intent(in)::comm
    real(real64),intent(in)::weights(:),local_potential(:)
    complex(real64),intent(in)::basis_values(:,:)
    complex(real64),intent(out)::local_rows(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::i,j,p,ierr,local_bad,global_bad
    local_rows=(0d0,0d0);ok=.false.;message='';local_bad=0
    if(size(weights)/=size(local_potential).or.size(basis_values,2)/=size(weights).or.&
        any(shape(local_rows)/=[size(basis_values,1),size(basis_values,1)]).or.&
        .not.all(ieee_is_finite(weights)).or.any(weights<=0d0).or.&
        .not.all(ieee_is_finite(local_potential)).or.&
        .not.all(ieee_is_finite(real(basis_values))).or..not.all(ieee_is_finite(aimag(basis_values))))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid fragment local-potential projection contract';return
    endif
    do p=1,size(weights)
      do j=1,size(basis_values,1)
        do i=1,size(basis_values,1)
          local_rows(i,j)=local_rows(i,j)+weights(p)*conjg(basis_values(i,p))*&
            local_potential(p)*basis_values(j,p)
        enddo
      enddo
    enddo
    ok=all(ieee_is_finite(real(local_rows))).and.all(ieee_is_finite(aimag(local_rows)))
    if(ok)then;message='';else;message='projected fragment local rows are invalid';endif
#else
    ok=.false.;message='MPI is required for fragment local-potential projection'
#endif
  end subroutine project_dg_hybrid_fragment_local_rows
end module dg_hybrid_variational_potential
