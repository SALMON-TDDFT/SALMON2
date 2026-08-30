#include "config.h"
module rt_dg_hybrid_stationarity
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi,only:MPI_Allreduce,MPI_DOUBLE_COMPLEX,MPI_DOUBLE_PRECISION,MPI_INTEGER,&
    MPI_MAX,MPI_SUM,MPI_SUCCESS
#endif
  implicit none
  private
  type,public::s_rt_dg_hybrid_stationarity_reference
    logical::valid=.false.
    real(real64)::total_energy=0d0,electron_count=0d0
    integer(int64),allocatable::row_ids(:)
    real(real64),allocatable::density(:),occupations(:)
    complex(real64),allocatable::projector(:,:)
  end type
  type,public::s_rt_dg_hybrid_stationarity_receipt
    logical::accepted=.false.
    real(real64)::density_drift=huge(1d0),energy_drift=huge(1d0),&
      projector_drift=huge(1d0),electron_drift=huge(1d0),&
      hamiltonian_residual=huge(1d0)
  end type
  public::initialize_rt_dg_hybrid_stationarity,evaluate_rt_dg_hybrid_stationarity
contains
  subroutine initialize_rt_dg_hybrid_stationarity(comm,row_ids,density,total_energy,&
      coefficients,s_coefficients,occupations,electron_count,hamiltonian_residual,reference,ok,message)
    integer,intent(in)::comm
    integer(int64),intent(in)::row_ids(:)
    real(real64),intent(in)::density(:),total_energy,occupations(:),electron_count,hamiltonian_residual
    complex(real64),intent(in)::coefficients(:,:),s_coefficients(:,:)
    type(s_rt_dg_hybrid_stationarity_reference),intent(out)::reference
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::local_bad,global_bad,ierr
    local_bad=merge(0,1,size(row_ids)==size(coefficients,1).and.&
      all(shape(coefficients)==shape(s_coefficients)).and.size(occupations)==size(coefficients,2).and.&
      size(density)>0.and.finite_real(density).and.finite_real(occupations).and.&
      finite_complex(coefficients).and.finite_complex(s_coefficients).and.&
      ieee_is_finite(total_energy).and.ieee_is_finite(electron_count).and.&
      ieee_is_finite(hamiltonian_residual))
#ifdef USE_MPI
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      ok=.false.;message='invalid hybrid RT stationarity reference';return
    endif
#else
    global_bad=local_bad
#endif
    allocate(reference%row_ids,source=row_ids);allocate(reference%density,source=density)
    allocate(reference%occupations,source=occupations)
    allocate(reference%projector(size(coefficients,1),size(coefficients,1)))
    reference%projector=matmul(coefficients,conjg(transpose(s_coefficients)))
    reference%total_energy=total_energy;reference%electron_count=electron_count
    reference%valid=.true.;ok=.true.;message=''
  end subroutine

  subroutine evaluate_rt_dg_hybrid_stationarity(comm,reference,density,total_energy,&
      coefficients,s_coefficients,electron_count,hamiltonian_residual,tolerances,receipt,ok,message)
    integer,intent(in)::comm
    type(s_rt_dg_hybrid_stationarity_reference),intent(in)::reference
    real(real64),intent(in)::density(:),total_energy,electron_count,hamiltonian_residual,tolerances(5)
    complex(real64),intent(in)::coefficients(:,:),s_coefficients(:,:)
    type(s_rt_dg_hybrid_stationarity_receipt),intent(out)::receipt
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::projector(:,:)
    real(real64)::local_values(4),global_values(4),projector_scale
    integer::local_bad,global_bad,ierr
    receipt%accepted=.false.;ok=.false.;message=''
    local_bad=merge(0,1,reference%valid.and.size(density)==size(reference%density).and.&
      size(coefficients,1)==size(reference%projector,1).and.&
      all(shape(coefficients)==shape(s_coefficients)).and.finite_real(density).and.&
      finite_complex(coefficients).and.finite_complex(s_coefficients).and.&
      ieee_is_finite(total_energy).and.ieee_is_finite(electron_count).and.&
      ieee_is_finite(hamiltonian_residual).and.finite_real(tolerances).and.all(tolerances>=0d0))
#ifdef USE_MPI
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid hybrid RT stationarity sample';return;endif
#else
    global_bad=local_bad
#endif
    allocate(projector(size(coefficients,1),size(coefficients,1)))
    projector=matmul(coefficients,conjg(transpose(s_coefficients)))
    local_values=[sum((density-reference%density)**2),sum(reference%density**2),&
      sum(abs(projector-reference%projector)**2),sum(abs(reference%projector)**2)]
#ifdef USE_MPI
    call MPI_Allreduce(local_values,global_values,4,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='hybrid RT stationarity reduction failed';return;endif
#else
    global_values=local_values
#endif
    receipt%density_drift=sqrt(global_values(1))/max(1d0,sqrt(global_values(2)))
    projector_scale=max(1d0,sqrt(global_values(4)))
    receipt%projector_drift=sqrt(global_values(3))/projector_scale
    receipt%energy_drift=abs(total_energy-reference%total_energy)/max(1d0,abs(reference%total_energy))
    receipt%electron_drift=abs(electron_count-reference%electron_count)/max(1d0,abs(reference%electron_count))
    receipt%hamiltonian_residual=abs(hamiltonian_residual)
    receipt%accepted=receipt%density_drift<=tolerances(1).and.receipt%energy_drift<=tolerances(2).and.&
      receipt%projector_drift<=tolerances(3).and.receipt%electron_drift<=tolerances(4).and.&
      receipt%hamiltonian_residual<=tolerances(5)
    ok=receipt%accepted
    if(ok)then;message='';else;message='hybrid RT zero-field stationarity gate failed';endif
  end subroutine

  logical function finite_real(values) result(ok)
    real(real64),intent(in)::values(:);ok=all(ieee_is_finite(values))
  end function
  logical function finite_complex(values) result(ok)
    complex(real64),intent(in)::values(:,:)
    ok=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function
end module rt_dg_hybrid_stationarity
