#include "config.h"
module rt_dg_hybrid_stationarity
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi,only:MPI_Allreduce,MPI_DOUBLE_COMPLEX,MPI_DOUBLE_PRECISION,MPI_INTEGER,&
    MPI_MAX,MPI_MIN,MPI_SUM,MPI_SUCCESS
#endif
  implicit none
  private
  type,public::s_rt_dg_hybrid_stationarity_reference
    logical::valid=.false.
    integer::certified_rank=0
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
  subroutine initialize_rt_dg_hybrid_stationarity(comm,certified_rank,row_ids,density,total_energy,&
      coefficients,s_coefficients,occupations,electron_count,hamiltonian_residual,reference,ok,message)
    integer,intent(in)::comm,certified_rank
    integer(int64),intent(in)::row_ids(:)
    real(real64),intent(in)::density(:),total_energy,occupations(:),electron_count,hamiltonian_residual
    complex(real64),intent(in)::coefficients(:,:),s_coefficients(:,:)
    type(s_rt_dg_hybrid_stationarity_reference),intent(out)::reference
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::local_s_coefficients(:,:),global_s_coefficients(:,:),weighted_coefficients(:,:)
    integer::local_bad,global_bad,ierr,i,j,local_grid_count,global_grid_count,&
      dimensions(4),minimum_dimensions(4),maximum_dimensions(4)
    reference%valid=.false.;reference%certified_rank=0;ok=.false.;message=''
    dimensions=[certified_rank,size(coefficients,2),size(s_coefficients,2),size(occupations)]
#ifdef USE_MPI
    call MPI_Allreduce(dimensions,minimum_dimensions,4,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(dimensions,maximum_dimensions,4,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(minimum_dimensions/=maximum_dimensions))then
      message='inconsistent hybrid RT stationarity dimensions';return
    endif
#else
    minimum_dimensions=dimensions;maximum_dimensions=dimensions
#endif
    local_bad=merge(0,1,size(row_ids)==size(coefficients,1).and.&
      all(shape(coefficients)==shape(s_coefficients)).and.size(occupations)==size(coefficients,2).and.&
      certified_rank>0.and.size(coefficients,2)>0.and.size(coefficients,2)<=certified_rank.and.&
      finite_real(density).and.finite_real(occupations).and.all(occupations>=0d0).and.&
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
    if(global_bad/=0)then;message='invalid hybrid RT stationarity reference';return;endif
#endif
    local_grid_count=size(density)
#ifdef USE_MPI
    call MPI_Allreduce(local_grid_count,global_grid_count,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_grid_count<=0)then
      ok=.false.;message='empty global hybrid RT stationarity density';return
    endif
#else
    global_grid_count=local_grid_count
    if(global_grid_count<=0)then;message='empty hybrid RT stationarity density';return;endif
#endif
    call validate_certified_row_catalog(comm,certified_rank,row_ids,ok,message)
    if(.not.ok)return
    allocate(local_s_coefficients(certified_rank,size(s_coefficients,2)),&
      global_s_coefficients(certified_rank,size(s_coefficients,2)))
    local_s_coefficients=(0d0,0d0)
    do i=1,size(row_ids)
      local_s_coefficients(int(row_ids(i)),:)=s_coefficients(i,:)
    enddo
#ifdef USE_MPI
    call MPI_Allreduce(local_s_coefficients,global_s_coefficients,size(local_s_coefficients),&
      MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='hybrid RT stationarity coefficient reduction failed';return;endif
#else
    global_s_coefficients=local_s_coefficients
#endif
    allocate(reference%row_ids,source=row_ids);allocate(reference%density,source=density)
    allocate(reference%occupations,source=occupations)
    allocate(weighted_coefficients,source=coefficients)
    do j=1,size(occupations);weighted_coefficients(:,j)=weighted_coefficients(:,j)*occupations(j);enddo
    allocate(reference%projector(size(coefficients,1),certified_rank))
    reference%projector=matmul(weighted_coefficients,conjg(transpose(global_s_coefficients)))
    reference%total_energy=total_energy;reference%electron_count=electron_count
    reference%certified_rank=certified_rank;reference%valid=.true.;ok=.true.;message=''
  end subroutine

  subroutine evaluate_rt_dg_hybrid_stationarity(comm,certified_rank,reference,density,total_energy,&
      coefficients,s_coefficients,electron_count,hamiltonian_residual,tolerances,receipt,ok,message)
    integer,intent(in)::comm,certified_rank
    type(s_rt_dg_hybrid_stationarity_reference),intent(in)::reference
    real(real64),intent(in)::density(:),total_energy,electron_count,hamiltonian_residual,tolerances(5)
    complex(real64),intent(in)::coefficients(:,:),s_coefficients(:,:)
    type(s_rt_dg_hybrid_stationarity_receipt),intent(out)::receipt
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::projector(:,:),local_s_coefficients(:,:),global_s_coefficients(:,:),&
      weighted_coefficients(:,:)
    real(real64)::local_values(4),global_values(4),projector_scale
    integer::local_bad,global_bad,ierr,i,j,dimensions(4),minimum_dimensions(4),maximum_dimensions(4)
    receipt%accepted=.false.;ok=.false.;message=''
    local_bad=merge(0,1,reference%valid.and.allocated(reference%row_ids).and.&
      allocated(reference%density).and.allocated(reference%occupations).and.allocated(reference%projector))
#ifdef USE_MPI
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid hybrid RT stationarity reference';return;endif
#else
    global_bad=local_bad
    if(global_bad/=0)then;message='invalid hybrid RT stationarity reference';return;endif
#endif
    dimensions=[certified_rank,reference%certified_rank,size(coefficients,2),size(s_coefficients,2)]
#ifdef USE_MPI
    call MPI_Allreduce(dimensions,minimum_dimensions,4,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(dimensions,maximum_dimensions,4,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(minimum_dimensions/=maximum_dimensions))then
      message='inconsistent hybrid RT stationarity dimensions';return
    endif
#else
    minimum_dimensions=dimensions;maximum_dimensions=dimensions
#endif
    local_bad=merge(0,1,certified_rank>0.and.reference%certified_rank==certified_rank.and.&
      size(density)==size(reference%density).and.size(coefficients,1)==size(reference%row_ids).and.&
      size(reference%projector,1)==size(reference%row_ids).and.size(reference%projector,2)==certified_rank.and.&
      size(coefficients,2)==size(reference%occupations).and.size(coefficients,2)<=certified_rank.and.&
      all(shape(coefficients)==shape(s_coefficients)).and.finite_real(density).and.&
      finite_complex(coefficients).and.finite_complex(s_coefficients).and.&
      ieee_is_finite(total_energy).and.ieee_is_finite(electron_count).and.&
      ieee_is_finite(hamiltonian_residual).and.finite_real(tolerances).and.all(tolerances>=0d0))
#ifdef USE_MPI
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid hybrid RT stationarity sample';return;endif
#else
    global_bad=local_bad
    if(global_bad/=0)then;message='invalid hybrid RT stationarity sample';return;endif
#endif
    call validate_certified_row_catalog(comm,certified_rank,reference%row_ids,ok,message)
    if(.not.ok)return
    allocate(local_s_coefficients(certified_rank,size(s_coefficients,2)),&
      global_s_coefficients(certified_rank,size(s_coefficients,2)))
    local_s_coefficients=(0d0,0d0)
    do i=1,size(reference%row_ids)
      local_s_coefficients(int(reference%row_ids(i)),:)=s_coefficients(i,:)
    enddo
#ifdef USE_MPI
    call MPI_Allreduce(local_s_coefficients,global_s_coefficients,size(local_s_coefficients),&
      MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='hybrid RT stationarity coefficient reduction failed';return;endif
#else
    global_s_coefficients=local_s_coefficients
#endif
    allocate(weighted_coefficients,source=coefficients)
    do j=1,size(reference%occupations)
      weighted_coefficients(:,j)=weighted_coefficients(:,j)*reference%occupations(j)
    enddo
    allocate(projector(size(coefficients,1),certified_rank))
    projector=matmul(weighted_coefficients,conjg(transpose(global_s_coefficients)))
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

  subroutine validate_certified_row_catalog(comm,certified_rank,row_ids,ok,message)
    integer,intent(in)::comm,certified_rank
    integer(int64),intent(in)::row_ids(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,allocatable::local_counts(:),global_counts(:)
    integer::local_rank,minimum_rank,maximum_rank,local_bad,global_bad,ierr,i
    ok=.false.;message='';local_rank=certified_rank
#ifdef USE_MPI
    call MPI_Allreduce(local_rank,minimum_rank,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_rank,maximum_rank,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_rank/=maximum_rank.or.minimum_rank<=0)then
      message='inconsistent certified hybrid RT rank';return
    endif
#else
    minimum_rank=local_rank;maximum_rank=local_rank
    if(certified_rank<=0)then;message='invalid certified hybrid RT rank';return;endif
#endif
    local_bad=merge(0,1,all(row_ids>=1_int64).and.all(row_ids<=int(certified_rank,int64)))
#ifdef USE_MPI
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid certified hybrid RT row IDs';return;endif
#else
    global_bad=local_bad
    if(global_bad/=0)then;message='invalid certified hybrid RT row IDs';return;endif
#endif
    allocate(local_counts(certified_rank),global_counts(certified_rank));local_counts=0
    do i=1,size(row_ids);local_counts(int(row_ids(i)))=local_counts(int(row_ids(i)))+1;enddo
#ifdef USE_MPI
    call MPI_Allreduce(local_counts,global_counts,certified_rank,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='certified hybrid RT row ownership reduction failed';return;endif
#else
    global_counts=local_counts
#endif
    if(any(global_counts/=1))then;message='certified hybrid RT rows are not owned exactly once';return;endif
    ok=.true.;message=''
  end subroutine validate_certified_row_catalog

  logical function finite_real(values) result(ok)
    real(real64),intent(in)::values(:);ok=all(ieee_is_finite(values))
  end function
  logical function finite_complex(values) result(ok)
    complex(real64),intent(in)::values(:,:)
    ok=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function
end module rt_dg_hybrid_stationarity
