#include "config.h"
module dg_overlapping_wannier_density
  use,intrinsic::iso_fortran_env,only:int64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::reconstruct_dg_overlapping_wannier_density
contains
  subroutine reconstruct_dg_overlapping_wannier_density(comm,physical_ids,weights,values,tail_generation,&
      expected_generation,coefficients,occupations,expected_core_count,overlap,tolerance,density,&
      integrated_charge,trace_charge,ok,message)
    integer,intent(in)::comm,expected_generation
    integer(int64),intent(in)::physical_ids(:),expected_core_count
    real(8),intent(in)::weights(:),occupations(:),tolerance
    complex(8),intent(in)::values(:,:),coefficients(:,:),overlap(:,:)
    integer,intent(in)::tail_generation(:,:)
    real(8),intent(out)::density(:),integrated_charge,trace_charge
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nlocal,nwann,nstate,p,i,j,ierr,local_bad,global_bad,scalar_local(3),scalar_min(3),scalar_max(3)
    integer,allocatable::owners(:)
    complex(8),allocatable::psi(:,:),density_matrix(:,:),sp(:,:),reference_coefficients(:,:),&
      reference_overlap(:,:)
    real(8),allocatable::reference_occupations(:)
    real(8)::local_charge,trace_imag,tol_min,tol_max,payload_defect,global_payload_defect
    complex(8)::trace_value
    integer(int64)::expected_min,expected_max

    ok=.false.;message='';density=0d0;integrated_charge=huge(1d0);trace_charge=huge(1d0)
    nlocal=size(physical_ids);nwann=size(values,1);nstate=size(coefficients,2);local_bad=0
    scalar_local=[nwann,nstate,expected_generation]
    call MPI_Allreduce(scalar_local,scalar_min,3,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(scalar_local,scalar_max,3,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(expected_core_count,expected_min,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(expected_core_count,expected_max,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    call MPI_Allreduce(tolerance,tol_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    call MPI_Allreduce(tolerance,tol_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(any(scalar_min/=scalar_max).or.expected_min/=expected_max.or.tol_min/=tol_max)then
      message='inconsistent replicated overlapping-Wannier density contract';return
    endif
    if(expected_core_count<=0_int64.or.expected_core_count>int(huge(1),int64))local_bad=1
    if(nwann<1.or.nstate<1.or.size(coefficients,1)/=nwann)local_bad=1
    if(int(nwann,int64)>int(huge(1),int64)/int(max(1,nwann),int64))local_bad=1
    if(int(nwann,int64)>int(huge(1),int64)/int(max(1,nstate),int64))local_bad=1
    if(expected_generation<=0)local_bad=1
    if(size(weights)/=nlocal.or.size(values,2)/=nlocal.or.any(shape(tail_generation)/=shape(values)))local_bad=1
    if(size(density)/=nlocal.or.size(occupations)/=nstate.or.any(shape(overlap)/=[nwann,nwann]))local_bad=1
    if(any(physical_ids<1_int64).or.any(physical_ids>expected_core_count))local_bad=1
    if(any(weights<=0d0).or.any(occupations<0d0).or.tolerance<=0d0)local_bad=1
    if(.not.all(ieee_is_finite(weights)).or..not.all(ieee_is_finite(occupations)))local_bad=1
    if(.not.finite_matrix(values).or..not.finite_matrix(coefficients).or..not.finite_matrix(overlap))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid overlapping-Wannier density contract';return;endif
    allocate(reference_coefficients,source=coefficients);allocate(reference_occupations,source=occupations)
    allocate(reference_overlap,source=overlap)
    call MPI_Bcast(reference_coefficients,nwann*nstate,MPI_DOUBLE_COMPLEX,0,comm,ierr)
    call MPI_Bcast(reference_occupations,nstate,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(reference_overlap,nwann*nwann,MPI_DOUBLE_COMPLEX,0,comm,ierr)
    payload_defect=max(maxval(abs(coefficients-reference_coefficients)),&
      max(maxval(abs(occupations-reference_occupations)),maxval(abs(overlap-reference_overlap))))
    call MPI_Allreduce(payload_defect,global_payload_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(global_payload_defect>tolerance)then
      message='inconsistent replicated density coefficients, occupations, or overlap';return
    endif
    allocate(owners(int(expected_core_count)));owners=0
    do p=1,nlocal;owners(int(physical_ids(p)))=owners(int(physical_ids(p)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,owners,size(owners),MPI_INTEGER,MPI_SUM,comm,ierr)
    if(any(owners/=1))then;message='duplicate or missing unique-core density ownership';return;endif
    allocate(density_matrix(nwann,nwann))
    density_matrix=(0d0,0d0)
    do i=1,nstate
      density_matrix=density_matrix+occupations(i)*&
        matmul(reshape(coefficients(:,i),[nwann,1]),reshape(conjg(coefficients(:,i)),[1,nwann]))
    enddo
    local_bad=0
    do p=1,nlocal
      do i=1,nwann
        if(tail_generation(i,p)/=expected_generation)local_bad=1
      enddo
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='missing Wannier tail pair at unique-core density point';return;endif
    allocate(psi(nlocal,nstate));psi=matmul(transpose(values),coefficients)
    do p=1,nlocal
      density(p)=sum(occupations*abs(psi(p,:))**2)
    enddo
    if(.not.all(ieee_is_finite(density)))then;message='nonfinite reconstructed core density';return;endif
    local_charge=sum(weights*density)
    call MPI_Allreduce(local_charge,integrated_charge,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    sp=matmul(overlap,density_matrix);trace_value=(0d0,0d0)
    do i=1,nwann;trace_value=trace_value+sp(i,i);enddo
    trace_charge=real(trace_value,8);trace_imag=abs(aimag(trace_value))
    if(.not.all(ieee_is_finite([integrated_charge,trace_charge,trace_imag])).or.&
        trace_imag>tolerance*max(1d0,abs(trace_charge)))then
      message='nonfinite or complex overlapping-Wannier electron count';return
    endif
    if(abs(integrated_charge-trace_charge)>tolerance*max(1d0,abs(trace_charge)))then
      message='integrated density and Tr(P S) disagree';return
    endif
    ok=.true.;message=''
#else
    ok=.false.;message='overlapping-Wannier density reconstruction requires MPI'
#endif
  end subroutine

  logical function finite_matrix(matrix)
    complex(8),intent(in)::matrix(:,:)
    finite_matrix=all(ieee_is_finite(real(matrix,8))).and.all(ieee_is_finite(aimag(matrix)))
  end function
end module
