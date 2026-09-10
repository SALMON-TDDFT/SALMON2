#include "config.h"
module dg_hybrid_ground_state_types
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  type,public::s_dg_hybrid_spectral_certification
    logical::valid=.false.,compatibility_dynamic_rank=.false.,proof_state_present=.false.
    integer::full_rank=0,occupied_rank=0,requested_rank=0,boundary_cluster_rank=0,certified_rank=0
    integer::extension_states=0,worst_operation=0
    real(real64)::energy_window=0d0,e_homo=0d0,requested_cutoff=0d0,certified_cutoff=0d0
    real(real64)::extension_energy=0d0,proof_energy=0d0
    real(real64)::occupied_subspace_defect=0d0,occupied_projector_defect=0d0
    real(real64)::target_subspace_defect=0d0,target_energy_defect=0d0,density_defect=0d0
    real(real64)::worst_operation_defect=0d0,maximum_physical_defect=0d0
    integer(int64)::fingerprint=0_int64
  end type s_dg_hybrid_spectral_certification
  type,public::s_dg_hybrid_ground_state
    logical::valid=.false.,converged=.false.
    logical::refinement_converged=.false.,refinement_exhausted=.false.
    integer::global_count=0,noccupied=0,final_eigensolve_count=0,additional_refinement_count=0
    integer(int64)::hybrid_basis_fingerprint=0_int64,metric_fingerprint=0_int64
    integer(int64)::operator_fingerprint=0_int64,position_fingerprint=0_int64
    integer(int64)::fingerprint=0_int64,workspace_peak_bytes=0_int64
    real(real64)::e_homo=0d0,terminal_density_change=huge(0d0),terminal_energy_change=huge(0d0)
    type(s_dg_hybrid_spectral_certification)::spectral_certification
    integer(int64),allocatable::owned_row_ids(:)
    complex(real64),allocatable::coefficients(:,:)
    real(real64),allocatable::occupations(:),eigenvalues(:)
  end type s_dg_hybrid_ground_state
  public::validate_dg_hybrid_ground_state
contains
  subroutine validate_dg_hybrid_ground_state(comm,global_count,noccupied,row_ids,coefficients,occupations,eigenvalues,&
      expected_electron_count,hybrid_basis_fingerprint,metric_fingerprint,operator_fingerprint,position_fingerprint,&
      tolerance,state,workspace_peak_bytes,fingerprint,ok,message)
    integer,intent(in)::comm,global_count,noccupied
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::coefficients(:,:)
    real(real64),intent(in)::occupations(:),eigenvalues(:),expected_electron_count,tolerance
    integer(int64),intent(in)::hybrid_basis_fingerprint,metric_fingerprint,operator_fingerprint,position_fingerprint
    type(s_dg_hybrid_ground_state),intent(out)::state
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::i,j,row,nowned,ierr,rank,local_bad,global_bad,allocation_status,minimum_integer,maximum_integer
    integer,allocatable::ownership_count(:)
    integer(int64)::minimum_receipt,maximum_receipt,local_hash,global_hash,entry_hash,bits,complex_elements,real_elements
    real(real64)::electron_count,minimum_real,maximum_real
    ok=.false.;message='';workspace_peak_bytes=0_int64;fingerprint=0_int64;state%valid=.false.
    nowned=size(row_ids);local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)then;message='hybrid state communicator failed';return;endif
    call agree_integer(global_count,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='rank-disagreeing hybrid state extent';return;endif
    call agree_integer(noccupied,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='rank-disagreeing occupied-state count';return;endif
    call agree_real_bits(tolerance,minimum_receipt,maximum_receipt,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_receipt/=maximum_receipt)then;message='rank-disagreeing hybrid state tolerance';return;endif
    call agree_real_bits(expected_electron_count,minimum_receipt,maximum_receipt,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_receipt/=maximum_receipt)then;message='rank-disagreeing electron count';return;endif
    call agree_receipt(hybrid_basis_fingerprint,minimum_receipt,maximum_receipt,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_receipt/=maximum_receipt)then;message='rank-disagreeing hybrid basis receipt';return;endif
    call agree_receipt(metric_fingerprint,minimum_receipt,maximum_receipt,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_receipt/=maximum_receipt)then;message='rank-disagreeing hybrid metric receipt';return;endif
    call agree_receipt(operator_fingerprint,minimum_receipt,maximum_receipt,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_receipt/=maximum_receipt)then;message='rank-disagreeing hybrid operator receipt';return;endif
    call agree_receipt(position_fingerprint,minimum_receipt,maximum_receipt,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_receipt/=maximum_receipt)then;message='rank-disagreeing hybrid position receipt';return;endif
    if(global_count<1.or.noccupied<1.or.noccupied>global_count)local_bad=1
    if(size(coefficients,1)/=nowned.or.size(coefficients,2)/=noccupied)local_bad=1
    if(size(occupations)/=noccupied.or.size(eigenvalues)/=noccupied)local_bad=1
    if(any(row_ids<1_int64).or.any(row_ids>int(max(0,global_count),int64)))local_bad=1
    if(.not.ieee_is_finite(tolerance).or.tolerance<1d-15.or.tolerance>1d-2)local_bad=1
    if(.not.ieee_is_finite(expected_electron_count).or.expected_electron_count<0d0)local_bad=1
    if(hybrid_basis_fingerprint==0_int64.or.metric_fingerprint==0_int64.or.&
      operator_fingerprint==0_int64.or.position_fingerprint==0_int64)local_bad=1
    if(.not.finite_complex(coefficients).or..not.finite_real(occupations).or..not.finite_real(eigenvalues))local_bad=1
    if(finite_real(occupations))then;if(any(occupations<0d0))local_bad=1;endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid hybrid occupied-state contract';return;endif
    local_bad=merge(0,1,all(occupations>64d0*epsilon(1d0)))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='hybrid state occupation threshold excludes a published column';return
    endif
    if(noccupied>1)then
      local_bad=merge(0,1,all(eigenvalues(2:)>=eigenvalues(:noccupied-1)))
    else
      local_bad=0
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='hybrid occupied eigenvalues must be ascending';return
    endif
    do i=1,noccupied
      call agree_real_bits(occupations(i),minimum_receipt,maximum_receipt,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_receipt/=maximum_receipt)then;message='rank-disagreeing occupations';return;endif
      call agree_real_bits(eigenvalues(i),minimum_receipt,maximum_receipt,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_receipt/=maximum_receipt)then;message='rank-disagreeing eigenvalues';return;endif
    enddo
    electron_count=sum(occupations)
    call MPI_Allreduce(electron_count,minimum_real,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='electron-count minimum failed';return;endif
    call MPI_Allreduce(electron_count,maximum_real,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_real/=maximum_real.or.&
      abs(electron_count-expected_electron_count)>tolerance*max(1d0,expected_electron_count))then
      message='hybrid occupied-state electron count is invalid';return
    endif
    if(int(nowned,int64)>huge(complex_elements)/int(noccupied,int64))local_bad=1
    if(local_bad==0)then
      complex_elements=int(nowned,int64)*int(noccupied,int64)
      real_elements=2_int64*int(noccupied,int64)
      if(complex_elements>huge(workspace_peak_bytes)/16_int64)local_bad=1
      if(real_elements>huge(workspace_peak_bytes)/8_int64)local_bad=1
      if(int(nowned,int64)>huge(workspace_peak_bytes)/8_int64)local_bad=1
      if(int(global_count,int64)>huge(workspace_peak_bytes)/4_int64)local_bad=1
      if(local_bad==0)then
        workspace_peak_bytes=16_int64*complex_elements+8_int64*real_elements+8_int64*int(nowned,int64)
        if(workspace_peak_bytes>huge(workspace_peak_bytes)-4_int64*int(global_count,int64))local_bad=1
        if(local_bad==0)workspace_peak_bytes=workspace_peak_bytes+4_int64*int(global_count,int64)
      endif
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='hybrid state workspace receipt overflow';return;endif
    allocate(ownership_count(global_count),state%owned_row_ids(nowned),state%coefficients(nowned,noccupied),&
      state%occupations(noccupied),state%eigenvalues(noccupied),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='cannot allocate hybrid occupied state';return;endif
    ownership_count=0
    do i=1,nowned;row=int(row_ids(i));ownership_count(row)=ownership_count(row)+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership_count,global_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid state ownership reduction failed';return;endif
    if(any(ownership_count/=1))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='hybrid state rows are not owned exactly once';return;endif
    local_hash=0_int64
    do i=1,nowned
      row=int(row_ids(i))
      do j=1,noccupied
        entry_hash=ieor(int(row,int64),ishftc(int(j,int64),11))
        bits=transfer(real(coefficients(i,j),real64),bits);entry_hash=ieor(entry_hash,ishftc(bits,19))
        bits=transfer(aimag(coefficients(i,j)),bits);entry_hash=ieor(entry_hash,ishftc(bits,37))
        local_hash=ieor(local_hash,entry_hash)
      enddo
    enddo
    call MPI_Allreduce(local_hash,global_hash,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid state fingerprint reduction failed';return;endif
    fingerprint=ieor(global_hash,hybrid_basis_fingerprint)
    fingerprint=ieor(fingerprint,ishftc(metric_fingerprint,7))
    fingerprint=ieor(fingerprint,ishftc(operator_fingerprint,13))
    fingerprint=ieor(fingerprint,ishftc(position_fingerprint,23))
    do i=1,noccupied
      bits=transfer(occupations(i),bits);fingerprint=ieor(fingerprint,ishftc(bits,mod(5*i,63)))
      bits=transfer(eigenvalues(i),bits);fingerprint=ieor(fingerprint,ishftc(bits,mod(9*i,63)))
    enddo
    if(fingerprint==0_int64)fingerprint=ieor(global_hash,719_int64)
    state%valid=.true.;state%converged=.false.;state%global_count=global_count;state%noccupied=noccupied
    state%final_eigensolve_count=0
    state%hybrid_basis_fingerprint=hybrid_basis_fingerprint;state%metric_fingerprint=metric_fingerprint
    state%operator_fingerprint=operator_fingerprint;state%position_fingerprint=position_fingerprint
    state%fingerprint=fingerprint;state%workspace_peak_bytes=workspace_peak_bytes
    state%e_homo=eigenvalues(noccupied)
    state%owned_row_ids=row_ids;state%coefficients=coefficients
    state%occupations=occupations;state%eigenvalues=eigenvalues
    ok=.true.;message=''
  contains
    subroutine agree_integer(value,minimum_value,maximum_value,status)
      integer,intent(in)::value
      integer,intent(out)::minimum_value,maximum_value,status
      call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER,MPI_MIN,comm,status);if(status/=MPI_SUCCESS)return
      call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER,MPI_MAX,comm,status)
    end subroutine agree_integer
    subroutine agree_receipt(value,minimum_value,maximum_value,status)
      integer(int64),intent(in)::value
      integer(int64),intent(out)::minimum_value,maximum_value
      integer,intent(out)::status
      call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER8,MPI_MIN,comm,status);if(status/=MPI_SUCCESS)return
      call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER8,MPI_MAX,comm,status)
    end subroutine agree_receipt
    subroutine agree_real_bits(value,minimum_value,maximum_value,status)
      real(real64),intent(in)::value
      integer(int64),intent(out)::minimum_value,maximum_value
      integer,intent(out)::status
      integer(int64)::value_bits
      value_bits=transfer(value,value_bits);call agree_receipt(value_bits,minimum_value,maximum_value,status)
    end subroutine agree_real_bits
    logical function finite_complex(values)
      complex(real64),intent(in)::values(:,:)
      finite_complex=all(ieee_is_finite(real(values,real64))).and.all(ieee_is_finite(aimag(values)))
    end function finite_complex
    logical function finite_real(values)
      real(real64),intent(in)::values(:)
      finite_real=all(ieee_is_finite(values))
    end function finite_real
    subroutine cleanup()
      if(allocated(ownership_count))deallocate(ownership_count)
      if(allocated(state%owned_row_ids))deallocate(state%owned_row_ids)
      if(allocated(state%coefficients))deallocate(state%coefficients)
      if(allocated(state%occupations))deallocate(state%occupations)
      if(allocated(state%eigenvalues))deallocate(state%eigenvalues)
      state%valid=.false.;state%converged=.false.;state%final_eigensolve_count=0
      state%e_homo=0d0
      state%spectral_certification=s_dg_hybrid_spectral_certification()
    end subroutine cleanup
#else
    ok=.false.;message='MPI is required for hybrid ground-state validation'
    workspace_peak_bytes=0_int64;fingerprint=0_int64
#endif
  end subroutine validate_dg_hybrid_ground_state
end module dg_hybrid_ground_state_types
