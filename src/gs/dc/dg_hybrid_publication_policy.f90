#include "config.h"
module dg_hybrid_publication_policy
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi, only: MPI_Allreduce, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_INTEGER8, MPI_MAX, MPI_MIN, MPI_SUCCESS
#endif
  implicit none
  private
  type,public::s_dg_hybrid_candidate_acceptance
    logical::valid=.false.,complete_lcfo=.false.,occupation_policy=.false.
    logical::occupied_gate=.false.,density_gate=.false.,spectral_certified=.false.
    logical::certified_basis_ready=.false.,publication_authorized=.false.
    logical::legacy_dynamic_rank=.false.,legacy_warning_required=.false.,legacy_warning_observed=.false.
    integer::phase=0,construction_rank=0,solved_rank=0,occupied_rank=0
    integer::requested_rank=0,certified_rank=0,rt_basis_rank=0,published_rt_rank=0
    integer::checkpoint_version=0
    real(real64)::energy_window=0d0,occupied_gate_defect=0d0,density_gate_defect=0d0,&
      gate_tolerance=0d0
    integer(int64)::eigensystem_fingerprint=0_int64,occupation_fingerprint=0_int64
    integer(int64)::certification_fingerprint=0_int64,basis_fingerprint=0_int64
    integer(int64)::operator_fingerprint=0_int64
  end type s_dg_hybrid_candidate_acceptance
  public::validate_dg_hybrid_v5_publication_rank_policy,commit_candidate_transition
contains
  subroutine validate_dg_hybrid_v5_publication_rank_policy(icomm,receipt,energy_window,&
      requested_rank,certified_rank,construction_rank,ok,message)
    integer,intent(in)::icomm,requested_rank,certified_rank,construction_rank
    type(s_dg_hybrid_candidate_acceptance),intent(in)::receipt
    real(real64),intent(in)::energy_window
    logical,intent(out)::ok
    character(*),intent(out)::message
    logical::valid,full_rank_legacy
    type(s_dg_hybrid_candidate_acceptance)::candidate,scratch
    integer::ierr,signature(3),signature_min(3),signature_max(3)
    integer(int64)::signature64(2),signature64_min(2),signature64_max(2)
    candidate=receipt;scratch=receipt
    signature=[requested_rank,certified_rank,construction_rank]
    signature64=[transfer(energy_window,signature64(1)),candidate_acceptance_fingerprint(receipt)]
    call MPI_Allreduce(signature,signature_min,3,MPI_INTEGER,MPI_MIN,icomm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(signature,signature_max,3,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(signature64,signature64_min,2,MPI_INTEGER8,MPI_MIN,icomm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(signature64,signature64_max,2,MPI_INTEGER8,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(signature_min/=signature_max).or.any(signature64_min/=signature64_max))then
      ok=.false.;message='Hybrid v5 publication rank policy arguments disagree across MPI ranks';return
    endif
    full_rank_legacy=energy_window==-1d0.and.receipt%valid.and.receipt%publication_authorized.and.&
      receipt%legacy_dynamic_rank.and.receipt%legacy_warning_observed.and.receipt%energy_window==-1d0.and.&
      receipt%certified_rank==construction_rank.and.receipt%construction_rank==construction_rank.and.&
      receipt%published_rt_rank==construction_rank
    valid=construction_rank>0.and.requested_rank>0.and.requested_rank<=certified_rank.and.&
      certified_rank<=construction_rank.and.ieee_is_finite(energy_window).and.&
      (energy_window>=0d0.or.energy_window==-1d0).and.&
      (certified_rank<construction_rank.or.full_rank_legacy)
    call commit_candidate_transition(icomm,valid,candidate,scratch,&
      'Hybrid v5 full-rank certification requires authenticated energy_window=-1 authorization',ok,message)
  end subroutine validate_dg_hybrid_v5_publication_rank_policy

  subroutine commit_candidate_transition(icomm,local_valid,candidate,receipt,failure_message,ok,message)
    integer,intent(in)::icomm
    logical,intent(in)::local_valid
    type(s_dg_hybrid_candidate_acceptance),intent(in)::candidate
    type(s_dg_hybrid_candidate_acceptance),intent(inout)::receipt
    character(*),intent(in)::failure_message
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::local_bad,global_bad,ierr
    integer(int64)::local_hash,minimum_hash,maximum_hash
    local_bad=merge(0,1,local_valid)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      ok=.false.;message=failure_message;return
    endif
    local_hash=candidate_acceptance_fingerprint(candidate)
    call MPI_Allreduce(local_hash,minimum_hash,1,MPI_INTEGER8,MPI_MIN,icomm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_hash,maximum_hash,1,MPI_INTEGER8,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash)then
      ok=.false.;message='rank-disagreeing Hybrid final-candidate acceptance';return
    endif
    receipt=candidate;ok=.true.;message=''
#else
    ok=.false.;message='Hybrid final-candidate acceptance requires MPI'
#endif
  end subroutine commit_candidate_transition

  integer(int64) function candidate_acceptance_fingerprint(receipt) result(hash)
    type(s_dg_hybrid_candidate_acceptance),intent(in)::receipt
    hash=int(z'6A09E667F3BCC909',int64)
    call mix(int(receipt%phase,int64));call mix(int(receipt%construction_rank,int64))
    call mix(int(receipt%solved_rank,int64));call mix(int(receipt%occupied_rank,int64))
    call mix(int(receipt%requested_rank,int64));call mix(int(receipt%certified_rank,int64))
    call mix(int(receipt%rt_basis_rank,int64));call mix(int(receipt%published_rt_rank,int64))
    call mix(int(receipt%checkpoint_version,int64));call mix(transfer(receipt%energy_window,hash))
    call mix(transfer(receipt%occupied_gate_defect,hash))
    call mix(transfer(receipt%density_gate_defect,hash));call mix(transfer(receipt%gate_tolerance,hash))
    call mix(receipt%eigensystem_fingerprint);call mix(receipt%occupation_fingerprint)
    call mix(receipt%certification_fingerprint);call mix(receipt%basis_fingerprint)
    call mix(receipt%operator_fingerprint)
    call mix(int(merge(1,0,receipt%valid),int64));call mix(int(merge(1,0,receipt%complete_lcfo),int64))
    call mix(int(merge(1,0,receipt%occupation_policy),int64))
    call mix(int(merge(1,0,receipt%occupied_gate),int64));call mix(int(merge(1,0,receipt%density_gate),int64))
    call mix(int(merge(1,0,receipt%spectral_certified),int64))
    call mix(int(merge(1,0,receipt%certified_basis_ready),int64))
    call mix(int(merge(1,0,receipt%publication_authorized),int64))
    call mix(int(merge(1,0,receipt%legacy_dynamic_rank),int64))
    call mix(int(merge(1,0,receipt%legacy_warning_required),int64))
    call mix(int(merge(1,0,receipt%legacy_warning_observed),int64))
  contains
    subroutine mix(value)
      integer(int64),intent(in)::value
      hash=ieor(ishftc(hash,11),value)
      if(hash==0_int64)hash=int(z'510E527FADE682D1',int64)
    end subroutine mix
  end function candidate_acceptance_fingerprint
end module dg_hybrid_publication_policy
