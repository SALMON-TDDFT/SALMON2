module dg_hybrid_interface_continuation
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private

  type,public::s_dg_hybrid_interface_continuation
    logical::valid=.false.
    integer::basis_generation=0
    integer(int64)::mapping_fingerprint=0_int64
    real(real64)::rate=0d0
    real(real64)::lambda=0d0
    integer::step_index=0
    integer::accepted_steps=0
    logical::finished=.false.
    integer(int64)::fingerprint=0_int64
  end type s_dg_hybrid_interface_continuation

  public::initialize_dg_hybrid_interface_continuation,accept_dg_hybrid_interface_point
contains
  subroutine initialize_dg_hybrid_interface_continuation(comm,basis_generation,mapping_fingerprint,&
      rate,state,ok,message,full_from_start)
    integer,intent(in)::comm,basis_generation
    integer(int64),intent(in)::mapping_fingerprint
    real(real64),intent(in)::rate
    type(s_dg_hybrid_interface_continuation),intent(inout)::state
    logical,intent(out)::ok
    character(*),intent(out)::message
    logical,optional,intent(in)::full_from_start
    type(s_dg_hybrid_interface_continuation)::candidate
    integer::controls(2),minimum_controls(2),maximum_controls(2),ierr,local_bad,global_bad
    integer(int64)::fingerprints(1),minimum_fingerprints(1),maximum_fingerprints(1)
    real(real64)::rates(1),minimum_rates(1),maximum_rates(1)
    logical::full

    ok=.false.;message='';full=.false.;if(present(full_from_start))full=full_from_start
    controls=[basis_generation,merge(1,0,full)];fingerprints=[mapping_fingerprint];rates=[rate]
    local_bad=merge(0,1,ieee_is_finite(rate))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DG interface continuation validation reduction failed';return;endif
    if(global_bad/=0)then;message='nonfinite DG interface continuation rate';return;endif
    local_bad=merge(0,1,basis_generation>0.and.mapping_fingerprint/=0_int64.and.rate>0d0.and.rate<=1d0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DG interface continuation validation reduction failed';return;endif
    if(global_bad/=0)then;message='invalid DG interface continuation controls';return;endif
    call MPI_Allreduce(controls,minimum_controls,2,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DG interface continuation control reduction failed';return;endif
    call MPI_Allreduce(controls,maximum_controls,2,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DG interface continuation control reduction failed';return;endif
    call MPI_Allreduce(fingerprints,minimum_fingerprints,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DG interface continuation fingerprint reduction failed';return;endif
    call MPI_Allreduce(fingerprints,maximum_fingerprints,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DG interface continuation fingerprint reduction failed';return;endif
    call MPI_Allreduce(rates,minimum_rates,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DG interface continuation rate reduction failed';return;endif
    call MPI_Allreduce(rates,maximum_rates,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DG interface continuation rate reduction failed';return;endif
    if(any(controls/=minimum_controls).or.any(controls/=maximum_controls).or.&
       any(fingerprints/=minimum_fingerprints).or.any(fingerprints/=maximum_fingerprints).or.&
       any(rates/=minimum_rates).or.any(rates/=maximum_rates))then
      message='DG interface continuation controls disagree across ranks';return
    endif
    candidate%valid=.true.
    candidate%basis_generation=basis_generation
    candidate%mapping_fingerprint=mapping_fingerprint
    candidate%rate=rate
    candidate%lambda=merge(1d0,0d0,full)
    candidate%step_index=0
    candidate%accepted_steps=0
    candidate%finished=.false.
    candidate%fingerprint=state_fingerprint(candidate)
    state=candidate;ok=.true.;message=''
  end subroutine initialize_dg_hybrid_interface_continuation

  subroutine accept_dg_hybrid_interface_point(comm,basis_generation,mapping_fingerprint,local_accept,&
      state,ok,message)
    integer,intent(in)::comm,basis_generation
    integer(int64),intent(in)::mapping_fingerprint
    logical,intent(in)::local_accept
    type(s_dg_hybrid_interface_continuation),intent(inout)::state
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(s_dg_hybrid_interface_continuation)::candidate
    integer::controls(6),minimum_controls(6),maximum_controls(6),ierr,local_bad,global_bad
    integer::local_vote,global_vote
    integer(int64)::fingerprints(3),minimum_fingerprints(3),maximum_fingerprints(3)
    real(real64)::reals(2),minimum_reals(2),maximum_reals(2)

    ok=.false.;message=''
    controls=[basis_generation,state%basis_generation,state%step_index,state%accepted_steps,&
      merge(1,0,state%valid),merge(1,0,state%finished)]
    fingerprints=[mapping_fingerprint,state%mapping_fingerprint,state%fingerprint]
    reals=[state%rate,state%lambda]
    local_bad=merge(0,1,ieee_is_finite(state%rate).and.ieee_is_finite(state%lambda))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DG interface point validation reduction failed';return;endif
    if(global_bad/=0)then;message='nonfinite DG interface continuation state';return;endif
    local_bad=merge(0,1,state%valid.and.state%basis_generation>0.and.&
      state%mapping_fingerprint/=0_int64.and.state%rate>0d0.and.state%rate<=1d0.and.&
      state%lambda>=0d0.and.state%lambda<=1d0.and.state%step_index>=0.and.&
      state%accepted_steps>=0.and.state%fingerprint==state_fingerprint(state))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DG interface point validation reduction failed';return;endif
    if(global_bad/=0)then;message='invalid or stale DG interface continuation state';return;endif
    call MPI_Allreduce(controls,minimum_controls,6,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DG interface point control reduction failed';return;endif
    call MPI_Allreduce(controls,maximum_controls,6,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DG interface point control reduction failed';return;endif
    call MPI_Allreduce(fingerprints,minimum_fingerprints,3,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DG interface point fingerprint reduction failed';return;endif
    call MPI_Allreduce(fingerprints,maximum_fingerprints,3,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DG interface point fingerprint reduction failed';return;endif
    call MPI_Allreduce(reals,minimum_reals,2,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DG interface point real-state reduction failed';return;endif
    call MPI_Allreduce(reals,maximum_reals,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DG interface point real-state reduction failed';return;endif
    if(any(controls/=minimum_controls).or.any(controls/=maximum_controls).or.&
       any(fingerprints/=minimum_fingerprints).or.any(fingerprints/=maximum_fingerprints).or.&
       any(reals/=minimum_reals).or.any(reals/=maximum_reals))then
      message='DG interface continuation controls disagree across ranks';return
    endif
    if(basis_generation/=state%basis_generation.or.mapping_fingerprint/=state%mapping_fingerprint)then
      message='invalid or stale DG interface continuation state';return
    endif
    if(state%finished)then;message='DG interface continuation already finished';return;endif
    local_vote=merge(1,0,local_accept)
    call MPI_Allreduce(local_vote,global_vote,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DG interface point acceptance reduction failed';return;endif
    if(global_vote==0)then;message='collective DG interface point rejected';return;endif

    candidate=state
    candidate%accepted_steps=candidate%accepted_steps+1
    if(candidate%lambda==1d0)then
      candidate%finished=.true.
    else
      candidate%lambda=min(1d0,candidate%lambda+candidate%rate)
      candidate%step_index=candidate%step_index+1
    endif
    candidate%fingerprint=state_fingerprint(candidate)
    state=candidate;ok=.true.;message=''
  end subroutine accept_dg_hybrid_interface_point

  pure integer(int64) function state_fingerprint(state) result(hash)
    type(s_dg_hybrid_interface_continuation),intent(in)::state
    integer(int64)::rate_bits,lambda_bits
    rate_bits=transfer(state%rate,rate_bits);lambda_bits=transfer(state%lambda,lambda_bits)
    hash=667182210696128091_int64
    hash=ieor(ishftc(hash,9),merge(1_int64,0_int64,state%valid))
    hash=ieor(ishftc(hash,9),int(state%basis_generation,int64))
    hash=ieor(ishftc(hash,9),state%mapping_fingerprint)
    hash=ieor(ishftc(hash,9),rate_bits)
    hash=ieor(ishftc(hash,9),lambda_bits)
    hash=ieor(ishftc(hash,9),int(state%step_index,int64))
    hash=ieor(ishftc(hash,9),int(state%accepted_steps,int64))
    hash=ieor(ishftc(hash,9),merge(1_int64,0_int64,state%finished))
    if(hash==0_int64)hash=1_int64
  end function state_fingerprint
end module dg_hybrid_interface_continuation
