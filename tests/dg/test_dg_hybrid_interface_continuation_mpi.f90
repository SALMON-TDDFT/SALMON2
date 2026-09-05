program test_dg_hybrid_interface_continuation_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_positive_inf,ieee_quiet_nan,ieee_value
  use dg_hybrid_interface_continuation,only:s_dg_hybrid_interface_continuation,&
    initialize_dg_hybrid_interface_continuation,accept_dg_hybrid_interface_point
  implicit none
  integer::ierr,rank,nproc

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  call check_schedule(0.2_real64,6)
  call check_schedule(0.3_real64,5)
  call check_initialization_rejections
  call check_acceptance_rejections
  call check_nonfinite_rejections
  call check_finished_rejection
  if(rank==0)write(*,'(a,i0,a)')'PASS DG interface continuation on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine check_schedule(rate,expected_acceptances)
    real(real64),intent(in)::rate
    integer,intent(in)::expected_acceptances
    type(s_dg_hybrid_interface_continuation)::state
    integer::lo,hi
    integer(int64)::fp_lo,fp_hi
    real(real64)::previous
    logical::ok
    character(512)::message

    call initialize_dg_hybrid_interface_continuation(MPI_COMM_WORLD,7,73491_int64,rate,state,ok,message)
    call require(ok,'valid schedule rejected: '//trim(message))
    call require(state%valid.and.state%lambda==0d0.and.state%step_index==0.and.&
      state%accepted_steps==0.and..not.state%finished,'initial continuation point is not exact zero')
    call require(state%fingerprint/=0_int64,'initial continuation fingerprint is zero')
    previous=-1d0
    do while(.not.state%finished)
      call require(state%lambda>=previous,'accepted continuation points are not monotone')
      previous=state%lambda
      call accept_dg_hybrid_interface_point(MPI_COMM_WORLD,7,73491_int64,.true.,state,ok,message)
      call require(ok,'valid continuation point rejected: '//trim(message))
    enddo
    call require(state%lambda==1d0,'terminal continuation point is not exactly one')
    call require(state%accepted_steps==expected_acceptances,'accepted point count is incorrect')
    call MPI_Allreduce(state%step_index,lo,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(state%step_index,hi,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    call require(lo==hi,'continuation step index differs across ranks')
    call MPI_Allreduce(state%fingerprint,fp_lo,1,MPI_INTEGER8,MPI_MIN,MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(state%fingerprint,fp_hi,1,MPI_INTEGER8,MPI_MAX,MPI_COMM_WORLD,ierr)
    call require(fp_lo==fp_hi,'continuation fingerprint differs across ranks')
  end subroutine check_schedule

  subroutine check_initialization_rejections
    type(s_dg_hybrid_interface_continuation)::state,snapshot
    logical::ok
    character(512)::message

    call initialize_dg_hybrid_interface_continuation(MPI_COMM_WORLD,7,73491_int64,0d0,state,ok,message)
    call require(.not.ok.and..not.state%valid,'zero continuation rate was accepted')
    call initialize_dg_hybrid_interface_continuation(MPI_COMM_WORLD,7,73491_int64,-0.1d0,state,ok,message)
    call require(.not.ok.and..not.state%valid,'negative continuation rate was accepted')
    call initialize_dg_hybrid_interface_continuation(MPI_COMM_WORLD,7,73491_int64,1.1d0,state,ok,message)
    call require(.not.ok.and..not.state%valid,'continuation rate above one was accepted')
    call initialize_dg_hybrid_interface_continuation(MPI_COMM_WORLD,7,73491_int64,0.2d0,state,ok,message)
    call require(ok,'reinitialization fixture setup failed: '//trim(message))
    snapshot=state
    call initialize_dg_hybrid_interface_continuation(MPI_COMM_WORLD,0,73491_int64,0.2d0,state,ok,message)
    call require(.not.ok.and.unchanged(state,snapshot),'failed reinitialization changed valid state')
    if(nproc>1)then
      call initialize_dg_hybrid_interface_continuation(MPI_COMM_WORLD,7,73491_int64,&
        merge(0.2d0,0.3d0,rank==0),state,ok,message)
      call require(.not.ok.and.index(message,'disagree')>0,'rank-disagreeing continuation rate was accepted')
      call initialize_dg_hybrid_interface_continuation(MPI_COMM_WORLD,7+merge(1,0,rank==0),73491_int64,&
        0.2d0,state,ok,message)
      call require(.not.ok.and.index(message,'disagree')>0,'rank-disagreeing basis generation was accepted')
      call initialize_dg_hybrid_interface_continuation(MPI_COMM_WORLD,7,&
        73491_int64+merge(1_int64,0_int64,rank==0),0.2d0,state,ok,message)
      call require(.not.ok.and.index(message,'disagree')>0,'rank-disagreeing mapping fingerprint was accepted')
    endif
  end subroutine check_initialization_rejections

  subroutine check_acceptance_rejections
    type(s_dg_hybrid_interface_continuation)::state,snapshot
    logical::ok
    character(512)::message

    call initialize_dg_hybrid_interface_continuation(MPI_COMM_WORLD,7,73491_int64,0.3d0,state,ok,message)
    call require(ok,'rejection fixture initialization failed: '//trim(message))
    snapshot=state
    call accept_dg_hybrid_interface_point(MPI_COMM_WORLD,7,73491_int64,&
      merge(.false.,.true.,rank==0),state,ok,message)
    call require(.not.ok.and.index(message,'rejected')>0,'rank-local point rejection was accepted')
    call require(unchanged(state,snapshot),'rank-local rejection advanced continuation state')
    if(nproc>1)then
      call accept_dg_hybrid_interface_point(MPI_COMM_WORLD,7+merge(1,0,rank==0),73491_int64,.true.,&
        state,ok,message)
      call require(.not.ok.and.index(message,'disagree')>0,'rank-disagreeing acceptance basis was accepted')
      call require(unchanged(state,snapshot),'basis disagreement advanced continuation state')
      call accept_dg_hybrid_interface_point(MPI_COMM_WORLD,7,&
        73491_int64+merge(1_int64,0_int64,rank==0),.true.,state,ok,message)
      call require(.not.ok.and.index(message,'disagree')>0,'rank-disagreeing acceptance mapping was accepted')
      call require(unchanged(state,snapshot),'mapping disagreement advanced continuation state')
    endif
  end subroutine check_acceptance_rejections

  subroutine check_nonfinite_rejections
    type(s_dg_hybrid_interface_continuation)::state,snapshot
    real(real64)::bad_rate
    integer(int64)::rate_bits,lambda_bits
    logical::ok
    character(512)::message

    bad_rate=ieee_value(0d0,ieee_quiet_nan)
    if(nproc>1.and.rank/=0)bad_rate=0.2d0
    call initialize_dg_hybrid_interface_continuation(MPI_COMM_WORLD,7,73491_int64,bad_rate,state,ok,message)
    call require(.not.ok,'rank-local NaN continuation rate was accepted')
    bad_rate=ieee_value(0d0,ieee_positive_inf)
    if(nproc>1.and.rank/=0)bad_rate=0.2d0
    call initialize_dg_hybrid_interface_continuation(MPI_COMM_WORLD,7,73491_int64,bad_rate,state,ok,message)
    call require(.not.ok,'rank-local infinite continuation rate was accepted')

    call initialize_dg_hybrid_interface_continuation(MPI_COMM_WORLD,7,73491_int64,0.2d0,state,ok,message)
    call require(ok,'nonfinite-state fixture initialization failed: '//trim(message))
    if(nproc==1.or.rank==0)state%lambda=ieee_value(0d0,ieee_quiet_nan)
    snapshot=state;rate_bits=transfer(state%rate,rate_bits);lambda_bits=transfer(state%lambda,lambda_bits)
    call accept_dg_hybrid_interface_point(MPI_COMM_WORLD,7,73491_int64,.true.,state,ok,message)
    call require(.not.ok.and.transfer(state%rate,rate_bits)==rate_bits.and.&
      transfer(state%lambda,lambda_bits)==lambda_bits.and.state%accepted_steps==snapshot%accepted_steps.and.&
      state%fingerprint==snapshot%fingerprint,'nonfinite internal state was accepted or mutated')

    call initialize_dg_hybrid_interface_continuation(MPI_COMM_WORLD,7,73491_int64,0.2d0,state,ok,message)
    call require(ok,'fingerprint-corruption fixture initialization failed: '//trim(message))
    state%fingerprint=state%fingerprint+1_int64;snapshot=state
    call accept_dg_hybrid_interface_point(MPI_COMM_WORLD,7,73491_int64,.true.,state,ok,message)
    call require(.not.ok.and.unchanged(state,snapshot),'corrupt fingerprint was accepted or repaired')
  end subroutine check_nonfinite_rejections

  subroutine check_finished_rejection
    type(s_dg_hybrid_interface_continuation)::state,snapshot
    logical::ok
    character(512)::message

    call initialize_dg_hybrid_interface_continuation(MPI_COMM_WORLD,7,73491_int64,1d0,state,ok,message)
    call require(ok,'finished-state fixture initialization failed: '//trim(message))
    call accept_dg_hybrid_interface_point(MPI_COMM_WORLD,7,73491_int64,.true.,state,ok,message)
    call require(ok.and.state%lambda==1d0.and..not.state%finished,'lambda-one point was not staged')
    call accept_dg_hybrid_interface_point(MPI_COMM_WORLD,7,73491_int64,.true.,state,ok,message)
    call require(ok.and.state%finished,'terminal lambda-one point was not accepted')
    snapshot=state
    call accept_dg_hybrid_interface_point(MPI_COMM_WORLD,7,73491_int64,.true.,state,ok,message)
    call require(.not.ok.and.unchanged(state,snapshot),'finished continuation was accepted again')
  end subroutine check_finished_rejection

  logical function unchanged(actual,expected)
    type(s_dg_hybrid_interface_continuation),intent(in)::actual,expected
    unchanged=actual%valid.eqv.expected%valid
    unchanged=unchanged.and.actual%basis_generation==expected%basis_generation.and.&
      actual%mapping_fingerprint==expected%mapping_fingerprint.and.actual%rate==expected%rate.and.&
      actual%lambda==expected%lambda.and.actual%step_index==expected%step_index.and.&
      actual%accepted_steps==expected%accepted_steps.and.(actual%finished.eqv.expected%finished).and.&
      actual%fingerprint==expected%fingerprint
  end function unchanged

  subroutine require(condition,detail)
    logical,intent(in)::condition
    character(*),intent(in)::detail
    integer::local_failure,global_failure,code
    local_failure=merge(0,1,condition)
    call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,code)
    if(global_failure/=0)then
      if(.not.condition)write(0,'(a,i0,2a)')'rank ',rank,': ',trim(detail)
      call MPI_Abort(MPI_COMM_WORLD,1,code)
    endif
  end subroutine require
end program test_dg_hybrid_interface_continuation_mpi
