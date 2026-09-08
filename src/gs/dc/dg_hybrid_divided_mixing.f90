module dg_hybrid_divided_mixing
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  implicit none
  private

  type,public::s_dg_hybrid_divided_mixing_state
    logical::valid=.false.
    character(16)::method=''
    real(real64)::mixrate=0d0
    integer::basis_generation=0
    integer(int64)::inventory_fingerprint=0_int64
    integer::history_length=0
    integer::reset_count=0
    character(64)::last_reset_reason=''
    integer(int64)::fingerprint=0_int64
  end type s_dg_hybrid_divided_mixing_state

  public::prepare_dg_hybrid_divided_mixing,accept_dg_hybrid_divided_mixing
contains
  subroutine prepare_dg_hybrid_divided_mixing(comm,requested_method,inherited_method,configured_mixrate,&
      basis_generation,inventory_fingerprint,collective_rollback,state,selected_method,effective_mixrate,&
      mixing_iteration,reset_required,reset_reason,ok,message)
    integer,intent(in)::comm,basis_generation
    character(*),intent(in)::requested_method,inherited_method
    real(real64),intent(in)::configured_mixrate
    integer(int64),intent(in)::inventory_fingerprint
    logical,intent(in)::collective_rollback
    type(s_dg_hybrid_divided_mixing_state),intent(inout)::state
    character(*),intent(out)::selected_method,reset_reason,message
    real(real64),intent(out)::effective_mixrate
    integer,intent(out)::mixing_iteration
    logical,intent(out)::reset_required,ok
    type(s_dg_hybrid_divided_mixing_state)::candidate
    character(16)::resolved_method
    integer::method_code,minimum_code,maximum_code,minimum_generation,maximum_generation
    integer::rollback_local,rollback_global,ierr,local_bad,global_bad
    integer(int64)::minimum_inventory,maximum_inventory
    real(real64)::minimum_mixrate,maximum_mixrate

    ok=.false.;message='';selected_method='';reset_reason='';effective_mixrate=0d0
    mixing_iteration=0;reset_required=.false.
    resolved_method=resolve_method(requested_method,inherited_method)
    method_code=encode_method(resolved_method)
    local_bad=merge(0,1,method_code>0.and.configured_mixrate>=0d0.and.configured_mixrate<=1d0.and.&
      basis_generation>0.and.inventory_fingerprint/=0_int64)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      message='divided Hybrid mixing validation reduction failed';return
    endif
    if(global_bad/=0)then
      if(method_code==0)then
        message='unsupported divided Hybrid mixing method'
      else
        message='invalid divided Hybrid mixing lifecycle input'
      endif
      return
    endif

    call MPI_Allreduce(method_code,minimum_code,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='divided Hybrid mixing method reduction failed';return;endif
    call MPI_Allreduce(method_code,maximum_code,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='divided Hybrid mixing method reduction failed';return;endif
    call MPI_Allreduce(configured_mixrate,minimum_mixrate,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='divided Hybrid mix rate reduction failed';return;endif
    call MPI_Allreduce(configured_mixrate,maximum_mixrate,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='divided Hybrid mix rate reduction failed';return;endif
    call MPI_Allreduce(basis_generation,minimum_generation,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='divided Hybrid basis generation reduction failed';return;endif
    call MPI_Allreduce(basis_generation,maximum_generation,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='divided Hybrid basis generation reduction failed';return;endif
    call MPI_Allreduce(inventory_fingerprint,minimum_inventory,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='divided Hybrid inventory reduction failed';return;endif
    call MPI_Allreduce(inventory_fingerprint,maximum_inventory,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='divided Hybrid inventory reduction failed';return;endif
    rollback_local=merge(1,0,collective_rollback)
    call MPI_Allreduce(rollback_local,rollback_global,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='divided Hybrid rollback reduction failed';return;endif
    if(minimum_code/=maximum_code.or.minimum_mixrate/=maximum_mixrate.or.&
       minimum_generation/=maximum_generation.or.minimum_inventory/=maximum_inventory)then
      message='divided Hybrid mixing controls disagree across ranks';return
    endif

    candidate=state
    if(.not.candidate%valid)then
      candidate%valid=.true.
      candidate%method=resolved_method
      candidate%mixrate=configured_mixrate
      candidate%basis_generation=basis_generation
      candidate%inventory_fingerprint=inventory_fingerprint
      candidate%history_length=0
      candidate%reset_count=0
      candidate%last_reset_reason=''
    else
      if(trim(candidate%method)/=trim(resolved_method).or.candidate%mixrate/=configured_mixrate)then
        message='divided Hybrid mixing method and mix rate must persist across density epochs';return
      endif
      if(rollback_global/=0)then
        reset_required=.true.;reset_reason='collective_rollback'
      elseif(candidate%basis_generation/=basis_generation)then
        reset_required=.true.;reset_reason='basis_generation'
      elseif(candidate%inventory_fingerprint/=inventory_fingerprint)then
        reset_required=.true.;reset_reason='common_inventory'
      endif
      if(reset_required)then
        candidate%history_length=0
        candidate%reset_count=candidate%reset_count+1
        candidate%last_reset_reason=reset_reason
      endif
      candidate%basis_generation=basis_generation
      candidate%inventory_fingerprint=inventory_fingerprint
    endif
    candidate%fingerprint=fingerprint_state(candidate)
    state=candidate
    selected_method=state%method
    effective_mixrate=state%mixrate
    mixing_iteration=state%history_length+1
    ok=.true.
  end subroutine prepare_dg_hybrid_divided_mixing

  subroutine accept_dg_hybrid_divided_mixing(comm,mixing_iteration,state,ok,message,local_accept_ok)
    integer,intent(in)::comm,mixing_iteration
    type(s_dg_hybrid_divided_mixing_state),intent(inout)::state
    logical,intent(out)::ok
    character(*),intent(out)::message
    logical,intent(in),optional::local_accept_ok
    type(s_dg_hybrid_divided_mixing_state)::candidate
    integer::minimum_iteration,maximum_iteration,local_bad,global_bad,ierr

    ok=.false.;message=''
    local_bad=merge(0,1,state%valid.and.mixing_iteration==state%history_length+1)
    if(present(local_accept_ok))then
      if(.not.local_accept_ok)local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='divided Hybrid mixing acceptance reduction failed';return;endif
    call MPI_Allreduce(mixing_iteration,minimum_iteration,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='divided Hybrid mixing iteration reduction failed';return;endif
    call MPI_Allreduce(mixing_iteration,maximum_iteration,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='divided Hybrid mixing iteration reduction failed';return;endif
    if(global_bad/=0.or.minimum_iteration/=maximum_iteration)then
      message='collective divided Hybrid mixing epoch was not accepted';return
    endif
    candidate=state
    candidate%history_length=mixing_iteration
    candidate%fingerprint=fingerprint_state(candidate)
    state=candidate
    ok=.true.
  end subroutine accept_dg_hybrid_divided_mixing

  pure function resolve_method(requested_method,inherited_method) result(method)
    character(*),intent(in)::requested_method,inherited_method
    character(16)::method
    method=adjustl(requested_method)
    if(trim(method)=='inherit')method=adjustl(inherited_method)
  end function resolve_method

  pure integer function encode_method(method)
    character(*),intent(in)::method
    select case(trim(method))
    case('simple');encode_method=1
    case('pulay');encode_method=2
    case('broyden');encode_method=3
    case default;encode_method=0
    end select
  end function encode_method

  pure integer(int64) function fingerprint_state(state)
    type(s_dg_hybrid_divided_mixing_state),intent(in)::state
    integer(int64)::rate_bits
    rate_bits=transfer(state%mixrate,rate_bits)
    fingerprint_state=ieor(state%inventory_fingerprint,ishft(int(state%basis_generation,int64),7))
    fingerprint_state=ieor(fingerprint_state,ishft(int(state%history_length,int64),19))
    fingerprint_state=ieor(fingerprint_state,ishft(int(state%reset_count,int64),31))
    fingerprint_state=ieor(fingerprint_state,ishft(int(encode_method(state%method),int64),3))
    fingerprint_state=ieor(fingerprint_state,rate_bits)
    if(fingerprint_state==0_int64)fingerprint_state=1_int64
  end function fingerprint_state
end module dg_hybrid_divided_mixing
