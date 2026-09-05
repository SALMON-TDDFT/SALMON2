program test_dg_hybrid_divided_mixing_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_divided_mixing,only:s_dg_hybrid_divided_mixing_state,&
    prepare_dg_hybrid_divided_mixing,accept_dg_hybrid_divided_mixing
  implicit none
  integer::ierr,rank,nproc,mixing_iteration
  type(s_dg_hybrid_divided_mixing_state)::state,snapshot
  character(16)::method
  character(64)::reset_reason
  character(512)::message
  real(real64)::effective_mixrate
  logical::ok,reset_required

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  call require(any(nproc==[2,4,8]),'test requires 2, 4, or 8 ranks')

  call prepare('inherit','pulay',0.25d0,4,8173_int64,.false.)
  call require(ok,'initial inherited Pulay preparation failed: '//trim(message))
  call require(method=='pulay'.and.abs(effective_mixrate-0.25d0)<1d-15,&
    'inherited method or mix rate was not retained')
  call require(mixing_iteration==1.and..not.reset_required,&
    'fresh history did not start at one without a diagnosed reset')
  call accept_dg_hybrid_divided_mixing(MPI_COMM_WORLD,mixing_iteration,state,ok,message)
  call require(ok.and.state%history_length==1,'first accepted epoch was not recorded')

  call prepare('inherit','pulay',0.25d0,4,8173_int64,.false.)
  call require(ok.and.mixing_iteration==2.and..not.reset_required,&
    'ordinary accepted epoch reset the Pulay history')
  call accept_dg_hybrid_divided_mixing(MPI_COMM_WORLD,mixing_iteration,state,ok,message)
  call require(ok.and.state%history_length==2,'second accepted epoch did not grow history')

  call prepare('pulay','simple',0.25d0,4,8173_int64,.false.)
  call require(ok.and.method=='pulay'.and.mixing_iteration==3,&
    'explicit Pulay dispatch did not preserve history')
  call accept_dg_hybrid_divided_mixing(MPI_COMM_WORLD,mixing_iteration,state,ok,message)
  call require(ok.and.state%history_length==3,'third accepted epoch did not grow history')

  call prepare('pulay','simple',0.25d0,5,8173_int64,.false.)
  call require(ok.and.reset_required.and.mixing_iteration==1,&
    'basis generation change did not reset history')
  call require(trim(reset_reason)=='basis_generation','wrong basis reset reason')
  call accept_dg_hybrid_divided_mixing(MPI_COMM_WORLD,mixing_iteration,state,ok,message)
  call require(ok.and.state%history_length==1,'basis-reset epoch was not accepted')

  call prepare('pulay','simple',0.25d0,5,9917_int64,.false.)
  call require(ok.and.reset_required.and.trim(reset_reason)=='common_inventory',&
    'common inventory extension did not diagnose a reset')
  call accept_dg_hybrid_divided_mixing(MPI_COMM_WORLD,mixing_iteration,state,ok,message)
  call require(ok,'inventory-reset epoch was not accepted')

  call prepare('pulay','simple',0.25d0,5,9917_int64,.true.)
  call require(ok.and.reset_required.and.trim(reset_reason)=='collective_rollback',&
    'collective rollback did not diagnose a reset')
  call accept_dg_hybrid_divided_mixing(MPI_COMM_WORLD,mixing_iteration,state,ok,message)
  call require(ok.and.state%reset_count==3,'diagnosed reset count is wrong')

  call exercise_dispatch('simple')
  call exercise_dispatch('broyden')

  snapshot=state
  call prepare('simple','pulay',0.25d0,5,9917_int64,.false.)
  call require(.not.ok,'mixer changed during persistent history')
  call require(state%fingerprint==snapshot%fingerprint,'method change changed accepted state')

  snapshot=state
  call prepare('pulay','pulay',0.20d0,5,9917_int64,.false.)
  call require(.not.ok,'mix rate changed during persistent history')
  call require(state%fingerprint==snapshot%fingerprint,'mix-rate change changed accepted state')

  snapshot=state
  call prepare('simple_potential','pulay',0.25d0,5,9917_int64,.false.)
  call require(.not.ok,'unsupported divided mixer was accepted')
  call require(state%fingerprint==snapshot%fingerprint,'unsupported method changed accepted state')

  snapshot=state
  call prepare('pulay','pulay',0.25d0+merge(0.01d0,0d0,rank==0),5,9917_int64,.false.)
  call require(.not.ok,'rank-disagreeing mix rate was accepted')
  call require(state%fingerprint==snapshot%fingerprint,'collective failure changed accepted state')

  if(rank==0)write(*,'(a,i0,a)')'PASS divided Hybrid mixing on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine prepare(requested,inherited,rate,generation,inventory,rollback)
    character(*),intent(in)::requested,inherited
    real(real64),intent(in)::rate
    integer,intent(in)::generation
    integer(int64),intent(in)::inventory
    logical,intent(in)::rollback
    call prepare_dg_hybrid_divided_mixing(MPI_COMM_WORLD,requested,inherited,rate,generation,&
      inventory,rollback,state,method,effective_mixrate,mixing_iteration,reset_required,&
      reset_reason,ok,message)
  end subroutine prepare

  subroutine exercise_dispatch(requested)
    character(*),intent(in)::requested
    type(s_dg_hybrid_divided_mixing_state)::fresh
    call prepare_dg_hybrid_divided_mixing(MPI_COMM_WORLD,requested,'pulay',0.15d0,1,77_int64,&
      .false.,fresh,method,effective_mixrate,mixing_iteration,reset_required,reset_reason,ok,message)
    call require(ok.and.trim(method)==trim(requested),'valid mixer dispatch failed: '//trim(requested))
    call require(abs(effective_mixrate-0.15d0)<1d-15,'valid mixer changed mix rate')
  end subroutine exercise_dispatch

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
end program test_dg_hybrid_divided_mixing_mpi
