! Copyright 2026 SALMONpy contributors. Apache-2.0.
! RT ownership extracted from fast-PR exx_native at c9e8ce431f43eefd4a88998134adb4796f8ec837.
! ACE/history algorithms remain in their original modules. No opaque handles
! are cloned or released here; ACE's metric communicator identifiers are borrowed.
module exx_rt_state
 use iso_fortran_env,only:int64
 use exx_k_exchange,only:exx_k_kernel,exx_k_kernel_destroy
 use exx_symmetry,only:s_exx_symmetry_map
 use exx_ace,only:s_exx_ace
 use exx_factor_history,only:s_factor_history
 use exx_spatial,only:spatial_exx_state
 implicit none
 private
 public::s_exx_rt_state,rt_taylor_stage,rt_ptcn_begin,rt_state_clear,rt_state_bytes,rt_values_copy
 type s_exx_rt_state
  type(spatial_exx_state)::spatial
  logical::recover_errors=.false.
  integer::step_error=0
  ! FFT owners are never intrinsically copied. Close before releasing this state.
  logical :: native_resources_enabled=.false.
  type(exx_k_kernel),allocatable :: native_kernel
  type(s_exx_symmetry_map) :: native_symmetry
  logical::exx_freeze=.false.,exx_force_full_action=.false.
  logical::exx_adaptive_ready=.false.,exx_support_changed=.false.
  real(8)::exx_exchange_energy=0d0
  real(8),allocatable::cached_occupation(:,:)
  complex(8),allocatable::cached_action(:,:,:)
  logical::finite_support_localized=.true.,radius_warning_reported=.false.
  logical::adaptive_active=.false.,cached_adaptive_ready=.false.
  type(s_exx_ace)::scf_damped_ace,ace,initial_ace,midpoint_ace
  logical::scf_operator_damped=.false.
  type(s_factor_history),allocatable::factor_history(:),ptcn_history(:)
  logical::factor_initialized=.false.,factor_configured=.false.,factor_experiment=.false.,factor_learned=.false.
  logical::factor_prediction_installed=.false.,factor_all_ready=.false.
  real(8),allocatable::factor_occupation(:,:,:)
  integer::factor_step=0
  character(512)::factor_restart_directory=''
  integer::factor_restart_step=-1
  integer::ptcn_accepted_step=0
  logical::ptcn_attempted=.false.,ptcn_predicted=.false.
  real(8)::ptcn_history_dt=0d0
  complex(8),allocatable::full_source(:,:,:),initial_source(:,:,:)
  logical::taylor_active=.false.,taylor_midpoint=.false.
  complex(8),allocatable::cached_source(:,:,:),target_work(:,:,:),action_work(:,:,:)
  logical::reported_team=.false.
 end type
contains
 subroutine rt_values_copy(state,values,restore)
  ! Same-owner snapshots copy numerical arrays/flags only. The live FFT/backend
  ! owner is temporarily moved aside, never cloned, destroyed or rebuilt here.
  type(s_exx_rt_state),intent(inout)::state,values
  logical,intent(in)::restore
  type(exx_k_kernel),allocatable::kernel
  call move_alloc(state%native_kernel,kernel)
  if(restore)then
   state=values
  else
   values=state
  endif
  call move_alloc(kernel,state%native_kernel)
 end subroutine
 subroutine rt_taylor_stage(state,stage,full_exact,status)
  use exx_ace,only:exx_ace_clear
  implicit none
  type(s_exx_rt_state),intent(inout)::state
  integer,intent(in)::stage
  logical,intent(in)::full_exact
  integer,intent(out)::status
  status=1
  select case(stage)
  case(0)
   ! Validate the source before changing trial flags or previous accepted data.
   if(full_exact.and..not.allocated(state%full_source))return
   state%taylor_active=.true.;state%taylor_midpoint=.false.
   if(full_exact)then
    state%initial_source=state%full_source
   else
    state%initial_ace=state%ace
   endif
  case(1)
   state%taylor_midpoint=.true.
  case(2)
   state%taylor_active=.false.;state%taylor_midpoint=.false.
   call exx_ace_clear(state%initial_ace)
   call exx_ace_clear(state%midpoint_ace)
   if(allocated(state%initial_source))deallocate(state%initial_source)
  case default
   return
  end select
  status=0
 end subroutine
 subroutine rt_ptcn_begin(state,delta_t,learned,status)
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  type(s_exx_rt_state),intent(inout)::state
  real(8),intent(in)::delta_t
  logical,intent(in)::learned
  integer,intent(out)::status
  status=1
  if(learned)then
   if(.not.ieee_is_finite(delta_t))return
   if(delta_t<=0d0)return
   if(.not.ieee_is_finite(state%ptcn_history_dt))return
   if(state%ptcn_history_dt<0d0)return
   if(state%ptcn_history_dt/=0d0)then
    if(abs(delta_t-state%ptcn_history_dt)>1d-14*state%ptcn_history_dt)then
     status=2;return
    endif
   endif
   if(state%ptcn_history_dt==0d0)state%ptcn_history_dt=delta_t
  endif
  state%ptcn_attempted=.false.;state%ptcn_predicted=.false.;status=0
 end subroutine
 subroutine rt_state_clear(state)
  implicit none
  type(s_exx_rt_state),intent(inout)::state
  type(s_exx_rt_state)::empty
  if(allocated(state%native_kernel))then
   call exx_k_kernel_destroy(state%native_kernel)
   deallocate(state%native_kernel)
  endif
  state=empty
 end subroutine
 subroutine rt_state_bytes(state,bytes,status)
  implicit none
  type(s_exx_rt_state),intent(in)::state
  integer(int64),intent(out)::bytes
  integer,intent(out)::status
  integer::i
  bytes=0;status=0
  ! FFT plans do not expose complete allocation sizes; do not report a false total.
  if(state%native_resources_enabled.or.allocated(state%native_kernel))then
   bytes=-1;status=2;return
  endif
  call ace_bytes(state%ace);call ace_bytes(state%scf_damped_ace)
  call ace_bytes(state%initial_ace);call ace_bytes(state%midpoint_ace)
  if(allocated(state%cached_occupation))call add(size(state%cached_occupation,kind=int64),8_int64)
  if(allocated(state%cached_action))call add(size(state%cached_action,kind=int64),16_int64)
  if(allocated(state%factor_occupation))call add(size(state%factor_occupation,kind=int64),8_int64)
  if(allocated(state%full_source))call add(size(state%full_source,kind=int64),16_int64)
  if(allocated(state%initial_source))call add(size(state%initial_source,kind=int64),16_int64)
  if(allocated(state%cached_source))call add(size(state%cached_source,kind=int64),16_int64)
  if(allocated(state%target_work))call add(size(state%target_work,kind=int64),16_int64)
  if(allocated(state%action_work))call add(size(state%action_work,kind=int64),16_int64)
  if(allocated(state%factor_history))then
   call add(size(state%factor_history,kind=int64),int(storage_size(state%factor_history)/8,int64))
   do i=1,size(state%factor_history)
    if(allocated(state%factor_history(i)%x))call add(size(state%factor_history(i)%x,kind=int64),16_int64)
   enddo
  endif
  if(allocated(state%ptcn_history))then
   call add(size(state%ptcn_history,kind=int64),int(storage_size(state%ptcn_history)/8,int64))
   do i=1,size(state%ptcn_history)
    if(allocated(state%ptcn_history(i)%x))call add(size(state%ptcn_history(i)%x,kind=int64),16_int64)
   enddo
  endif
 contains
  subroutine add(n,width)
   implicit none
   integer(int64),intent(in)::n,width
   if(status/=0)return
   if(n>huge(bytes)/width)then
    status=1;return
   endif
   if(bytes>huge(bytes)-n*width)then
    status=1;return
   endif
   bytes=bytes+n*width
  end subroutine
  subroutine ace_bytes(ace)
   implicit none
   type(s_exx_ace),intent(in)::ace
   if(allocated(ace%factors))call add(size(ace%factors,kind=int64),16_int64)
   if(allocated(ace%values))call add(size(ace%values,kind=int64),16_int64)
   if(allocated(ace%metric_factor))call add(size(ace%metric_factor,kind=int64),16_int64)
   if(allocated(ace%metric_rows))call add(size(ace%metric_rows,kind=int64),int(storage_size(0)/8,int64))
   if(allocated(ace%metric_cols))call add(size(ace%metric_cols,kind=int64),int(storage_size(0)/8,int64))
   if(allocated(ace%offset))call add(size(ace%offset,kind=int64),int(storage_size(0)/8,int64))
   if(allocated(ace%row))call add(size(ace%row,kind=int64),int(storage_size(0)/8,int64))
  end subroutine
 end subroutine
end module
