#include "config.h"
! SALMON adapter: legacy full-grid/orbital layout distributes k points.
! Gamma hybrid SCF and native mesh RT preserve three-dimensional Cartesian spatial blocks.
! Its source and ACE factors retain only local grid rows; overlaps are reduced.
module exx_native
  use rt_step_recovery,only:recovery_agree
#ifdef USE_MPI
  use exx_worker_service, only: service_enabled,service_exchange
#endif
  use exx_blas_threads, only: exx_blas_thread_control
  use exx_sparse_orbitals, only: sparse_clear
  use exx_functional, only: exchange_fraction,exchange_screening
  use exx_cufft, only: exx_cufft_create
  use exx_k_backend, only: k_backend_factory
  use exx_k_cufft, only: exx_k_cufft_create
  use iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use structures
  use occupation, only: get_occupation_input
  use plusU_global, only: PLUS_U_ON
  use exx_k_exchange, only: exx_k_kernel,exx_k_kernel_init,exx_k_kernel_apply_distributed
  use exx_spinor, only: exx_spinor_apply_distributed,exx_spinor_weight_sources
  use exx_k_orbits, only: apply_k_orbits,apply_spinor_k_orbits
  use exx_ace
  use exx_rt_state,only:s_exx_rt_state,rt_taylor_stage,rt_ptcn_begin,rt_state_clear
  use exx_factor_history,only:s_factor_history,history_accept,history_predict,history_write,history_read
  use exx_distributed_gauge, only: gauge_tiles_rotate
  use exx_orbitals, only: orbital_ace_build,orbital_ace_apply,orbital_layout,orbital_rotate,orbital_hermitian_action
  use exx_adaptive_support, only: adaptive_source_mask
  use exx_spatial
  use exx_wannier
  use exx_symmetry
  use communication, only: comm_summation,comm_alltoall,comm_get_max,comm_bcast
  use salmon_global, only: xc,yn_periodic,yn_spinorbit,yn_jm,yn_dc,yn_md,yn_opt,yn_predictor_corrector, &
    yn_symmetrized_stencil,propagator,hse_omega, &
    pbeh_coulomb_radius,theory,yn_conventional_from_dcdft,num_rgrid, &
    yn_hse_wannier,exx_mlwf_interval,exx_mlwf_maxiter,exx_mlwf_tolerance,exx_mlwf_radius,exx_mlwf_norm_fraction,exx_local_fft, &
    exx_factor_exact_interval,exx_factor_warmup_steps,exx_factor_n_teach, &
    exx_factor_mode,yn_exx_fft_orbits, &
    yn_exx_dc_mlwf,yn_hse_sr_pair_symmetry,exx_pre_scf_active,exx_ace_support,exx_pair_screening,exx_pair_tolerance, &
    hse_sr_tolerance,hse_block_rows,hse_q_reduction, &
    exx_local_backend,exx_gpu_batch_size,exx_kpoint_backend, &
    hse_fft_layout
  implicit none
  private
  public :: exx_validate_vector_potential
  public :: exx_set_phase,exx_pre_scf
  public :: exx_pre_scf_initialize,exx_pre_scf_begin,exx_pre_scf_not_ready,exx_pre_scf_count_ready
  public :: exx_bind_system,exx_unbind_system,exx_context,exx_system_energy
  public :: exx_check_localization
  public :: exx_enabled,exx_refresh,exx_add_action,exx_exchange_energy,exx_freeze
  public :: exx_pack,exx_unpack
  public :: exx_taylor_stage,exx_core_exchange,exx_force_full_action
  type(s_exx_wannier),save :: wannier
  public :: exx_scf_prepare,exx_scf_residual
  public :: exx_bind_gs_system
  public :: exx_set_scf_flags,exx_adaptive_is_ready,exx_is_frozen,exx_support_has_changed
  public::exx_factor_begin,exx_factor_endpoint,exx_factor_checkpoint,exx_factor_restart_pending
  public::exx_factor_mode_metadata
  public::exx_ptcn_begin,exx_ptcn_predict_inner,exx_ptcn_accept
  ! Compatibility aliases refer only to the default legacy owner.
  type(s_exx_rt_state),target,save :: rt_owned
  logical,pointer,save,public :: exx_adaptive_ready=>rt_owned%exx_adaptive_ready
  logical,pointer,save,public :: exx_support_changed=>rt_owned%exx_support_changed
  logical,pointer,save :: exx_force_full_action=>rt_owned%exx_force_full_action
  real(8),pointer,save :: exx_exchange_energy=>rt_owned%exx_exchange_energy
  logical,pointer,save :: exx_freeze=>rt_owned%exx_freeze
  ! Rank-specific loops keep IEEE inquiries scalar on Fujitsu compilers.
  private :: salmon_all_finite,finite_real_3d
  interface salmon_all_finite
    module procedure finite_real_3d
  end interface
public::exx_recovery_status
contains
  integer(c_int) function sp_fast_hse_symmetry_history_protocol() bind(C)
    use iso_c_binding,only:c_int
    sp_fast_hse_symmetry_history_protocol=1_c_int
  end function
  integer function exx_recovery_status(system) result(code)
    type(s_dft_system),intent(in)::system
    type(s_exx_rt_state),pointer::owner
    owner=>resolve_rt_context(system=system)
    code=0
    if(owner%recover_errors)code=owner%step_error
  end function

  ! Value copies isolate HSE decisions; uninitialized legacy systems retain fallback inputs.
  subroutine hse_occupation_input(system,nelec,temperature,nstate)
    type(s_dft_system),intent(in)::system
    integer,intent(out)::nelec,nstate
    real(8),intent(out)::temperature
    type(s_occupation_config)::value
    call get_occupation_input(system,value)
    nelec=value%nelec
    temperature=value%temperature
    nstate=value%nstate
  end subroutine hse_occupation_input
  function exx_context(system) result(owner)
    type(s_dft_system),intent(in)::system
    type(s_exx_rt_state),pointer::owner
    owner=>resolve_rt_context(system=system)
  end function
  subroutine exx_pre_scf_initialize(system,threshold,steps,status)
    type(s_dft_system),intent(inout) :: system
    real(8),intent(in) :: threshold
    integer,intent(in) :: steps
    integer,intent(out) :: status
    status=1
    if(.not.ieee_is_finite(threshold))return
    if(threshold<0d0.or.steps<1)return
    call exx_set_phase(system,threshold>0d0)
    system%pre_scf=s_pre_scf_state()
    system%pre_scf%threshold=threshold;system%pre_scf%required_steps=steps
    system%pre_scf%initialized=.true.
    status=0
  end subroutine
  subroutine exx_pre_scf_begin(system,completed_iterations)
    type(s_dft_system),intent(inout) :: system
    integer,intent(in) :: completed_iterations
    if(.not.system%pre_scf%initialized)error stop 'Uninitialized pre-SCF state'
    if(completed_iterations<0.or.completed_iterations==huge(0))error stop 'Invalid SCF iteration'
    system%pre_scf%count=0
    system%pre_scf%first_iteration=completed_iterations+1
  end subroutine
  integer function exx_pre_scf_not_ready(system,iteration,residual) result(not_ready)
    type(s_dft_system),intent(in) :: system
    integer,intent(in) :: iteration
    real(8),intent(in) :: residual
    not_ready=0
    if(iteration==system%pre_scf%first_iteration.or..not.(residual<system%pre_scf%threshold))not_ready=1
  end function
  logical function exx_pre_scf_count_ready(system,not_ready_total) result(ready)
    type(s_dft_system),intent(inout) :: system
    integer,intent(in) :: not_ready_total
    if(not_ready_total==0)then
      system%pre_scf%count=system%pre_scf%count+1
    else
      system%pre_scf%count=0
    endif
    ready=system%pre_scf%count>=system%pre_scf%required_steps
  end function
  subroutine exx_set_phase(system,pre_scf)
    type(s_dft_system),intent(inout) :: system
    logical,intent(in) :: pre_scf
    logical :: binding_enabled,orbital_profile
    ! Phase transitions invalidate only this explicitly bound owner's operator/history.
    ! Unbound systems can select a phase before binding; they must not reuse legacy caches.
    if(exx_pre_scf(system).neqv.pre_scf)then
      if(associated(system%salmonpy_hse_owner))then
        select type(owner=>system%salmonpy_hse_owner)
        type is(s_exx_rt_state)
          binding_enabled=owner%native_resources_enabled
          orbital_profile=owner%orbital_scf
          call rt_state_clear(owner)
          owner%native_resources_enabled=binding_enabled
          owner%orbital_scf=orbital_profile
        class default
          error stop 'Invalid HSE owner for phase transition'
        end select
      endif
    endif
    ! The caller must refresh exchange after returning to production.
    system%xc_pre_scf_active=pre_scf
    system%xc_phase_owned=.true.
  end subroutine
  logical function exx_pre_scf(system) result(active)
    type(s_dft_system),intent(in) :: system
    if(system%xc_phase_owned)then
      active=system%xc_pre_scf_active
    else
      ! Compatibility boundary for unextracted GS/DC phase state machines.
      active=exx_pre_scf_active
    endif
  end function
  function system_functional(system) result(name)
    type(s_dft_system),intent(in) :: system
    character(256) :: name
    if(system%xc_input_ready)then
      name=system%xc_input%xcname
    else
      name=xc
    endif
  end function
  real(8) function system_coulomb_radius(system) result(value)
    type(s_dft_system),intent(in) :: system
    if(system%xc_input_ready)then
      value=system%xc_input%pbeh_coulomb_radius
    else
      value=pbeh_coulomb_radius
    endif
  end function
  real(8) function system_exchange_fraction(system) result(value)
    type(s_dft_system),intent(in) :: system
    if(system%xc_input_ready)then
      value=system%xc_input%exchange_fraction
    else
      value=exchange_fraction()
    endif
  end function
  real(8) function system_exchange_screening(system) result(value)
    type(s_dft_system),intent(in) :: system
    if(system%xc_input_ready)then
      value=0d0
      if(system%xc_input%screened_exchange)value=system%xc_input%hse_omega
    else
      value=exchange_screening()
    endif
  end function
  real(8) function exx_system_energy(system) result(value)
    type(s_dft_system),intent(in)::system
    type(s_exx_rt_state),pointer::owner
    owner=>exx_context(system)
    value=owner%exx_exchange_energy
  end function
  subroutine exx_bind_system(system,owner,status)
    type(s_dft_system),intent(inout)::system
    type(s_exx_rt_state),target,intent(inout)::owner
    integer,intent(out)::status
    status=1
    if(associated(system%salmonpy_hse_owner).or.owner%native_resources_enabled)return
    owner%native_resources_enabled=.true.
    system%salmonpy_hse_owner=>owner
    status=0
  end subroutine
  subroutine exx_bind_gs_system(system,info,owner,status,allow_symmetry)
    type(s_dft_system),intent(inout)::system
    type(s_parallel_info),intent(in)::info
    type(s_exx_rt_state),target,intent(inout)::owner
    integer,intent(out)::status
    logical,optional,intent(in)::allow_symmetry
    ! 2 means the GS entry retains the legacy owner; 1 is a binding conflict.
    status=2
    if(theory/='dft'.or.yn_opt/='n')return
    ! Native retained GS opts into owned symmetry; CLI legacy policy stays unchanged.
    if(system%kpoint_input%use_symmetry)then
      if(.not.present(allow_symmetry))return
      if(.not.allow_symmetry)return
    endif
    if(info%isize_r>1)then
      if(.not.present(allow_symmetry))return
      if(.not.allow_symmetry)return
    endif
    if(.not.native_owner_supported(system,info))return
    call exx_bind_system(system,owner,status)
    if(status==0)owner%orbital_scf=product(system%kmesh)==1.and. &
      (info%isize_o>1.or.yn_hse_wannier=='y')
  end subroutine
  subroutine exx_unbind_system(system)
    type(s_dft_system),intent(inout)::system
    nullify(system%salmonpy_hse_owner)
  end subroutine
  logical function native_owner_supported(system,info) result(ok)
    type(s_dft_system),intent(in)::system
    type(s_parallel_info),intent(in)::info
    integer :: nelec,nstate,j
    real(8) :: temperature,offdiag(3,3)
    logical :: orbital_gs,kspace_gs
    call hse_occupation_input(system,nelec,temperature,nstate)
    orbital_gs=.false.
    offdiag=system%primitive_a
    do j=1,3
      offdiag(j,j)=0d0
    enddo
    if(theory=='dft'.and.(info%isize_o>1.or.info%isize_r>1.or.yn_hse_wannier=='y'))then
      offdiag=system%primitive_a
      do j=1,3
        offdiag(j,j)=0d0
      enddo
      orbital_gs=info%isize_k==1.and.info%isize_o<=system%no.and. &
        all(system%kmesh==1).and..not.system%kpoint_input%use_symmetry.and. &
        maxval(abs(offdiag))<=1d-12.and.maxval(abs(system%vec_k))<=1d-12.and. &
        exx_local_backend=='cpu'.and.exx_pair_screening=='off'
    endif
    kspace_gs=theory=='dft'.and.product(system%kmesh)>1.and.info%isize_o==1.and. &
      .not.system%kpoint_input%use_symmetry.and.yn_hse_wannier=='n'.and. &
      maxval(abs(offdiag))<=1d-12
    ok=system_functional(system)=='hse06'.and.yn_periodic=='y'.and.yn_spinorbit=='n'.and.system%nspin==1.and. &
      yn_dc=='n'.and.(yn_hse_wannier=='n'.or.orbital_gs).and.temperature<0d0.and. &
      (info%isize_r==1.or.orbital_gs.or.kspace_gs).and.(info%isize_o==1.or.orbital_gs).and. &
      info%isize_k>=1.and.info%numm==1.and. &
      exx_mlwf_radius==0d0.and.exx_mlwf_norm_fraction==0d0.and.exx_kpoint_backend=='cpu'
    if(.not.ok)return
    if(system%kpoint_input%use_symmetry)then
      ok=.false.
      ! Explicit-owner symmetry scope: one/two k ranks, CPU full-q scalar exchange.
      ! Histories belong to local representatives and use the same full-q action.
      if(theory/='dft')then
        if(exx_factor_mode/='off'.and.exx_factor_mode/='strict'.and.exx_factor_mode/='learned')return
        if(propagator/='hse_taylor4'.and.propagator/='hse_taylor4_full'.and.propagator/='hse_ptcn')return
        if(propagator=='hse_taylor4_full'.and.exx_factor_mode/='off')return
      endif
      if(info%isize_k>2.or.info%isize_k>system%nk.or.any(hse_q_reduction/=1).or.any(system%kmesh<1))return
      if(.not.allocated(system%kpoint_input%symmat_a))return
      if(.not.allocated(system%kpoint_input%symmat_b))return
      if(any(shape(system%kpoint_input%symmat_a)/=shape(system%kpoint_input%symmat_b)))return
      if(size(system%kpoint_input%symmat_a,1)/=3.or.size(system%kpoint_input%symmat_a,2)/=4)return
      if(size(system%kpoint_input%symmat_a,3)<1)return
      if(any(lbound(system%kpoint_input%symmat_a)/=1).or.any(lbound(system%kpoint_input%symmat_b)/=1))return
      ok=.true.
    endif
  end function
  subroutine exx_validate_vector_potential(system,ac,status)
    type(s_dft_system),intent(in)::system
    real(8),intent(in)::ac(3)
    integer,intent(out)::status
    type(s_rgrid)::grid_view
    real(8)::cartesian(3,3)
    integer::field_status,i,j,k
    integer(int64)::points
    status=1
    do i=1,3
      if(.not.ieee_is_finite(ac(i)))return
    enddo
    ! No fallback to the legacy process-global HSE owner or sym_sub settings.
    if(.not.associated(system%salmonpy_hse_owner))return
    select type(owner=>system%salmonpy_hse_owner)
    type is(s_exx_rt_state)
      if(.not.owner%native_resources_enabled)return
      if(.not.system%kpoint_input%use_symmetry)then
        if(owner%native_symmetry%nsym/=0.or.owner%native_symmetry%nfull/=0.or.owner%native_symmetry%nrep/=0)return
        status=0;return
      endif
      call symmetry_check_field(owner%native_symmetry,ac,field_status)
      if(field_status==1)return
      if(owner%native_symmetry%nrep/=system%nk)return
      grid_view%num=owner%native_symmetry%n
      if(.not.owned_symmetry_grid_matches(system,grid_view))return
      do j=1,3
        if(.not.ieee_is_finite(owner%native_symmetry%h(j)))return
        do i=1,3
          if(.not.ieee_is_finite(owner%native_symmetry%cell(i,j)))return
          if(.not.ieee_is_finite(owner%native_symmetry%inverse_cell(i,j)))return
        enddo
      enddo
      if(any(owner%native_symmetry%cell/=system%primitive_a).or.any(owner%native_symmetry%h/=system%hgs))return
      points=1_int64
      do i=1,3
        if(system%kmesh(i)<1)return
        if(int(system%kmesh(i),int64)>int(owner%native_symmetry%nfull,int64)/points)return
        points=points*int(system%kmesh(i),int64)
      enddo
      if(points/=int(owner%native_symmetry%nfull,int64))return
      if(.not.allocated(system%kpoint_input%symmat_a).or..not.allocated(system%kpoint_input%symmat_b))return
      if(any(shape(system%kpoint_input%symmat_a)/=[3,4,owner%native_symmetry%nsym]))return
      if(any(shape(system%kpoint_input%symmat_b)/=[3,4,owner%native_symmetry%nsym]))return
      if(any(lbound(system%kpoint_input%symmat_a)/=1).or.any(lbound(system%kpoint_input%symmat_b)/=1))return
      if(.not.allocated(owner%native_symmetry%rotation).or..not.allocated(owner%native_symmetry%translation))return
      if(any(shape(owner%native_symmetry%rotation)/=[3,3,owner%native_symmetry%nsym]))return
      if(any(shape(owner%native_symmetry%translation)/=[3,owner%native_symmetry%nsym]))return
      if(any(lbound(owner%native_symmetry%rotation)/=1).or.any(lbound(owner%native_symmetry%translation)/=1))return
      do k=1,owner%native_symmetry%nsym
        do j=1,4;do i=1,3
          if(.not.ieee_is_finite(system%kpoint_input%symmat_a(i,j,k)))return
          if(.not.ieee_is_finite(system%kpoint_input%symmat_b(i,j,k)))return
        enddo;enddo
        do j=1,3;do i=1,3
          if(.not.ieee_is_finite(owner%native_symmetry%rotation(i,j,k)))return
        enddo;enddo
        do i=1,3
          if(.not.ieee_is_finite(owner%native_symmetry%translation(i,k)))return
        enddo
        ! The immutable fractional A and reciprocal B tables are computed by
        ! different matrix products. Match the same tolerance used at map init;
        ! exact A==B rejects valid sheared cells at roundoff (observed 2.9e-18).
        if(maxval(abs(system%kpoint_input%symmat_a(:,1:3,k)-owner%native_symmetry%rotation(:,:,k)))>1d-10)return
        if(any(system%kpoint_input%symmat_b(:,1:3,k)/=owner%native_symmetry%rotation(:,:,k)))return
        if(any(system%kpoint_input%symmat_a(:,4,k)/=owner%native_symmetry%translation(:,k)))return
        cartesian=matmul(system%primitive_a, &
          matmul(owner%native_symmetry%rotation(:,:,k),owner%native_symmetry%inverse_cell))
        do j=1,3;do i=1,3
          if(.not.ieee_is_finite(cartesian(i,j)))return
        enddo;enddo
        if(maxval(abs(cartesian-owner%native_symmetry%cart_rotation(:,:,k)))>1d-10)return
      enddo
      status=field_status
    class default
      return
    end select
  end subroutine exx_validate_vector_potential
  logical function owned_symmetry_grid_matches(system,mg) result(ok)
    type(s_dft_system),intent(in)::system
    type(s_rgrid),intent(in)::mg
    integer(int64)::points
    real(8)::scale,ratio
    integer::i,j
    ok=.false.
    if(system%ngrid<1.or.any(mg%num<1))return
    points=1_int64
    do i=1,3
      if(int(mg%num(i),int64)>int(system%ngrid,int64)/points)return
      points=points*int(mg%num(i),int64)
    enddo
    if(points/=int(system%ngrid,int64))return
    do i=1,3
      if(.not.ieee_is_finite(system%hgs(i)))return
      if(system%hgs(i)<=0d0)return
      do j=1,3
        if(.not.ieee_is_finite(system%primitive_a(j,i)))return
      enddo
      scale=maxval(abs(system%primitive_a(:,i)))
      if(scale<=0d0)return
      ratio=(scale/system%hgs(i))*sqrt(sum((system%primitive_a(:,i)/scale)**2))
      if(.not.ieee_is_finite(ratio))return
      if(abs(ratio-dble(mg%num(i)))>1d-4)return
    enddo
    ok=.true.
  end function owned_symmetry_grid_matches
  function resolve_rt_context(rt_context,system) result(owner)
    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::owner
    type(s_dft_system),optional,intent(in)::system
    if(present(system))then
      if(associated(system%salmonpy_hse_owner))then
        select type(bound=>system%salmonpy_hse_owner)
        type is(s_exx_rt_state)
          owner=>bound
        class default
          error stop "Invalid HSE owner type"
        end select
        if(present(rt_context))then
          if(.not.associated(owner,rt_context))error stop "Conflicting HSE owners"
        endif
        return
      endif
    endif
    if(present(rt_context))then
      owner=>rt_context
    else
      owner=>rt_owned
    endif
  end function resolve_rt_context
  subroutine exx_set_scf_flags(system,adaptive_ready,freeze,support_changed)
    type(s_dft_system),intent(in)::system
    logical,optional,intent(in)::adaptive_ready,freeze,support_changed
    type(s_exx_rt_state),pointer::owner
    owner=>resolve_rt_context(system=system)
    if(present(adaptive_ready))owner%exx_adaptive_ready=adaptive_ready
    if(present(freeze))owner%exx_freeze=freeze
    if(present(support_changed))owner%exx_support_changed=support_changed
  end subroutine

  logical function exx_adaptive_is_ready(system)
    type(s_dft_system),intent(in)::system
    type(s_exx_rt_state),pointer::owner
    owner=>resolve_rt_context(system=system)
    exx_adaptive_is_ready=owner%exx_adaptive_ready
  end function

  logical function exx_is_frozen(system)
    type(s_dft_system),intent(in)::system
    type(s_exx_rt_state),pointer::owner
    owner=>resolve_rt_context(system=system)
    exx_is_frozen=owner%exx_freeze
  end function

  logical function exx_support_has_changed(system)
    type(s_dft_system),intent(in)::system
    type(s_exx_rt_state),pointer::owner
    owner=>resolve_rt_context(system=system)
    exx_support_has_changed=owner%exx_support_changed
  end function

  logical function thermal_scf_damping(system)
    type(s_dft_system),intent(in) :: system
    integer :: nelec,nstate
    real(8) :: temperature
    call hse_occupation_input(system,nelec,temperature,nstate)
    thermal_scf_damping=exx_enabled(system).and.theory=='dft'.and.yn_spinorbit=='y'.and. &
      yn_dc=='n'.and.temperature>=0d0.and.propagator/='hse_taylor4_full'
  end function

  subroutine exx_scf_prepare(system,mg,info,psi,alpha,rt_context)

    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    type(s_dft_system),intent(in)::system
    type(s_rgrid),intent(in)::mg
    type(s_parallel_info),intent(in)::info
    type(s_orbital),intent(in)::psi
    real(8),intent(in)::alpha
    type(s_exx_ace)::mixed
    complex(8),allocatable::targets(:,:,:)
    real(8)::residual
    integer::ierr,total_error
    exx_state_owner=>resolve_rt_context(rt_context,system)
    if(.not.thermal_scf_damping(system))return
    allocate(targets(2*product(mg%num),system%no,info%numk))
    call exx_pack(psi,mg,info,targets)
    call exx_ace_scf_mix(exx_state_owner%scf_damped_ace,exx_state_owner%ace,targets,alpha,mixed,residual,ierr)
    call comm_summation(ierr,total_error,info%icomm_k)
    if(total_error/=0)error stop 'Thermal HSE SCF operator damping failed'
    exx_state_owner%ace=mixed;exx_state_owner%scf_damped_ace=mixed;exx_state_owner%scf_operator_damped=.true.
  end subroutine

  subroutine exx_scf_residual(system,mg,info,psi,residual,rt_context)

    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    type(s_dft_system),intent(in)::system
    type(s_rgrid),intent(in)::mg
    type(s_parallel_info),intent(in)::info
    type(s_orbital),intent(in)::psi
    real(8),intent(out)::residual
    complex(8),allocatable::targets(:,:,:)
    real(8)::local,maximum(1)
    integer::ierr,total_error
    exx_state_owner=>resolve_rt_context(rt_context,system)
    residual=0d0
    if(.not.thermal_scf_damping(system))return
    if(.not.exx_ace_ready(exx_state_owner%scf_damped_ace))then
      residual=huge(1d0);return
    endif
    if(exx_state_owner%scf_operator_damped)error stop 'HSE SCF residual requires restored physical exchange'
    allocate(targets(2*product(mg%num),system%no,info%numk))
    call exx_pack(psi,mg,info,targets)
    call exx_ace_action_residual(exx_state_owner%ace,exx_state_owner%scf_damped_ace,targets,local,ierr)
    call comm_summation(ierr,total_error,info%icomm_k)
    if(total_error/=0)error stop 'Thermal HSE SCF residual failed'
    call comm_get_max([local],maximum,1,info%icomm_k)
    residual=maximum(1)
  end subroutine



  ! Experimental canonical factor history: Gamma/grid or full k-only MPI.
  logical function native_integer_occupation(system) result(valid)
    implicit none
    type(s_dft_system),intent(in) :: system
    integer :: occupied
    integer :: nelec,nstate
    real(8) :: temperature
    call hse_occupation_input(system,nelec,temperature,nstate)
    occupied=nelec/2
    valid=.false.
    if(yn_spinorbit=='y')then
      if(temperature>=0d0.or.system%no/=nelec.or.nelec<1)return
      if(maxval(abs(system%rocc-1d0))>1d-12)return
      valid=.true.;return
    endif
    if(temperature>=0d0)return
    if(mod(nelec,2)/=0.or.occupied<1.or.occupied>system%no)return
    if(maxval(abs(system%rocc(:occupied,:,:)-2d0))>1d-12)return
    if(system%no>occupied)then
      if(yn_dc/='n'.or.product(system%kmesh)==1)return
      if(maxval(abs(system%rocc(occupied+1:,:,:)))>1d-12)return
    endif
    valid=.true.
  end function native_integer_occupation

  ! Thermal spinor histories retain all target states; only their source weights
  ! are fractional. The scalar integer-prefix contract remains separate.
  logical function thermal_factor_supported(system,info,rt_context) result(valid)

    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    type(s_dft_system),intent(in)::system
    type(s_parallel_info),intent(in)::info
    integer :: nelec,nstate
    real(8) :: temperature
    call hse_occupation_input(system,nelec,temperature,nstate)
    exx_state_owner=>resolve_rt_context(rt_context,system)
    valid=yn_spinorbit=='y'.and.temperature>=0d0.and.exx_factor_mode=='learned'.and. &
      product(system%kmesh)>1.and.exx_kpoint_backend=='cpu'.and. &
      info%isize_r==1.and.info%isize_o==1.and.info%numm==1.and. &
      all(hse_q_reduction==1).and.exx_state_owner%factor_restart_step<0
  end function

  subroutine factor_check_occupation(system,info,initialize,rt_context)

    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    type(s_dft_system),intent(in)::system
    type(s_parallel_info),intent(in)::info
    logical,intent(in)::initialize
    integer::i,j,k,bad,total
    integer :: nelec,nstate
    real(8) :: temperature
    call hse_occupation_input(system,nelec,temperature,nstate)
    exx_state_owner=>resolve_rt_context(rt_context,system)
    if(.not.exx_state_owner%factor_experiment.or.yn_spinorbit/='y'.or.temperature<0d0)return
    if(.not.allocated(exx_state_owner%factor_occupation))then
      if(.not.initialize)return
      exx_state_owner%factor_occupation=system%rocc
    endif
    bad=0
    if(any(shape(exx_state_owner%factor_occupation)/=shape(system%rocc)))then
      bad=1
    else
      do k=1,size(system%rocc,3);do j=1,size(system%rocc,2);do i=1,size(system%rocc,1)
        if(.not.ieee_is_finite(system%rocc(i,j,k)))then
          bad=1
        else if(abs(system%rocc(i,j,k)-exx_state_owner%factor_occupation(i,j,k))>1d-14)then
          bad=1
        endif
      enddo;enddo;enddo
    endif
    call comm_summation(bad,total,info%icomm_rko)
    if(total/=0)error stop 'Thermal SOC factor history requires fixed occupations'
  end subroutine

  subroutine factor_configure(system,info,rt_context)
    use exx_distributed_metric,only:distributed_metric_available
    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    type(s_dft_system),intent(in)::system
    type(s_parallel_info),intent(in)::info
    integer::ik
    exx_state_owner=>resolve_rt_context(rt_context,system)
    ! Dedicated CN history never uses Taylor's frozen-expectation hook.
    if(propagator=='hse_ptcn')return
    if(.not.exx_state_owner%factor_configured)then
      exx_state_owner%factor_experiment=exx_factor_mode/='off'
      if(exx_state_owner%factor_experiment)then
        if(exx_factor_mode/='strict'.and.exx_factor_mode/='learned')error stop 'Invalid factor history mode'
        exx_state_owner%factor_learned=exx_factor_mode=='learned'

        if(system_functional(system)/='hse06'.or.yn_dc/='n'.or.yn_md/='n'.or. &
          propagator/='hse_taylor4'.or. &
          product(system%kmesh)<2.or.hse_sr_tolerance/=0d0.or. &
          exx_mlwf_radius/=0d0.or. &
          exx_mlwf_norm_fraction/=0d0.or.exx_ace_support/='occupied'.or. &
          exx_pair_screening/='off'.or. &
          (.not.native_integer_occupation(system).and..not.thermal_factor_supported(system,info, &
      & rt_context=exx_state_owner))) &
          error stop 'Factor history requires static full-support dense-k HSE and compatible occupations'
        if(info%isize_o>1.and..not.distributed_metric_available(info%icomm_ro)) &
          error stop 'Orbital-distributed factor history requires ScaLAPACK'
        if(product(system%kmesh)>1)then
          if(info%isize_r/=1.or.info%isize_o/=1.or.info%numm/=1) &
            error stop 'Factor k history requires k-only MPI'
          if(info%numk<1.or.info%ik_s<1.or.info%ik_e>system%nk.or. &
             info%ik_e-info%ik_s+1/=info%numk)error stop 'Invalid factor representative k layout'
          if(.not.system%kpoint_input%use_symmetry.and.product(system%kmesh)/=system%nk) &
            error stop 'Factor k history requires full mesh without symmetry'
          if(any(hse_q_reduction/=1))error stop 'Factor k history requires all q points'
          do ik=1,system%nk
            if(.not.ieee_is_finite(system%wtk(ik)))error stop 'Nonfinite factor k weight'
          enddo
          ! Symmetry retains zero-padded storage beyond active representatives.
          if(any(system%wtk(:system%nk)<=0d0).or. &
             abs(sum(system%wtk(:system%nk))-1d0)>1d-12)error stop 'Invalid factor k weights'
          if(info%id_rko==0)write(*,'(a,l1,2i10)') &
            'FACTOR_HISTORY symmetry/full/representatives: ',system%kpoint_input%use_symmetry,product(system%kmesh),system%nk
        endif
        if(exx_factor_exact_interval<1)error stop 'Factor interval must be positive'
        allocate(exx_state_owner%factor_history(info%numk))
        exx_state_owner%factor_history%interval=exx_factor_exact_interval
        exx_state_owner%factor_history%n_teach=exx_factor_n_teach
      endif
      exx_state_owner%factor_configured=.true.
    endif
  end subroutine

  subroutine exx_factor_begin(step,system,mg,info,psi,rt_context)
    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    integer,intent(in)::step
    type(s_dft_system),intent(in)::system
    type(s_rgrid),intent(in)::mg
    type(s_parallel_info),intent(in)::info
    type(s_orbital),intent(in)::psi
    exx_state_owner=>resolve_rt_context(rt_context,system)
    if(.not.exx_state_owner%factor_initialized)then
      call factor_configure(system,info, &
      & rt_context=exx_state_owner)
      if(exx_state_owner%factor_experiment)then
        call factor_check_occupation(system,info,.true., &
      & rt_context=exx_state_owner)
        call factor_materialize(info, &
      & rt_context=exx_state_owner)
        if(exx_state_owner%factor_restart_step>=0.and.exx_state_owner%factor_learned)then
          call exx_factor_restore(system,mg,info, &
      & rt_context=exx_state_owner)
        else
          call exx_factor_endpoint(system,info,0, &
      & rt_context=exx_state_owner)
        endif
        if(info%id_rko==0)write(*,'(a,l1)')'FACTOR_HISTORY learned: ',exx_state_owner%factor_learned
      endif
      exx_state_owner%factor_initialized=.true.
    endif
    if(.not.exx_state_owner%factor_experiment)return
    call factor_check_occupation(system,info,.false., &
      & rt_context=exx_state_owner)
    exx_state_owner%factor_step=step;exx_state_owner%factor_prediction_installed=.false.
    exx_state_owner%exx_freeze=exx_state_owner%factor_learned.and.step>exx_factor_warmup_steps.and. &
      mod(step,exx_factor_exact_interval)/=0.and.exx_state_owner%factor_all_ready
    ! At strict correction steps refresh the true starting operator before Taylor stage zero.
    if(.not.exx_state_owner%exx_freeze.and.step>exx_factor_warmup_steps.and.exx_state_owner%factor_learned)call exx_refresh(system,mg,info,psi, &
      & rt_context=exx_state_owner)
    if(info%id_rko==0)write(*,'(a,i8,l2)')'FACTOR_HISTORY step/predicted: ',step,exx_state_owner%exx_freeze
  end subroutine

  subroutine exx_factor_mode_metadata(directory,step,info,writing)
    implicit none
    character(*),intent(in)::directory
    integer,intent(in)::step
    type(s_parallel_info),intent(in)::info
    logical,intent(in)::writing
    integer::mode,saved_mode,unit,status,total,ios,version,saved_step
    logical::exists,has_history
    mode=0;status=0
    select case(exx_factor_mode)
    case('off');mode=0
    case('strict');mode=1
    case('learned');mode=2
    case default;status=1
    end select
    call comm_summation(status,total,info%icomm_rko)
    if(total/=0)error stop 'Invalid factor mode at checkpoint/restart'
    if(info%id_rko==0)then
      inquire(file=trim(directory)//'factor-mode.bin',exist=exists)
      if(writing.or.exists)then
        if(writing)then
          open(newunit=unit,file=trim(directory)//'factor-mode.bin', &
            form='unformatted',status='replace',iostat=status)
        else
          open(newunit=unit,file=trim(directory)//'factor-mode.bin', &
            form='unformatted',status='old',iostat=status)
        endif
        if(status==0)then
          if(writing)then
            write(unit,iostat=status)1,step,mode
          else
            read(unit,iostat=status)version,saved_step,saved_mode
            if(status==0)then
              if(version/=1.or.saved_step/=step.or.saved_mode/=mode)status=1
            endif
          endif
          close(unit,iostat=ios)
          if(ios/=0)status=1
        endif
      else
        inquire(file=trim(directory)//'factor-history-00000000.bin',exist=has_history)
        if(mode/=0.or.has_history)status=1
      endif
    endif
    call comm_bcast(status,info%icomm_rko)
    if(status/=0)error stop 'Factor mode checkpoint missing or mismatched'
  end subroutine

  ! Deferred restore follows initialization of the exchange backend.
  subroutine exx_factor_restart_pending(directory,step,rt_context)
    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    character(*),intent(in)::directory
    integer,intent(in)::step
    exx_state_owner=>resolve_rt_context(rt_context)
    exx_state_owner%factor_restart_directory=directory
    exx_state_owner%factor_restart_step=step
  end subroutine

  ! Persist the route selection as well as the physical SR radius. Transient FFT
  ! plans and wave/source caches are reconstructed, never serialized.
  subroutine factor_checkpoint_config(config)
    implicit none
    character(64),intent(out)::config(6)
    character(40),parameter::names(6)=[character(40):: &
      'SALMON_HSE_SR_RADIUS_AU','SALMON_HSE_SR_GROUP_FFT','SALMON_HSE_SR_WAVE_FIRST', &
      'SALMON_HSE_SR_RETURN_SUM','SALMON_HSE_SR_WAVE_BATCH','SALMON_HSE_SR_PAIR_SYMMETRY']
    integer::i,stat
    config=''
    do i=1,6
      call get_environment_variable(trim(names(i)),config(i),status=stat)
      if(stat==1)config(i)=''
      if(stat/=0.and.stat/=1)error stop 'Invalid factor checkpoint route setting'
    enddo
  end subroutine

  subroutine exx_factor_checkpoint(directory,step,system,mg,info,rt_context)
    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    character(*),intent(in)::directory
    integer,intent(in)::step
    type(s_dft_system),intent(in)::system
    type(s_rgrid),intent(in)::mg
    type(s_parallel_info),intent(in)::info
    integer::unit,status,total,ik,ios
    character(64)::checkpoint_config(6)
    character(600)::filename
    exx_state_owner=>resolve_rt_context(rt_context,system)
    if(.not.exx_state_owner%factor_learned)return
    status=0
    if(.not.exx_state_owner%factor_initialized.or.exx_state_owner%ace%packed.or..not.allocated(exx_state_owner%ace%factors))status=1
    call comm_summation(status,total,info%icomm_rko)
    if(total/=0)error stop 'Learned ACE checkpoint requires initialized dense factors'
    write(filename,'(a,a,i8.8,a)')trim(directory),'factor-history-',info%id_rko,'.bin'
    open(newunit=unit,file=trim(filename),form='unformatted',access='stream',status='replace',iostat=status)
    if(status==0)then
      write(unit,iostat=status)3,step,mg%is,mg%ie,system%no,system%nk, &
        info%isize_r,info%isize_k,info%isize_o,info%ik_s,info%numk, &
        exx_factor_exact_interval,exx_factor_warmup_steps,shape(exx_state_owner%ace%factors), &
        exx_state_owner%ace%dv,exx_state_owner%ace%condition,exx_state_owner%ace%grid_rows,exx_state_owner%factor_all_ready, &
        hse_sr_tolerance,exx_mlwf_norm_fraction,exx_mlwf_radius
      call factor_checkpoint_config(checkpoint_config)
      if(status==0)write(unit,iostat=status)info%io_s,info%numo,info%id_rko,info%id_r,info%id_o
      if(status==0)write(unit,iostat=status)checkpoint_config
      if(status==0)write(unit,iostat=status)hse_q_reduction
      if(status==0)write(unit,iostat=status)exx_state_owner%ace%factors
      do ik=1,size(exx_state_owner%factor_history)
        if(status/=0)exit
        call history_write(unit,exx_state_owner%factor_history(ik),status)
      enddo
      close(unit,iostat=ios)
      if(ios/=0)status=1
    endif
    status=merge(1,0,status/=0)
    call comm_summation(status,total,info%icomm_rko)
    if(total/=0)error stop 'Learned ACE checkpoint write failed'
  end subroutine

  subroutine exx_factor_restore(system,mg,info,rt_context)
    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    type(s_dft_system),intent(in)::system
    type(s_rgrid),intent(in)::mg
    type(s_parallel_info),intent(in)::info
    type(s_factor_history),allocatable::loaded(:)
    complex(8),allocatable::factors(:,:,:)
    integer::unit,status,total,ios,ik,version,step,header(15),expected(15),dims(3),rows,saved_q_reduction(3)
    integer::saved_layout(5),expected_layout(5)
    character(64)::saved_config(6),expected_config(6)
    real(8)::dv,condition,settings(3)
    logical::ready
    character(600)::filename
    exx_state_owner=>resolve_rt_context(rt_context,system)
    expected=[mg%is,mg%ie,system%no,system%nk,info%isize_r,info%isize_k, &
      info%isize_o,info%ik_s,info%numk,exx_factor_exact_interval,exx_factor_warmup_steps]
    write(filename,'(a,a,i8.8,a)')trim(exx_state_owner%factor_restart_directory),'factor-history-',info%id_rko,'.bin'
    open(newunit=unit,file=trim(filename),form='unformatted',access='stream',status='old',iostat=status)
    if(status==0)then
      read(unit,iostat=status)version,step,header,dims,dv,condition,rows,ready,settings
      saved_q_reduction=1
      expected_layout=[info%io_s,info%numo,info%id_rko,info%id_r,info%id_o]
      call factor_checkpoint_config(expected_config)
      if(status==0.and.version==3)then
        read(unit,iostat=status)saved_layout
        if(status==0)read(unit,iostat=status)saved_config
        if(status==0)then
          if(any(saved_layout/=expected_layout).or.any(saved_config/=expected_config))status=1
        endif
      else if(status==0.and.info%isize_o/=1)then
        status=1
      endif
      if(status==0.and.(version==2.or.version==3))read(unit,iostat=status)saved_q_reduction
      if(status==0)then
        if((version/=1.and.version/=2.and.version/=3).or.step/=exx_state_owner%factor_restart_step.or.any(header/=expected))status=1
        if(any(saved_q_reduction/=hse_q_reduction))status=1
        if(any(dims/=shape(exx_state_owner%ace%factors)))status=1
        if(.not.ieee_is_finite(dv).or..not.ieee_is_finite(condition))status=1
        if(dv/=system%hvol.or.rows/=exx_state_owner%ace%grid_rows)status=1
        if(any(settings/=[hse_sr_tolerance,exx_mlwf_norm_fraction,exx_mlwf_radius]))status=1
        do ik=1,3
          if(.not.ieee_is_finite(settings(ik)))status=1
        enddo
      endif
      if(status==0)then
        allocate(factors(dims(1),dims(2),dims(3)),loaded(info%numk))
        read(unit,iostat=status)factors
        if(status==0)then
          if(.not.finite_complex_3d_restart(factors))status=1
        endif
        do ik=1,info%numk
          if(status/=0)exit
          loaded(ik)%n_teach=exx_factor_n_teach
          call history_read(unit,loaded(ik),dims(1),dims(2),status)
          if(status/=0)exit
          if(loaded(ik)%interval/=exx_factor_exact_interval)status=1
          if(loaded(ik)%last/=step-mod(step,exx_factor_exact_interval))status=1

        enddo
      endif
      close(unit,iostat=ios)
      if(ios/=0)status=1
    endif
    status=merge(1,0,status/=0)
    call comm_summation(status,total,info%icomm_rko)
    if(total/=0)error stop 'Learned ACE restart state missing, incompatible or corrupt'
    status=count(.not.loaded%ready)
    call comm_summation(status,total,info%icomm_rko)
    if(ready.neqv.(total==0))error stop "Learned ACE restart readiness mismatch"
    exx_state_owner%ace%factors=factors;exx_state_owner%ace%dv=dv;exx_state_owner%ace%condition=condition
    exx_state_owner%factor_history=loaded;exx_state_owner%factor_all_ready=ready
    if(allocated(exx_state_owner%cached_source))deallocate(exx_state_owner%cached_source)
    if(allocated(exx_state_owner%cached_occupation))deallocate(exx_state_owner%cached_occupation)
    exx_state_owner%factor_restart_step=-1
  end subroutine

  logical function finite_complex_3d_restart(a) result(ok)
    implicit none
    complex(8),intent(in)::a(:,:,:)
    integer::i,j,k
    ok=.false.
    do k=1,size(a,3)
      do j=1,size(a,2)
        do i=1,size(a,1)
          if(.not.ieee_is_finite(real(a(i,j,k),8)))return
          if(.not.ieee_is_finite(aimag(a(i,j,k))))return
        enddo
      enddo
    enddo
    ok=.true.
  end function

  subroutine factor_materialize(info,rt_context)
    use exx_distributed_metric,only:distributed_metric_rotate
    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    type(s_parallel_info),intent(in)::info
    complex(8),allocatable::rowbuf(:),rotated(:,:)
    integer,allocatable::counts(:)
    integer::i,j,k,n,ng,first,stat
    exx_state_owner=>resolve_rt_context(rt_context)
    if(.not.exx_state_owner%ace%packed)return
    n=size(exx_state_owner%ace%offset)-1;ng=exx_state_owner%ace%grid_rows
    allocate(exx_state_owner%ace%factors(ng,n,1));exx_state_owner%ace%factors=0
    do j=1,n
      do k=exx_state_owner%ace%offset(j),exx_state_owner%ace%offset(j+1)-1
        exx_state_owner%ace%factors(exx_state_owner%ace%row(k),j,1)=exx_state_owner%ace%values(k)
      enddo
    enddo
    if(exx_state_owner%ace%metric_distributed)then
      call orbital_layout(n,info%icomm_r,info%icomm_o,counts,first,stat)
      if(stat/=0)error stop 'Factor history orbital layout failed'
      allocate(rotated(ng,n))
      call distributed_metric_rotate(exx_state_owner%ace,exx_state_owner%ace%factors(:,:,1),info%icomm_o,counts,first,rotated)
      exx_state_owner%ace%factors(:,:,1)=rotated
      deallocate(exx_state_owner%ace%metric_rows,exx_state_owner%ace%metric_cols)
      exx_state_owner%ace%metric_distributed=.false.
    else
!$omp parallel private(rowbuf,i,j)
    allocate(rowbuf(n))
!$omp do
    do i=1,ng
      rowbuf=exx_state_owner%ace%factors(i,:,1)
      do j=1,n
        exx_state_owner%ace%factors(i,j,1)=sum(rowbuf*exx_state_owner%ace%metric_factor(:,j))
      enddo
    enddo
!$omp end do
    deallocate(rowbuf)
!$omp end parallel
    endif
    deallocate(exx_state_owner%ace%offset,exx_state_owner%ace%row,exx_state_owner%ace%values,exx_state_owner%ace%metric_factor)
    exx_state_owner%ace%packed=.false.
  end subroutine

  subroutine exx_factor_endpoint(system,info,step,rt_context)
    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    type(s_dft_system),intent(in)::system
    type(s_parallel_info),intent(in)::info
    integer,intent(in)::step
    integer::stat,ik,not_ready,total_not_ready
    exx_state_owner=>resolve_rt_context(rt_context,system)
    if(.not.exx_state_owner%factor_learned)return
    if(mod(step,exx_factor_exact_interval)/=0)return
    if(exx_state_owner%exx_freeze)error stop 'Predicted endpoint cannot enter exact factor history'
    call factor_materialize(info, &
      & rt_context=exx_state_owner)
    if(.not.allocated(exx_state_owner%ace%factors))error stop 'Factor endpoint missing local factors'
    if(size(exx_state_owner%ace%factors,3)/=info%numk.or.size(exx_state_owner%factor_history)/=info%numk) &
      error stop 'Factor endpoint representative/history mismatch'
    do ik=1,info%numk
      call history_accept(exx_state_owner%factor_history(ik),exx_state_owner%ace%factors(:,:,ik),system%hvol,step,sumgrid,stat,aligngrid)
        if(stat/=0)then
          if(.not.exx_state_owner%recover_errors)error stop 'Factor k history alignment/teacher failed'
          exx_state_owner%step_error=33
        endif
      if(info%id_r==0)write(*,'(a,4i8,2es18.9,l2)')'FACTOR_HISTORY k/teacher/step/count/coeff/ready: ', &
        info%ik_s+ik-1,exx_state_owner%factor_history(ik)%teachers,step,exx_state_owner%factor_history(ik)%count, &
        exx_state_owner%factor_history(ik)%coeff,exx_state_owner%factor_history(ik)%ready
    enddo
    if(exx_state_owner%recover_errors)then
      call recovery_agree(info,33,exx_state_owner%step_error/=0,exx_state_owner%step_error)
      if(exx_state_owner%step_error/=0)return
    endif
    ! Readiness can change only after exact teachers; synchronize here, not every RT step.
    not_ready=count(.not.exx_state_owner%factor_history%ready)
    call comm_summation(not_ready,total_not_ready,info%icomm_rko)
    exx_state_owner%factor_all_ready=total_not_ready==0
  contains
    subroutine aligngrid(x,ref,status)
      use exx_distributed_alignment,only:distributed_factor_align
      use exx_distributed_metric,only:distributed_metric_available
      use exx_factor_history,only:factor_align
      implicit none
      complex(8),intent(inout)::x(:,:)
      complex(8),intent(in)::ref(:,:)
      integer,intent(out)::status
      if(distributed_metric_available(info%icomm_ro))then
        call distributed_factor_align(x,ref,info%icomm_r,info%icomm_o,info%icomm_ro,status)
      else
        call factor_align(x,ref,sumgrid,status)
      endif
    end subroutine
    subroutine sumgrid(a)
    implicit none
      complex(8),intent(inout)::a(:,:)
      complex(8)::total(size(a,1),size(a,2))
      call comm_summation(a,total,size(a),info%icomm_ro);a=total
    end subroutine
  end subroutine

  subroutine factor_predicted_expectation(system,mg,info,psi,rt_context)
    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    type(s_dft_system),intent(in)::system
    type(s_rgrid),intent(in)::mg
    type(s_parallel_info),intent(in)::info
    type(s_orbital),intent(in)::psi
    real(8)::ex
    integer::j,io,stat,ik,bad
    exx_state_owner=>resolve_rt_context(rt_context,system)
    if(.not.exx_state_owner%factor_prediction_installed)then
      bad=0
      do ik=1,info%numk
        call history_predict(exx_state_owner%factor_history(ik),exx_state_owner%factor_step-exx_state_owner%factor_history(ik)%last,exx_state_owner%ace%factors(:,:,ik),stat)
        bad=max(bad,abs(stat))
        if(stat/=0.and..not.exx_state_owner%recover_errors)error stop 'Factor k predictor failed'
      enddo
      if(exx_state_owner%recover_errors)then
        call recovery_agree(info,33,bad/=0,exx_state_owner%step_error)
        if(exx_state_owner%step_error/=0)return
      endif
      exx_state_owner%factor_prediction_installed=.true.
    endif
    if(.not.allocated(exx_state_owner%target_work))allocate(exx_state_owner%target_work(product(mg%num)*size(psi%zwf,4),info%numo,info%numk))
    if(.not.allocated(exx_state_owner%action_work))allocate(exx_state_owner%action_work(product(mg%num)*size(psi%zwf,4),info%numo,info%numk))
    call exx_pack(psi,mg,info,exx_state_owner%target_work)
    if(info%isize_r>1.or.info%isize_o>1.or.exx_state_owner%ace%packed)then
      call orbital_ace_apply(exx_state_owner%ace,exx_state_owner%target_work,exx_state_owner%action_work,info%icomm_r,info%icomm_o,stat)
    else
      call exx_ace_apply(exx_state_owner%ace,exx_state_owner%target_work,exx_state_owner%action_work,stat,thread_control=exx_blas_thread_control)
    endif
    if(exx_state_owner%recover_errors)then
      call recovery_agree(info,30,stat/=0,exx_state_owner%step_error)
      if(exx_state_owner%step_error/=0)return
    endif
    if(stat/=0)error stop 'Predicted ACE expectation failed'
    ex=0d0
    do ik=1,info%numk
    do j=1,info%numo
      io=info%io_s+j-1
      ex=ex+.5d0*system_exchange_fraction(system)*system%hvol*system%rocc(io,info%ik_s+ik-1,1)*system%wtk(info%ik_s+ik-1) &
        *real(sum(conjg(exx_state_owner%target_work(:,j,ik))*exx_state_owner%action_work(:,j,ik)),8)
    enddo
    enddo
    call comm_summation(ex,exx_state_owner%exx_exchange_energy,info%icomm_rko)
  end subroutine

  subroutine exx_ptcn_begin(dt,rt_context)
    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    real(8),intent(in)::dt
    integer::status
    exx_state_owner=>resolve_rt_context(rt_context)
    call rt_ptcn_begin(exx_state_owner,dt,exx_factor_mode=='learned',status)
    if(status==1)error stop 'Learned CN requires finite positive dt'
    if(status==2)error stop 'Learned CN requires constant dt'
  end subroutine

  subroutine ptcn_collect_exact(system,info,step,rt_context)
    use exx_distributed_alignment,only:distributed_factor_align
    use exx_distributed_metric,only:distributed_metric_available
    use exx_factor_history,only:factor_align

    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    type(s_dft_system),intent(in)::system
    type(s_parallel_info),intent(in)::info
    integer,intent(in)::step
    integer::k,status,bad
    complex(8),allocatable::aligned(:,:,:)
    exx_state_owner=>resolve_rt_context(rt_context,system)
    if(exx_state_owner%exx_freeze.or.exx_state_owner%scf_operator_damped)error stop 'CN teacher must be physical exact exchange'
    call factor_materialize(info, &
      & rt_context=exx_state_owner)
    if(.not.allocated(exx_state_owner%ace%factors))error stop 'CN teacher missing dense ACE'
    if(.not.allocated(exx_state_owner%ptcn_history))then
      allocate(exx_state_owner%ptcn_history(info%numk));exx_state_owner%ptcn_history%interval=1
      exx_state_owner%ptcn_history%n_teach=exx_factor_n_teach
    endif
    aligned=exx_state_owner%ace%factors
    bad=0
    do k=1,info%numk
      call history_accept(exx_state_owner%ptcn_history(k),aligned(:,:,k),system%hvol,step,sumgrid,status,aligngrid)
      bad=max(bad,abs(status))
      if(.not.exx_state_owner%recover_errors)then
        call comm_get_max(bad,info%icomm_rko)
        if(bad/=0)error stop 'CN exact teacher alignment failed'
      endif
    enddo
    if(exx_state_owner%recover_errors)then
      call recovery_agree(info,33,bad/=0,exx_state_owner%step_error)
      if(exx_state_owner%step_error/=0)return
    endif
    exx_state_owner%ace%factors=aligned
  contains
    subroutine sumgrid(a)
      complex(8),intent(inout)::a(:,:)
      complex(8)::total(size(a,1),size(a,2))
      call comm_summation(a,total,size(a),info%icomm_ro);a=total
    end subroutine
    subroutine aligngrid(x,ref,status)
      complex(8),intent(inout)::x(:,:)
      complex(8),intent(in)::ref(:,:)
      integer,intent(out)::status
      if(distributed_metric_available(info%icomm_ro))then
        call distributed_factor_align(x,ref,info%icomm_r,info%icomm_o,info%icomm_ro,status)
      else
        call factor_align(x,ref,sumgrid,status)
      endif
    end subroutine
  end subroutine

  subroutine exx_ptcn_predict_inner(system,info,rt_context)

    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    type(s_dft_system),intent(in)::system
    type(s_parallel_info),intent(in)::info
    complex(8),allocatable::prediction(:,:,:)
    integer::k,status,bad
    exx_state_owner=>resolve_rt_context(rt_context,system)
    if(exx_factor_mode/='learned'.or.exx_state_owner%ptcn_attempted)return
    exx_state_owner%ptcn_attempted=.true.
    ! The first exact RHS action already installed H[u_n]; record only step zero.
    if(.not.allocated(exx_state_owner%ptcn_history))call ptcn_collect_exact(system,info,0, &
      & rt_context=exx_state_owner)
    if(exx_state_owner%recover_errors.and.exx_state_owner%step_error/=0)return
    bad=count(.not.exx_state_owner%ptcn_history%ready);call comm_get_max(bad,info%icomm_rko)
    if(bad/=0)return
    prediction=exx_state_owner%ace%factors
    bad=0
    do k=1,info%numk
      call history_predict(exx_state_owner%ptcn_history(k),1,prediction(:,:,k),status)
      bad=max(bad,abs(status))
    enddo
    call comm_get_max(bad,info%icomm_rko)
    if(bad/=0)return ! Collective exact-ACE fallback; never install a partial prediction.
    exx_state_owner%ace%factors=prediction;exx_state_owner%ptcn_predicted=.true.
    ! Invalidate the physical-source cache even if endpoint psi equals input psi.
    ! The next fresh action MUST reconstruct physical exchange, not return a cache hit.
    exx_state_owner%scf_operator_damped=.true.
  end subroutine

  subroutine exx_ptcn_accept(system,info,rt_context)

    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    type(s_dft_system),intent(in)::system
    type(s_parallel_info),intent(in)::info
    exx_state_owner=>resolve_rt_context(rt_context,system)
    if(exx_factor_mode/='learned')return
    ! Caller invokes only after fresh residual, Gram and electron gates succeeded.
    call ptcn_collect_exact(system,info,exx_state_owner%ptcn_accepted_step+1, &
      & rt_context=exx_state_owner)
    exx_state_owner%ptcn_accepted_step=exx_state_owner%ptcn_accepted_step+1
  end subroutine

  subroutine exx_check_localization(rt_context,system)
    implicit none
    type(s_dft_system),optional,intent(in)::system
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    exx_state_owner=>resolve_rt_context(rt_context,system)
    if(dc_canonical())return
    ! A density criterion cannot certify a gauge-dependent truncated operator.
    if(exx_mlwf_norm_fraction>0d0.and.exx_mlwf_radius==0d0.and..not.exx_state_owner%adaptive_active) &
      error stop 'Adaptive EXX: localized support not established; SCF result rejected'
    if(exx_mlwf_radius>0d0.and..not.exx_state_owner%finite_support_localized) &
      error stop 'EXX MLWF finite support: localization not converged; SCF result rejected'
  end subroutine exx_check_localization

  subroutine exx_taylor_stage(stage,rt_context)
    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    integer,intent(in)::stage
    integer::status
    exx_state_owner=>resolve_rt_context(rt_context)
    call rt_taylor_stage(exx_state_owner,stage,propagator=='hse_taylor4_full',status)
    if(status/=0)error stop 'Invalid HSE Taylor stage or missing initial source'
  end subroutine


  subroutine warn_fixed_radius(max_loss,rt_context)
    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    real(8),intent(in) :: max_loss
    real(8) :: target
    exx_state_owner=>resolve_rt_context(rt_context)
    if(exx_mlwf_radius<=0d0.or.exx_state_owner%radius_warning_reported)return
    target=merge(exx_mlwf_norm_fraction,.999d0,exx_mlwf_norm_fraction>0d0)
    if(1d0-max_loss>=target-64d0*epsilon(1d0))return
    write(*,'(a,3es18.9)') 'WARNING EXX fixed radius retains less than target; R(bohr)/min retained/target: ', &
      exx_mlwf_radius,1d0-max_loss,target
    exx_state_owner%radius_warning_reported=.true.
  end subroutine warn_fixed_radius

  logical function dc_gamma_pairs_requested(system)
    implicit none
    type(s_dft_system),intent(in) :: system
    dc_gamma_pairs_requested=dc_canonical().and.product(system%kmesh)==1.and. &
      yn_hse_sr_pair_symmetry=='y'.and.hse_sr_tolerance==0d0
  end function

  logical function dc_canonical()
    implicit none
    dc_canonical=yn_dc=='y'.and.yn_exx_dc_mlwf=='n'
  end function

  logical function exx_enabled(system)
    use exx_functional, only: is_hybrid
    implicit none
    type(s_dft_system),optional,intent(in) :: system
    if(present(system))then
      exx_enabled=is_hybrid(system_functional(system)).and..not.exx_pre_scf(system)
    else
      exx_enabled=is_hybrid(xc).and..not.exx_pre_scf_active
    endif
  end function

  subroutine exx_pack(psi,mg,info,a)
    implicit none
    type(s_orbital),intent(in) :: psi
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    complex(8),intent(out) :: a(:,:,:)
    integer :: ik,io,ix,iy,iz,g,is,ng
    ng=product(mg%num)
!$omp parallel do collapse(3) default(shared) private(ik,io,is,ix,iy,iz,g) schedule(static)
    do ik=info%ik_s,info%ik_e;do io=info%io_s,info%io_e;do is=1,size(psi%zwf,4)
      do iz=mg%is(3),mg%ie(3);do iy=mg%is(2),mg%ie(2);do ix=mg%is(1),mg%ie(1)
        g=1+(ix-mg%is(1))+mg%num(1)*((iy-mg%is(2))+mg%num(2)*(iz-mg%is(3)))+(is-1)*ng
        a(g,io-info%io_s+1,ik-info%ik_s+1)=psi%zwf(ix,iy,iz,is,io,ik,1)
      enddo;enddo;enddo
    enddo;enddo;enddo
!$omp end parallel do
  end subroutine

  subroutine exx_unpack(a,psi,mg,info)
    implicit none
    complex(8),intent(in) :: a(:,:,:)
    type(s_orbital),intent(inout) :: psi
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    integer :: ik,io,ix,iy,iz,g,is,ng
    ng=product(mg%num)
!$omp parallel do collapse(3) default(shared) private(ik,io,is,ix,iy,iz,g) schedule(static)
    do ik=info%ik_s,info%ik_e;do io=info%io_s,info%io_e;do is=1,size(psi%zwf,4)
      do iz=mg%is(3),mg%ie(3);do iy=mg%is(2),mg%ie(2);do ix=mg%is(1),mg%ie(1)
        g=1+(ix-mg%is(1))+mg%num(1)*((iy-mg%is(2))+mg%num(2)*(iz-mg%is(3)))+(is-1)*ng
        psi%zwf(ix,iy,iz,is,io,ik,1)=a(g,io-info%io_s+1,ik-info%ik_s+1)
      enddo;enddo;enddo
    enddo;enddo;enddo
!$omp end parallel do
    psi%update_zwf_overlap=.false.
  end subroutine

  recursive subroutine exx_refresh(system,mg,info,psi,rt_context,rt_status)
    use exx_functional, only: is_global_hybrid
    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    integer,optional,intent(out)::rt_status
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: psi
    complex(8),allocatable :: w(:,:,:),local(:,:,:),weighted(:,:,:)
    real(8) :: ex,electrons
    logical :: thermal_spinor
    integer :: ierr,total_error,ng,nk,no,n,mesh,j,first_full,count_full
    procedure(k_backend_factory),pointer :: create_backend
    integer :: nelec,nstate
    real(8) :: temperature
    call hse_occupation_input(system,nelec,temperature,nstate)
    exx_state_owner=>resolve_rt_context(rt_context,system)
    if(exx_state_owner%recover_errors.and.exx_state_owner%step_error/=0)return
    if(present(rt_status))rt_status=0
    if(.not.associated(exx_state_owner,rt_owned))then
      if(.not.exx_state_owner%native_resources_enabled)then
        if(present(rt_status))then
          rt_status=1;return
        endif
        error stop 'Explicit HSE RT context supports serial dense ACE/lifecycle only'
      endif
    endif
    if(exx_state_owner%native_resources_enabled)then
      if(.not.native_owner_supported(system,info))then
        if(present(rt_status))then
          rt_status=2;return
        endif
        error stop 'Independent HSE owner requires serial full-support non-spinor configuration'
      endif
    endif
    if(exx_state_owner%native_resources_enabled)then
      if(system%kpoint_input%use_symmetry)then
        if(.not.owned_symmetry_grid_matches(system,mg))then
          if(present(rt_status))then
            rt_status=2;return
          endif
          error stop 'Independent HSE symmetry requires matching owned full grid'
        endif
      endif
    endif
    if(.not.exx_enabled(system))return
    if(theory=='tddft_response'.or.theory=='tddft_pulse')call factor_configure(system,info, &
      & rt_context=exx_state_owner)
    call factor_check_occupation(system,info,.false., &
      & rt_context=exx_state_owner)
    if(exx_state_owner%exx_freeze)then
      if(exx_state_owner%factor_experiment)call factor_predicted_expectation(system,mg,info,psi, &
      & rt_context=exx_state_owner)
      return
    endif
    if(yn_periodic/='y'.or.(system%nspin/=1.and.yn_spinorbit/='y').or..not.allocated(psi%zwf)) &
      error stop 'HSE06: periodic complex unpolarized orbitals required'
    if(yn_md=='y'.and.theory/='dft_md')then
      if((.not.is_global_hybrid(system_functional(system))).or. &
         (theory/='tddft_response'.and.theory/='tddft_pulse').or.yn_conventional_from_dcdft/='y') &
        error stop 'Hybrid: unsupported real-time ionic extension'
    endif
    if(yn_jm/='n'.or.yn_symmetrized_stencil=='y') &
      error stop 'HSE06: unsupported Hamiltonian/ionic extension'
    if(PLUS_U_ON)error stop 'HSE06: DFT+U combination unsupported'
    if(allocated(system%Ac_micro%v))error stop 'HSE06: microscopic vector potential unsupported'
    ! Select occupation semantics before the integer-prefix native route.
    ! Nonnegative temperature retains the existing occupation-weighted source.
    if(temperature>=0d0.and.yn_spinorbit/="y")then
      if(info%isize_r>1.or.info%isize_o>1.or.dc_gamma_pairs_requested(system).or. &
         (exx_mlwf_norm_fraction>0d0.and.exx_mlwf_radius==0d0))then
        call refresh_spatial(system,mg,info,psi, &
      & rt_context=exx_state_owner)
      else
        call refresh_wannier(system,mg,info,psi, &
      & rt_context=exx_state_owner)
      endif
      return
    endif
    if(theory=='dft'.and.system%nk>1.and.info%isize_r>1.and. &
       exx_state_owner%native_resources_enabled)then
      call refresh_kspace(system,mg,info,psi,exx_state_owner)
      return
    endif
    if(theory=='dft'.and.exx_state_owner%orbital_scf)then
      call refresh_spatial(system,mg,info,psi,rt_context=exx_state_owner)
      return
    endif
    if((.not.exx_state_owner%factor_experiment.or.product(system%kmesh)==1).and. &
       (info%isize_r>1.or.info%isize_o>1.or.dc_gamma_pairs_requested(system).or. &
        (exx_mlwf_norm_fraction>0d0.and.exx_mlwf_radius==0d0)).and. &
       ((theory=='dft').or. &
        (yn_dc=='n'.and.(yn_conventional_from_dcdft=='y'.or.is_global_hybrid(system_functional(system))).and. &
         (theory=='tddft_response'.or.theory=='tddft_pulse'))))then
      call refresh_spatial(system,mg,info,psi, &
      & rt_context=exx_state_owner)
      return
    endif
    if(info%isize_r/=1.or.info%isize_o/=1.or.info%numm/=1) &
      error stop 'HSE06: initial native support requires k-only MPI distribution'

    if(use_wannier_exchange(system).and.(.not.exx_state_owner%factor_experiment.or.product(system%kmesh)==1))then
      call refresh_wannier(system,mg,info,psi, &
      & rt_context=exx_state_owner)
      return
    endif
    ng=product(mg%num);nk=system%nk;no=system%no;n=mg%num(1);mesh=nint(real(nk,8)**(1d0/3d0))
    if(spinor_pair_fft_requested().and.(exx_kpoint_backend/='cpu'.or.any(hse_q_reduction/=1))) &
      error stop 'Spinor pair FFT requires CPU full-q native exchange'
    if(system%kpoint_input%use_symmetry.and.any(hse_q_reduction/=1))error stop 'HSE q reduction requires an unreduced k mesh'
    if(system%kpoint_input%use_symmetry)mesh=system%kmesh(1)
    if(.not.system%kpoint_input%use_symmetry.and.product(system%kmesh)/=nk) &
      error stop 'EXX: full uniform k mesh required'
    if(system%kpoint_input%use_symmetry)then
      if(.not.exx_state_owner%native_resources_enabled)then
        if(any(mg%num/=num_rgrid))error stop 'EXX symmetry requires the full real-space grid'
      endif
    endif
    thermal_spinor=yn_spinorbit=='y'.and.temperature>=0d0
    if(.not.system%kpoint_input%use_symmetry.and.maxval(abs(system%wtk-1d0/nk))>1d-12) &
      error stop 'HSE06: uniform k weights required'
    if(.not.thermal_spinor.and..not.native_integer_occupation(system)) &
      error stop 'HSE06: occupied-prefix spin pairs required'
    if(thermal_spinor)then
      if(any(hse_q_reduction/=1).or.info%isize_r/=1.or.info%isize_o/=1) &
        error stop 'Thermal SOC HSE: full uniform k-only mesh required'
    endif
    if(info%io_s/=1.or.info%io_e/=no.or.info%numk<1)error stop 'HSE06: unsupported orbital layout'
    if(.not.allocated(exx_state_owner%native_kernel))allocate(exx_state_owner%native_kernel)
    if(all(exx_state_owner%native_kernel%n==0))then
      nullify(create_backend)
      if(exx_kpoint_backend=='cufft')create_backend=>exx_k_cufft_create
      ! Keep representatives persistent; expand only rank-local stars for EXX.
      if(system%kpoint_input%use_symmetry)then
        call symmetry_init(exx_state_owner%native_symmetry,mg%num,system%hgs,system%vec_k(:,:nk),system%wtk(:nk), &
          system%kpoint_input%symmat_a,system%kpoint_input%symmat_b,system%kmesh,ierr,cell=system%primitive_a)
        call comm_summation(ierr,total_error,info%icomm_rko)
        if(total_error/=0)error stop 'HSE symmetry: invalid grid, stars or weights'
        call symmetry_validate_atoms(exx_state_owner%native_symmetry,system%Rion,system%kion,ierr)
        call comm_summation(ierr,total_error,info%icomm_rko)
        if(total_error/=0)error stop 'HSE symmetry: operation does not preserve atoms'
        first_full=exx_state_owner%native_symmetry%first(info%ik_s)
        count_full=exx_state_owner%native_symmetry%first(info%ik_e+1)-first_full
        call exx_k_kernel_init(exx_state_owner%native_kernel,mg%num,system%kmesh,system%hgs, &
          exx_state_owner%native_symmetry%full_k,system_exchange_screening(system), &
          max(1,min(16,64/info%isize_k)),ierr,first_full,count_full, &
          block_rows=hse_block_rows,fft_layout=hse_fft_layout, &
          create_backend=create_backend,coulomb_radius=system_coulomb_radius(system), &
          cell=system%primitive_a,reciprocal=system%primitive_b,q_reduction=hse_q_reduction, &
          defer_fft=spinor_pair_fft_requested().or.(yn_spinorbit=='y'.and.yn_exx_fft_orbits=='y'), &
          fft_layout_sync=synchronize_fft_layout)
      else
        call exx_k_kernel_init(exx_state_owner%native_kernel,mg%num,system%kmesh,system%hgs,system%vec_k, &
          system_exchange_screening(system), &
          max(1,min(16,64/info%isize_k)),ierr,info%ik_s,info%numk, &
          block_rows=hse_block_rows,fft_layout=hse_fft_layout, &
          create_backend=create_backend,coulomb_radius=system_coulomb_radius(system), &
          cell=system%primitive_a,reciprocal=system%primitive_b,q_reduction=hse_q_reduction, &
          defer_fft=spinor_pair_fft_requested().or.(yn_spinorbit=='y'.and.yn_exx_fft_orbits=='y'), &
          fft_layout_sync=synchronize_fft_layout)
      endif
      call comm_summation(ierr,total_error,info%icomm_rko)
      if(total_error/=0)error stop 'HSE06: kernel initialization failed'
      if(allocated(exx_state_owner%native_kernel%accelerator).and.info%id_k==0) &
        write(*,'(a)')'EXX_KPOINT_BACKEND=cufft (experimental; GPU stage includes packing, FFT, kernel and transfers)'
      call check_k_algorithm(exx_state_owner,system)
      if(info%id_k==0)then
        call report_k_algorithm(exx_state_owner,system)
        write(*,'(a)')'EXX_DISTRIBUTED_K: full-support density tiles'
        write(*,'(a,3i6,a,i10)')'EXX_Q_REDUCTION: ',hse_q_reduction,' q points=',exx_state_owner%native_kernel%nq
        if(k_pairs_enabled(exx_state_owner,system))write(*,'(a,2i18)')'EXX_K_PAIR_SYMMETRY unique/full pairs: ', &
          int(exx_state_owner%native_kernel%ng,int64)*(int(exx_state_owner%native_kernel%ng,int64)+1)/2,int(exx_state_owner%native_kernel%ng,int64)**2
        if(propagator=='hse_taylor4_full')write(*,'(a)')'EXX_DIRECT: ACE construction disabled'
      endif
    endif
    if(yn_spinorbit=='y')ng=2*ng
    allocate(local(ng,no,info%numk))
    call exx_pack(psi,mg,info,local)
    if(thermal_spinor)then
      allocate(weighted(ng,no,info%numk))
      call exx_spinor_weight_sources(local,system%rocc(:,info%ik_s:info%ik_e,1),weighted,ierr)
      call comm_summation(ierr,total_error,info%icomm_k)
      if(total_error/=0)error stop 'Thermal SOC HSE: invalid source occupations/orbitals'
      electrons=0d0
      do j=1,system%nk
        if(.not.ieee_is_finite(system%wtk(j)))error stop 'Thermal SOC HSE: nonfinite k weight'
        electrons=electrons+sum(system%rocc(:,j,1))*system%wtk(j)
      enddo
      if(.not.ieee_is_finite(electrons))error stop 'Thermal SOC HSE: nonfinite electron count'
      if(abs(electrons-nelec)>1d-8*max(1,nelec))error stop 'Thermal SOC HSE: electron count mismatch'
    endif
    ierr=1
    if(allocated(exx_state_owner%cached_source))then
      if(all(shape(exx_state_owner%cached_source)==shape(local)))then
        if(all(exx_state_owner%cached_source==local))ierr=0
      endif
    endif
    if(thermal_spinor)then
      if(.not.allocated(exx_state_owner%cached_occupation))then
        ierr=1
      else if(any(shape(exx_state_owner%cached_occupation)/=shape(system%rocc(:,:,1))))then
        ierr=1
      else if(any(exx_state_owner%cached_occupation/=system%rocc(:,:,1)))then
        ierr=1
      endif
    endif
    ! A damped SCF ACE must be replaced by physical exchange even if psi
    ! and occupations are unchanged (e.g. a zero-CG diagnostic).
    if(exx_state_owner%scf_operator_damped)ierr=1
    call comm_summation(ierr,total_error,info%icomm_k)
    if(total_error==0)return
    exx_state_owner%scf_operator_damped=.false.
    allocate(w(ng,no,info%numk))
    if(propagator=='hse_taylor4_full')then
      if(thermal_spinor)then
        exx_state_owner%full_source=weighted
      else
        exx_state_owner%full_source=local
      endif
    endif
    if(thermal_spinor)then
      call apply_distributed(weighted,local,w,info,ierr,exx_state_owner,system)
    else
      call apply_distributed(local,local,w,info,ierr,exx_state_owner,system)
    endif
    call comm_summation(ierr,total_error,info%icomm_k)
    if(exx_state_owner%recover_errors)then
      call recovery_agree(info,30,total_error/=0,exx_state_owner%step_error)
      if(exx_state_owner%step_error/=0)return
    endif
    if(total_error/=0)error stop 'HSE06: distributed exchange action failed'
    if(.not.exx_state_owner%reported_team.and.info%id_k==0)write(*,'(a,i0)')'HSE_OPENMP threads=',exx_state_owner%native_kernel%threads_used
    if(.not.exx_state_owner%reported_team.and.info%id_k==0)write(*,'(a,l1)')'HSE_FFT_CONTIGUOUS=',exx_state_owner%native_kernel%contiguous_fft
    exx_state_owner%reported_team=.true.
    if(propagator/='hse_taylor4_full')then
      call exx_ace_build(exx_state_owner%ace,local,w,system%hvol,ierr,thread_control=exx_blas_thread_control)
      call comm_summation(ierr,total_error,info%icomm_k)
      if(exx_state_owner%recover_errors)then
      call recovery_agree(info,30,total_error/=0,exx_state_owner%step_error)
      if(exx_state_owner%step_error/=0)return
    endif
    if(total_error/=0)error stop 'HSE06: ACE metric failed'
    endif
    exx_state_owner%cached_source=local
    if(thermal_spinor)exx_state_owner%cached_occupation=system%rocc(:,:,1)
    ex=0d0
    do j=1,info%numk
      ex=ex+.5d0*system_exchange_fraction(system)*sum(sum(real(conjg(local(:,:,j))*w(:,:,j),8),dim=1) &
        *system%rocc(:,info%ik_s+j-1,1))*system%hvol*system%wtk(info%ik_s+j-1)
    enddo
    call comm_summation(ex,exx_state_owner%exx_exchange_energy,info%icomm_k)
  contains
    subroutine synchronize_fft_layout(times,contiguous,status,automatic)
      use exx_fft_consensus,only:agree_fft_layout
      implicit none
      real(8),intent(inout)::times(2)
      logical,intent(inout)::contiguous
      logical,intent(in)::automatic
      integer,intent(inout)::status
      call agree_fft_layout(times,contiguous,status,info%icomm_k,automatic)
    end subroutine
  end subroutine

  ! Explicit SCF bridge: distributed accepted rows, replicated full-grid FFT
  ! scratch per spatial rank. This is the existing fast k kernel, not a second
  ! exchange approximation. k ranks retain their normal distributed ring.
  subroutine refresh_kspace(system,mg,info,psi,owner)
    type(s_dft_system),intent(in)::system
    type(s_rgrid),intent(in)::mg
    type(s_parallel_info),intent(in)::info
    type(s_orbital),intent(in)::psi
    type(s_exx_rt_state),target,intent(inout)::owner
    type(s_rgrid)::whole
    type(s_parallel_info)::ki
    type(s_orbital)::fullpsi
    complex(8),allocatable::local(:,:,:)
    integer::ng,j,k,x,y,z,g,f,status,total
    whole=mg;whole%is=1;whole%ie=num_rgrid;whole%num=num_rgrid
    ki=info;ki%isize_r=1;ki%id_r=0;ki%if_divide_rspace=.false.
    ng=product(mg%num)
    allocate(local(ng,info%numo,info%numk))
    call exx_pack(psi,mg,info,local)
    call gather(local)
    allocate(fullpsi%zwf(1:num_rgrid(1),1:num_rgrid(2),1:num_rgrid(3),1, &
      info%io_s:info%io_e,info%ik_s:info%ik_e,1))
    call exx_unpack(local,fullpsi,whole,ki)
    call gather(owner%cached_source)
    call gather(owner%cached_action)
    call gather(owner%full_source)
    call gather(owner%initial_source)
    call gather(owner%ace%factors)
    call gather(owner%scf_damped_ace%factors)
    call gather(owner%initial_ace%factors)
    call gather(owner%midpoint_ace%factors)
    call exx_refresh(system,whole,ki,fullpsi,rt_context=owner,rt_status=status)
    call comm_summation(status,total,info%icomm_rko)
    ! Always restore row-local ownership, including a handled failure.
    call scatter(owner%cached_source)
    call scatter(owner%cached_action)
    call scatter(owner%full_source)
    call scatter(owner%initial_source)
    call scatter(owner%ace%factors)
    call scatter(owner%scf_damped_ace%factors)
    call scatter(owner%initial_ace%factors)
    call scatter(owner%midpoint_ace%factors)
    if(total/=0)error stop 'Native k/spatial SCF exchange refresh failed'
  contains
    subroutine gather(a)
      complex(8),allocatable,intent(inout)::a(:,:,:)
      complex(8),allocatable::part(:,:,:),full(:,:,:)
      if(.not.allocated(a))return
      allocate(part(product(num_rgrid),size(a,2),size(a,3)),source=(0d0,0d0))
      do z=mg%is(3),mg%ie(3)
      do y=mg%is(2),mg%ie(2)
      do x=mg%is(1),mg%ie(1)
        g=1+x-mg%is(1)+mg%num(1)*(y-mg%is(2)+mg%num(2)*(z-mg%is(3)))
        f=x+num_rgrid(1)*(y-1+num_rgrid(2)*(z-1))
        part(f,:,:)=a(g,:,:)
      enddo
      enddo
      enddo
      allocate(full,mold=part)
      call comm_summation(part,full,size(part),info%icomm_r)
      call move_alloc(full,a)
    end subroutine
    subroutine scatter(a)
      complex(8),allocatable,intent(inout)::a(:,:,:)
      complex(8),allocatable::part(:,:,:)
      if(.not.allocated(a))return
      allocate(part(ng,size(a,2),size(a,3)))
      do z=mg%is(3),mg%ie(3)
      do y=mg%is(2),mg%ie(2)
      do x=mg%is(1),mg%ie(1)
        g=1+x-mg%is(1)+mg%num(1)*(y-mg%is(2)+mg%num(2)*(z-mg%is(3)))
        f=x+num_rgrid(1)*(y-1+num_rgrid(2)*(z-1))
        part(g,:,:)=a(f,:,:)
      enddo
      enddo
      enddo
      call move_alloc(part,a)
    end subroutine
  end subroutine

  ! Crystal reuse takes precedence. Unsupported backends retain their reference route.
  logical function k_pairs_enabled(exx_state_owner,system)
    use salmon_global, only: yn_exx_k_pair_symmetry
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_exx_rt_state),intent(in)::exx_state_owner
    k_pairs_enabled=yn_spinorbit=='n'.and.yn_exx_k_pair_symmetry=='y'.and. &
          .not.system%kpoint_input%use_symmetry.and.yn_dc=='n'.and. &
      exx_state_owner%native_kernel%nk>1.and.exx_state_owner%native_kernel%nq==exx_state_owner%native_kernel%nk.and..not.allocated(exx_state_owner%native_kernel%accelerator).and. &
      (theory=='dft'.or.theory=='tddft_response'.or.theory=='tddft_pulse')
  end function k_pairs_enabled

  subroutine apply_distributed(source,target,action,info,ierr,exx_state_owner,system)
    use exx_k_pairs, only: apply_k_pairs
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_exx_rt_state),intent(inout)::exx_state_owner
    complex(8),intent(in) :: source(:,:,:),target(:,:,:)
    complex(8),intent(out) :: action(:,:,:)
    type(s_parallel_info),intent(in) :: info
    integer,intent(out) :: ierr
    integer :: local_layout(2*info%isize_k),layout(2*info%isize_k),np,j,index,total_error
    integer :: band,ng
    complex(8),allocatable :: expanded_target(:,:,:),expanded_action(:,:,:),rotated(:,:),expanded_source(:,:,:)
    integer :: nelec,nstate
    real(8) :: temperature
    call hse_occupation_input(system,nelec,temperature,nstate)
    np=info%isize_k;local_layout=0
    local_layout(info%id_k+1)=info%ik_s
    local_layout(np+info%id_k+1)=info%numk
    if(system%kpoint_input%use_symmetry)then
      local_layout(info%id_k+1)=exx_state_owner%native_symmetry%first(info%ik_s)
      local_layout(np+info%id_k+1)=exx_state_owner%native_symmetry%first(info%ik_e+1)-exx_state_owner%native_symmetry%first(info%ik_s)
    endif
    call comm_summation(local_layout,layout,size(layout),info%icomm_k)
    if(yn_spinorbit=='y')then
      if(system%kpoint_input%use_symmetry)then
        call symmetry_expand(exx_state_owner%native_symmetry,info%ik_s,source,target,expanded_source,expanded_target,ierr)
        call comm_summation(ierr,total_error,info%icomm_k)
        if(total_error/=0)then
          ierr=1;return
        endif
        if(yn_exx_fft_orbits=='y')then
          ng=exx_state_owner%native_kernel%ng
          !$omp parallel do collapse(2) schedule(static)
          do j=1,size(expanded_source,3)
            do band=1,size(expanded_source,2)
              expanded_source(:ng,band,j)=expanded_source(:ng,band,j)*exx_state_owner%native_kernel%phase(:,j)
              expanded_source(ng+1:,band,j)=expanded_source(ng+1:,band,j)*exx_state_owner%native_kernel%phase(:,j)
            enddo
          enddo
          !$omp end parallel do
          call apply_spinor_k_orbits(exx_state_owner%native_kernel,exx_state_owner%native_symmetry,expanded_source,target,action,layout(:np),layout(np+1:), &
            info%ik_s,info%id_k,info%icomm_k,transpose_orbit_tiles,ierr)
          return
        endif
#ifdef USE_MPI
        if(service_enabled())then
          call service_exchange(exx_state_owner%native_kernel,exx_state_owner%native_symmetry,expanded_source,target,action,layout(:np), &
            layout(np+1:),info%ik_s,info%id_k,info%icomm_k,ierr)
          return
        endif
#endif
        allocate(expanded_action(size(expanded_target,1),size(expanded_target,2),size(expanded_target,3)))
        call exx_spinor_apply_distributed(exx_state_owner%native_kernel,expanded_source,expanded_target,expanded_action,layout(:np), &
          layout(np+1:),info%id_k,info%icomm_k,transpose_k_tiles,ierr,pair_tiles=spinor_pair_fft_requested())
        if(ierr==0)then
          do j=1,info%numk
            index=exx_state_owner%native_symmetry%first(info%ik_s+j-1)-exx_state_owner%native_symmetry%first(info%ik_s)+1
            action(:,:,j)=expanded_action(:,:,index)
          enddo
        endif
      else
        call exx_spinor_apply_distributed(exx_state_owner%native_kernel,source,target,action,layout(:np),layout(np+1:), &
          info%id_k,info%icomm_k,transpose_k_tiles,ierr,pair_tiles=spinor_pair_fft_requested())
      endif
      return
    endif
    if(system%kpoint_input%use_symmetry)then
      ierr=0
      if(.not.salmon_all_finite(real(source)).or..not.salmon_all_finite(aimag(source)))ierr=1
      call comm_summation(ierr,total_error,info%icomm_k)
      if(total_error/=0)then
        ierr=1;return
      endif
#ifdef USE_MPI
      if(service_enabled())then
        ! Keep the established one-spin occupied projector and little-group average.
        ! Only exact exchange uses world workers; SCF and learned ACE stay on k owners.
        index=size(source,2)
        if(temperature<0d0.and.yn_dc=='n')index=nelec/2
        call symmetry_expand(exx_state_owner%native_symmetry,info%ik_s,source(:,:index,:),target, &
          expanded_source,expanded_target,ierr)
        call comm_summation(ierr,total_error,info%icomm_k)
        if(total_error/=0)then
          ierr=1;return
        endif
        call service_exchange(exx_state_owner%native_kernel,exx_state_owner%native_symmetry,expanded_source,target,action,layout(:np), &
          layout(np+1:),info%ik_s,info%id_k,info%icomm_k,ierr)
        return
      endif
#endif
      ! Opt-in route. Keep the established full route as the default.
      if(yn_exx_fft_orbits=='y')then
        index=size(source,2)
        if(temperature<0d0.and.yn_dc=='n')index=nelec/2
        call apply_k_orbits(exx_state_owner%native_kernel,exx_state_owner%native_symmetry,source(:,:index,:),target,action,layout(:np),layout(np+1:), &
          info%ik_s,info%id_k,info%icomm_k,transpose_orbit_tiles,ierr)
        return
      endif
      allocate(rotated(size(source,1),size(source,2)))
      allocate(expanded_target(size(target,1),size(target,2),layout(np+info%id_k+1)))
      do j=1,size(expanded_target,3)
        index=layout(info%id_k+1)+j-1
        call symmetry_transform(exx_state_owner%native_symmetry,index,1,target(:,:,exx_state_owner%native_symmetry%owner(index)-info%ik_s+1), &
          expanded_target(:,:,j),ierr)
        if(ierr/=0)error stop 'HSE symmetry target transformation failed'
      enddo
      allocate(expanded_action(size(expanded_target,1),size(expanded_target,2),size(expanded_target,3)))
      ! Source argument is only a valid layout placeholder when the density callback is supplied.
      call exx_k_kernel_apply_distributed(exx_state_owner%native_kernel,expanded_target,expanded_target,expanded_action,layout(:np), &
        layout(np+1:),info%id_k,info%icomm_k,transpose_k_tiles,ierr,fill_density_contiguous)
      if(ierr==0)then
        do j=1,info%numk
          index=exx_state_owner%native_symmetry%first(info%ik_s+j-1)-exx_state_owner%native_symmetry%first(info%ik_s)+1
          action(:,:,j)=expanded_action(:,:,index)
        enddo
      endif
      return
    endif
    if(k_pairs_enabled(exx_state_owner,system))then
      if(temperature<0d0)then
        call apply_k_pairs(exx_state_owner%native_kernel,source(:,:nelec/2,:),target,action,layout(:np),layout(np+1:), &
          info%id_k,info%icomm_k,transpose_k_tiles,ierr)
      else
        call apply_k_pairs(exx_state_owner%native_kernel,source,target,action,layout(:np),layout(np+1:), &
          info%id_k,info%icomm_k,transpose_k_tiles,ierr)
      endif
    else if(temperature<0d0.and.yn_dc=='n')then
      call exx_k_kernel_apply_distributed(exx_state_owner%native_kernel,source(:,:nelec/2,:),target,action,layout(:np),layout(np+1:), &
        info%id_k,info%icomm_k,transpose_k_tiles,ierr)
    else
      call exx_k_kernel_apply_distributed(exx_state_owner%native_kernel,source,target,action,layout(:np),layout(np+1:), &
        info%id_k,info%icomm_k,transpose_k_tiles,ierr)
    endif
  contains
    subroutine fill_density_contiguous(j,lo,rows,density)
      implicit none
      integer,intent(in) :: j,lo,rows
      complex(8),intent(out) :: density(:,:)
      call fill_density(j,lo,rows,density)
    end subroutine fill_density_contiguous

    subroutine fill_density(j,lo,rows,density)
      implicit none
      integer,intent(in) :: j,lo,rows
      complex(8),intent(out) :: density(:,:)
      integer :: full_index,rep,op,g,stat,ng,no,blas_state(3)
      complex(8) :: coefficient
  ! BLAS explicit interface: sequence association preserves strided row offsets.
  interface
    subroutine zgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
      implicit none
      character(1),intent(in) :: transa,transb
      integer,intent(in) :: m,n,k,lda,ldb,ldc
      complex(8),intent(in) :: alpha,beta,a(lda,*),b(ldb,*)
      complex(8),intent(inout) :: c(ldc,*)
    end subroutine zgemm
  end interface

      ng=size(source,1);no=size(source,2)
      if(temperature<0d0.and.yn_dc=='n')no=nelec/2
      full_index=layout(info%id_k+1)+j-1
      rep=exx_state_owner%native_symmetry%owner(full_index)-info%ik_s+1
      coefficient=cmplx(1d0/exx_state_owner%native_symmetry%multiplicity(full_index),0d0,8)
      density=0d0
      blas_state=0
      call exx_blas_thread_control(1,blas_state)
      ! Average little-group projectors one orbital block at a time, not copies of all orbitals.
      do op=1,exx_state_owner%native_symmetry%multiplicity(full_index)
        call symmetry_transform(exx_state_owner%native_symmetry,full_index,op,source(:,:no,rep),rotated(:,:no),stat)
        if(stat/=0)error stop 'HSE symmetry source transformation failed'
        do g=1,no
          rotated(:,g)=rotated(:,g)*exx_state_owner%native_kernel%phase(:,j)
        enddo
        call zgemm('N','C',rows,ng,no,coefficient,rotated(lo,1),ng,rotated(1,1),ng, &
          (1d0,0d0),density,size(density,1))
      enddo
      call exx_blas_thread_control(0,blas_state)
    end subroutine
  end subroutine

  subroutine transpose_k_tiles(send,recv,count,comm)
    ! Explicit context avoids a callback capturing an enclosing info argument.
    implicit none
    complex(8),intent(in) :: send(:)
    complex(8),intent(out) :: recv(:)
    integer,intent(in) :: count,comm
    call comm_alltoall(send,recv,comm,count)
  end subroutine transpose_k_tiles

  subroutine transpose_orbit_tiles(send,recv,sendcounts,recvcounts,displs,comm)
#ifdef USE_MPI
    use mpi
#endif
    implicit none
    complex(8),intent(in) :: send(:)
    complex(8),intent(out) :: recv(:)
    integer,intent(in) :: sendcounts(:),recvcounts(:),displs(:),comm
#ifdef USE_MPI
    integer :: ierr
#endif
#ifdef USE_MPI
    call MPI_Alltoallv(send,sendcounts,displs,MPI_DOUBLE_COMPLEX,recv,recvcounts,displs,MPI_DOUBLE_COMPLEX,comm,ierr)
    if(ierr/=MPI_SUCCESS)error stop 'HSE orbit transpose failed'
#else
    recv(:recvcounts(1))=send(:sendcounts(1))
#endif
  end subroutine transpose_orbit_tiles

  subroutine check_k_algorithm(exx_state_owner,system)
    use salmon_global, only: exx_k_algorithm
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_exx_rt_state),intent(in)::exx_state_owner
    if(exx_k_algorithm=='pair')then
      if(yn_spinorbit=='y')then
        if(.not.spinor_pair_fft_requested())error stop 'Manual spinor pair algorithm was not selected'
      else
        if(.not.k_pairs_enabled(exx_state_owner,system))error stop 'Manual scalar pair algorithm is unavailable for this layout'
      endif
    else if(exx_k_algorithm=='orbits')then
      if(.not.system%kpoint_input%use_symmetry) &
        error stop 'Manual orbit algorithm is unavailable for this layout'
    endif
  end subroutine check_k_algorithm

  subroutine report_k_algorithm(exx_state_owner,system)
    use salmon_global, only: exx_k_algorithm,exx_k_algorithm_reason
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_exx_rt_state),intent(in)::exx_state_owner
    character(16) :: selected
    selected='reference'
    if(yn_spinorbit=='y')then
      if(spinor_pair_fft_requested())selected='spinor_pair'
      if(system%kpoint_input%use_symmetry.and.yn_exx_fft_orbits=='y')selected='spinor_orbits'
    else if(system%kpoint_input%use_symmetry.and.yn_exx_fft_orbits=='y')then
      selected='orbits'
    else if(k_pairs_enabled(exx_state_owner,system))then
      selected='scalar_pair'
    endif
    write(*,'(4a)')'EXX_K_ALGORITHM requested=',trim(exx_k_algorithm),' selected=',trim(selected)
    write(*,'(2a)')'EXX_K_ALGORITHM reason=',trim(exx_k_algorithm_reason)
  end subroutine report_k_algorithm

  logical function spinor_pair_fft_requested()
    use salmon_global, only: yn_exx_spinor_pair_fft
    implicit none
    spinor_pair_fft_requested=yn_spinorbit=='y'.and.yn_exx_spinor_pair_fft=='y'
  end function

  subroutine exx_add_action(psi,hpsi,system,mg,info,rt_context,rt_status)
    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    type(spatial_exx_state),pointer::spatial
    integer,optional,intent(out)::rt_status
    type(s_orbital),intent(in) :: psi
    type(s_orbital),intent(inout) :: hpsi
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    integer :: ierr,ng,total_error,info_error
    integer,allocatable :: orbital_comm
    real(8) :: action_scale
    exx_state_owner=>resolve_rt_context(rt_context,system)
    spatial=>exx_state_owner%spatial
    if(exx_state_owner%recover_errors.and.exx_state_owner%step_error/=0)return
    if(present(rt_status))rt_status=0
    if(.not.associated(exx_state_owner,rt_owned).and..not.exx_state_owner%native_resources_enabled)then
      if(exx_state_owner%exx_force_full_action.or.propagator=='hse_taylor4_full'.or. &
          exx_state_owner%ace%packed.or.exx_state_owner%ace%metric_distributed.or. &
          exx_state_owner%initial_ace%packed.or.exx_state_owner%initial_ace%metric_distributed.or. &
          info%isize_r/=1.or.info%isize_o/=1.or.info%isize_k/=1.or.system%nspin/=1.or.yn_spinorbit=='y')then
        if(present(rt_status))then
          rt_status=1;return
        endif
        error stop 'Explicit HSE RT context supports serial dense ACE/lifecycle only'
      endif
    endif
    if(exx_state_owner%native_resources_enabled)then
      if(.not.native_owner_supported(system,info))then
        if(present(rt_status))then
          rt_status=2;return
        endif
        error stop 'Independent HSE owner requires serial full-support non-spinor configuration'
      endif
    endif
    if(exx_state_owner%native_resources_enabled)then
      if(system%kpoint_input%use_symmetry)then
        if(.not.owned_symmetry_grid_matches(system,mg))then
          if(present(rt_status))then
            rt_status=2;return
          endif
          error stop 'Independent HSE symmetry requires matching owned full grid'
        endif
      endif
    endif
    if(.not.exx_enabled(system))return
    if(propagator=='hse_taylor4_full')then
      if(.not.allocated(exx_state_owner%full_source))error stop 'EXX: direct exchange source is not initialized'
    else
      if(.not.exx_ace_ready(exx_state_owner%ace))error stop 'HSE06: occupied exchange source is not initialized'
    endif
    if(info%isize_o>1)orbital_comm=info%icomm_o
    ng=product(mg%num)*size(psi%zwf,4)
    if(exx_state_owner%ace%packed.and.product(system%kmesh)==1.and..not.exx_state_owner%exx_force_full_action.and. &
       .not.(exx_state_owner%taylor_active.and.propagator=='hse_taylor4_full'))then
      if(.not.(exx_state_owner%taylor_active.and.exx_state_owner%taylor_midpoint).or.exx_state_owner%initial_ace%packed)then
        call apply_packed_blocks()
        return
      endif
    endif
    if(allocated(exx_state_owner%target_work))then
      if(any(shape(exx_state_owner%target_work)/=[ng,info%numo,info%numk]))deallocate(exx_state_owner%target_work)
    endif
    if(allocated(exx_state_owner%action_work))then
      if(any(shape(exx_state_owner%action_work)/=[ng,info%numo,info%numk]))deallocate(exx_state_owner%action_work)
    endif
    if(.not.allocated(exx_state_owner%target_work))allocate(exx_state_owner%target_work(ng,info%numo,info%numk))
    if(.not.allocated(exx_state_owner%action_work))allocate(exx_state_owner%action_work(ng,info%numo,info%numk))
    call exx_pack(psi,mg,info,exx_state_owner%target_work)
    action_scale=system_exchange_fraction(system)
    if(use_wannier_exchange(system).and.exx_state_owner%exx_force_full_action)then
      if(info%isize_r>1.or.info%isize_o>1.or.dc_gamma_pairs_requested(system).or. &
         (exx_mlwf_norm_fraction>0d0.and.exx_mlwf_radius==0d0))then
        call spatial_exx_apply(spatial,num_rgrid,system%hgs,[info%isize_x,info%isize_y,info%isize_z], &
          [info%id_x,info%id_y,info%id_z],[info%icomm_x,info%icomm_y,info%icomm_z],info%icomm_r, &
          system_coulomb_radius(system),exx_state_owner%target_work,exx_state_owner%action_work,info_error, &
          omega=system_exchange_screening(system),comm_o=orbital_comm)
        if(info_error/=0)error stop 'Spatial EXX: full action failed'
      else
        call apply_wannier_collective(exx_state_owner%target_work,exx_state_owner%action_work,info,system)
      endif
      ierr=0
    else if(exx_state_owner%exx_force_full_action)then
      call apply_distributed(exx_state_owner%cached_source,exx_state_owner%target_work,exx_state_owner%action_work, &
          info,ierr,exx_state_owner,system)
      call comm_summation(ierr,total_error,info%icomm_k)
      if(total_error/=0)error stop 'EXX: full trial action failed'
    else if(propagator=='hse_taylor4_full')then
      if(exx_state_owner%taylor_active.and.exx_state_owner%taylor_midpoint)then
        ! Keep occupied-prefix selection separate for the two endpoint sources.
        call apply_distributed(exx_state_owner%initial_source,exx_state_owner%target_work, &
          exx_state_owner%action_work,info,ierr,exx_state_owner,system)
        call comm_summation(ierr,total_error,info%icomm_k)
        if(total_error/=0)error stop 'HSE Taylor initial full target action failed'
        call add_mesh_action(.5d0*system_exchange_fraction(system))
        call apply_distributed(exx_state_owner%full_source,exx_state_owner%target_work,exx_state_owner%action_work, &
          info,ierr,exx_state_owner,system)
        action_scale=.5d0*system_exchange_fraction(system)
      else if(exx_state_owner%taylor_active)then
        call apply_distributed(exx_state_owner%initial_source,exx_state_owner%target_work, &
          exx_state_owner%action_work,info,ierr,exx_state_owner,system)
      else
        call apply_distributed(exx_state_owner%full_source,exx_state_owner%target_work,exx_state_owner%action_work, &
          info,ierr,exx_state_owner,system)
      endif
      call comm_summation(ierr,total_error,info%icomm_k)
      if(total_error/=0)error stop 'HSE Taylor full target action failed'
    else
      if(exx_state_owner%taylor_active.and.exx_state_owner%taylor_midpoint)then
        call apply_endpoint(exx_state_owner%initial_ace)
        if(exx_state_owner%recover_errors)then
          call recovery_agree(info,30,ierr/=0,exx_state_owner%step_error)
          if(exx_state_owner%step_error/=0)return
        endif
        if(ierr/=0)error stop 'HSE06: initial endpoint ACE application failed'
        call add_mesh_action(.5d0*system_exchange_fraction(system))
        action_scale=.5d0*system_exchange_fraction(system)
      endif
      call apply_endpoint(exx_state_owner%ace)
    endif
    if(exx_state_owner%recover_errors)then
      call recovery_agree(info,30,ierr/=0,exx_state_owner%step_error)
      if(exx_state_owner%step_error/=0)return
    endif
    if(ierr/=0)error stop 'HSE06: ACE application failed'
    call add_mesh_action(action_scale)
  contains
    subroutine apply_packed_blocks()
      implicit none
      integer,parameter :: block_size=32
      integer :: largest,first,count,ix,iy,iz,j,g,io
      real(8) :: block_weight
      largest=info%numo
      call comm_get_max(largest,info%icomm_ro)
      do first=1,largest,block_size
        count=max(0,min(block_size,info%numo-first+1))
        ! Zero-column orbital ranks still join every packed-ACE collective.
        if(allocated(exx_state_owner%target_work))then
          if(any(shape(exx_state_owner%target_work)/=[ng,count,1]))deallocate(exx_state_owner%target_work)
        endif
        if(allocated(exx_state_owner%action_work))then
          if(any(shape(exx_state_owner%action_work)/=[ng,count,1]))deallocate(exx_state_owner%action_work)
        endif
        if(.not.allocated(exx_state_owner%target_work))allocate(exx_state_owner%target_work(ng,count,1))
        if(.not.allocated(exx_state_owner%action_work))allocate(exx_state_owner%action_work(ng,count,1))
!$omp parallel do collapse(3) default(none) schedule(static) &
!$omp private(j,iz,iy,ix,g,io) shared(count,mg,info,first,psi,exx_state_owner)
        do j=1,count
          do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
              do ix=mg%is(1),mg%ie(1)
                g=1+(ix-mg%is(1))+mg%num(1)*((iy-mg%is(2))+mg%num(2)*(iz-mg%is(3)))
                io=info%io_s+first+j-2
                exx_state_owner%target_work(g,j,1)=psi%zwf(ix,iy,iz,1,io,info%ik_s,1)
              enddo
            enddo
          enddo
        enddo
!$omp end parallel do
        block_weight=system_exchange_fraction(system)
        if(exx_state_owner%taylor_active.and.exx_state_owner%taylor_midpoint)then
          call apply_endpoint(exx_state_owner%initial_ace)
          if(ierr/=0)error stop 'EXX: initial block ACE application failed'
          block_weight=.5d0*system_exchange_fraction(system)
          call add_mesh_action(block_weight,first,first+count-1)
        endif
        call apply_endpoint(exx_state_owner%ace)
        if(ierr/=0)error stop 'EXX: block ACE application failed'
        call add_mesh_action(block_weight,first,first+count-1)
      enddo
    end subroutine apply_packed_blocks
    subroutine apply_endpoint(state)
      implicit none
      type(s_exx_ace),intent(in) :: state
      if(info%isize_o>1.or.state%packed.or.exx_state_owner%orbital_scf)then
        call orbital_ace_apply(state,exx_state_owner%target_work,exx_state_owner%action_work,info%icomm_r,info%icomm_o,ierr)
      else if(info%isize_r>1)then
        call exx_ace_apply(state,exx_state_owner%target_work,exx_state_owner%action_work,ierr,sum_spatial,exx_blas_thread_control)
      else
        call exx_ace_apply(state,exx_state_owner%target_work,exx_state_owner%action_work,ierr,thread_control=exx_blas_thread_control)
      endif
    end subroutine
    subroutine add_mesh_action(weight,first,last)
      implicit none
      real(8),intent(in) :: weight
      integer,optional,intent(in) :: first,last
      integer :: ix,iy,iz,io,ik,g,is,io_first,io_last,ngrid
      io_first=info%io_s;io_last=info%io_e;ngrid=product(mg%num)
      if(present(first))io_first=info%io_s+first-1
      if(present(last))io_last=info%io_s+last-1
!$omp parallel do collapse(3) default(shared) private(ik,io,is,ix,iy,iz,g) schedule(static)
      do ik=info%ik_s,info%ik_e;do io=io_first,io_last;do is=1,size(hpsi%zwf,4)
        do iz=mg%is(3),mg%ie(3);do iy=mg%is(2),mg%ie(2);do ix=mg%is(1),mg%ie(1)
          g=1+(ix-mg%is(1))+mg%num(1)*((iy-mg%is(2))+mg%num(2)*(iz-mg%is(3)))+(is-1)*ngrid
          hpsi%zwf(ix,iy,iz,is,io,ik,1)=hpsi%zwf(ix,iy,iz,is,io,ik,1) &
            +weight*exx_state_owner%action_work(g,io-io_first+1,ik-info%ik_s+1)
        enddo;enddo;enddo
      enddo;enddo;enddo
!$omp end parallel do
    end subroutine
    subroutine sum_spatial(a)
      implicit none
      complex(8),intent(inout) :: a(:,:)
      complex(8) :: total(size(a,1),size(a,2))
      call comm_summation(a,total,size(a),info%icomm_r)
      a=total
    end subroutine
  end subroutine
  logical function use_wannier_exchange(system)
    implicit none
    type(s_dft_system),intent(in) :: system
    integer :: nelec,nstate
    real(8) :: temperature
    call hse_occupation_input(system,nelec,temperature,nstate)
    if(yn_spinorbit=='y')then
      use_wannier_exchange=.false.
      return
    endif
    ! Full-support ordinary multi-k exchange uses density-tile MPI for every
    ! hybrid. Retain Wannier for DC, fractional/empty states and localized support.
    use_wannier_exchange=temperature>=0d0.or.yn_dc=='y'.or. &
      exx_mlwf_radius>0d0.or.exx_mlwf_norm_fraction>0d0.or. &
      (yn_hse_wannier=='y'.and.(product(system%kmesh)==1.or.temperature>=0d0.or. &
        (nstate>0.and.nstate*2/=nelec)))
  end function

  subroutine refresh_spatial(system,mg,info,psi,rt_context)
    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    type(spatial_exx_state),pointer::spatial
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: psi
    complex(8),allocatable :: local(:,:,:),w(:,:,:)
    real(8) :: ex,offdiag(3,3)
    integer :: status,total,changed,maxiter,j,io
    integer,allocatable :: orbital_comm
    real(8),allocatable :: radii(:),loss(:)
    logical,allocatable :: protected(:)
    integer(int64) :: fft_work(3)
    logical :: was_active,screen_fallback,radius_covers_cell
    integer :: requested_screen_mode
    real(8) :: correction_norm,accepted_bound
    integer :: adaptive_bad
    real(8) :: mask_diagnostic(2),mask_maximum(2)
    exx_state_owner=>resolve_rt_context(rt_context,system)
    spatial=>exx_state_owner%spatial
    if((exx_mlwf_norm_fraction>0d0.or.exx_mlwf_radius>0d0).and.theory/='dft')exx_state_owner%exx_adaptive_ready=.true.
    if(info%isize_k/=1.or.info%numm/=1.or.system%nk/=1) &
      error stop 'Spatial EXX: Gamma spatial domains required'
    if(info%isize_o>system%no)error stop 'Spatial EXX: each orbital group must own at least one state'
    if(any(system%kmesh/=1).or.maxval(abs(system%vec_k))>1d-12.or.system%kpoint_input%use_symmetry) &
      error stop 'Spatial EXX: unshifted Gamma required'
    if(theory/='dft'.and.maxval(abs(system%rocc-2d0))>1d-12) &
      error stop 'Spatial EXX: occupied spin pairs required'
    if(theory/='dft'.and.propagator/='hse_taylor4')error stop 'Spatial EXX: Taylor4 ACE required'
    if(any(mg%num/=num_rgrid/[info%isize_x,info%isize_y,info%isize_z]).or. &
       any(mg%is/=1+[info%id_x,info%id_y,info%id_z]*mg%num)) &
      error stop 'Spatial EXX: Cartesian mesh layout mismatch'
    offdiag=system%primitive_a
    do j=1,3
      offdiag(j,j)=0d0
    enddo
    if(maxval(abs(offdiag))>1d-12)error stop 'Spatial EXX: orthogonal cell required'
    if(info%isize_o>1.or.exx_state_owner%orbital_scf)orbital_comm=info%icomm_o
    ! Apply and refresh never run concurrently. Reuse their mesh workspaces.
    if(allocated(exx_state_owner%target_work))then
      if(all(shape(exx_state_owner%target_work)==[product(mg%num),info%numo,1]))then
        call move_alloc(exx_state_owner%target_work,local)
      else
        deallocate(exx_state_owner%target_work)
      endif
    endif
    if(.not.allocated(local))allocate(local(product(mg%num),info%numo,1))
    call exx_pack(psi,mg,info,local)
    changed=1
    if(allocated(exx_state_owner%cached_source).and.allocated(exx_state_owner%cached_occupation))then
      if(all(shape(exx_state_owner%cached_source)==shape(local)).and.all(shape(exx_state_owner%cached_occupation)==shape(system%rocc(:,:,1))))then
        if(all(exx_state_owner%cached_source==local).and.all(exx_state_owner%cached_occupation==system%rocc(:,:,1)))changed=0
      endif
    endif
    if((exx_mlwf_norm_fraction>0d0.or.exx_mlwf_radius>0d0).and.(exx_state_owner%exx_adaptive_ready.neqv.exx_state_owner%cached_adaptive_ready))changed=1
    call comm_summation(changed,total,info%icomm_ro)
    if(total==0)then
      call move_alloc(local,exx_state_owner%target_work)
      return
    endif
    ! The old source has served its equality check; the new source replaces it.
    if(allocated(exx_state_owner%cached_source))deallocate(exx_state_owner%cached_source)
    maxiter=0
    if(mod(spatial%updates,exx_mlwf_interval)==0)maxiter=exx_mlwf_maxiter
    if((exx_mlwf_norm_fraction>0d0.or.exx_mlwf_radius>0d0).and.exx_state_owner%exx_adaptive_ready.and.spatial%last_localization_status/=0) &
      maxiter=exx_mlwf_maxiter
    spatial%seed_localized=exx_mlwf_norm_fraction>0d0.or.exx_mlwf_radius>0d0
    spatial%retain_accepted_gauge=theory/='dft'
    if(dc_canonical())then
      maxiter=0
      call spatial_exx_canonical_source(spatial,local,system%rocc(info%io_s:info%io_e,:,1), &
        info%icomm_r,status,comm_o=orbital_comm)
      if(info%id_ro==0.and.spatial%updates==1)then
        write(*,'(a)')'EXX_DC canonical full-fragment source (spatial)'
      endif
    else
    call spatial_exx_refresh(spatial,num_rgrid,system%hgs,[info%isize_x,info%isize_y,info%isize_z], &
      [info%id_x,info%id_y,info%id_z],[info%icomm_x,info%icomm_y,info%icomm_z],info%icomm_r,local, &
      maxiter,exx_mlwf_tolerance,status,occupation=system%rocc(info%io_s:info%io_e,:,1),comm_o=orbital_comm)
    endif
    if(status/=0)error stop 'Spatial EXX: source refresh failed'
    call sparse_clear(spatial%sparse_source)
    if((exx_mlwf_norm_fraction>0d0.or.exx_mlwf_radius>0d0).and.exx_mlwf_norm_fraction<1d0.and.theory/='dft')then
      if(spatial%last_localization_status/=0) &
        error stop 'Adaptive RT requires an accepted transported MLWF gauge; refine initial localization'
    endif
    if(spatial%retained_gauge.and.info%id_ro==0)write(*,'(a,i10)') &
      'EXX_SPATIAL retained accepted transported gauge at refresh: ',spatial%updates
    was_active=exx_state_owner%adaptive_active
    ! A sphere this large cannot truncate any point under the periodic metric;
    ! as for fraction=1, its exact action does not require a converged gauge.
    radius_covers_cell=exx_mlwf_radius>0d0.and. &
      exx_mlwf_radius>=sqrt(sum((.5d0*num_rgrid*system%hgs)**2))
    exx_state_owner%adaptive_active=.not.dc_canonical().and. &
      (exx_mlwf_norm_fraction>0d0.or.exx_mlwf_radius>0d0).and.exx_state_owner%exx_adaptive_ready.and. &
      (spatial%last_localization_status==0.or.radius_covers_cell.or. &
       (exx_mlwf_norm_fraction==1d0.and.exx_mlwf_radius==0d0))
    if(exx_state_owner%adaptive_active)then
      allocate(radii(info%numo),loss(info%numo),protected(info%numo))
      call adaptive_source_mask(num_rgrid,system%hgs,mg%is-1,mg%num,info%icomm_r,spatial%source, &
        merge(exx_mlwf_norm_fraction,.999d0,exx_mlwf_norm_fraction>0d0),radii,loss,protected,status, &
        fixed_radius=exx_mlwf_radius)
      call comm_summation(status,adaptive_bad,info%icomm_ro)
      if(adaptive_bad/=0)error stop 'Adaptive EXX: source mask failed'
      mask_diagnostic=[maxval(loss),maxval(radii)]
      call comm_get_max(mask_diagnostic,mask_maximum,2,info%icomm_ro)
      if(info%id_ro==0)then
        if(exx_mlwf_radius>0d0)then
          write(*,'(a,3es18.9)') 'EXX_FIXED radius/min retained norm/max loss: ', &
            exx_mlwf_radius,1d0-mask_maximum(1),mask_maximum(1)
          call warn_fixed_radius(mask_maximum(1), &
      & rt_context=exx_state_owner)
        else
          write(*,'(a,3es18.9)')'EXX_ADAPTIVE fraction/max radius/max norm loss: ', &
            exx_mlwf_norm_fraction,mask_maximum(2),mask_maximum(1)
        endif
      endif
    endif
    if(exx_mlwf_radius>0d0)exx_state_owner%finite_support_localized=exx_state_owner%adaptive_active
    if(exx_state_owner%adaptive_active.neqv.was_active)exx_state_owner%exx_support_changed=.true.
    exx_state_owner%cached_adaptive_ready=exx_state_owner%exx_adaptive_ready
    spatial%compact=exx_state_owner%adaptive_active.and.exx_local_fft=='auto'
    nullify(spatial%local_batch,spatial%create_local_backend)
    if(exx_local_backend=='cufft')then
      spatial%create_local_backend=>exx_cufft_create
      spatial%local_batch_size=exx_gpu_batch_size
      if(spatial%updates==1.and.info%id_ro==0)write(*,'(a,i0)') &
        'EXX experimental cuFFT compact backend; CPU global fallback; per-rank batch=',exx_gpu_batch_size
    endif
    requested_screen_mode=0
    if(exx_pair_screening=='diagnose')requested_screen_mode=1
    if(exx_pair_screening=='on')requested_screen_mode=2
    if(dc_canonical())requested_screen_mode=0
    spatial%screen_mode=requested_screen_mode;spatial%screen_tolerance=exx_pair_tolerance/2d0
    fft_work=0_int64
    if(allocated(exx_state_owner%action_work))then
      if(all(shape(exx_state_owner%action_work)==[product(mg%num),info%numo,1]))then
        call move_alloc(exx_state_owner%action_work,w)
      else
        deallocate(exx_state_owner%action_work)
      endif
    endif
    if(.not.allocated(w))allocate(w(product(mg%num),info%numo,1))
    call apply_exchange_action()
    if(status/=0)error stop 'Spatial EXX: exchange action failed'
    if(requested_screen_mode/=0.and.info%id_ro==0)write(*,'(a,i2,2i18,2es18.9)') &
      'EXX_PAIR mode/candidates/skipped/action bound/max rank CPU seconds: ',requested_screen_mode, &
      spatial%screen_candidates,spatial%screen_skipped,spatial%screen_bound,spatial%screen_cpu_seconds
    if(requested_screen_mode/=0.and.info%id_ro==0)write(*,'(a,2i18)') &
      'EXX_PAIR generated grid products/catalogue entries: ',spatial%pair_products,spatial%pair_catalog_entries
    if(requested_screen_mode/=0.and.info%id_ro==0)write(*,'(a,i18)') &
      'EXX_PAIR evaluated product points: ',spatial%pair_product_points
    correction_norm=0d0;accepted_bound=0d0
    if(requested_screen_mode==2.and.spatial%screen_skipped>0)then
      call orbital_hermitian_action(local,w,system%hvol,info%icomm_r,info%icomm_o, &
        exx_pair_tolerance-spatial%screen_bound,correction_norm,status)
      accepted_bound=spatial%screen_bound+correction_norm
    endif
    if(status==0)call build_exchange_ace()
    screen_fallback=status/=0.and.requested_screen_mode==2.and.spatial%screen_skipped>0
    if(screen_fallback)then
      ! Pair-dependent omissions need not define a Hermitian input metric.
      ! Never relax the existing Hermitian/positive ACE validation to accept them.
      spatial%screen_mode=0;accepted_bound=0d0
      call apply_exchange_action()
      if(status/=0)error stop 'Spatial EXX: unscreened fallback action failed'
      call build_exchange_ace()
    endif
    if(status/=0)error stop 'Spatial EXX: ACE build failed'
    if(requested_screen_mode/=0.and.info%id_ro==0)write(*,'(a,l1,a,es18.9)') &
      'EXX_PAIR unscreened ACE fallback: ',screen_fallback,' accepted action bound: ', &
      accepted_bound
    if(exx_state_owner%adaptive_active.and.info%id_ro==0)write(*,'(a,3i18)') &
      'EXX_ADAPTIVE local/global pairs/local FFT points (orbital group 0): ', &
      fft_work
    exx_state_owner%cached_occupation=system%rocc(:,:,1)
    ex=0d0
    do j=1,info%numo
      io=info%io_s+j-1
      ex=ex+.5d0*system_exchange_fraction(system)*system%hvol*system%rocc(io,1,1)*system%wtk(1) &
        *real(sum(conjg(local(:,j,1))*w(:,j,1)),8)
    enddo
    call comm_summation(ex,exx_state_owner%exx_exchange_energy,info%icomm_ro)
    call move_alloc(local,exx_state_owner%cached_source)
    if(yn_dc=='y')then
      call move_alloc(w,exx_state_owner%cached_action)
    else
      if(allocated(exx_state_owner%cached_action))deallocate(exx_state_owner%cached_action)
      call move_alloc(w,exx_state_owner%action_work)
    endif
    if(.not.dc_canonical().and.info%id_ro==0.and.(spatial%updates==1.or.maxiter>0)) &
      write(*,'(a,3i8,3es16.7)')'EXX_SPATIAL refresh/iterations/status/spread/gradient/overlap: ', &
      spatial%updates,spatial%iterations,spatial%localization_status,spatial%spread,spatial%gradient,spatial%min_singular
  contains
    subroutine apply_exchange_action()
!$    use omp_lib, only: omp_get_max_threads,omp_in_parallel
      implicit none
      complex(8),allocatable :: localized_action(:,:,:),adjoint(:,:)
      integer,allocatable :: counts(:)
      integer :: first,workers,thread_state(3)
      if(requested_screen_mode==0.and.dc_gamma_pairs_requested(system))then
        call spatial_exx_apply(spatial,num_rgrid,system%hgs,[info%isize_x,info%isize_y,info%isize_z], &
          [info%id_x,info%id_y,info%id_z],[info%icomm_x,info%icomm_y,info%icomm_z],info%icomm_r, &
          system_coulomb_radius(system),local,w,status,omega=system_exchange_screening(system),comm_o=orbital_comm, &
          pair_weight=.5d0*system%rocc(info%io_s:info%io_e,1,1))
      else if(requested_screen_mode==0)then
        call spatial_exx_apply(spatial,num_rgrid,system%hgs,[info%isize_x,info%isize_y,info%isize_z], &
          [info%id_x,info%id_y,info%id_z],[info%icomm_x,info%icomm_y,info%icomm_z],info%icomm_r, &
          system_coulomb_radius(system),local,w,status,omega=system_exchange_screening(system),comm_o=orbital_comm)
      else
        allocate(localized_action(size(w,1),size(w,2),1))
        call spatial_exx_apply(spatial,num_rgrid,system%hgs,[info%isize_x,info%isize_y,info%isize_z], &
          [info%id_x,info%id_y,info%id_z],[info%icomm_x,info%icomm_y,info%icomm_z],info%icomm_r, &
          system_coulomb_radius(system),spatial%previous,localized_action,status, &
          omega=system_exchange_screening(system),comm_o=orbital_comm)
        if(status/=0)return
        if(allocated(spatial%gauge_tiles%matrix))then
          call gauge_tiles_rotate(spatial%gauge_tiles,localized_action(:,:,1),info%icomm_r,info%icomm_o, &
            w(:,:,1),status,adjoint=.true.)
          if(status/=0)return
        else
          adjoint=conjg(transpose(spatial%gauge(:,:,1)))
          if(info%isize_o>1)then
            call orbital_layout(info%numo,info%icomm_r,info%icomm_o,counts,first,status)
            if(status/=0)return
            call orbital_rotate(localized_action(:,:,1),adjoint,info%icomm_o,counts,first,w(:,:,1))
          else
            workers=1;thread_state=0
!$          workers=omp_get_max_threads()
!$          if(omp_in_parallel())workers=1
            if(workers>1)call exx_blas_thread_control(workers,thread_state)
            w(:,:,1)=matmul(localized_action(:,:,1),adjoint)
            call exx_blas_thread_control(0,thread_state)
          endif
        endif
      endif
      call record_fft_work()
    end subroutine
    subroutine record_fft_work()
      ! Count rejected ACE attempts too: each exchange call resets its counters.
      implicit none
      fft_work=fft_work+[spatial%local_pairs,spatial%global_pairs,spatial%local_points]
    end subroutine
    subroutine build_exchange_ace()
      implicit none
      if(info%isize_o>1.or.info%isize_r>1.or.exx_state_owner%orbital_scf)then
        call orbital_ace_build(exx_state_owner%ace,local,w,system%hvol,info%icomm_r,info%icomm_o,status,comm_matrix=info%icomm_ro)
      else
        call exx_ace_build(exx_state_owner%ace,local,w,system%hvol,status,sum_spatial,exx_blas_thread_control)
      endif
      call comm_summation(status,adaptive_bad,info%icomm_ro)
      status=adaptive_bad
    end subroutine
    subroutine sum_spatial(a)
      implicit none
      complex(8),intent(inout) :: a(:,:)
      complex(8) :: total(size(a,1),size(a,2))
      call comm_summation(a,total,size(a),info%icomm_r)
      a=total
    end subroutine
  end subroutine

  subroutine refresh_wannier(system,mg,info,psi,rt_context)
    use communication, only: comm_bcast
    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: psi
    complex(8),allocatable :: local(:,:,:),send(:,:,:),allpsi(:,:,:),allw(:,:,:),w(:,:,:)
    real(8) :: ex,offdiag(3,3)
    integer :: ng,no,nk,ik,j,status,changed,total_changed,maxiter
    integer :: nelec,nstate
    real(8) :: temperature
    call hse_occupation_input(system,nelec,temperature,nstate)
    exx_state_owner=>resolve_rt_context(rt_context,system)
    ng=product(mg%num);no=system%no;nk=system%nk
    if(system%kpoint_input%use_symmetry.or.nk/=product(system%kmesh))error stop 'HSE Wannier: full k mesh required'
    if(maxval(abs(system%wtk-1d0/nk))>1d-12)error stop 'HSE Wannier: uniform k weights required'
    if(temperature>=0d0.and.yn_dc=='n')then
      do ik=1,nk
        do j=1,no
          if(.not.ieee_is_finite(system%rocc(j,ik,1)).or.system%rocc(j,ik,1)<0d0.or. &
             system%rocc(j,ik,1)>2d0)error stop 'Thermal HSE: invalid source occupations'
        enddo
      enddo
      ex=sum(sum(system%rocc(:,:,1),dim=1)*system%wtk)
      if(abs(ex-real(nelec,8))>1d-8)error stop 'Thermal HSE: source electron count mismatch'
      if(wannier%updates==0.and.info%id_k==0) &
        write(*,'(a,3es24.15)')'EXX_THERMAL electrons/min/max occupation: ', &
          ex,minval(system%rocc(:,:,1)),maxval(system%rocc(:,:,1))
    endif
    offdiag=system%primitive_a
    do j=1,3
      offdiag(j,j)=0d0
    enddo
    if(maxval(abs(offdiag))>1d-12)error stop 'HSE Wannier: orthogonal cell required'
    if(info%io_s/=1.or.info%io_e/=no.or.info%numk<1)error stop 'HSE Wannier: invalid orbital layout'
    allocate(local(ng,no,info%numk))
    call exx_pack(psi,mg,info,local)
    changed=1
    if(allocated(exx_state_owner%cached_source).and.allocated(exx_state_owner%cached_occupation))then
      if(all(shape(exx_state_owner%cached_source)==shape(local)).and.all(shape(exx_state_owner%cached_occupation)==shape(system%rocc(:,:,1))))then
        if(all(exx_state_owner%cached_source==local).and.all(exx_state_owner%cached_occupation==system%rocc(:,:,1)))changed=0
      endif
    endif
    call comm_summation(changed,total_changed,info%icomm_k)
    if(total_changed==0)return
    allocate(send(ng,no,nk),allpsi(ng,no,nk),allw(ng,no,nk),w(ng,no,info%numk))
    send=0d0;send(:,:,info%ik_s:info%ik_e)=local
    call comm_summation(send,allpsi,size(send),info%icomm_k)
    status=0
    if(info%id_k==0)then
      if(wannier%ng==0)then
        if(system_functional(system)=='hse06')then
          call wannier_init(wannier,mg%num,system%kmesh,system%hgs,system%vec_k,system_exchange_screening(system),status)
        else
          call wannier_init(wannier,mg%num,system%kmesh,system%hgs,system%vec_k,0d0,status,system_coulomb_radius(system))
          write(*,*)'PBEh40 Coulomb radius (bohr): ', &
            merge(system_coulomb_radius(system),.5d0*minval(mg%num*system%kmesh*system%hgs),system_coulomb_radius(system)>0d0)
        endif
      endif
      if(status==0)then
        maxiter=0
        if(mod(wannier%updates,exx_mlwf_interval)==0)maxiter=exx_mlwf_maxiter
        if(dc_canonical())maxiter=0
        call wannier_refresh_source(wannier,allpsi,system%rocc(:,:,1),maxiter,exx_mlwf_tolerance,status, &
          localize=.not.dc_canonical().and..not.(temperature>=0d0.and.yn_dc=='n'))
        if(dc_canonical().and.wannier%updates==1)write(*,'(a)') 'EXX_DC canonical full-fragment source (k mesh)'
      endif
      if(status==0)then
        wannier%use_local_fft=exx_local_fft=='auto'
        call wannier_truncate_source(wannier,exx_mlwf_radius,status)
        call warn_fixed_radius(wannier%max_discarded_norm_fraction, &
      & rt_context=exx_state_owner)
        exx_state_owner%finite_support_localized=wannier%discarded_norm_fraction==0d0.or.wannier%last_localization_status==0
      endif
      if(status==0)call wannier_apply(wannier,allpsi,allw,status)
      if(.not.dc_canonical().and.status==0.and.(wannier%updates==1.or.maxiter>0))then
        write(*,'(a,3i7,3es16.7)')'EXX_WANNIER refresh/iterations/status/spread/gradient/overlap: ', &
        wannier%updates,wannier%localization_iterations,wannier%localization_status, &
        wannier%spread,wannier%gradient,wannier%min_singular
        if(wannier%use_local_fft)then
          write(*,'(a,4i18)')'EXX_FFT local/global pairs/actual/full pair grid points: ', &
            wannier%local_fft_pairs_executed,wannier%fft_pairs_executed-wannier%local_fft_pairs_executed, &
            wannier%fft_pair_grid_points,wannier%fft_pairs_executed*int(wannier%ngs,int64)
        endif
        if(exx_mlwf_radius>0d0)then
          write(*,'(a,es16.7,i8,2es16.7)')'EXX_MLWF radius (bohr)/protected/total loss/max factor loss: ', &
            exx_mlwf_radius,wannier%protected_sources,wannier%discarded_norm_fraction,wannier%max_discarded_norm_fraction
          if(wannier%localization_status/=0) &
            write(*,'(a)')'EXX_MLWF: localization not converged; finite-radius source approximation is gauge dependent.'
        else if(wannier%localization_status/=0)then
          write(*,'(a)')'EXX_WANNIER: localization not converged; retaining full-support exact exchange.'
        endif
      endif
    endif
    call comm_bcast(status,info%icomm_k,0)
    if(status/=0)error stop 'HSE Wannier: collective exchange refresh failed'
    call comm_bcast(exx_state_owner%finite_support_localized,info%icomm_k,0)
    call comm_bcast(allw,info%icomm_k,0)
    w=allw(:,:,info%ik_s:info%ik_e)
    call exx_ace_build(exx_state_owner%ace,local,w,system%hvol,status,thread_control=exx_blas_thread_control)
    call comm_summation(status,total_changed,info%icomm_k)
    if(total_changed/=0)error stop 'HSE Wannier: ACE construction metric failed'
    exx_state_owner%cached_source=local;exx_state_owner%cached_occupation=system%rocc(:,:,1)
    if(yn_dc=='y')then
      exx_state_owner%cached_action=w
    else if(allocated(exx_state_owner%cached_action))then
      deallocate(exx_state_owner%cached_action)
    endif
    ex=0d0
    do ik=1,info%numk
      do j=1,no
        ex=ex+.5d0*system_exchange_fraction(system)*system%rocc(j,info%ik_s+ik-1,1)*system%wtk(info%ik_s+ik-1)*system%hvol &
          *real(sum(conjg(local(:,j,ik))*w(:,j,ik)),8)
      enddo
    enddo
    call comm_summation(ex,exx_state_owner%exx_exchange_energy,info%icomm_k)
  end subroutine

  subroutine apply_wannier_collective(target,action,info,system)
    use communication, only: comm_bcast
    implicit none
    type(s_dft_system),intent(in) :: system
    complex(8),intent(in) :: target(:,:,:)
    complex(8),intent(out) :: action(:,:,:)
    type(s_parallel_info),intent(in) :: info
    complex(8),allocatable :: send(:,:,:),alltarget(:,:,:),allaction(:,:,:)
    integer :: ng,nt,nk,status
    ng=size(target,1);nt=size(target,2);nk=product(system%kmesh)
    allocate(send(ng,nt,nk),alltarget(ng,nt,nk),allaction(ng,nt,nk))
    send=0d0;send(:,:,info%ik_s:info%ik_e)=target
    call comm_summation(send,alltarget,size(send),info%icomm_k)
    status=0
    if(info%id_k==0)call wannier_apply(wannier,alltarget,allaction,status)
    call comm_bcast(status,info%icomm_k,0)
    if(status/=0)error stop 'HSE Wannier: full trial action failed'
    call comm_bcast(allaction,info%icomm_k,0)
    action=allaction(:,:,info%ik_s:info%ik_e)
  end subroutine

  subroutine exx_core_exchange(system,mg,info,psi,core,energy,rt_context)
    implicit none
    type(s_exx_rt_state),target,optional,intent(inout)::rt_context
    type(s_exx_rt_state),pointer::exx_state_owner
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: psi
    integer,intent(in) :: core(3)
    real(8),intent(out) :: energy
    real(8) :: local,prefactor
    integer :: ix,iy,iz,io,ik,g
    exx_state_owner=>resolve_rt_context(rt_context,system)
    if(.not.allocated(exx_state_owner%cached_action))error stop 'DC HSE: refreshed exchange action missing'
    prefactor=.5d0*system_exchange_fraction(system)*system%hvol
    local=0d0
!$omp parallel do collapse(4) default(none) schedule(static) &
!$omp private(ik,io,iz,iy,ix,g) shared(info,mg,core,system,psi,exx_state_owner,prefactor) reduction(+:local)
    do ik=info%ik_s,info%ik_e
      do io=info%io_s,info%io_e
        do iz=mg%is(3),min(mg%ie(3),core(3))
          do iy=mg%is(2),min(mg%ie(2),core(2))
            do ix=mg%is(1),min(mg%ie(1),core(1))
              g=1+(ix-mg%is(1))+mg%num(1)*((iy-mg%is(2))+mg%num(2)*(iz-mg%is(3)))
              local=local+prefactor*system%rocc(io,ik,1)*system%wtk(ik) &
                *real(conjg(psi%zwf(ix,iy,iz,1,io,ik,1))*exx_state_owner%cached_action(g,io-info%io_s+1,ik-info%ik_s+1),8)
            enddo
          enddo
        enddo
      enddo
    enddo
!$omp end parallel do
    call comm_summation(local,energy,info%icomm_rko)
  end subroutine

  pure logical function finite_real_3d(values) result(finite)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in) :: values(:,:,:)
    real(8) :: value
    integer :: i,j,k
    finite=.false.
    do k=1,size(values,3)
      do j=1,size(values,2)
        do i=1,size(values,1)
          value=values(i,j,k)
          if(.not.ieee_is_finite(value))return
        enddo
      enddo
    enddo
    finite=.true.
  end function
end module
