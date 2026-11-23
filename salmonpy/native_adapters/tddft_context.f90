! Copyright 2026 SALMONpy contributors. Apache-2.0.
! Internal serial/k-only MPI development context; not a stable or failure-atomic public ABI.
! Global input configuration must remain common and immutable. Never copy contexts.
#include "config.h"
module salmonpy_tddft_context
 use filesystem,only:resolve_path
 use structures
 use salmonpy_gs_seed,only:s_gs_seed
 use salmonpy_step_output
 use hse_checkpoint,only:numeric_payload,hse_checkpoint_install,hse_checkpoint_capture
 use exx_rt_state,only:s_exx_rt_state,rt_state_clear
 use exx_native,only:exx_bind_system,exx_unbind_system
 use initialization_rt_sub,only:initialization_rt
 use salmon_xc,only:finalize_xc
 use salmon_global
 use parallelization,only:nproc_size_global,nproc_group_global,nproc_id_global
 use communication,only:comm_summation,comm_bcast
#ifdef USE_MPI
 use communication,only:comm_begin_children,comm_end_children,comm_release_children
#endif
 use exx_native,only:exx_validate_vector_potential
 use sym_sub,only:use_symmetry,SymMatA,SymMatB,Amat,Bmat,Ainv,Binv
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 public::rt_model_settings
 public::rt_history_info,rt_step_transaction,rt_input_config,tddft_context,rt_open,rt_step,rt_close,rt_probe_write,rt_checkpoint,rt_accepted_step,rt_current,rt_vector_potential
 type fixed_settings
  real(8)::omega,electronic_temperature
  integer::grid(3),mesh(3),bands,energy_interval,electrons,hamil,exact_interval,warmup,teachers
  character(32)::functional,propagation,factor
  character(1)::predictor
 end type
 type,public::s_rt_model_settings
  private
  type(fixed_settings)::values
  logical::pinned=.false.
 end type
type :: checkpoint_payload
  type(s_orbital) :: spsi_in,spsi_out,tpsi
  type(s_scalar) :: rho,rho_jm,Vh,Vh_stock1,Vh_stock2,Vbox,Vpsl
  type(s_scalar),allocatable :: rho_s(:),V_local(:),Vxc(:)
  type(s_rt) :: rt
  type(s_poisson) :: poisson
  type(s_dft_energy) :: energy
  real(8) :: fields(3,5)
  integer :: Mit,accepted
end type
 type tddft_context
 private
 integer::mpi_child_owner=0
type(s_rgrid) :: lg
type(s_rgrid) :: mg
type(s_dft_system)  :: system
type(s_rt) :: rt
type(s_parallel_info) :: info
type(s_poisson) :: poisson
type(s_stencil) :: stencil
type(s_xc_functional) :: xc_func
type(s_reciprocal_grid) :: fg
type(s_ewald_ion_ion) :: ewald
type(s_dft_energy) :: energy
type(s_md) :: md
type(s_ofile) :: ofl
type(s_scalar) :: Vpsl
type(s_scalar) :: rho,rho_jm,Vh,Vh_stock1,Vh_stock2,Vbox
type(s_scalar),allocatable :: rho_s(:),V_local(:),Vxc(:)
type(s_orbital) :: spsi_in,spsi_out
type(s_orbital) :: tpsi ! temporary wavefunctions
type(s_sendrecv_grid) :: srg,srg_scalar
type(s_pp_info) :: pp
type(s_pp_grid) :: ppg
type(s_pp_nlcc) :: ppn
type(s_singlescale) :: singlescale
type(s_unfold) :: unfold

 type(s_exx_rt_state)::hse
 type(s_rt_run_config)::run_settings
 logical::checkpoint_failed=.false.
 type(fixed_settings)::settings
 ! Protected numerical buffers and staged restore candidates belong to this
 ! context. Never copy owning FFT, MPI, XC-library or output handles here.
 type(checkpoint_payload),allocatable::step_saved,pending
 type(numeric_payload),allocatable::step_exchange,pending_hse
 type(s_xc_operator_payload),allocatable::step_xc
 type(step_output_state)::step_outputs
 real(8)::step_cval=0d0
 logical::step_restored=.false.
 integer::Mit=0,accepted=0
 end type
contains
 subroutine rt_model_settings(reset,owner)
 logical,intent(in)::reset
 type(s_rt_model_settings),intent(inout)::owner
 owner%pinned=.not.reset
 if(.not.reset)owner%values=capture_settings()
 end subroutine
 subroutine rt_step_transaction(ctx,mode,status)
 use communication,only:comm_get_max
 type(tddft_context),pointer,intent(inout)::ctx
 integer,intent(in)::mode
 integer,intent(out)::status
 integer::reference
 status=0;reference=mode
 call comm_bcast(reference,nproc_group_global)
 if(mode/=reference.or.mode<0.or.mode>3)status=11
 if(.not.associated(ctx))status=12
 call comm_get_max(status,nproc_group_global)
 if(status/=0)return
 if(mode==0)then
   if(allocated(ctx%step_saved).or.allocated(ctx%pending))status=14
   if(yn_out_dns_rt/='n'.or.yn_out_dns_ac_je/='n'.or.yn_out_micro_je/='n')status=15
    if(yn_out_elf_rt/='n'.or.yn_out_estatic_rt/='n'.or.yn_out_current_decomposed/='n')status=15
    if(projection_option/='no'.or.dm_unfold_option/='no'.or.yn_out_rvf_rt/='n')status=15
 else if(mode/=2)then
   if(.not.allocated(ctx%step_saved))status=14
   if(mode==1.and.ctx%step_restored)status=14
 endif
 call comm_get_max(status,nproc_group_global)
 if(status/=0)return
 select case(mode)
 case(0)
   allocate(ctx%step_saved,ctx%step_exchange,ctx%step_xc)
   call hse_checkpoint_capture(ctx%hse,ctx%step_exchange,status)
   call comm_get_max(status,nproc_group_global)
   if(status==0)then
     call checkpoint_copy(ctx,ctx%step_saved,.false.)
     ctx%step_xc=ctx%system%xc_payload;ctx%step_cval=ctx%xc_func%cval
     if(nproc_id_global==0)call step_output_capture(ctx%step_outputs,ctx%ofl,status)
     call comm_get_max(status,nproc_group_global)
   endif
   if(status/=0)then
     deallocate(ctx%step_saved,ctx%step_exchange,ctx%step_xc)
   else
     ctx%step_restored=.false.
   endif
 case(1)
   call hse_checkpoint_install(ctx%hse,ctx%step_exchange)
   call checkpoint_copy(ctx,ctx%step_saved,.true.)
   ctx%system%xc_payload=ctx%step_xc;ctx%xc_func%cval=ctx%step_cval
   block
   use hamiltonian,only:update_kvector_nonlocalpt
   call update_kvector_nonlocalpt(ctx%info%ik_s,ctx%info%ik_e,ctx%system,ctx%ppg)
   end block
   ctx%step_restored=.true.
   if(nproc_id_global==0)call step_output_close(ctx%step_outputs,status)
 case(2)
   if(allocated(ctx%step_saved))then
     if(nproc_id_global==0)then
       call step_output_reopen(ctx%step_outputs,status,ctx%ofl)
       if(status==0)call step_output_flush(ctx%step_outputs,status)
     endif
     call comm_get_max(status,nproc_group_global)
     ! Retain ownership metadata if reconnect/flush failed. Old closed units
     ! may belong to another context; retry or close must still skip them.
     if(status==0)then
       deallocate(ctx%step_saved,ctx%step_exchange,ctx%step_xc)
       ctx%step_restored=.false.
     endif
   endif
 case(3)
   if(nproc_id_global==0)call step_output_reopen(ctx%step_outputs,status,ctx%ofl)
 end select
 call comm_get_max(status,nproc_group_global)
 end subroutine
 subroutine rt_input_config(value)
 type(s_rt_run_config),intent(out)::value
 value%delta_t=dt;value%max_steps=nt
 value%stencil_mode=stencil_openmp_mode
 value%current_mode=current_openmp_mode
 value%force_mode=force_openmp_mode
 end subroutine
 logical function valid_mode(mode)
 character(*),intent(in)::mode
 valid_mode=mode=='auto'.or.mode=='orbital'.or.mode=='rgrid'
 end function
 function capture_settings() result(value)
 type(fixed_settings)::value
 value%grid=num_rgrid;value%mesh=num_kgrid
 value%omega=hse_omega;value%electronic_temperature=temperature
 value%electrons=nelec;value%hamil=n_hamil;value%factor=exx_factor_mode
 value%exact_interval=exx_factor_exact_interval;value%warmup=exx_factor_warmup_steps;value%teachers=exx_factor_n_teach
 value%predictor=yn_predictor_corrector
 value%energy_interval=out_rt_energy_step
 value%bands=nstate;value%functional=xc;value%propagation=propagator
 end function
 logical function settings_match(value) result(same)
 type(fixed_settings),intent(in)::value
 same=hse_omega==value%omega.and.temperature==value%electronic_temperature.and.nelec==value%electrons.and. &
      n_hamil==value%hamil.and.exx_factor_mode==value%factor.and. &
      exx_factor_exact_interval==value%exact_interval.and.exx_factor_warmup_steps==value%warmup.and. &
      exx_factor_n_teach==value%teachers.and.yn_predictor_corrector==value%predictor.and. &
      out_rt_energy_step==value%energy_interval.and.all(num_rgrid==value%grid).and. &
      all(num_kgrid==value%mesh).and.nstate==value%bands.and.xc==value%functional.and.propagator==value%propagation
 end function
 logical function context_symmetry_matches(ctx) result(ok)
 type(tddft_context),intent(in)::ctx
 ok=.false.
 ! Density/current symmetry still uses immutable process-wide tables. Until
 ! those are extracted, refuse drift rather than combining different groups.
 if(use_symmetry.neqv.ctx%system%kpoint_input%use_symmetry)return
 if(.not.use_symmetry)then
  ok=.true.;return
 endif
 if(.not.allocated(SymMatA).or..not.allocated(SymMatB))return
 if(.not.allocated(ctx%system%kpoint_input%symmat_a).or..not.allocated(ctx%system%kpoint_input%symmat_b))return
 if(any(shape(SymMatA)/=shape(ctx%system%kpoint_input%symmat_a)))return
 if(any(shape(SymMatB)/=shape(ctx%system%kpoint_input%symmat_b)))return
 if(any(lbound(SymMatA)/=1).or.any(lbound(SymMatB)/=1))return
 if(any(SymMatA/=ctx%system%kpoint_input%symmat_a))return
 if(any(SymMatB/=ctx%system%kpoint_input%symmat_b))return
 if(any(Amat/=ctx%system%primitive_a).or.any(Bmat/=ctx%system%primitive_b))return
 if(any(Ainv/=transpose(ctx%system%primitive_b)/(2d0*acos(-1d0))))return
 if(any(Binv/=transpose(ctx%system%primitive_a)/(2d0*acos(-1d0))))return
 ok=.true.
 end function
 subroutine rt_open(ctx,name,status,delta_t,openmp_mode,config,allow_symmetry,memory_gs,model_settings)
 type(tddft_context),pointer,intent(inout)::ctx
 character(*),intent(in)::name
 integer,intent(out)::status
 character(*),optional,intent(in)::openmp_mode
 real(8),optional,intent(in)::delta_t
 type(s_rt_run_config),optional,intent(in)::config
 type(s_gs_seed),optional,intent(in)::memory_gs
 logical,optional,intent(in)::allow_symmetry
 type(s_rt_model_settings),intent(inout)::model_settings
 type(s_rt_model_settings)::owner_settings
 type(s_rt_run_config)::selected
 character(len(sysname))::saved_name
 logical::exists
 integer::i,child_owner
 owner_settings=model_settings
 status=1
#if defined(USE_OPENACC) || defined(USE_FFTW)
 goto 900
#endif
 if(theory/='tddft_response'.or.out_rt_energy_step/=1)goto 900
 if(spin/='unpolarized'.or.temperature>=0d0.or.any(hse_q_reduction/=1))goto 900
 if(associated(ctx))goto 900
 if(present(config))then
  if(present(delta_t).or.present(openmp_mode))goto 900
  selected=config
 else
  call rt_input_config(selected)
  if(present(delta_t))selected%delta_t=delta_t
  if(present(openmp_mode))then
   if(.not.valid_mode(openmp_mode))goto 900
   selected%stencil_mode=openmp_mode;selected%current_mode=openmp_mode;selected%force_mode=openmp_mode
  endif
 endif
 if(.not.ieee_is_finite(selected%delta_t))goto 900
 if(selected%delta_t<=0d0)goto 900
 if(selected%max_steps<1.or.selected%max_steps>=huge(0))goto 900
 if(.not.valid_mode(selected%stencil_mode).or..not.valid_mode(selected%current_mode).or. &
    .not.valid_mode(selected%force_mode))goto 900
 if(owner_settings%pinned)then
  if(.not.settings_match(owner_settings%values))then
   status=2;goto 900
  endif
 endif
 if(nproc_size_global>1)then
  if(nproc_ob/=1.or.any(nproc_rgrid/=1).or.nproc_k/=nproc_size_global)goto 900
 endif
 if(xc/='hse06'.or.yn_periodic/='y'.or.yn_spinorbit/='n'.or. &
    yn_hse_wannier/='n'.or.yn_md/='n'.or.yn_dc/='n'.or.trans_longi/='tr')goto 900
 if(use_symmetry)then
  ! Internal opt-in only. Existing worker/C/Python callers keep the default prohibition.
  if(.not.present(allow_symmetry))goto 900
  if(.not.allow_symmetry)goto 900
  if(nproc_size_global<1.or.nproc_size_global>2)goto 900
  if(exx_factor_mode/='off'.and.exx_factor_mode/='strict'.and.exx_factor_mode/='learned')goto 900
  if(propagator/='hse_taylor4'.and.propagator/='hse_taylor4_full'.and.propagator/='hse_ptcn')goto 900
  if(propagator=='hse_taylor4_full'.and.exx_factor_mode/='off')goto 900
 endif
 if(yn_restart/='n'.or.e_impulse/=0d0.or.ae_shape1/='impulse'.or.ae_shape2/='none')goto 900
 if(yn_out_current_decomposed/='n'.or.yn_out_intraband_current/='n'.or.yn_out_dns_ac_je/='n'.or. &
    yn_out_mag_decomposed_rt/='n'.or.yn_out_mag_micro_rt/='n'.or.yn_out_spin_current_decomposed/='n'.or. &
    yn_out_spin_current_micro/='n'.or.yn_out_mom_distr_rt/='n'.or.yn_out_rt_energy_components/='n')goto 900
 if(projection_option/='no'.or.dm_unfold_option/='no'.or.yn_out_rvf_rt/='n')goto 900
 if(len_trim(name)<1.or.len_trim(name)>len(sysname)-32)goto 900
 do i=1,len_trim(name)
  if(index('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-',name(i:i))==0)goto 900
 enddo
 inquire(file=resolve_path(trim(base_directory)//trim(name)//'_rt.data'),exist=exists)
 if(exists)goto 900
 inquire(file=resolve_path(trim(base_directory)//trim(name)//'_rt_energy.data'),exist=exists)
 if(exists)goto 900
 if(yn_out_dns_rt/='n'.or.yn_out_micro_je/='n'.or.yn_out_elf_rt/='n'.or.yn_out_estatic_rt/='n')goto 900
 status=0
900 call checkpoint_collective_status(status)
 if(status/=0)return
 child_owner=0
#ifdef USE_MPI
 call comm_begin_children(child_owner,status)
#endif
 call checkpoint_collective_status(status)
 if(status/=0)then
#ifdef USE_MPI
  call comm_end_children(child_owner)
#endif
#ifdef USE_MPI
  call comm_release_children(child_owner)
#endif
  return
 endif
 allocate(ctx);ctx%mpi_child_owner=child_owner
 call exx_bind_system(ctx%system,ctx%hse,status)
 call checkpoint_collective_status(status)
 if(status/=0)then
#ifdef USE_MPI
  call comm_end_children(child_owner)
#endif
#ifdef USE_MPI
  call comm_release_children(child_owner)
#endif
  deallocate(ctx);return
 endif
 ! SYSname is a construction-only output label; no numerical setting is switched.
 saved_name=sysname;sysname=name
call initialization_rt( ctx%Mit, ctx%system, ctx%energy, ctx%ewald, ctx%rt, ctx%md, &
                        ctx%singlescale,  &
                        ctx%stencil, ctx%fg, ctx%poisson,  &
                        ctx%lg, ctx%mg,   &
                        ctx%info,  &
                        ctx%xc_func, ctx%ofl,  &
                        ctx%srg, ctx%srg_scalar,  &
                        ctx%spsi_in, ctx%spsi_out, ctx%tpsi, ctx%rho, ctx%rho_jm, ctx%rho_s,  &
                        ctx%V_local, ctx%Vbox, ctx%Vh, ctx%Vh_stock1, ctx%Vh_stock2, ctx%Vxc, ctx%Vpsl,&
                        ctx%pp, ctx%ppg, ctx%ppn, ctx%unfold, run_config=selected, memory_gs=memory_gs )
 sysname=saved_name
#ifdef USE_MPI
 call comm_end_children(child_owner)
#endif
 ctx%run_settings=selected
 ctx%settings=capture_settings()
 owner_settings%values=ctx%settings;owner_settings%pinned=.true.
 model_settings=owner_settings
 ctx%accepted=ctx%Mit
 status=0
 end subroutine
 subroutine rt_step(ctx,ac,current,total_energy,status,checked)
#ifdef USE_MPI
 use communication,only:comm_resume_children,comm_end_children,comm_get_max
#endif
 type(tddft_context),pointer,intent(inout)::ctx
 real(8),intent(in)::ac(3)
 real(8),intent(out)::current(3),total_energy
 integer,intent(out)::status
 logical,optional,intent(in)::checked
 logical::verify
#ifdef USE_MPI
 logical::entered
#endif
 verify=.true.;if(present(checked))verify=checked
 current=0d0;total_energy=0d0;status=0
 if(.not.associated(ctx))status=1
#ifdef USE_MPI
 if(verify)call comm_get_max(status,nproc_group_global)
#endif
 if(status/=0)return
#ifdef USE_MPI
 call comm_resume_children(ctx%mpi_child_owner,status)
 entered=status==0
 if(verify)call comm_get_max(status,nproc_group_global)
 if(status/=0)then
  if(entered)call comm_end_children(ctx%mpi_child_owner)
  return
 endif
#endif
 call rt_step_impl(ctx,ac,current,total_energy,status,verify)
#ifdef USE_MPI
 if(entered)call comm_end_children(ctx%mpi_child_owner)
#endif
 end subroutine
 subroutine rt_step_impl(ctx,ac,current,total_energy,status,checked)
 type(tddft_context),pointer,intent(inout)::ctx
 real(8),intent(in)::ac(3)
 real(8),intent(out)::current(3),total_energy
 integer,intent(out)::status
 logical,intent(in)::checked
 integer::itt,i,guard_status,field_status
 type(checkpoint_payload)::accepted_payload
 type(numeric_payload)::accepted_exchange
 type(s_xc_operator_payload)::accepted_xc
 real(8)::accepted_cval
 status=1;current=0d0;total_energy=0d0
 if(.not.associated(ctx))return
 guard_status=0
 if(checked)then
  if(.not.settings_match(ctx%settings).or..not.context_symmetry_matches(ctx))guard_status=2
 endif
 if(allocated(ctx%pending))guard_status=2
 if(ctx%checkpoint_failed)guard_status=2
 if(ctx%accepted>=ctx%rt%max_steps)guard_status=2
 if(checked)then
  do i=1,3
   if(.not.ieee_is_finite(ac(i)))guard_status=2
  enddo
 endif
 if(guard_status==0)then
  call exx_validate_vector_potential(ctx%system,ac,field_status)
  if(field_status==1)guard_status=2
  if(field_status==2)guard_status=34
 endif
 block
 use communication,only:comm_get_max
 status=guard_status
 if(checked)call comm_get_max(status,nproc_group_global)
 end block
 if(status/=0)return
 if(checked)then
 call hse_checkpoint_capture(ctx%hse,accepted_exchange,status)
 call checkpoint_collective_status(status)
 if(status/=0)return
 call checkpoint_copy(ctx,accepted_payload,.false.)
 accepted_xc=ctx%system%xc_payload;accepted_cval=ctx%xc_func%cval
 endif
 ctx%rt%recover_errors=checked;ctx%rt%step_error=0
 ctx%hse%recover_errors=checked;ctx%hse%step_error=0
 ctx%rt%factor_endpoint_before_output=.true.
 itt=ctx%accepted+1
 ctx%rt%Ac_ext(:,itt)=ac;ctx%rt%Ac_tot(:,itt)=ac
 if(mod(itt-ctx%Mit,2)==1)then
    call time_evolution_step(ctx%Mit,ctx%rt%max_steps,itt,ctx%lg,ctx%mg,ctx%system,ctx%rt,ctx%info,ctx%stencil,ctx%xc_func &
     & ,ctx%srg,ctx%srg_scalar,ctx%pp,ctx%ppg,ctx%ppn,ctx%spsi_in,ctx%spsi_out,ctx%tpsi,ctx%rho,ctx%rho_jm,ctx%rho_s,ctx%V_local,ctx%Vbox,ctx%Vh,ctx%Vh_stock1,ctx%Vh_stock2,ctx%Vxc &
     & ,ctx%Vpsl,ctx%fg,ctx%energy,ctx%ewald,ctx%md,ctx%ofl,ctx%poisson,ctx%singlescale,ctx%unfold)
 else
    call time_evolution_step(ctx%Mit,ctx%rt%max_steps,itt,ctx%lg,ctx%mg,ctx%system,ctx%rt,ctx%info,ctx%stencil,ctx%xc_func &
     & ,ctx%srg,ctx%srg_scalar,ctx%pp,ctx%ppg,ctx%ppn,ctx%spsi_out,ctx%spsi_in,ctx%tpsi,ctx%rho,ctx%rho_jm,ctx%rho_s,ctx%V_local,ctx%Vbox,ctx%Vh,ctx%Vh_stock1,ctx%Vh_stock2,ctx%Vxc &
     & ,ctx%Vpsl,ctx%fg,ctx%energy,ctx%ewald,ctx%md,ctx%ofl,ctx%poisson,ctx%singlescale,ctx%unfold)
 endif
 guard_status=ctx%hse%step_error
 if(guard_status==0)guard_status=ctx%rt%step_error
 block
 use communication,only:comm_get_max
 if(checked)call comm_get_max(guard_status,nproc_group_global)
 end block
 ctx%rt%factor_endpoint_before_output=.false.
 ctx%rt%recover_errors=.false.;ctx%rt%step_error=0
 ctx%hse%recover_errors=.false.;ctx%hse%step_error=0
 if(guard_status/=0)then
   if(.not.checked)then
     status=guard_status;return
   endif
   call hse_checkpoint_install(ctx%hse,accepted_exchange)
   call checkpoint_copy(ctx,accepted_payload,.true.)
   ctx%system%xc_payload=accepted_xc;ctx%xc_func%cval=accepted_cval
   block
   use hamiltonian,only:update_kvector_nonlocalpt
   call update_kvector_nonlocalpt(ctx%info%ik_s,ctx%info%ik_e,ctx%system,ctx%ppg)
   end block
   status=guard_status
   return
 endif
 ctx%accepted=itt
 current=-ctx%rt%curr(:,itt)
 total_energy=ctx%energy%E_tot
 status=0
 end subroutine
 ! Physical charge current at the accepted endpoint, without propagation.
 subroutine rt_current(ctx,current,status)
 use density_matrix,only:calc_current
 use hamiltonian,only:update_kvector_nonlocalpt
 type(tddft_context),pointer,intent(inout)::ctx
 real(8),intent(out)::current(3)
 integer,intent(out)::status
 real(8)::matter(3,2)
 integer::finite_i
 status=1;current=0d0
 if(.not.associated(ctx))return
 if(ctx%checkpoint_failed.or..not.settings_match(ctx%settings).or..not.context_symmetry_matches(ctx))return
 if(ctx%accepted>ctx%Mit)then
  current=-ctx%rt%curr(:,ctx%accepted)
 else
  ! Initialization has not evaluated the t=0 observable; query the initial GS.
  ! Nonlocal phases follow the same endpoint convention as rt_step.
  call update_kvector_nonlocalpt(ctx%info%ik_s,ctx%info%ik_e,ctx%system,ctx%ppg)
  call calc_current(ctx%system,ctx%mg,ctx%stencil,ctx%info,ctx%srg,ctx%spsi_in,ctx%ppg, &
                    matter(:,1:ctx%system%nspin))
  ctx%spsi_in%update_zwf_overlap=.true.
  current=-sum(matter(:,1:ctx%system%nspin),dim=2)
 endif
 do finite_i=1,3
  if(.not.ieee_is_finite(current(finite_i)))return
 enddo
 status=0
 end subroutine rt_current
 subroutine rt_vector_potential(ctx,A,status)
 type(tddft_context),pointer,intent(in)::ctx
 real(8),intent(out)::A(3)
 integer,intent(out)::status
 status=1;A=0d0
 if(.not.associated(ctx))return
 if(ctx%checkpoint_failed)return
 A=ctx%rt%Ac_tot(:,ctx%accepted)
 status=0
 end subroutine
 subroutine rt_close(ctx)
 use sendrecv_grid,only:dealloc_cache
 type(tddft_context),pointer,intent(inout)::ctx
 if(.not.associated(ctx))return
 ! A rolled-back transaction may have closed these units already; another
 ! context can now own their old numbers. Only close our still-connected streams.
 if(nproc_id_global==0)then
  if(allocated(ctx%step_saved))then
   block
   integer::output_status
   call step_output_close(ctx%step_outputs,output_status)
   end block
  else
   close(ctx%ofl%fh_rt);close(ctx%ofl%fh_rt_energy)
  endif
 endif
 call dealloc_cache(ctx%srg)
 call dealloc_cache(ctx%srg_scalar)
 call finalize_xc(ctx%xc_func)
 call exx_unbind_system(ctx%system)
 call rt_state_clear(ctx%hse)
#ifdef USE_MPI
 call comm_release_children(ctx%mpi_child_owner)
#endif
 deallocate(ctx)
 nullify(ctx)
 end subroutine
 subroutine rt_history_info(ctx,values,status)
 ! Rank-local read-only diagnostics; no propagation, MPI or checkpoint mutation.
 type(tddft_context),pointer,intent(in)::ctx
 integer,intent(out)::values(8),status
 values=0;status=1
 if(.not.associated(ctx))return
 if(ctx%checkpoint_failed)return
 if(allocated(ctx%hse%factor_history))then
   values(1)=size(ctx%hse%factor_history)
   values(2)=count(ctx%hse%factor_history%ready)
 endif
 if(allocated(ctx%hse%ptcn_history))then
   values(3)=size(ctx%hse%ptcn_history)
   values(4)=count(ctx%hse%ptcn_history%ready)
 endif
 values(5)=merge(1,0,ctx%hse%factor_prediction_installed)
 values(6)=merge(1,0,ctx%hse%ptcn_predicted)
 values(7)=ctx%hse%factor_step;values(8)=ctx%hse%ptcn_accepted_step
 status=0
 end subroutine
 subroutine rt_probe_write(ctx,unit)
 ! Synthetic verification payload only: not a restart or a user GS archive.
 type(tddft_context),pointer,intent(in)::ctx
 integer,intent(in)::unit
 integer::k
 write(unit)ctx%accepted,ctx%rt%curr(:,ctx%accepted),ctx%energy%E_tot,ctx%energy%E_xc
 write(unit)ctx%spsi_in%zwf,ctx%spsi_out%zwf
 write(unit)ctx%rho%f,ctx%Vh%f,ctx%Vh_stock1%f,ctx%Vh_stock2%f
 do k=1,size(ctx%Vxc)
  write(unit)ctx%Vxc(k)%f,ctx%V_local(k)%f
 enddo
 write(unit)ctx%hse%exx_exchange_energy,ctx%hse%factor_step,ctx%hse%ptcn_accepted_step
 if(allocated(ctx%hse%ace%factors))write(unit)ctx%hse%ace%factors
 if(allocated(ctx%hse%cached_source))write(unit)ctx%hse%cached_source
 if(allocated(ctx%hse%factor_history))then
  do k=1,size(ctx%hse%factor_history)
   if(allocated(ctx%hse%factor_history(k)%x))write(unit)ctx%hse%factor_history(k)%x
  enddo
 endif
 if(allocated(ctx%hse%ptcn_history))then
  do k=1,size(ctx%hse%ptcn_history)
   if(allocated(ctx%hse%ptcn_history(k)%x))write(unit)ctx%hse%ptcn_history(k)%x
  enddo
 endif
 end subroutine

 ! Same-build, same-input rank-local native archive. Caller verifies the binary/GS
 ! provenance and authenticates the file before load. Opaque context, FFT, XC,
 ! MPI and output handles are never serialized or intrinsically copied.
 ! Restore only into a fresh rt_open; no propagation is performed. Output files
 ! are fresh streams, not a resumed text log. On a late restore error discard
 ! the context (HSE owner may have been partly loaded).
 subroutine rt_checkpoint(ctx,path,restore,status,stage)
 use hse_checkpoint,only:hse_checkpoint_write,hse_checkpoint_read
 use hamiltonian,only:update_kvector_nonlocalpt
 use exx_native,only:exx_refresh
 type(tddft_context),pointer,intent(inout)::ctx
 character(*),intent(in)::path
 logical,intent(in)::restore
 integer,intent(out)::status
 integer,optional,intent(in)::stage
 type(checkpoint_payload),allocatable::copy
 type(fixed_settings)::saved_settings
 type(s_rt_run_config)::saved_run
 integer::u,ios,closing,allocation_status,header(32),expected(32)
 integer::hse_position,need_kernel,finite_i,finite_j
 integer(8)::checkpoint_array_limit
 logical::opened
 real(8)::cell(3,3),spacing(3),saved_initial_energy
 character(32)::magic
 character(1)::extra
 logical::hse_loaded
 status=1;hse_loaded=.false.
 ios=0
 if(.not.associated(ctx))ios=1
 call checkpoint_collective_status(ios)
 if(ios/=0)return
 if(allocated(ctx%step_saved).and.(restore.or.present(stage)))ios=14
 call checkpoint_collective_status(ios)
 if(ios/=0)return
 if(.not.present(stage).and.allocated(ctx%pending))return
 if(present(stage))then
  if(stage==4)then
   if(allocated(ctx%pending))deallocate(ctx%pending)
   if(allocated(ctx%pending_hse))deallocate(ctx%pending_hse)
   status=0;return
  endif
  if(stage==3)then
   if(.not.allocated(ctx%pending).or..not.allocated(ctx%pending_hse))return
   ! Recheck the immutable process contract at commit, not just while staging.
   ios=0
   if(.not.settings_match(ctx%settings).or..not.context_symmetry_matches(ctx))ios=1
   call checkpoint_collective_status(ios)
   if(ios/=0)return
   call hse_checkpoint_install(ctx%hse,ctx%pending_hse)
   call checkpoint_copy(ctx,ctx%pending,.true.)
   call update_kvector_nonlocalpt(ctx%info%ik_s,ctx%info%ik_e,ctx%system,ctx%ppg)
   ctx%spsi_in%update_zwf_overlap=.true.
   ctx%spsi_out%update_zwf_overlap=.true.
   ctx%tpsi%update_zwf_overlap=.true.
   deallocate(ctx%pending,ctx%pending_hse);status=0;return
  endif
  if(stage/=2.or.allocated(ctx%pending).or.allocated(ctx%pending_hse))return
 endif
 ios=0
 if(.not.associated(ctx))ios=1
 call checkpoint_collective_status(ios)
 if(ios/=0)return
 if(ctx%checkpoint_failed.or..not.settings_match(ctx%settings).or..not.context_symmetry_matches(ctx))ios=1
 if(len_trim(path)==0)ios=1
 if(restore.and..not.present(stage).and.ctx%accepted/=ctx%Mit)ios=1
 if(xc/='hse06'.or.yn_spinorbit/='n'.or.yn_dc/='n'.or.yn_md/='n')ios=1
 if(projection_option/='no'.or.dm_unfold_option/='no'.or.ctx%singlescale%flag_use)ios=1
 if(ctx%info%isize_r/=1.or.ctx%info%isize_o/=1.or.ctx%info%isize_k/=nproc_size_global)ios=1
 call checkpoint_collective_status(ios)
 if(ios/=0)return
 expected=[2,ctx%rt%max_steps,ctx%system%nspin,ctx%system%no,ctx%system%nk, &
    ctx%lg%num,ctx%mg%is,ctx%mg%ie,ctx%mg%is_array,ctx%mg%ie_array, &
    ctx%info%io_s,ctx%info%io_e,ctx%info%ik_s,ctx%info%ik_e, &
    storage_size(1d0),storage_size(1),storage_size((0d0,0d0)), &
    nproc_size_global,nproc_id_global,ctx%info%isize_r,ctx%info%isize_o,ctx%info%isize_k]
 checkpoint_array_limit=max(64_8*(int(ctx%rt%max_steps,8)+1_8),64_8*int(ctx%system%ngrid,8))
 if(allocated(ctx%spsi_in%zwf))checkpoint_array_limit=max(checkpoint_array_limit,size(ctx%spsi_in%zwf,kind=8))
 if(allocated(ctx%spsi_in%rwf))checkpoint_array_limit=max(checkpoint_array_limit,size(ctx%spsi_in%rwf,kind=8))
 allocate(copy,stat=allocation_status)
 call checkpoint_collective_status(allocation_status)
 if(allocation_status/=0)return
 call checkpoint_copy(ctx,copy,.false.)
 if(.not.restore.and.ctx%system%kpoint_input%use_symmetry)then
  ! Overlap caches are not serialized. Record their invalid state explicitly in
  ! the symmetry profile; never mutate the live flags or change legacy V1 bytes.
  copy%spsi_in%update_zwf_overlap=.true.
  copy%spsi_out%update_zwf_overlap=.true.
  copy%tpsi%update_zwf_overlap=.true.
 endif
 if(restore)then
  open(newunit=u,file=resolve_path(trim(path)),access='stream',form='unformatted',status='old',action='read',iostat=ios)
 else
  open(newunit=u,file=resolve_path(trim(path)),access='stream',form='unformatted',status='new',action='write',iostat=ios)
 endif
 opened=ios==0
 call checkpoint_collective_status(ios)
 if(ios/=0)then
  if(opened)close(u)
  return
 endif
 magic='SALMONPY_FAST_HSE_RT_V1';header=expected
 saved_settings=ctx%settings;saved_run=ctx%run_settings
 cell=ctx%system%primitive_a;spacing=ctx%system%hgs;saved_initial_energy=ctx%energy%E_tot0
 if(restore)then
  read(u,iostat=ios)magic,header,saved_settings,saved_run,cell,spacing,saved_initial_energy
  if(ios==0)then
   if(magic/='SALMONPY_FAST_HSE_RT_V1'.or.any(header/=expected))ios=1
   if(.not.settings_match(saved_settings))ios=1
   if(saved_run%delta_t/=ctx%run_settings%delta_t.or.saved_run%max_steps/=ctx%run_settings%max_steps)ios=1
   if(saved_run%stencil_mode/=ctx%run_settings%stencil_mode.or. &
      saved_run%current_mode/=ctx%run_settings%current_mode.or.saved_run%force_mode/=ctx%run_settings%force_mode)ios=1
   if(any(cell/=ctx%system%primitive_a).or.any(spacing/=ctx%system%hgs))ios=1
   if(.not.ieee_is_finite(saved_initial_energy).or.saved_initial_energy/=ctx%energy%E_tot0)ios=1
  endif
 else
  write(u,iostat=ios)magic,header,saved_settings,saved_run,cell,spacing,saved_initial_energy
 endif
 if(ios==0)call checkpoint_io(u,restore,copy,ios,checkpoint_array_limit)
 if(ios==0.and.restore)then
  if(copy%Mit/=ctx%Mit.or.copy%accepted<copy%Mit.or.copy%accepted>ctx%rt%max_steps)ios=1
  do finite_j=1,4
  do finite_i=1,3
   if(.not.ieee_is_finite(copy%fields(finite_i,finite_j)))ios=1
  enddo
  enddo
  if(copy%energy%E_tot0/=saved_initial_energy)ios=1
 endif
 call checkpoint_collective_status(ios)
 if(ios==0)then
  if(restore)then
   ! Rebuild the fresh CPU kernel only if initialization did not need exchange.
   ! This prepares a static operator from the initial GS, never an RT step.
   need_kernel=0
   if(.not.allocated(ctx%hse%native_kernel))need_kernel=1
   call checkpoint_collective_status(need_kernel)
   if(need_kernel/=0.and.present(stage))ios=1
   if(need_kernel/=0.and..not.present(stage))then
    call exx_refresh(ctx%system,ctx%mg,ctx%info,ctx%spsi_in,rt_context=ctx%hse,rt_status=ios)
   endif
   call checkpoint_collective_status(ios)
   if(ios==0)then
    inquire(unit=u,pos=hse_position)
    call hse_checkpoint_read(u,ctx%hse,ios,staged=ctx%pending_hse)
   endif
  else
   call hse_checkpoint_write(u,ctx%hse,ios)
  endif
 endif
 if(ios==0)then
  magic='SALMONPY_FAST_HSE_RT_COMPLETE_1'
  if(restore)then
   read(u,iostat=ios)magic
   if(ios==0.and.magic/='SALMONPY_FAST_HSE_RT_COMPLETE_1')ios=1
   if(ios==0)then
    read(u,iostat=ios)extra
    if(is_iostat_end(ios))then
     ios=0
    else
     ios=1
    endif
   endif
  else
   write(u,iostat=ios)magic
  endif
 endif
 ! Every shard and both archive footers must validate before any owner commits.
 call checkpoint_collective_status(ios)

 close(u,iostat=closing)
 if(closing/=0)ios=closing
 call checkpoint_collective_status(ios)
 if(ios/=0)then
  if(allocated(ctx%pending_hse))deallocate(ctx%pending_hse)
  return
 endif
 if(restore.and.present(stage))then
  call move_alloc(copy,ctx%pending)
 else if(restore)then
  call hse_checkpoint_install(ctx%hse,ctx%pending_hse)
  deallocate(ctx%pending_hse)
  call checkpoint_copy(ctx,copy,.true.)
  call update_kvector_nonlocalpt(ctx%info%ik_s,ctx%info%ik_e,ctx%system,ctx%ppg)
  ctx%spsi_in%update_zwf_overlap=.true.
  ctx%spsi_out%update_zwf_overlap=.true.
  ctx%tpsi%update_zwf_overlap=.true.
 endif
 status=0
 end subroutine rt_checkpoint

 subroutine checkpoint_collective_status(status)
 integer,intent(inout)::status
 integer::local_failure,total_failure
 local_failure=0
 if(status/=0)local_failure=1
 call comm_summation(local_failure,total_failure,nproc_group_global)
 status=total_failure
 end subroutine

 integer function rt_accepted_step(ctx) result(step)
 type(tddft_context),pointer,intent(in)::ctx
 step=-1
 if(associated(ctx))then
  if(.not.ctx%checkpoint_failed)step=ctx%accepted
 endif
 end function rt_accepted_step

subroutine checkpoint_copy(ctx,copy,restore)
  type(tddft_context),intent(inout)::ctx
  type(checkpoint_payload),intent(inout)::copy
  logical,intent(in)::restore
  if(restore)then
    ctx%spsi_in=copy%spsi_in
    ctx%spsi_out=copy%spsi_out
    ctx%tpsi=copy%tpsi
    ctx%rho=copy%rho
    ctx%rho_jm=copy%rho_jm
    ctx%Vh=copy%Vh
    ctx%Vh_stock1=copy%Vh_stock1
    ctx%Vh_stock2=copy%Vh_stock2
    ctx%Vbox=copy%Vbox
    ctx%Vpsl=copy%Vpsl
    if(allocated(copy%rho_s))then
      ctx%rho_s=copy%rho_s
    else
      if(allocated(ctx%rho_s))deallocate(ctx%rho_s)
    endif
    if(allocated(copy%V_local))then
      ctx%V_local=copy%V_local
    else
      if(allocated(ctx%V_local))deallocate(ctx%V_local)
    endif
    if(allocated(copy%Vxc))then
      ctx%Vxc=copy%Vxc
    else
      if(allocated(ctx%Vxc))deallocate(ctx%Vxc)
    endif
    ctx%energy=copy%energy
    ctx%poisson=copy%poisson
    ctx%Mit=copy%Mit
    ctx%accepted=copy%accepted
    ctx%rt%Dp0_e=copy%rt%Dp0_e
    if(allocated(copy%rt%dDp_e))then
      ctx%rt%dDp_e=copy%rt%dDp_e
    else
      if(allocated(ctx%rt%dDp_e))deallocate(ctx%rt%dDp_e)
    endif
    if(allocated(copy%rt%Dp_e))then
      ctx%rt%Dp_e=copy%rt%Dp_e
    else
      if(allocated(ctx%rt%Dp_e))deallocate(ctx%rt%Dp_e)
    endif
    if(allocated(copy%rt%Dp_i))then
      ctx%rt%Dp_i=copy%rt%Dp_i
    else
      if(allocated(ctx%rt%Dp_i))deallocate(ctx%rt%Dp_i)
    endif
    if(allocated(copy%rt%Qp_e))then
      ctx%rt%Qp_e=copy%rt%Qp_e
    else
      if(allocated(ctx%rt%Qp_e))deallocate(ctx%rt%Qp_e)
    endif
    if(allocated(copy%rt%rIe))then
      ctx%rt%rIe=copy%rt%rIe
    else
      if(allocated(ctx%rt%rIe))deallocate(ctx%rt%rIe)
    endif
    if(allocated(copy%rt%curr))then
      ctx%rt%curr=copy%rt%curr
    else
      if(allocated(ctx%rt%curr))deallocate(ctx%rt%curr)
    endif
    if(allocated(copy%rt%E_ext))then
      ctx%rt%E_ext=copy%rt%E_ext
    else
      if(allocated(ctx%rt%E_ext))deallocate(ctx%rt%E_ext)
    endif
    if(allocated(copy%rt%E_ind))then
      ctx%rt%E_ind=copy%rt%E_ind
    else
      if(allocated(ctx%rt%E_ind))deallocate(ctx%rt%E_ind)
    endif
    if(allocated(copy%rt%E_tot))then
      ctx%rt%E_tot=copy%rt%E_tot
    else
      if(allocated(ctx%rt%E_tot))deallocate(ctx%rt%E_tot)
    endif
    if(allocated(copy%rt%Ac_ext))then
      ctx%rt%Ac_ext=copy%rt%Ac_ext
    else
      if(allocated(ctx%rt%Ac_ext))deallocate(ctx%rt%Ac_ext)
    endif
    if(allocated(copy%rt%Ac_ind))then
      ctx%rt%Ac_ind=copy%rt%Ac_ind
    else
      if(allocated(ctx%rt%Ac_ind))deallocate(ctx%rt%Ac_ind)
    endif
    if(allocated(copy%rt%Ac_tot))then
      ctx%rt%Ac_tot=copy%rt%Ac_tot
    else
      if(allocated(ctx%rt%Ac_tot))deallocate(ctx%rt%Ac_tot)
    endif
    if(allocated(copy%rt%zc))then
      ctx%rt%zc=copy%rt%zc
    else
      if(allocated(ctx%rt%zc))deallocate(ctx%rt%zc)
    endif
    if(allocated(copy%rt%vloc_t))then
      ctx%rt%vloc_t=copy%rt%vloc_t
    else
      if(allocated(ctx%rt%vloc_t))deallocate(ctx%rt%vloc_t)
    endif
    if(allocated(copy%rt%vloc_new))then
      ctx%rt%vloc_new=copy%rt%vloc_new
    else
      if(allocated(ctx%rt%vloc_new))deallocate(ctx%rt%vloc_new)
    endif
    if(allocated(copy%rt%vloc_old))then
      ctx%rt%vloc_old=copy%rt%vloc_old
    else
      if(allocated(ctx%rt%vloc_old))deallocate(ctx%rt%vloc_old)
    endif
    if(allocated(copy%rt%rho0_s))then
      ctx%rt%rho0_s=copy%rt%rho0_s
    else
      if(allocated(ctx%rt%rho0_s))deallocate(ctx%rt%rho0_s)
    endif
    ctx%rt%vonf=copy%rt%vonf
    ctx%rt%j_e=copy%rt%j_e
    ctx%rt%E_old=copy%rt%E_old
    ctx%rt%rho0_loaded=copy%rt%rho0_loaded
    ctx%system%vec_Ac=copy%fields(:,1)
    ctx%system%vec_Ac_ext=copy%fields(:,2)
    ctx%system%vec_E=copy%fields(:,3)
    ctx%system%vec_E_ext=copy%fields(:,4)
  else
    copy%spsi_in=ctx%spsi_in
    copy%spsi_out=ctx%spsi_out
    copy%tpsi=ctx%tpsi
    copy%rho=ctx%rho
    copy%rho_jm=ctx%rho_jm
    copy%Vh=ctx%Vh
    copy%Vh_stock1=ctx%Vh_stock1
    copy%Vh_stock2=ctx%Vh_stock2
    copy%Vbox=ctx%Vbox
    copy%Vpsl=ctx%Vpsl
    if(allocated(ctx%rho_s))then
      copy%rho_s=ctx%rho_s
    else
      if(allocated(copy%rho_s))deallocate(copy%rho_s)
    endif
    if(allocated(ctx%V_local))then
      copy%V_local=ctx%V_local
    else
      if(allocated(copy%V_local))deallocate(copy%V_local)
    endif
    if(allocated(ctx%Vxc))then
      copy%Vxc=ctx%Vxc
    else
      if(allocated(copy%Vxc))deallocate(copy%Vxc)
    endif
    copy%energy=ctx%energy
    copy%poisson=ctx%poisson
    copy%Mit=ctx%Mit
    copy%accepted=ctx%accepted
    copy%rt%Dp0_e=ctx%rt%Dp0_e
    if(allocated(ctx%rt%dDp_e))then
      copy%rt%dDp_e=ctx%rt%dDp_e
    else
      if(allocated(copy%rt%dDp_e))deallocate(copy%rt%dDp_e)
    endif
    if(allocated(ctx%rt%Dp_e))then
      copy%rt%Dp_e=ctx%rt%Dp_e
    else
      if(allocated(copy%rt%Dp_e))deallocate(copy%rt%Dp_e)
    endif
    if(allocated(ctx%rt%Dp_i))then
      copy%rt%Dp_i=ctx%rt%Dp_i
    else
      if(allocated(copy%rt%Dp_i))deallocate(copy%rt%Dp_i)
    endif
    if(allocated(ctx%rt%Qp_e))then
      copy%rt%Qp_e=ctx%rt%Qp_e
    else
      if(allocated(copy%rt%Qp_e))deallocate(copy%rt%Qp_e)
    endif
    if(allocated(ctx%rt%rIe))then
      copy%rt%rIe=ctx%rt%rIe
    else
      if(allocated(copy%rt%rIe))deallocate(copy%rt%rIe)
    endif
    if(allocated(ctx%rt%curr))then
      copy%rt%curr=ctx%rt%curr
    else
      if(allocated(copy%rt%curr))deallocate(copy%rt%curr)
    endif
    if(allocated(ctx%rt%E_ext))then
      copy%rt%E_ext=ctx%rt%E_ext
    else
      if(allocated(copy%rt%E_ext))deallocate(copy%rt%E_ext)
    endif
    if(allocated(ctx%rt%E_ind))then
      copy%rt%E_ind=ctx%rt%E_ind
    else
      if(allocated(copy%rt%E_ind))deallocate(copy%rt%E_ind)
    endif
    if(allocated(ctx%rt%E_tot))then
      copy%rt%E_tot=ctx%rt%E_tot
    else
      if(allocated(copy%rt%E_tot))deallocate(copy%rt%E_tot)
    endif
    if(allocated(ctx%rt%Ac_ext))then
      copy%rt%Ac_ext=ctx%rt%Ac_ext
    else
      if(allocated(copy%rt%Ac_ext))deallocate(copy%rt%Ac_ext)
    endif
    if(allocated(ctx%rt%Ac_ind))then
      copy%rt%Ac_ind=ctx%rt%Ac_ind
    else
      if(allocated(copy%rt%Ac_ind))deallocate(copy%rt%Ac_ind)
    endif
    if(allocated(ctx%rt%Ac_tot))then
      copy%rt%Ac_tot=ctx%rt%Ac_tot
    else
      if(allocated(copy%rt%Ac_tot))deallocate(copy%rt%Ac_tot)
    endif
    if(allocated(ctx%rt%zc))then
      copy%rt%zc=ctx%rt%zc
    else
      if(allocated(copy%rt%zc))deallocate(copy%rt%zc)
    endif
    if(allocated(ctx%rt%vloc_t))then
      copy%rt%vloc_t=ctx%rt%vloc_t
    else
      if(allocated(copy%rt%vloc_t))deallocate(copy%rt%vloc_t)
    endif
    if(allocated(ctx%rt%vloc_new))then
      copy%rt%vloc_new=ctx%rt%vloc_new
    else
      if(allocated(copy%rt%vloc_new))deallocate(copy%rt%vloc_new)
    endif
    if(allocated(ctx%rt%vloc_old))then
      copy%rt%vloc_old=ctx%rt%vloc_old
    else
      if(allocated(copy%rt%vloc_old))deallocate(copy%rt%vloc_old)
    endif
    if(allocated(ctx%rt%rho0_s))then
      copy%rt%rho0_s=ctx%rt%rho0_s
    else
      if(allocated(copy%rt%rho0_s))deallocate(copy%rt%rho0_s)
    endif
    copy%rt%vonf=ctx%rt%vonf
    copy%rt%j_e=ctx%rt%j_e
    copy%rt%E_old=ctx%rt%E_old
    copy%rt%rho0_loaded=ctx%rt%rho0_loaded
    copy%fields=0d0
    copy%fields(:,1)=ctx%system%vec_Ac
    copy%fields(:,2)=ctx%system%vec_Ac_ext
    copy%fields(:,3)=ctx%system%vec_E
    copy%fields(:,4)=ctx%system%vec_E_ext
  endif
end subroutine
subroutine checkpoint_io(unit,reading,value,ios,array_limit)
  integer,intent(in)::unit
  logical,intent(in)::reading
  type(checkpoint_payload),intent(inout)::value
  integer,intent(inout)::ios
  integer(8),intent(in)::array_limit
  logical::present_array
  integer::lo(7),hi(7),i1,i2,k,allocation_status
  integer(8)::elements,extent
  if(ios/=0)return
  call checkpoint_io_s_orbital(unit,reading,value%spsi_in,ios,array_limit)
  if(ios/=0)return
  call checkpoint_io_s_orbital(unit,reading,value%spsi_out,ios,array_limit)
  if(ios/=0)return
  call checkpoint_io_s_orbital(unit,reading,value%tpsi,ios,array_limit)
  if(ios/=0)return
  call checkpoint_io_s_scalar(unit,reading,value%rho,ios,array_limit)
  if(ios/=0)return
  call checkpoint_io_s_scalar(unit,reading,value%rho_jm,ios,array_limit)
  if(ios/=0)return
  call checkpoint_io_s_scalar(unit,reading,value%Vh,ios,array_limit)
  if(ios/=0)return
  call checkpoint_io_s_scalar(unit,reading,value%Vh_stock1,ios,array_limit)
  if(ios/=0)return
  call checkpoint_io_s_scalar(unit,reading,value%Vh_stock2,ios,array_limit)
  if(ios/=0)return
  call checkpoint_io_s_scalar(unit,reading,value%Vbox,ios,array_limit)
  if(ios/=0)return
  call checkpoint_io_s_scalar(unit,reading,value%Vpsl,ios,array_limit)
  if(ios/=0)return
  call checkpoint_io_s_rt(unit,reading,value%rt,ios,array_limit)
  if(ios/=0)return
  call checkpoint_io_s_poisson(unit,reading,value%poisson,ios,array_limit)
  if(ios/=0)return
  call checkpoint_io_s_dft_energy(unit,reading,value%energy,ios,array_limit)
  if(ios/=0)return
  present_array=allocated(value%rho_s)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%rho_s))deallocate(value%rho_s)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%rho_s))then
        if(any(lo(1:1)/=lbound(value%rho_s)).or.any(hi(1:1)/=ubound(value%rho_s)))then
          ios=1;return
        endif
      else
        allocate(value%rho_s(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%rho_s);hi(1:1)=ubound(value%rho_s)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    do i1=lbound(value%rho_s,1),ubound(value%rho_s,1)
    call checkpoint_io_s_scalar(unit,reading,value%rho_s(i1),ios,array_limit)
    if(ios/=0)return
    enddo
  endif
  present_array=allocated(value%V_local)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%V_local))deallocate(value%V_local)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%V_local))then
        if(any(lo(1:1)/=lbound(value%V_local)).or.any(hi(1:1)/=ubound(value%V_local)))then
          ios=1;return
        endif
      else
        allocate(value%V_local(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%V_local);hi(1:1)=ubound(value%V_local)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    do i1=lbound(value%V_local,1),ubound(value%V_local,1)
    call checkpoint_io_s_scalar(unit,reading,value%V_local(i1),ios,array_limit)
    if(ios/=0)return
    enddo
  endif
  present_array=allocated(value%Vxc)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%Vxc))deallocate(value%Vxc)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%Vxc))then
        if(any(lo(1:1)/=lbound(value%Vxc)).or.any(hi(1:1)/=ubound(value%Vxc)))then
          ios=1;return
        endif
      else
        allocate(value%Vxc(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%Vxc);hi(1:1)=ubound(value%Vxc)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    do i1=lbound(value%Vxc,1),ubound(value%Vxc,1)
    call checkpoint_io_s_scalar(unit,reading,value%Vxc(i1),ios,array_limit)
    if(ios/=0)return
    enddo
  endif
  if(reading)then
    read(unit,iostat=ios) value%fields
  else
    write(unit,iostat=ios) value%fields
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%Mit
  else
    write(unit,iostat=ios) value%Mit
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%accepted
  else
    write(unit,iostat=ios) value%accepted
  endif
  if(ios/=0)return
end subroutine
subroutine checkpoint_io_s_scalar(unit,reading,value,ios,array_limit)
  integer,intent(in)::unit
  logical,intent(in)::reading
  type(s_scalar),intent(inout)::value
  integer,intent(inout)::ios
  integer(8),intent(in)::array_limit
  logical::present_array
  integer::lo(7),hi(7),i1,i2,k,allocation_status
  integer(8)::elements,extent
  if(ios/=0)return
  present_array=allocated(value%f)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%f))deallocate(value%f)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%f))then
        if(any(lo(1:3)/=lbound(value%f)).or.any(hi(1:3)/=ubound(value%f)))then
          ios=1;return
        endif
      else
        allocate(value%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%f);hi(1:3)=ubound(value%f)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%f
    else
      write(unit,iostat=ios) value%f
    endif
    if(ios/=0)return
  endif
end subroutine
subroutine checkpoint_io_s_vector(unit,reading,value,ios,array_limit)
  integer,intent(in)::unit
  logical,intent(in)::reading
  type(s_vector),intent(inout)::value
  integer,intent(inout)::ios
  integer(8),intent(in)::array_limit
  logical::present_array
  integer::lo(7),hi(7),i1,i2,k,allocation_status
  integer(8)::elements,extent
  if(ios/=0)return
  present_array=allocated(value%v)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%v))deallocate(value%v)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:4),hi(1:4)
      if(ios/=0)return
      elements=1_8
      do k=1,4
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%v))then
        if(any(lo(1:4)/=lbound(value%v)).or.any(hi(1:4)/=ubound(value%v)))then
          ios=1;return
        endif
      else
        allocate(value%v(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),lo(4):hi(4)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:4)=lbound(value%v);hi(1:4)=ubound(value%v)
      write(unit,iostat=ios)lo(1:4),hi(1:4)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%v
    else
      write(unit,iostat=ios) value%v
    endif
    if(ios/=0)return
  endif
end subroutine
subroutine checkpoint_io_s_orbital(unit,reading,value,ios,array_limit)
  integer,intent(in)::unit
  logical,intent(in)::reading
  type(s_orbital),intent(inout)::value
  integer,intent(inout)::ios
  integer(8),intent(in)::array_limit
  logical::present_array
  integer::lo(7),hi(7),i1,i2,k,allocation_status
  integer(8)::elements,extent
  if(ios/=0)return
  present_array=allocated(value%rwf)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%rwf))deallocate(value%rwf)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:7),hi(1:7)
      if(ios/=0)return
      elements=1_8
      do k=1,7
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%rwf))then
        if(any(lo(1:7)/=lbound(value%rwf)).or.any(hi(1:7)/=ubound(value%rwf)))then
          ios=1;return
        endif
      else
        allocate(value%rwf(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),lo(4):hi(4), &
          lo(5):hi(5),lo(6):hi(6),lo(7):hi(7)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:7)=lbound(value%rwf);hi(1:7)=ubound(value%rwf)
      write(unit,iostat=ios)lo(1:7),hi(1:7)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%rwf
    else
      write(unit,iostat=ios) value%rwf
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%zwf)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%zwf))deallocate(value%zwf)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:7),hi(1:7)
      if(ios/=0)return
      elements=1_8
      do k=1,7
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%zwf))then
        if(any(lo(1:7)/=lbound(value%zwf)).or.any(hi(1:7)/=ubound(value%zwf)))then
          ios=1;return
        endif
      else
        allocate(value%zwf(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),lo(4):hi(4), &
          lo(5):hi(5),lo(6):hi(6),lo(7):hi(7)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:7)=lbound(value%zwf);hi(1:7)=ubound(value%zwf)
      write(unit,iostat=ios)lo(1:7),hi(1:7)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%zwf
    else
      write(unit,iostat=ios) value%zwf
    endif
    if(ios/=0)return
  endif
  if(reading)then
    read(unit,iostat=ios) value%update_zwf_overlap
  else
    write(unit,iostat=ios) value%update_zwf_overlap
  endif
  if(ios/=0)return
end subroutine
subroutine checkpoint_io_s_dft_energy(unit,reading,value,ios,array_limit)
  integer,intent(in)::unit
  logical,intent(in)::reading
  type(s_dft_energy),intent(inout)::value
  integer,intent(inout)::ios
  integer(8),intent(in)::array_limit
  logical::present_array
  integer::lo(7),hi(7),i1,i2,k,allocation_status
  integer(8)::elements,extent
  if(ios/=0)return
  present_array=allocated(value%esp)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%esp))deallocate(value%esp)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%esp))then
        if(any(lo(1:3)/=lbound(value%esp)).or.any(hi(1:3)/=ubound(value%esp)))then
          ios=1;return
        endif
      else
        allocate(value%esp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%esp);hi(1:3)=ubound(value%esp)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%esp
    else
      write(unit,iostat=ios) value%esp
    endif
    if(ios/=0)return
  endif
  if(reading)then
    read(unit,iostat=ios) value%E_tot
  else
    write(unit,iostat=ios) value%E_tot
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%E_kin
  else
    write(unit,iostat=ios) value%E_kin
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%E_h
  else
    write(unit,iostat=ios) value%E_h
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%E_xc
  else
    write(unit,iostat=ios) value%E_xc
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%E_ion_ion
  else
    write(unit,iostat=ios) value%E_ion_ion
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%E_ion_loc
  else
    write(unit,iostat=ios) value%E_ion_loc
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%E_ion_nloc
  else
    write(unit,iostat=ios) value%E_ion_nloc
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%E_U
  else
    write(unit,iostat=ios) value%E_U
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%E_tot0
  else
    write(unit,iostat=ios) value%E_tot0
  endif
  if(ios/=0)return
end subroutine
subroutine checkpoint_io_s_poisson(unit,reading,value,ios,array_limit)
  integer,intent(in)::unit
  logical,intent(in)::reading
  type(s_poisson),intent(inout)::value
  integer,intent(inout)::ios
  integer(8),intent(in)::array_limit
  logical::present_array
  integer::lo(7),hi(7),i1,i2,k,allocation_status
  integer(8)::elements,extent
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%iterVh
  else
    write(unit,iostat=ios) value%iterVh
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%npole_partial
  else
    write(unit,iostat=ios) value%npole_partial
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%npole_total
  else
    write(unit,iostat=ios) value%npole_total
  endif
  if(ios/=0)return
  present_array=allocated(value%ipole_tbl)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ipole_tbl))deallocate(value%ipole_tbl)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ipole_tbl))then
        if(any(lo(1:1)/=lbound(value%ipole_tbl)).or.any(hi(1:1)/=ubound(value%ipole_tbl)))then
          ios=1;return
        endif
      else
        allocate(value%ipole_tbl(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%ipole_tbl);hi(1:1)=ubound(value%ipole_tbl)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ipole_tbl
    else
      write(unit,iostat=ios) value%ipole_tbl
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ig_num)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ig_num))deallocate(value%ig_num)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ig_num))then
        if(any(lo(1:1)/=lbound(value%ig_num)).or.any(hi(1:1)/=ubound(value%ig_num)))then
          ios=1;return
        endif
      else
        allocate(value%ig_num(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%ig_num);hi(1:1)=ubound(value%ig_num)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ig_num
    else
      write(unit,iostat=ios) value%ig_num
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ig)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ig))deallocate(value%ig)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ig))then
        if(any(lo(1:3)/=lbound(value%ig)).or.any(hi(1:3)/=ubound(value%ig)))then
          ios=1;return
        endif
      else
        allocate(value%ig(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ig);hi(1:3)=ubound(value%ig)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ig
    else
      write(unit,iostat=ios) value%ig
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ig_bound)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ig_bound))deallocate(value%ig_bound)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ig_bound))then
        if(any(lo(1:3)/=lbound(value%ig_bound)).or.any(hi(1:3)/=ubound(value%ig_bound)))then
          ios=1;return
        endif
      else
        allocate(value%ig_bound(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ig_bound);hi(1:3)=ubound(value%ig_bound)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ig_bound
    else
      write(unit,iostat=ios) value%ig_bound
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%wkbound)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%wkbound))deallocate(value%wkbound)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%wkbound))then
        if(any(lo(1:1)/=lbound(value%wkbound)).or.any(hi(1:1)/=ubound(value%wkbound)))then
          ios=1;return
        endif
      else
        allocate(value%wkbound(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%wkbound);hi(1:1)=ubound(value%wkbound)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%wkbound
    else
      write(unit,iostat=ios) value%wkbound
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%wkbound2)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%wkbound2))deallocate(value%wkbound2)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%wkbound2))then
        if(any(lo(1:1)/=lbound(value%wkbound2)).or.any(hi(1:1)/=ubound(value%wkbound2)))then
          ios=1;return
        endif
      else
        allocate(value%wkbound2(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%wkbound2);hi(1:1)=ubound(value%wkbound2)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%wkbound2
    else
      write(unit,iostat=ios) value%wkbound2
    endif
    if(ios/=0)return
  endif
  if(reading)then
    read(unit,iostat=ios) value%n_multipole_xyz
  else
    write(unit,iostat=ios) value%n_multipole_xyz
  endif
  if(ios/=0)return
  present_array=allocated(value%zrhoG_ele)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%zrhoG_ele))deallocate(value%zrhoG_ele)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%zrhoG_ele))then
        if(any(lo(1:3)/=lbound(value%zrhoG_ele)).or.any(hi(1:3)/=ubound(value%zrhoG_ele)))then
          ios=1;return
        endif
      else
        allocate(value%zrhoG_ele(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%zrhoG_ele);hi(1:3)=ubound(value%zrhoG_ele)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%zrhoG_ele
    else
      write(unit,iostat=ios) value%zrhoG_ele
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff1x)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff1x))deallocate(value%ff1x)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff1x))then
        if(any(lo(1:3)/=lbound(value%ff1x)).or.any(hi(1:3)/=ubound(value%ff1x)))then
          ios=1;return
        endif
      else
        allocate(value%ff1x(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff1x);hi(1:3)=ubound(value%ff1x)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff1x
    else
      write(unit,iostat=ios) value%ff1x
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff1y)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff1y))deallocate(value%ff1y)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff1y))then
        if(any(lo(1:3)/=lbound(value%ff1y)).or.any(hi(1:3)/=ubound(value%ff1y)))then
          ios=1;return
        endif
      else
        allocate(value%ff1y(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff1y);hi(1:3)=ubound(value%ff1y)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff1y
    else
      write(unit,iostat=ios) value%ff1y
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff1z)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff1z))deallocate(value%ff1z)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff1z))then
        if(any(lo(1:3)/=lbound(value%ff1z)).or.any(hi(1:3)/=ubound(value%ff1z)))then
          ios=1;return
        endif
      else
        allocate(value%ff1z(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff1z);hi(1:3)=ubound(value%ff1z)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff1z
    else
      write(unit,iostat=ios) value%ff1z
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff2x)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff2x))deallocate(value%ff2x)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff2x))then
        if(any(lo(1:3)/=lbound(value%ff2x)).or.any(hi(1:3)/=ubound(value%ff2x)))then
          ios=1;return
        endif
      else
        allocate(value%ff2x(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff2x);hi(1:3)=ubound(value%ff2x)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff2x
    else
      write(unit,iostat=ios) value%ff2x
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff2y)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff2y))deallocate(value%ff2y)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff2y))then
        if(any(lo(1:3)/=lbound(value%ff2y)).or.any(hi(1:3)/=ubound(value%ff2y)))then
          ios=1;return
        endif
      else
        allocate(value%ff2y(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff2y);hi(1:3)=ubound(value%ff2y)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff2y
    else
      write(unit,iostat=ios) value%ff2y
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff2z)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff2z))deallocate(value%ff2z)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff2z))then
        if(any(lo(1:3)/=lbound(value%ff2z)).or.any(hi(1:3)/=ubound(value%ff2z)))then
          ios=1;return
        endif
      else
        allocate(value%ff2z(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff2z);hi(1:3)=ubound(value%ff2z)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff2z
    else
      write(unit,iostat=ios) value%ff2z
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff1)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff1))deallocate(value%ff1)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff1))then
        if(any(lo(1:3)/=lbound(value%ff1)).or.any(hi(1:3)/=ubound(value%ff1)))then
          ios=1;return
        endif
      else
        allocate(value%ff1(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff1);hi(1:3)=ubound(value%ff1)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff1
    else
      write(unit,iostat=ios) value%ff1
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff2)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff2))deallocate(value%ff2)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff2))then
        if(any(lo(1:3)/=lbound(value%ff2)).or.any(hi(1:3)/=ubound(value%ff2)))then
          ios=1;return
        endif
      else
        allocate(value%ff2(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff2);hi(1:3)=ubound(value%ff2)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff2
    else
      write(unit,iostat=ios) value%ff2
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff3x)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff3x))deallocate(value%ff3x)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff3x))then
        if(any(lo(1:3)/=lbound(value%ff3x)).or.any(hi(1:3)/=ubound(value%ff3x)))then
          ios=1;return
        endif
      else
        allocate(value%ff3x(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff3x);hi(1:3)=ubound(value%ff3x)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff3x
    else
      write(unit,iostat=ios) value%ff3x
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff3y)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff3y))deallocate(value%ff3y)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff3y))then
        if(any(lo(1:3)/=lbound(value%ff3y)).or.any(hi(1:3)/=ubound(value%ff3y)))then
          ios=1;return
        endif
      else
        allocate(value%ff3y(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff3y);hi(1:3)=ubound(value%ff3y)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff3y
    else
      write(unit,iostat=ios) value%ff3y
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff3z)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff3z))deallocate(value%ff3z)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff3z))then
        if(any(lo(1:3)/=lbound(value%ff3z)).or.any(hi(1:3)/=ubound(value%ff3z)))then
          ios=1;return
        endif
      else
        allocate(value%ff3z(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff3z);hi(1:3)=ubound(value%ff3z)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff3z
    else
      write(unit,iostat=ios) value%ff3z
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff4x)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff4x))deallocate(value%ff4x)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff4x))then
        if(any(lo(1:3)/=lbound(value%ff4x)).or.any(hi(1:3)/=ubound(value%ff4x)))then
          ios=1;return
        endif
      else
        allocate(value%ff4x(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff4x);hi(1:3)=ubound(value%ff4x)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff4x
    else
      write(unit,iostat=ios) value%ff4x
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff4y)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff4y))deallocate(value%ff4y)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff4y))then
        if(any(lo(1:3)/=lbound(value%ff4y)).or.any(hi(1:3)/=ubound(value%ff4y)))then
          ios=1;return
        endif
      else
        allocate(value%ff4y(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff4y);hi(1:3)=ubound(value%ff4y)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff4y
    else
      write(unit,iostat=ios) value%ff4y
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff4z)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff4z))deallocate(value%ff4z)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff4z))then
        if(any(lo(1:3)/=lbound(value%ff4z)).or.any(hi(1:3)/=ubound(value%ff4z)))then
          ios=1;return
        endif
      else
        allocate(value%ff4z(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff4z);hi(1:3)=ubound(value%ff4z)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff4z
    else
      write(unit,iostat=ios) value%ff4z
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%dgf)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%dgf))deallocate(value%dgf)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%dgf))then
        if(any(lo(1:3)/=lbound(value%dgf)).or.any(hi(1:3)/=ubound(value%dgf)))then
          ios=1;return
        endif
      else
        allocate(value%dgf(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%dgf);hi(1:3)=ubound(value%dgf)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%dgf
    else
      write(unit,iostat=ios) value%dgf
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%a_ffte)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%a_ffte))deallocate(value%a_ffte)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%a_ffte))then
        if(any(lo(1:3)/=lbound(value%a_ffte)).or.any(hi(1:3)/=ubound(value%a_ffte)))then
          ios=1;return
        endif
      else
        allocate(value%a_ffte(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%a_ffte);hi(1:3)=ubound(value%a_ffte)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%a_ffte
    else
      write(unit,iostat=ios) value%a_ffte
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%b_ffte)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%b_ffte))deallocate(value%b_ffte)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%b_ffte))then
        if(any(lo(1:3)/=lbound(value%b_ffte)).or.any(hi(1:3)/=ubound(value%b_ffte)))then
          ios=1;return
        endif
      else
        allocate(value%b_ffte(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%b_ffte);hi(1:3)=ubound(value%b_ffte)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%b_ffte
    else
      write(unit,iostat=ios) value%b_ffte
    endif
    if(ios/=0)return
  endif
end subroutine
subroutine checkpoint_io_s_rt(unit,reading,value,ios,array_limit)
  integer,intent(in)::unit
  logical,intent(in)::reading
  type(s_rt),intent(inout)::value
  integer,intent(inout)::ios
  integer(8),intent(in)::array_limit
  logical::present_array
  integer::lo(7),hi(7),i1,i2,k,allocation_status
  integer(8)::elements,extent
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%Dp0_e
  else
    write(unit,iostat=ios) value%Dp0_e
  endif
  if(ios/=0)return
  present_array=allocated(value%dDp_e)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%dDp_e))deallocate(value%dDp_e)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
      elements=1_8
      do k=1,2
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%dDp_e))then
        if(any(lo(1:2)/=lbound(value%dDp_e)).or.any(hi(1:2)/=ubound(value%dDp_e)))then
          ios=1;return
        endif
      else
        allocate(value%dDp_e(lo(1):hi(1),lo(2):hi(2)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:2)=lbound(value%dDp_e);hi(1:2)=ubound(value%dDp_e)
      write(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%dDp_e
    else
      write(unit,iostat=ios) value%dDp_e
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%Dp_e)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%Dp_e))deallocate(value%Dp_e)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
      elements=1_8
      do k=1,2
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%Dp_e))then
        if(any(lo(1:2)/=lbound(value%Dp_e)).or.any(hi(1:2)/=ubound(value%Dp_e)))then
          ios=1;return
        endif
      else
        allocate(value%Dp_e(lo(1):hi(1),lo(2):hi(2)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:2)=lbound(value%Dp_e);hi(1:2)=ubound(value%Dp_e)
      write(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%Dp_e
    else
      write(unit,iostat=ios) value%Dp_e
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%Dp_i)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%Dp_i))deallocate(value%Dp_i)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
      elements=1_8
      do k=1,2
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%Dp_i))then
        if(any(lo(1:2)/=lbound(value%Dp_i)).or.any(hi(1:2)/=ubound(value%Dp_i)))then
          ios=1;return
        endif
      else
        allocate(value%Dp_i(lo(1):hi(1),lo(2):hi(2)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:2)=lbound(value%Dp_i);hi(1:2)=ubound(value%Dp_i)
      write(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%Dp_i
    else
      write(unit,iostat=ios) value%Dp_i
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%Qp_e)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%Qp_e))deallocate(value%Qp_e)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%Qp_e))then
        if(any(lo(1:3)/=lbound(value%Qp_e)).or.any(hi(1:3)/=ubound(value%Qp_e)))then
          ios=1;return
        endif
      else
        allocate(value%Qp_e(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%Qp_e);hi(1:3)=ubound(value%Qp_e)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%Qp_e
    else
      write(unit,iostat=ios) value%Qp_e
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%rIe)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%rIe))deallocate(value%rIe)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%rIe))then
        if(any(lo(1:1)/=lbound(value%rIe)).or.any(hi(1:1)/=ubound(value%rIe)))then
          ios=1;return
        endif
      else
        allocate(value%rIe(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%rIe);hi(1:1)=ubound(value%rIe)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%rIe
    else
      write(unit,iostat=ios) value%rIe
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%curr)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%curr))deallocate(value%curr)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
      elements=1_8
      do k=1,2
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%curr))then
        if(any(lo(1:2)/=lbound(value%curr)).or.any(hi(1:2)/=ubound(value%curr)))then
          ios=1;return
        endif
      else
        allocate(value%curr(lo(1):hi(1),lo(2):hi(2)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:2)=lbound(value%curr);hi(1:2)=ubound(value%curr)
      write(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%curr
    else
      write(unit,iostat=ios) value%curr
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%E_ext)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%E_ext))deallocate(value%E_ext)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
      elements=1_8
      do k=1,2
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%E_ext))then
        if(any(lo(1:2)/=lbound(value%E_ext)).or.any(hi(1:2)/=ubound(value%E_ext)))then
          ios=1;return
        endif
      else
        allocate(value%E_ext(lo(1):hi(1),lo(2):hi(2)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:2)=lbound(value%E_ext);hi(1:2)=ubound(value%E_ext)
      write(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%E_ext
    else
      write(unit,iostat=ios) value%E_ext
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%E_ind)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%E_ind))deallocate(value%E_ind)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
      elements=1_8
      do k=1,2
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%E_ind))then
        if(any(lo(1:2)/=lbound(value%E_ind)).or.any(hi(1:2)/=ubound(value%E_ind)))then
          ios=1;return
        endif
      else
        allocate(value%E_ind(lo(1):hi(1),lo(2):hi(2)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:2)=lbound(value%E_ind);hi(1:2)=ubound(value%E_ind)
      write(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%E_ind
    else
      write(unit,iostat=ios) value%E_ind
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%E_tot)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%E_tot))deallocate(value%E_tot)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
      elements=1_8
      do k=1,2
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%E_tot))then
        if(any(lo(1:2)/=lbound(value%E_tot)).or.any(hi(1:2)/=ubound(value%E_tot)))then
          ios=1;return
        endif
      else
        allocate(value%E_tot(lo(1):hi(1),lo(2):hi(2)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:2)=lbound(value%E_tot);hi(1:2)=ubound(value%E_tot)
      write(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%E_tot
    else
      write(unit,iostat=ios) value%E_tot
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%Ac_ext)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%Ac_ext))deallocate(value%Ac_ext)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
      elements=1_8
      do k=1,2
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%Ac_ext))then
        if(any(lo(1:2)/=lbound(value%Ac_ext)).or.any(hi(1:2)/=ubound(value%Ac_ext)))then
          ios=1;return
        endif
      else
        allocate(value%Ac_ext(lo(1):hi(1),lo(2):hi(2)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:2)=lbound(value%Ac_ext);hi(1:2)=ubound(value%Ac_ext)
      write(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%Ac_ext
    else
      write(unit,iostat=ios) value%Ac_ext
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%Ac_ind)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%Ac_ind))deallocate(value%Ac_ind)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
      elements=1_8
      do k=1,2
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%Ac_ind))then
        if(any(lo(1:2)/=lbound(value%Ac_ind)).or.any(hi(1:2)/=ubound(value%Ac_ind)))then
          ios=1;return
        endif
      else
        allocate(value%Ac_ind(lo(1):hi(1),lo(2):hi(2)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:2)=lbound(value%Ac_ind);hi(1:2)=ubound(value%Ac_ind)
      write(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%Ac_ind
    else
      write(unit,iostat=ios) value%Ac_ind
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%Ac_tot)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%Ac_tot))deallocate(value%Ac_tot)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
      elements=1_8
      do k=1,2
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%Ac_tot))then
        if(any(lo(1:2)/=lbound(value%Ac_tot)).or.any(hi(1:2)/=ubound(value%Ac_tot)))then
          ios=1;return
        endif
      else
        allocate(value%Ac_tot(lo(1):hi(1),lo(2):hi(2)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:2)=lbound(value%Ac_tot);hi(1:2)=ubound(value%Ac_tot)
      write(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%Ac_tot
    else
      write(unit,iostat=ios) value%Ac_tot
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%zc)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%zc))deallocate(value%zc)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%zc))then
        if(any(lo(1:1)/=lbound(value%zc)).or.any(hi(1:1)/=ubound(value%zc)))then
          ios=1;return
        endif
      else
        allocate(value%zc(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%zc);hi(1:1)=ubound(value%zc)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%zc
    else
      write(unit,iostat=ios) value%zc
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%vloc_t)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%vloc_t))deallocate(value%vloc_t)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%vloc_t))then
        if(any(lo(1:1)/=lbound(value%vloc_t)).or.any(hi(1:1)/=ubound(value%vloc_t)))then
          ios=1;return
        endif
      else
        allocate(value%vloc_t(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%vloc_t);hi(1:1)=ubound(value%vloc_t)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    do i1=lbound(value%vloc_t,1),ubound(value%vloc_t,1)
    call checkpoint_io_s_scalar(unit,reading,value%vloc_t(i1),ios,array_limit)
    if(ios/=0)return
    enddo
  endif
  present_array=allocated(value%vloc_new)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%vloc_new))deallocate(value%vloc_new)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%vloc_new))then
        if(any(lo(1:1)/=lbound(value%vloc_new)).or.any(hi(1:1)/=ubound(value%vloc_new)))then
          ios=1;return
        endif
      else
        allocate(value%vloc_new(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%vloc_new);hi(1:1)=ubound(value%vloc_new)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    do i1=lbound(value%vloc_new,1),ubound(value%vloc_new,1)
    call checkpoint_io_s_scalar(unit,reading,value%vloc_new(i1),ios,array_limit)
    if(ios/=0)return
    enddo
  endif
  present_array=allocated(value%vloc_old)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%vloc_old))deallocate(value%vloc_old)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
      elements=1_8
      do k=1,2
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%vloc_old))then
        if(any(lo(1:2)/=lbound(value%vloc_old)).or.any(hi(1:2)/=ubound(value%vloc_old)))then
          ios=1;return
        endif
      else
        allocate(value%vloc_old(lo(1):hi(1),lo(2):hi(2)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:2)=lbound(value%vloc_old);hi(1:2)=ubound(value%vloc_old)
      write(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
    endif
    do i1=lbound(value%vloc_old,1),ubound(value%vloc_old,1)
    do i2=lbound(value%vloc_old,2),ubound(value%vloc_old,2)
    call checkpoint_io_s_scalar(unit,reading,value%vloc_old(i1,i2),ios,array_limit)
    if(ios/=0)return
    enddo
    enddo
  endif
  present_array=allocated(value%rho0_s)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%rho0_s))deallocate(value%rho0_s)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>array_limit)then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>array_limit/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%rho0_s))then
        if(any(lo(1:1)/=lbound(value%rho0_s)).or.any(hi(1:1)/=ubound(value%rho0_s)))then
          ios=1;return
        endif
      else
        allocate(value%rho0_s(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%rho0_s);hi(1:1)=ubound(value%rho0_s)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    do i1=lbound(value%rho0_s,1),ubound(value%rho0_s,1)
    call checkpoint_io_s_scalar(unit,reading,value%rho0_s(i1),ios,array_limit)
    if(ios/=0)return
    enddo
  endif
  call checkpoint_io_s_scalar(unit,reading,value%vonf,ios,array_limit)
  if(ios/=0)return
  call checkpoint_io_s_vector(unit,reading,value%j_e,ios,array_limit)
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%E_old
  else
    write(unit,iostat=ios) value%E_old
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%rho0_loaded
  else
    write(unit,iostat=ios) value%rho0_loaded
  endif
  if(ios/=0)return
end subroutine
end module
