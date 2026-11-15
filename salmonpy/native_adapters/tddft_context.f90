! Copyright 2026 SALMONpy contributors. Apache-2.0.
! Internal serial development context; not a stable or failure-atomic public ABI.
! Global input configuration must remain common and immutable. Never copy contexts.
#include "config.h"
module salmonpy_tddft_context
 use structures
 use exx_rt_state,only:s_exx_rt_state,rt_state_clear
 use exx_native,only:exx_bind_system,exx_unbind_system
 use initialization_rt_sub,only:initialization_rt
 use salmon_xc,only:finalize_xc
 use salmon_global
 use parallelization,only:nproc_size_global
 use sym_sub,only:use_symmetry
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 public::rt_input_config,tddft_context,rt_open,rt_step,rt_close,rt_probe_write
 type fixed_settings
  real(8)::omega,electronic_temperature
  integer::grid(3),mesh(3),bands,energy_interval,electrons,hamil,exact_interval,warmup,teachers
  character(32)::functional,propagation,factor
  character(1)::predictor
 end type
 type(fixed_settings),save::process_settings
 logical,save::settings_pinned=.false.
 type tddft_context
 private
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
 type(fixed_settings)::settings
 integer::Mit=0,accepted=0
 end type
contains
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
 subroutine rt_open(ctx,name,status,delta_t,openmp_mode,config)
 type(tddft_context),pointer,intent(inout)::ctx
 character(*),intent(in)::name
 integer,intent(out)::status
 character(*),optional,intent(in)::openmp_mode
 real(8),optional,intent(in)::delta_t
 type(s_rt_run_config),optional,intent(in)::config
 type(s_rt_run_config)::selected
 character(len(sysname))::saved_name
 logical::exists
 integer::i
 status=1
#if defined(USE_MPI) || defined(USE_OPENACC) || defined(USE_FFTW)
 return
#endif
 if(theory/='tddft_response'.or.out_rt_energy_step/=1)return
 if(spin/='unpolarized'.or.temperature>=0d0.or.any(hse_q_reduction/=1))return
 if(associated(ctx))return
 if(present(config))then
  if(present(delta_t).or.present(openmp_mode))return
  selected=config
 else
  call rt_input_config(selected)
  if(present(delta_t))selected%delta_t=delta_t
  if(present(openmp_mode))then
   if(.not.valid_mode(openmp_mode))return
   selected%stencil_mode=openmp_mode;selected%current_mode=openmp_mode;selected%force_mode=openmp_mode
  endif
 endif
 if(.not.ieee_is_finite(selected%delta_t))return
 if(selected%delta_t<=0d0)return
 if(selected%max_steps<1.or.selected%max_steps>=huge(0))return
 if(.not.valid_mode(selected%stencil_mode).or..not.valid_mode(selected%current_mode).or. &
    .not.valid_mode(selected%force_mode))return
 if(settings_pinned)then
  if(.not.settings_match(process_settings))then
   status=2;return
  endif
 endif
 if(nproc_size_global/=1.or.xc/='hse06'.or.yn_periodic/='y'.or.yn_spinorbit/='n'.or. &
    yn_hse_wannier/='n'.or.use_symmetry.or.yn_md/='n'.or.yn_dc/='n'.or.trans_longi/='tr')return
 if(yn_restart/='n'.or.e_impulse/=0d0.or.ae_shape1/='impulse'.or.ae_shape2/='none')return
 if(yn_out_current_decomposed/='n'.or.yn_out_intraband_current/='n'.or.yn_out_dns_ac_je/='n'.or. &
    yn_out_mag_decomposed_rt/='n'.or.yn_out_mag_micro_rt/='n'.or.yn_out_spin_current_decomposed/='n'.or. &
    yn_out_spin_current_micro/='n'.or.yn_out_mom_distr_rt/='n'.or.yn_out_rt_energy_components/='n')return
 if(projection_option/='no'.or.dm_unfold_option/='no'.or.yn_out_rvf_rt/='n')return
 if(len_trim(name)<1.or.len_trim(name)>len(sysname)-32)return
 do i=1,len_trim(name)
  if(index('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-',name(i:i))==0)return
 enddo
 inquire(file=trim(base_directory)//trim(name)//'_rt.data',exist=exists)
 if(exists)return
 inquire(file=trim(base_directory)//trim(name)//'_rt_energy.data',exist=exists)
 if(exists)return
 if(yn_out_dns_rt/='n'.or.yn_out_micro_je/='n'.or.yn_out_elf_rt/='n'.or.yn_out_estatic_rt/='n')return
 allocate(ctx)
 call exx_bind_system(ctx%system,ctx%hse,status)
 if(status/=0)then
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
                        ctx%pp, ctx%ppg, ctx%ppn, ctx%unfold, run_config=selected )
 sysname=saved_name
 ctx%settings=capture_settings()
 process_settings=ctx%settings;settings_pinned=.true.
 ctx%accepted=ctx%Mit
 status=0
 end subroutine
 subroutine rt_step(ctx,ac,current,total_energy,status)
 type(tddft_context),pointer,intent(inout)::ctx
 real(8),intent(in)::ac(3)
 real(8),intent(out)::current(3),total_energy
 integer,intent(out)::status
 integer::itt,i
 status=1;current=0d0;total_energy=0d0
 if(.not.associated(ctx))return
 if(.not.settings_match(ctx%settings))then
  status=2;return
 endif
 if(ctx%accepted>=ctx%rt%max_steps)return
 do i=1,3
  if(.not.ieee_is_finite(ac(i)))return
 enddo
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
 ctx%accepted=itt
 current=-ctx%rt%curr(:,itt)
 total_energy=ctx%energy%E_tot
 status=0
 end subroutine
 subroutine rt_close(ctx)
 type(tddft_context),pointer,intent(inout)::ctx
 if(.not.associated(ctx))return
 close(ctx%ofl%fh_rt);close(ctx%ofl%fh_rt_energy)
 call finalize_xc(ctx%xc_func)
 call exx_unbind_system(ctx%system)
 call rt_state_clear(ctx%hse)
 deallocate(ctx)
 nullify(ctx)
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
end module
