#include "config.h"
! Sequential GS lifecycle regression, not a concurrent GS context API.
program probe
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_is_finite
 use structures,only:s_gs_run_config,s_mixing,gs_density_normalization
 use mixing_sub,only:reset_mixing_rate,check_mixing_half
 use salmon_global
 use parallelization
 use inputoutput
 use math_constants
 use timer
 use salmonpy_gs_context,only:s_gs_context,gs_run,gs_close,gs_probe_write,gs_input_config
#ifdef SALMONPY_GS_IMPORT_PROBE
 use salmonpy_gs_context,only:gs_export_seed
 use salmonpy_gs_seed,only:s_gs_seed
#ifdef USE_MPI
 use mpi,only:MPI_COMM_WORLD,MPI_Finalize
 use communication,only:comm_attach,comm_child_counts
#endif
#endif
 use plusU_global,only:read_Hubbard_parameters
 use sym_sub,only:read_sw_symmetry
 implicit none
 type(s_mixing)::mix_a,mix_b
 type(s_gs_context),target::a,b
 integer::status,unit,saved_checkpoint
 character(32)::saved_metric,mode
 type(s_gs_run_config)::config_a,config_b
 call set_math_constants
 call setup_parallel
 call read_input
 calc_mode='GS'
 call read_Hubbard_parameters
 call read_sw_symmetry(yn_symmetry)
 call timer_initialize
 if(theory/='dft'.or.yn_dc/='n'.or.yn_opt/='n')stop 1
#ifdef SALMONPY_GS_IMPORT_PROBE
 call get_command_argument(1,mode)
 if(mode=='import')then
#ifdef USE_MPI
  call comm_attach(MPI_COMM_WORLD,status)
  if(status/=0)error stop 'GS probe attach failed'
#endif
  call import_control
  call end_parallel
#ifdef USE_MPI
  call MPI_Finalize(status)
#endif
  stop
 endif
#endif
 call gs_input_config(config_a);config_b=config_a
 ! Count selection is independent of all legacy module variables.
 config_a%normalization%nelec=4;config_a%normalization%nelec_spin=[2,3]
 if(gs_density_normalization(config_a%normalization,1)/=4d0)stop 19
 if(gs_density_normalization(config_a%normalization,2)/=5d0)stop 20
 config_a%normalization%nelec_spin=[0,0]
 if(gs_density_normalization(config_a%normalization,2)/=4d0)stop 21
 config_a=config_b;config_a%convergence_metric='rho_dne'
 config_a%normalization%nelec=0;config_a%normalization%nelec_spin=0
 call gs_run(a,status,config=config_a);if(status/=1)stop 22
 config_a=config_b
 config_a%max_iterations=0
 call gs_run(a,status,config=config_a);if(status/=1)stop 5
 config_a=config_b;config_a%max_iterations=huge(0)
 call gs_run(a,status,config=config_a);if(status/=1)stop 6
 config_a=config_b;config_a%threshold=ieee_value(0d0,ieee_quiet_nan)
 call gs_run(a,status,config=config_a);if(status/=1)stop 7
 config_a=config_b;config_a%pre_scf_steps=0
 call gs_run(a,status,config=config_a);if(status/=1)stop 8
 config_a=config_b;config_a%pre_scf_threshold=-1d0
 call gs_run(a,status,config=config_a);if(status/=1)stop 9
 config_a=config_b;config_a%pre_scf_threshold=1d8
 saved_checkpoint=checkpoint_interval;checkpoint_interval=1
 call gs_run(a,status,config=config_a);if(status/=1)stop 10
 checkpoint_interval=saved_checkpoint
 saved_metric=config_a%convergence_metric;config_a%convergence_metric='norm_pot'
 call gs_run(a,status,config=config_a);if(status/=1)stop 11
 config_a%convergence_metric=saved_metric
 yn_restart='y'
 call gs_run(a,status,config=config_a);if(status/=1)stop 12
 yn_restart='n'
 config_a=config_b;config_a%convergence_metric='unknown'
 call gs_run(a,status,config=config_a);if(status/=1)stop 13
 config_a=config_b;config_a%mixing%nmemory_p=22
 call gs_run(a,status,config=config_a);if(status/=1)stop 14
 config_a=config_b;config_a%mixing%mixrate=ieee_value(0d0,ieee_quiet_nan)
 call gs_run(a,status,config=config_a);if(status/=1)stop 15
 config_a=config_b;config_a%mixing%method='invalid'
 call gs_run(a,status,config=config_a);if(status/=1)stop 16
 config_a=config_b
 config_b%max_iterations=config_b%max_iterations+2
 call get_command_argument(1,mode)
 if(mode=='normalization')config_b%normalization%nelec=2*config_b%normalization%nelec
 if(mode=='distinct')then
  config_b%mixing%mixrate=0.175d0
  config_b%mixing%alpha_mb=0.325d0
  config_b%mixing%beta_p=0.275d0
  config_b%mixing%nmemory_p=1
 endif
 ! Every selected SCF control must come from the value, not these legacy globals.
 nscf=1;threshold=0d0;exx_pre_scf_threshold=0d0;exx_pre_scf_steps=99
 convergence='invalid-shared'
 method_mixing='invalid-shared';mixrate=-9d0;alpha_mb=-8d0;beta_p=-7d0
 nmemory_p=99;yn_auto_mixing='invalid';update_mixing_ratio=-1d0;nscf_init_mix_zero=99
 ! Distinct settings retain separate adaptive coefficients and reset baselines.
 mix_a%input=config_a%mixing;mix_b%input=config_b%mixing
 mix_a%input%method='simple';mix_a%input%mixrate=0.25d0;mix_a%input%update_ratio=0d0
 mix_b%input%method='simple';mix_b%input%mixrate=0.75d0
 call reset_mixing_rate(mix_a);call reset_mixing_rate(mix_b)
 call check_mixing_half(1,1d0,mix_a)
 if(mix_a%mixrate/=0.125d0.or.mix_b%mixrate/=0.75d0)stop 17
 call reset_mixing_rate(mix_a)
 if(mix_a%mixrate/=0.25d0.or.mix_b%mixrate/=0.75d0)stop 18
 base_directory='./first/'
 call gs_run(a,status,config=config_a);if(status/=0)stop 2
 call gs_run(a,status);if(status/=1)stop 3
 open(newunit=unit,file='a_before.bin',access='stream',form='unformatted',status='new')
 call gs_probe_write(a,unit);close(unit)
 config_a%normalization%nelec=-1;config_a%normalization%nelec_spin=-1
 config_a%mixing%method='invalid-caller';config_a%mixing%mixrate=-1d0
 config_a%max_iterations=1;config_a%threshold=0d0;config_a%convergence_metric='invalid-caller'
 base_directory='./second/'
 call gs_run(b,status,config=config_b);if(status/=0)stop 4
 open(newunit=unit,file='a_after.bin',access='stream',form='unformatted',status='new')
 call gs_probe_write(a,unit);close(unit)
 call gs_close(b);call gs_close(b)
 open(newunit=unit,file='a_after_close.bin',access='stream',form='unformatted',status='new')
 call gs_probe_write(a,unit);close(unit)
 call gs_close(a);call gs_close(a)
 call end_parallel
 print *,'Sequential GS lifecycle PASS'
#ifdef SALMONPY_GS_IMPORT_PROBE
contains
 subroutine import_control
  type(s_gs_context),target::empty
  type(s_gs_seed)::first,second,copy
  type(s_input_parser)::parser
  character(:),allocatable::before,after,before_b
  character(64)::first_xc
  real(8)::first_position(3)
  complex(8),allocatable::flat(:)
  integer::code,live,slots,j
  call gs_export_seed(empty,copy,code)
  if(code==0)error stop 'empty GS exported'
  call gs_run(a,code,import_only=.true.)
  if(code/=0)error stop 'import A failed'
  first_xc=xc
  before=snapshot(a)
  call gs_export_seed(a,first,code)
  if(code/=0)error stop 'export A failed'
  flat=reshape(first%psi%zwf,[size(first%psi%zwf)])
  do j=1,size(flat)
   if(.not.ieee_is_finite(real(flat(j),8)))error stop 'nonfinite imported real orbital'
   if(.not.ieee_is_finite(aimag(flat(j))))error stop 'nonfinite imported imaginary orbital'
  enddo
  call read_input('b.inp',parser)
  call read_sw_symmetry(yn_symmetry)
  call gs_run(b,code,import_only=.true.)
  if(code/=0)error stop 'import B failed'
  after=snapshot(a)
  if(before/=after)error stop 'import B changed A'
  before_b=snapshot(b)
  call gs_export_seed(b,second,code)
  if(code/=0)error stop 'export B failed'
  call first%inputs%activate(code)
  if(code/=0.or.xc/=first_xc)error stop 'first model inputs mixed'
  first_position=Rion_red(:,1)
  call second%inputs%activate(code)
  if(code/=0)error stop 'second model inputs unavailable'
  if(xc==first_xc.and.all(Rion_red(:,1)==first_position))error stop 'second model inputs mixed'
  copy=first;copy%psi%zwf=(9d0,2d0);copy%occupations=0d0
  call gs_export_seed(a,copy,code)
  if(code/=0.or.any(copy%psi%zwf/=first%psi%zwf))error stop 'detached orbitals alias'
  if(any(copy%occupations/=first%occupations))error stop 'detached occupations alias'
  call gs_close(a);call gs_close(a)
  after=snapshot(b)
  if(before_b/=after)error stop 'close A changed B'
  call read_input('a.inp',parser)
  call read_sw_symmetry(yn_symmetry)
  call gs_run(a,code,import_only=.true.)
  if(code/=0)error stop 'fresh import A failed'
  after=snapshot(a)
  if(before/=after)error stop 'fresh imported A differs'
  call gs_close(b);call gs_close(b)
  after=snapshot(a)
  if(before/=after)error stop 'close B changed A'
  call gs_close(a)
#ifdef USE_MPI
  call comm_child_counts(live,slots)
  if(live/=0)error stop 'GS import leaked communicators'
#endif
  print *,'GS_IMPORT independent/seed/fresh/close PASS'
 end subroutine
 function snapshot(value) result(bytes)
  type(s_gs_context),intent(in)::value
  character(:),allocatable::bytes
  integer::u,n
  open(newunit=u,status='scratch',access='stream',form='unformatted')
  call gs_probe_write(value,u)
  inquire(unit=u,size=n)
  allocate(character(n)::bytes)
  rewind(u);read(u)bytes;close(u)
 end function
#endif
end program
