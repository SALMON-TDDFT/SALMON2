! Sequential GS lifecycle regression, not a concurrent GS context API.
program probe
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
 use structures,only:s_gs_run_config
 use salmon_global
 use parallelization
 use inputoutput
 use math_constants
 use timer
 use salmonpy_gs_context,only:s_gs_context,gs_run,gs_close,gs_probe_write,gs_input_config
 use plusU_global,only:read_Hubbard_parameters
 use sym_sub,only:read_sw_symmetry
 implicit none
 type(s_gs_context),target::a,b
 integer::status,unit,saved_checkpoint
 character(32)::saved_metric
 type(s_gs_run_config)::config_a,config_b
 call set_math_constants
 call setup_parallel
 call read_input
 calc_mode='GS'
 call read_Hubbard_parameters
 call read_sw_symmetry(yn_symmetry)
 call timer_initialize
 if(theory/='dft'.or.yn_dc/='n'.or.yn_opt/='n')stop 1
 call gs_input_config(config_a);config_b=config_a
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
 saved_metric=convergence;convergence='norm_pot'
 call gs_run(a,status,config=config_a);if(status/=1)stop 11
 convergence=saved_metric
 yn_restart='y'
 call gs_run(a,status,config=config_a);if(status/=1)stop 12
 yn_restart='n'
 config_a=config_b
 config_b%max_iterations=config_b%max_iterations+2
 ! Every selected SCF control must come from the value, not these legacy globals.
 nscf=1;threshold=0d0;exx_pre_scf_threshold=0d0;exx_pre_scf_steps=99
 base_directory='./first/'
 call gs_run(a,status,config=config_a);if(status/=0)stop 2
 call gs_run(a,status);if(status/=1)stop 3
 open(newunit=unit,file='a_before.bin',access='stream',form='unformatted',status='new')
 call gs_probe_write(a,unit);close(unit)
 config_a%max_iterations=1;config_a%threshold=0d0
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
end program
