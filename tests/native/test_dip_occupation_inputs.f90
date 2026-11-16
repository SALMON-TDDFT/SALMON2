! Electron-count diagnostics use each system's immutable occupation inputs.
program probe
 use structures
 use parallelization,only:setup_parallel,end_parallel
 use salmon_global
 use communication,only:COMM_GROUP_SELF
 use dip,only:subdip
 use reference_dip,only:reference_subdip=>subdip
 implicit none
 type(s_dft_system)::system
 type(s_rgrid)::grid
 type(s_scalar)::rho(2)
 type(s_rt)::current,reference
 type(s_poisson)::poisson
 type(s_pp_info)::pp
 integer::i,case_id,selected
 real(8)::actual,expected,count
 character(16)::mode,case_arg
 call setup_parallel
 call get_command_argument(1,mode)
 call get_command_argument(2,case_arg)
 selected=0
 if(len_trim(case_arg)>0)read(case_arg,*)selected
 quiet=.true.;yn_jm='y';iperiodic=3;out_rt_energy_step=1;dt=.1d0
 grid%is=1;grid%ie=1;grid%num=1
 system%hvol=1d0;system%hgs=1d0
 allocate(current%Dp_i(3,0:1),current%Dp_e(3,0:1),current%dDp_e(3,0:1),current%rIe(0:1))
 current%Dp_i=0;current%Dp_e=0;current%dDp_e=0;current%rIe=0
 current%Dp0_e=0;current%time_step=dt
 do i=1,2
  allocate(rho(i)%f(1:1,1:1,1:1))
 enddo
 do case_id=1,4
  if(selected/=0.and.case_id/=selected)cycle
  system%nspin=1;count=4d0
  system%occupation_input%nelec=4;system%occupation_input%nelec_spin=0
  system%occupation_input_ready=.true.
  if(case_id==2)then
   system%nspin=2;count=3d0
   system%occupation_input%nelec=99;system%occupation_input%nelec_spin=[2,1]
  else if(case_id==3)then
   system%nspin=2;count=6d0;system%occupation_input%nelec=6
  else if(case_id==4)then
   system%occupation_input_ready=.false.
  endif
  nelec=system%occupation_input%nelec;nelec_spin=system%occupation_input%nelec_spin
  if(trim(mode)=='lower')count=count/10d0
  if(trim(mode)=='upper')count=count*10d0
  do i=1,2
   rho(i)%f=count/system%nspin
  enddo
  reference=current
  if(trim(mode)=='')then
   call reference_subdip(COMM_GROUP_SELF,1,reference,grid,grid,system,rho(:system%nspin),expected,poisson,0d0,pp)
  endif
  if(system%occupation_input_ready)then
   nelec=1000;nelec_spin=[1000,1000]
  endif
  call subdip(COMM_GROUP_SELF,1,current,grid,grid,system,rho(:system%nspin),actual,poisson,0d0,pp)
  if(trim(mode)/='')error stop 'boundary accepted'
  if(actual/=expected.or.any(current%Dp_e/=reference%Dp_e).or. &
     any(current%dDp_e/=reference%dDp_e).or.any(current%rIe/=reference%rIe))error stop 'output changed'
 enddo
 call end_parallel
 print *,'Dip occupation isolation PASS'
end program
