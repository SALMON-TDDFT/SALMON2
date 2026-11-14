! Real exchange adapter routing: explicit RT states must ignore legacy sentinels.
program probe
 use iso_fortran_env,only:int64
 use structures
 use exx_rt_state,only:s_exx_rt_state,rt_state_clear,rt_state_bytes
 use exx_native,only:exx_add_action,exx_refresh,exx_taylor_stage,exx_ptcn_begin,exx_freeze, &
   exx_force_full_action,exx_adaptive_ready,exx_support_changed,exx_exchange_energy
 use exx_ace,only:exx_ace_build,exx_ace_apply
 use salmon_global,only:xc,propagator,num_kgrid,yn_spinorbit,yn_hse_wannier,yn_dc,temperature, &
   exx_mlwf_radius,exx_mlwf_norm_fraction,nstate,nelec,exx_pre_scf_active,exx_factor_mode
 use communication,only:COMM_GROUP_SELF
 implicit none
 type(s_exx_rt_state),target::a,b
 type(s_rgrid)::mg
 type(s_parallel_info)::info
 type(s_dft_system)::system
 type(s_orbital)::psi,hpsi
 complex(8)::u(8,2,1),w(8,2,1),target(8,2,1),expected(8,2,1),initial(8,2,1),actual(8,2,1)
 integer::status,i,j
 integer(int64)::bytes_a,bytes_b,before_bytes
 xc='hse06';propagator='hse_taylor4';num_kgrid=1;yn_spinorbit='n';yn_hse_wannier='n';yn_dc='n'
 temperature=-1d0;exx_mlwf_radius=0d0;exx_mlwf_norm_fraction=0d0;nstate=2;nelec=4
 exx_pre_scf_active=.false.;exx_factor_mode='learned'
 mg%is=1;mg%ie=2;mg%num=2
 info%io_s=1;info%io_e=2;info%ik_s=1;info%ik_e=1;info%numo=2;info%numk=1
 info%isize_o=1;info%isize_r=1;info%isize_k=1;info%isize_ro=1;info%numm=1
 info%icomm_o=COMM_GROUP_SELF;info%icomm_r=COMM_GROUP_SELF;info%icomm_ro=COMM_GROUP_SELF
 info%icomm_k=COMM_GROUP_SELF;info%icomm_rko=COMM_GROUP_SELF
 system%nspin=1;system%no=2;system%nk=1;system%Hvol=1d0
 allocate(psi%zwf(2,2,2,1,2,1,1),hpsi%zwf(2,2,2,1,2,1,1))
 u=0;u(1,1,1)=1d0;u(2,2,1)=1d0
 do j=1,2;do i=1,8
  target(i,j,1)=cmplx(.01d0*i,.02d0*j,8)
 enddo;enddo
 psi%zwf=reshape(target,shape(psi%zwf))
 w=-u;call exx_ace_build(a%ace,u,w,1d0,status)
 if(status/=0)stop 1
 w=-2d0*u;call exx_ace_build(b%ace,u,w,1d0,status)
 if(status/=0)stop 2
 ! Deliberately incompatible default flags: explicit calls must not consult them.
 exx_freeze=.true.;exx_force_full_action=.true.;exx_adaptive_ready=.true.;exx_support_changed=.true.
 exx_exchange_energy=99d0
 a%exx_exchange_energy=1.25d0;b%exx_exchange_energy=2.5d0
 call exx_ptcn_begin(.001d0,rt_context=a)
 call exx_ptcn_begin(.002d0,rt_context=b)
 if(a%ptcn_history_dt/=.001d0.or.b%ptcn_history_dt/=.002d0)stop 3
 do i=1,4
  call check_action(a)
  call check_action(b)
 enddo
 ! Unsupported geometry routes return before touching accepted buffers.
 call rt_state_bytes(a,before_bytes,status)
 if(status/=0)stop 30
 hpsi%zwf=cmplx(.37d0,.19d0,8)
 a%exx_force_full_action=.true.
 call exx_add_action(psi,hpsi,system,mg,info,rt_context=a,rt_status=status)
 if(status==0.or.any(hpsi%zwf/=cmplx(.37d0,.19d0,8)))stop 31
 a%exx_force_full_action=.false.;a%ace%packed=.true.
 call exx_add_action(psi,hpsi,system,mg,info,rt_context=a,rt_status=status)
 if(status==0.or.any(hpsi%zwf/=cmplx(.37d0,.19d0,8)))stop 32
 a%ace%packed=.false.;a%ace%metric_distributed=.true.
 call exx_add_action(psi,hpsi,system,mg,info,rt_context=a,rt_status=status)
 if(status==0.or.any(hpsi%zwf/=cmplx(.37d0,.19d0,8)))stop 36
 a%ace%metric_distributed=.false.;info%isize_o=2
 call exx_add_action(psi,hpsi,system,mg,info,rt_context=a,rt_status=status)
 if(status==0.or.any(hpsi%zwf/=cmplx(.37d0,.19d0,8)))stop 37
 info%isize_o=1;system%nspin=2
 call exx_add_action(psi,hpsi,system,mg,info,rt_context=a,rt_status=status)
 if(status==0.or.any(hpsi%zwf/=cmplx(.37d0,.19d0,8)))stop 38
 system%nspin=1;propagator='hse_taylor4_full'
 call exx_add_action(psi,hpsi,system,mg,info,rt_context=a,rt_status=status)
 if(status==0.or.any(hpsi%zwf/=cmplx(.37d0,.19d0,8)))stop 33
 propagator='hse_taylor4'
 call exx_refresh(system,mg,info,psi,rt_context=a,rt_status=status)
 if(status==0)stop 34
 call rt_state_bytes(a,bytes_a,status)
 if(status/=0.or.bytes_a/=before_bytes)stop 35
 call exx_taylor_stage(0,rt_context=a)
 call exx_ace_apply(a%ace,target,initial,status)
 if(status/=0)stop 4
 w=-3d0*u;call exx_ace_build(a%ace,u,w,1d0,status)
 if(status/=0)stop 5
 call exx_taylor_stage(1,rt_context=a)
 call exx_ace_apply(a%ace,target,expected,status)
 if(status/=0)stop 6
 expected=.125d0*(initial+expected)
 hpsi%zwf=0d0
 call exx_add_action(psi,hpsi,system,mg,info,rt_context=a)
 actual=reshape(hpsi%zwf,shape(actual))
 if(any(abs(actual-expected)>1d-14))stop 7
 call check_action(b)
 call exx_taylor_stage(2,rt_context=a)
 if(a%taylor_active.or.b%taylor_active)stop 8
 if(.not.exx_freeze.or..not.exx_force_full_action.or..not.exx_adaptive_ready.or..not.exx_support_changed)stop 9
 if(exx_exchange_energy/=99d0.or.a%exx_exchange_energy/=1.25d0.or.b%exx_exchange_energy/=2.5d0)stop 10
 call rt_state_bytes(a,bytes_a,status)
 if(status/=0.or.bytes_a<=0)stop 11
 call rt_state_bytes(b,bytes_b,status)
 if(status/=0.or.bytes_b<=0)stop 12
 call rt_state_clear(a)
 call check_action(b)
 call rt_state_bytes(a,bytes_a,status)
 if(status/=0.or.bytes_a/=0)stop 13
 w=-4d0*u;call exx_ace_build(a%ace,u,w,1d0,status)
 if(status/=0)stop 14
 call rt_state_clear(b)
 call check_action(a)
 call rt_state_clear(a)
 call rt_state_bytes(a,bytes_a,status)
 if(status/=0.or.bytes_a/=0)stop 15
 call rt_state_bytes(b,bytes_b,status)
 if(status/=0.or.bytes_b/=0)stop 16
 ! Explicit begin calls must not have bound the legacy default's timestep.
 call exx_ptcn_begin(.003d0)
 print *,'Explicit HSE RT routing: interleaved exchange/trials/dt, independent workspaces/default flags, release PASS'
contains
 subroutine check_action(state)
  implicit none
  type(s_exx_rt_state),target,intent(inout)::state
  complex(8)::reference(8,2,1),found(8,2,1)
  integer::code
  call exx_ace_apply(state%ace,target,reference,code)
  if(code/=0)stop 20
  hpsi%zwf=0d0
  call exx_add_action(psi,hpsi,system,mg,info,rt_context=state)
  found=reshape(hpsi%zwf,shape(found))
  if(any(abs(found-.25d0*reference)>1d-14))stop 21
 end subroutine
end program
