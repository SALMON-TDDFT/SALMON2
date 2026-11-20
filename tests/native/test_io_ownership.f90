#ifdef SALMONPY_TBMBJ_INPUT_PROBE
#include "config.h"
program tbmbj_input_host
 use mpi
 use iso_c_binding
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 use inputoutput,only:read_input,s_input_parser
 use math_constants,only:set_math_constants
 use parallelization,only:setup_parallel
 use structures,only:s_xc_config,s_xc_functional,s_pp_info
 use salmon_xc,only:xc_input_config,init_xc_config,calc_xc,finalize_xc
 use salmon_global,only:global_xc=>xc,global_c=>cval
 implicit none
 interface
  integer(c_int) function attach(comm) bind(C,name='sp_attach')
   import c_int,c_int64_t
   integer(c_int64_t),value::comm
  end function
  integer(c_int) function detach() bind(C,name='sp_detach')
   import c_int
  end function
 end interface
 type(s_input_parser)::parser
 type(s_xc_config)::config
 type(s_xc_functional)::a,b
 type(s_pp_info)::pp
 real(8)::rho(2,2,2),grad(2,2,2,3),lap(2,2,2),tau(2,2,2),current(2,2,2,3)
 real(8)::ea(2,2,2),va(2,2,2),eb(2,2,2),vb(2,2,2),again(2,2,2),vagain(2,2,2)
 integer::status,ierr,i,j,k,u
 character(16)::mode
 call MPI_Init(ierr)
 status=attach(int(MPI_COMM_WORLD,c_int64_t))
 if(status/=0)stop 1
 call setup_parallel
 call set_math_constants
 call get_command_argument(1,mode)
 call read_input('inputfile',parser,status)
 if(status/=0.or.global_xc/='tbmbj')stop 2
 call xc_input_config(config)
 if(mode=='auto')then
  if(config%cval/=-1d0)stop 3
 else
  if(config%cval/=1.2d0)stop 4
 endif
 call init_xc_config(a,config,status)
 if(status/=0)stop 5
 if(mode=='auto')then
  config%cval=1.2d0
 else
  config%cval=-1d0
 endif
 call init_xc_config(b,config,status)
 if(status/=0)stop 6
 rho=.2d0;grad=.03d0;lap=.01d0;tau=.3d0;current=0d0;pp%flag_nlcc=.false.
 call calc_xc(a,pp,rho=rho,grho=grad,rlrho=lap,tau=tau,rj=current,eexc=ea,vxc=va)
 call calc_xc(b,pp,rho=rho,grho=grad,rlrho=lap,tau=tau,rj=current,eexc=eb,vxc=vb)
 if(maxval(abs(va-vb))<1d-10)stop 7
 global_xc='pz';global_c=9d0
 call calc_xc(a,pp,rho=rho,grho=grad,rlrho=lap,tau=tau,rj=current,eexc=again,vxc=vagain)
 if(any(ea/=again).or.any(va/=vagain))stop 8
 open(newunit=u,file='trace.txt',status='new')
 do k=1,2
 do j=1,2
 do i=1,2
  if(.not.ieee_is_finite(ea(i,j,k)).or..not.ieee_is_finite(va(i,j,k)))stop 9
  if(.not.ieee_is_finite(eb(i,j,k)).or..not.ieee_is_finite(vb(i,j,k)))stop 10
  write(u,'(4(es26.17e3,1x))')ea(i,j,k),va(i,j,k),eb(i,j,k),vb(i,j,k)
 enddo
 enddo
 enddo
 close(u)
 call finalize_xc(a);call finalize_xc(b)
 status=detach()
 if(status/=0)stop 11
 call MPI_Barrier(MPI_COMM_WORLD,ierr)
 call MPI_Finalize(ierr)
 print *, 'TB-mBJ input and independent XC evaluation PASS'
end program
#elif defined(SALMONPY_HSE_LEGACY_WARNING_PROBE)
program legacy_warning_host
 use parallelization,only:setup_parallel
 use, intrinsic :: iso_fortran_env,only:error_unit
 use mpi
 use iso_c_binding
 use salmonpy_tddft_context,only:rt_open,rt_model_settings,s_rt_model_settings,tddft_context
 use salmon_global,only:theory
 implicit none
 interface
  integer(c_int) function attach(comm) bind(C,name='sp_attach')
   import c_int,c_int64_t
   integer(c_int64_t),value::comm
  end function
  integer(c_int) function detach() bind(C,name='sp_detach')
   import c_int
  end function
 end interface
 type(tddft_context),pointer::ctx=>null()
 type(s_rt_model_settings)::owner
 integer::ierr,status
 call MPI_Init(ierr)
 status=attach(int(MPI_COMM_WORLD,c_int64_t))
 if(status/=0)call MPI_Abort(MPI_COMM_WORLD,91,ierr)
 call setup_parallel
 theory='warning_probe'
 call rt_model_settings(.true.,owner)
 call rt_open(ctx,'probe',status,model_settings=owner)
 if(status==0.or.associated(ctx))then
  write(error_unit,*)'explicit-owner status/context:',status,associated(ctx),'theory:',trim(theory)
  call MPI_Abort(MPI_COMM_WORLD,92,ierr)
 endif
 write(*,*)'explicit owner attempted before any fallback'
 write(error_unit,'(a)')'explicit-owner-complete'
 flush(error_unit)
 call rt_model_settings(.true.)
 call rt_model_settings(.true.)
 call rt_open(ctx,'probe',status)
 if(status==0.or.associated(ctx))call MPI_Abort(MPI_COMM_WORLD,93,ierr)
 status=detach()
 if(status/=0)call MPI_Abort(MPI_COMM_WORLD,94,ierr)
 call MPI_Finalize(ierr)
end program
#elif defined(SALMONPY_PARSER_HOST_PROBE)
program parser_status_only
 use mpi
 use iso_c_binding
 use inputoutput,only:read_input,s_input_parser
 use math_constants,only:set_math_constants
 implicit none
 interface
  integer(c_int) function attach(comm) bind(C,name='sp_attach')
   import c_int,c_int64_t
   integer(c_int64_t),value::comm
  end function
  integer(c_int) function detach() bind(C,name='sp_detach')
   import c_int
  end function
 end interface
 integer::ierr,status
 type(s_input_parser)::parser
 call MPI_Init(ierr)
 status=attach(int(MPI_COMM_WORLD,c_int64_t))
 if(status/=0)call MPI_Abort(MPI_COMM_WORLD,91,ierr)
 call set_math_constants
 call read_input('inputfile',parser,status)
 if(status/=1010)then
  write(*,*)'status-only Maxwell returned',status,'expected',1010
  call MPI_Abort(MPI_COMM_WORLD,92,ierr)
 endif
 call read_input('inputfile',parser,status,native_host=.false.)
 if(status/=1010)call MPI_Abort(MPI_COMM_WORLD,94,ierr)
 call read_input('inputfile',parser,status,native_host=.true.)
 if(status/=0)call MPI_Abort(MPI_COMM_WORLD,95,ierr)
 write(*,*)'parser explicit native-host contract PASS'
 status=detach()
 if(status/=0)call MPI_Abort(MPI_COMM_WORLD,93,ierr)
 call MPI_Finalize(ierr)
end program
#elif defined(SALMONPY_ABI_PROTOCOL_PROBE)
program test_native_protocol
 use iso_c_binding
 use mpi
 implicit none
 interface
#ifdef SALMONPY_ABI_FAST_HSE
  integer(c_int) function one() bind(C,name="sp_fast_hse_handle_protocol")
#else
  integer(c_int) function one() bind(C,name="sp_tddft_handle_protocol")
#endif
   import c_int
  end function
#ifdef SALMONPY_ABI_FAST_HSE
  integer(c_int) function two() bind(C,name="sp_fast_hse_model_protocol")
#else
  integer(c_int) function two() bind(C,name="sp_tddft_model_protocol")
#endif
   import c_int
  end function
#ifdef SALMONPY_ABI_FAST_HSE
  integer(c_int) function three() bind(C,name="sp_fast_hse_registry_protocol")
#else
  integer(c_int) function three() bind(C,name="sp_tddft_registry_protocol")
#endif
   import c_int
  end function
#ifdef SALMONPY_ABI_FAST_HSE
  integer(c_int) function four() bind(C,name="sp_fast_hse_checkpoint_protocol")
#else
  integer(c_int) function four() bind(C,name="sp_eh_checkpoint_protocol")
#endif
   import c_int
  end function
 end interface
 integer(c_int)::a,b,c,d
 integer::ierr
 logical::initialized
 call MPI_Initialized(initialized,ierr)
 if(ierr/=MPI_SUCCESS.or.initialized)stop 41
 a=one();b=two();c=three();d=four()
 if(a/=1.or.b/=1.or.c/=1.or.d/=1)stop 42
 call MPI_Initialized(initialized,ierr)
 if(ierr/=MPI_SUCCESS.or.initialized)stop 43
 write(*,'(4(i0,1x))')a,b,c,d
end program
#elif defined(SALMONPY_FIELD_CLOSE_PROBE)
program test_field_close
 use field_close_probe
 implicit none
 integer::status,complete,candidate,saved
 do complete=0,1
 do candidate=0,1
 do saved=0,1
  allocate(state);state%step=complete;state%capacity=1;state%fs%srg_ng=1
  allocate(state%fe%buffer(2),state%ld%buffer(3))
  if(candidate==1)then
   allocate(pending);allocate(pending%buffer(4))
  endif
  if(saved==1)then
   allocate(step_saved);allocate(step_saved%buffer(5))
  endif
  finalized=0
  status=sp_eh_close()
  if(status/=0)stop 11
  if(allocated(state).or.allocated(pending).or.allocated(step_saved))stop 12
  if(finalized/=complete)stop 13
  status=sp_eh_close()
  if(status/=1.or.allocated(pending).or.allocated(step_saved))stop 14
 enddo
 enddo
 enddo
 print *, 'field close early/full/candidate/transaction/repeated PASS'
end program
#else
#ifdef SALMONPY_SCALAR_IEEE_PROBE
program test_scalar_ieee
 use scalar_ieee_probe
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_positive_inf,ieee_negative_inf
 implicit none
 real(8)::a(3),b(2,3),c(2,2,2),d(2,2,2,2),bad(3),empty(0)
 type(s_factor_history)::h
 integer::i,j,k,l,n
 logical::valid
 bad=[ieee_value(0d0,ieee_quiet_nan),ieee_value(0d0,ieee_positive_inf),ieee_value(0d0,ieee_negative_inf)]
 a=1d0;b=2d0;c=3d0;d=4d0
 valid=.true.
 call accumulate_finite_1d(a,valid);call accumulate_finite_2d(b,valid)
 call accumulate_finite_3d(c,valid);call accumulate_finite_4d(d,valid)
 if(.not.valid)stop 11
 call accumulate_finite_1d(empty,valid)
 if(.not.valid)stop 12
 valid=.false.;call accumulate_finite_1d(a,valid)
 if(valid)stop 13
 do n=1,3
  do i=1,3
   a=1d0;a(i)=bad(n);valid=.true.;call accumulate_finite_1d(a,valid)
   if(valid)stop 14
  enddo
  do j=1,3
  do i=1,2
   b=2d0;b(i,j)=bad(n);valid=.true.;call accumulate_finite_2d(b,valid)
   if(valid)stop 15
  enddo
  enddo
  do k=1,2
  do j=1,2
  do i=1,2
   c=3d0;c(i,j,k)=bad(n);valid=.true.;call accumulate_finite_3d(c,valid)
   if(valid)stop 16
   do l=1,2
    d=4d0;d(i,j,k,l)=bad(n);valid=.true.;call accumulate_finite_4d(d,valid)
    if(valid)stop 17
   enddo
  enddo
  enddo
  enddo
  h=s_factor_history()
  if(.not.valid_history(h))stop 18
  do j=1,2
  do i=1,2
   h%g=0d0;h%g(i,j)=bad(n)
   if(valid_history(h))stop 19
  enddo
  enddo
  h%g=0d0
  do i=1,2
   h%b=0d0;h%b(i)=bad(n)
   if(valid_history(h))stop 20
  enddo
  h%b=0d0
  do i=1,2
   h%coeff=0d0;h%coeff(i)=bad(n)
   if(valid_history(h))stop 21
  enddo
 enddo
 h=s_factor_history();h%count=-1
 if(valid_history(h))stop 22
 print *, 'scalar IEEE finite/NaN/+Inf/-Inf/all positions/history PASS'
end program
#else
#ifdef SALMONPY_CHECKPOINT_SCOPE_PROBE
program test_checkpoint_scope
 use structures,only:s_scalar,s_vector
 use checkpoint_scope_probe
 implicit none
 type(s_scalar)::small,large,a,b,rejected
 type(s_vector)::vector,back
 integer::u,v,ios
 allocate(small%f(-1:0,2:3,0:1));small%f=0.125d0
 allocate(large%f(0:7,-2:1,1:3));large%f=-0.75d0
 allocate(vector%v(-1:0,0:1,1:2,1:3));vector%v=0.5d0
 open(newunit=u,file='small.bin',access='stream',form='unformatted',status='new')
 ios=0;call checkpoint_io_s_scalar(u,.false.,small,ios,8_8)
 if(ios/=0)stop 11
 rewind(u);ios=0;call checkpoint_io_s_scalar(u,.true.,a,ios,8_8)
 if(ios/=0.or.any(a%f/=small%f).or.any(lbound(a%f)/=lbound(small%f)))stop 12
 open(newunit=v,file='large.bin',access='stream',form='unformatted',status='new')
 ios=0;call checkpoint_io_s_scalar(v,.false.,large,ios,96_8)
 if(ios/=0)stop 13
 rewind(v);ios=0;call checkpoint_io_s_scalar(v,.true.,rejected,ios,8_8)
 if(ios==0.or.allocated(rejected%f))stop 14
 rewind(v);ios=0;call checkpoint_io_s_scalar(v,.true.,b,ios,96_8)
 if(ios/=0.or.any(b%f/=large%f).or.any(lbound(b%f)/=lbound(large%f)))stop 15
 rewind(u);ios=0;call checkpoint_io_s_scalar(u,.true.,a,ios,8_8)
 if(ios/=0.or.any(a%f/=small%f))stop 16
 close(u);close(v)
 open(newunit=u,file='vector.bin',access='stream',form='unformatted',status='new')
 ios=0;call checkpoint_io_s_vector(u,.false.,vector,ios,24_8)
 if(ios/=0)stop 17
 rewind(u);ios=0;call checkpoint_io_s_vector(u,.true.,back,ios,24_8)
 if(ios/=0.or.any(back%v/=vector%v).or.any(lbound(back%v)/=lbound(vector%v)))stop 18
 close(u)
 print *, 'checkpoint call-local bounds/small-large-small/reject/roundtrip PASS'
end program
#else
#ifdef SALMONPY_HOST_UNIT_PROBE
module salmonpy_test_host_unit
#ifdef SALMONPY_HOST_TIMER_PROBE
 use timer,only:timer_initialize,timer_set,timer_get,LOG_TOTAL,LOG_INIT_RT,LOG_CALC_RHO,LOG_CALC_HARTREE
#endif
#ifdef SALMONPY_HOST_MPI_PROBE
 use communication,only:comm_begin_children,comm_end_children,comm_release_children,comm_create_group,comm_get_globalinfo
#endif
 use iso_c_binding
 implicit none
#ifdef SALMONPY_HOST_MPI_PROBE
 integer::test_owner=0
#endif
contains
#ifdef SALMONPY_HOST_MPI_PROBE
 integer(c_int) function test_host_mpi(mode) bind(C)
  integer(c_int),value::mode
  integer::status,group,rank,n,child
  test_host_mpi=1
  if(mode==0)then
   call comm_begin_children(test_owner,status)
   if(status/=0)return
   call comm_get_globalinfo(group,rank,n)
   child=comm_create_group(group,0,rank)
   call comm_end_children(test_owner)
  else if(mode==1)then
   call comm_release_children(test_owner)
   test_owner=0
  else
   return
  endif
  test_host_mpi=0
 end function
#endif
#ifdef SALMONPY_HOST_TIMER_PROBE
 integer(c_int) function test_host_timers(mode) bind(C)
  integer(c_int),value::mode
  test_host_timers=1
  if(mode==0)then
   call timer_initialize
   call timer_set(LOG_TOTAL,7d0);call timer_set(LOG_INIT_RT,11d0)
   call timer_set(LOG_CALC_RHO,13d0);call timer_set(LOG_CALC_HARTREE,17d0)
  else if(mode/=1)then
   return
  endif
  if(timer_get(LOG_TOTAL)/=7d0.or.timer_get(LOG_INIT_RT)/=11d0)return
  if(timer_get(LOG_CALC_RHO)/=13d0.or.timer_get(LOG_CALC_HARTREE)/=17d0)return
  test_host_timers=0
 end function
#endif
 integer(c_int) function test_host_unit(path,mode) bind(C)
  character(c_char),intent(in)::path(*)
  integer(c_int),value::mode
  character(4096)::filename
  integer::i,ios
  logical::opened
  test_host_unit=1
  inquire(unit=30,opened=opened,iostat=ios)
  if(ios/=0)return
  select case(mode)
  case(0)
   if(opened)return
   filename=''
   do i=1,len(filename)
    if(path(i)==c_null_char)exit
    filename(i:i)=path(i)
   enddo
   if(i>len(filename))return
   open(unit=30,file=trim(filename),status='new',action='write',iostat=ios)
   if(ios/=0)return
   write(30,'(a)',iostat=ios)'host-before'
  case(1)
   if(.not.opened)return
   write(30,'(a)',iostat=ios)'host-after'
  case(2)
   if(.not.opened)return
   close(30,iostat=ios)
  case default
   return
  end select
  if(ios==0)test_host_unit=0
 end function
end module
#else
! Synthetic fixtures only: no GS wavefunction/density data or SCF calculation.
program test_io_ownership
 use structures
 use salmon_global
 use inputoutput,only:initialize_inputoutput_units,uenergy_from_au,ulength_from_au, &
   utime_from_au,t_unit_energy,au_energy_ev,au_length_aa,au_time_fs,nsplit_voxel_data
 use salmonpy_model_inputs,only:s_model_inputs
 use parallelization,only:setup_parallel,end_parallel,nproc_id_global,nproc_group_global
 use math_constants,only:set_math_constants
 use write_sub,only:write_dos,write_pdos
 use band_dft_sub,only:write_band,calc_band_write
 use write_file3d,only:write_avs
 implicit none
 type(s_model_inputs)::a,b,bad
 type(s_dft_system)::system
 type(s_dft_energy)::energy
 type(s_rgrid)::lg,mg
 type(s_parallel_info)::info
 type(s_pp_info)::pp
 type(s_orbital)::psi
 type(s_band_dft)::band
 integer::status,j,unit,i,host(4)=[1,77,100,101]
 character(32)::mode,variant
 character(80)::names(4),actual
 character(60)::suffix
 character(20)::header
 logical::opened,reference
 real(8)::field(2,1,1)
 call setup_parallel
 call set_math_constants
 call get_command_argument(1,mode)
 call get_command_argument(2,variant)
 reference=variant=='reference'
 unit_system='au';unit_time='au';unit_length='au';unit_energy='au';unit_charge='au'
 call initialize_inputoutput_units
 call a%capture()
 unit_time='fs';unit_length='AA';unit_energy='eV'
 call initialize_inputoutput_units
 call b%capture()
 if(mode=='units')then
   call a%activate(status)
   if(status/=0)stop 11
   if(uenergy_from_au/=1d0.or.ulength_from_au/=1d0.or.utime_from_au/=1d0)stop 12
   if(t_unit_energy%name/='a.u.')stop 13
   call b%activate(status)
   if(status/=0)stop 14
   if(uenergy_from_au/=au_energy_ev.or.ulength_from_au/=au_length_aa.or.utime_from_au/=au_time_fs)stop 15
   if(t_unit_energy%name/='eV')stop 16
   call a%activate(status)
   if(status/=0.or.uenergy_from_au/=1d0)stop 17
   unit_energy='invalid';call bad%capture()
   unit_energy='au'
   call bad%activate(status)
   if(status==0.or.unit_energy/='au'.or.uenergy_from_au/=1d0)stop 19
 else
   ! Exercise byte-preserving output, separately from model-switch conversion.
   unit_time='au';unit_length='au';unit_energy='au';unit_charge='au'
   call initialize_inputoutput_units
   if(.not.reference)then
     do j=1,4
       write(names(j),'(a,i0,a,i0)')'host-',nproc_id_global,'-',host(j)
       open(unit=host(j),file=trim(names(j)),status='replace')
       write(host(j),'(a)')'before'
     enddo
   endif
   system%nspin=1;system%nk=1;system%no=2;system%Hvol=.5d0
   allocate(system%wtk(1),system%vec_k(3,1),system%Rion(3,1))
   system%wtk=1d0;system%vec_k=0d0;system%Rion=0d0
   system%primitive_b=0d0
   do j=1,3
     system%primitive_b(j,j)=1d0
   enddo
   allocate(energy%esp(2,1,1));energy%esp(:,1,1)=[.25d0,.75d0]
   base_directory='./';sysname='synthetic';natom=1;nelec=2;nstate=2;yn_spinorbit='n'
   temperature=-1d0;yn_out_dos_set_fe_origin='n';out_dos_function='gaussian'
   out_dos_start=0d0;out_dos_end=1d0;out_dos_width=.1d0;out_dos_nenergy=5
   select case(trim(mode))
   case('dos')
     call write_dos(system,energy)
   case('band')
     allocate(band%band_kpt(3,1),band%check_conv_esp(2,1,1));band%band_kpt=0d0
     info%ik_s=1;info%ik_e=1
     call calc_band_write(1,system,band,info)
     call write_band(system,energy)
   case('avs','avs-unsplit')
     nsplit_voxel_data=2
     if(mode=='avs-unsplit')nsplit_voxel_data=1
     lg%is=[1,1,1];lg%ie=[2,1,1];lg%num=[2,1,1]
     field(:,1,1)=[.125d0,.25d0];suffix='synthetic-field';header='fixture'
     call write_avs(lg,77,suffix,header,field)
   case('pdos')
     allocate(kion(1));kion=1
     lg%is=1;lg%ie=1;lg%num=1;mg%is=1;mg%ie=1;mg%num=1
     allocate(lg%coordinate(1,3));lg%coordinate(1,:)=[1d0,0d0,0d0]
     pp%nrmax=3;allocate(pp%mlps(1),pp%upp_f(3,0:0,1),pp%rad(4,1))
     pp%mlps=0;pp%upp_f=1d0;pp%rad(:,1)=[0d0,2d0,4d0,6d0]
     allocate(psi%rwf(1,1,1,1,2,1,1));psi%rwf=1d0
     info%ik_s=1;info%ik_e=1;info%io_s=1;info%io_e=2
     info%icomm_r=nproc_group_global;info%icomm_ko=nproc_group_global
     call write_pdos(lg,mg,system,info,pp,energy,psi)
   case default
     stop 18
   end select
   if(.not.reference)then
     do j=1,4
       inquire(unit=host(j),opened=opened,name=actual)
       if(.not.opened)stop 21
       if(trim(actual)/=trim(names(j)))stop 22
       write(host(j),'(a)')'after'
       close(host(j))
     enddo
   endif
 endif
 call end_parallel
 print *, 'IO ownership PASS ',trim(mode)
end program

#endif

#endif

#endif

#endif
