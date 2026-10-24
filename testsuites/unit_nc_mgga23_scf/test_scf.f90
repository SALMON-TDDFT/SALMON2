program test_scf
 use structures
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_positive_inf
 use noncollinear_module,only:calc_dm_noncollinear,simple_mixing_so
 use salmon_global,only:xc
 use sym_sub,only:use_symmetry
 use nc_mgga23_scf,only:nc23_scf_prepare,nc23_scf_residual,nc23_spin_density_residual
 implicit none
 type(s_nc_mgga23_payload)::p,history
 type(s_dft_system)::system
 type(s_rgrid)::mg
 type(s_parallel_info)::info
 type(s_orbital)::psi
 type(s_scalar)::rho(2)
 type(s_mixing)::mixing
 complex(8)::before(1,1,1,2,2),after(1,1,1,2,2)
 real(8)::res,x
 character(32)::mode
 integer::i
 call get_command_argument(1,mode)
 p%active=.true.
 allocate(p%local(1,1,1,4),p%tau(1,1,1,4),p%current(1,1,1,4,3))
 p%local=0d0;p%tau=0d0;p%current=0d0
 p%band_expectation=7d0;p%core_deda=8d0
 if(index(mode,'--nan-')==1.or.index(mode,'--inf-')==1)then
  history=p
  select case(trim(mode))
  case('--nan-local');p%local(1,1,1,4)=ieee_value(0d0,ieee_quiet_nan)
  case('--nan-tau');p%tau(1,1,1,4)=ieee_value(0d0,ieee_quiet_nan)
  case('--nan-current');p%current(1,1,1,4,3)=ieee_value(0d0,ieee_quiet_nan)
  case('--inf-tau');p%tau(1,1,1,4)=ieee_value(0d0,ieee_positive_inf)
  case('--nan-history');history%current(1,1,1,4,3)=ieee_value(0d0,ieee_quiet_nan)
  endselect
  if(mode/='--nan-density'.and.mode/='--inf-density')then
   call nc23_scf_residual(p,history,0,res)
   error stop 'invalid coefficient accepted'
  endif
 endif
 call nc23_scf_prepare(p,history,.15d0)
 do i=1,100
  x=p%tau(1,1,1,4)
  p%tau(1,1,1,4)=1d0-3d0*x
  call nc23_scf_residual(p,history,0,res)
  call nc23_scf_prepare(p,history,.15d0)
 enddo
 if(abs(p%tau(1,1,1,4)-.25d0)>1d-12.or.res>1d-12) error stop 'unstable tau map not damped'
 p%current(1,1,1,4,3)=1d0
 call nc23_scf_residual(p,history,0,res)
 if(res<.9d0) error stop 'charge-independent spin-current residual missed'
 call nc23_scf_prepare(p,history,.15d0)
 if(abs(p%current(1,1,1,4,3)-.15d0)>1d-14) error stop 'spin current not mixed'
 if(p%band_expectation/=7d0.or.any(p%core_deda/=8d0)) error stop 'physical energy/current changed'
 before=(0d0,0d0);before(1,1,1,1,1)=1d0;before(1,1,1,2,2)=1d0
 after=before;after(1,1,1,1,2)=(0d0,.1d0);after(1,1,1,2,1)=(0d0,-.1d0)
 if(mode=='--inf-density')after(1,1,1,2,1)=cmplx(0d0,ieee_value(0d0,ieee_positive_inf),8)
 if(mode=='--nan-density')after(1,1,1,2,1)=cmplx(0d0,ieee_value(0d0,ieee_quiet_nan),8)
 call nc23_spin_density_residual(after,before,1d0,.15d0,0,res)
 if(abs(res-2d0/3d0)>1d-13) error stop 'charge-invisible LSDA spin change missed'
 mg%is=1;mg%ie=1;system%nspin=2;system%no=1;system%nk=1
 info%io_s=1;info%io_e=1;info%ik_s=1;info%ik_e=1;info%im_s=1;info%im_e=1;info%icomm_ko=0
 allocate(system%rocc(1,1,2),system%wtk(1),psi%zwf(1,1,1,2,1,1,1))
 allocate(rho(1)%f(1,1,1),rho(2)%f(1,1,1));system%rocc=1d0;system%wtk=1d0
 xc='pz';use_symmetry=.false.
 psi%zwf(1,1,1,:,1,1,1)=[cmplx(1d0,0d0,8),cmplx(0d0,1d0,8)]/sqrt(2d0)
 call calc_dm_noncollinear(psi,system,info,mg)
 call simple_mixing_so(mg,system,0d0,1d0,rho,mixing)
 if(abs(rho(1)%f(1,1,1)-1d0)>1d-13.or.abs(rho(2)%f(1,1,1))>1d-13) &
  error stop 'transverse LSDA eigenchannels replaced by lab diagonals'
 print *, 'PASS full gKS damping, tau fixed point, hidden spin-current residual, physical diagnostics preserved'
end program
