module structures
  type s_rgrid; integer::num(3)=[1,1,1]; end type
  type s_dft_system; integer::nspin=1;real(8)::hvol=1d0; end type
  type s_parallel_info; integer::icomm_r=0; end type
  type s_stencil; end type
  type s_xc_functional
    logical::use_gradient=.false.,use_laplacian=.false.,use_kinetic_energy=.false.,use_current=.false.
  end type
  type s_pp_info; end type
  type s_pp_nlcc; end type
  type s_orbital; end type
  type s_sendrecv_grid; end type
  type s_poisson; end type
  type s_reciprocal_grid; end type
  type s_scalar; real(8),allocatable::f(:,:,:); end type
end module
module communication
  integer::forced_integer_bad=0,integer_sum_calls=0
  interface comm_summation
    module procedure sum_i,sum_r
  end interface
contains
  subroutine sum_i(a,b,c);integer,intent(in)::a,c;integer,intent(out)::b;integer_sum_calls=integer_sum_calls+1;b=a+forced_integer_bad;end
  subroutine sum_r(a,b,c);real(8),intent(in)::a;integer,intent(in)::c;real(8),intent(out)::b;b=a;end
end module
module salmon_global
  character(1)::yn_hse='n'
end module
module hartree_sub
  integer::hartree_calls=0
contains
  subroutine hartree(lg,mg,info,system,fg,poisson,srg,stencil,rho,vh)
    use structures
    type(s_rgrid)::lg,mg;type(s_parallel_info)::info;type(s_dft_system)::system
    type(s_reciprocal_grid)::fg;type(s_poisson)::poisson;type(s_sendrecv_grid)::srg
    type(s_stencil)::stencil;type(s_scalar)::rho,vh
    hartree_calls=hartree_calls+1
    vh%f=.5d0*rho%f
  end subroutine
end module
module salmon_xc
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  integer::xc_calls=0
  logical::force_xc_nan=.false.
  real(8),allocatable::eexc_tmp(:,:,:)
contains
  subroutine exchange_correlation(system,xc,mg,srgs,srg,rhos,pp,ppn,info,psi,stencil,vxc,exc)
    use structures
    type(s_dft_system)::system;type(s_xc_functional)::xc;type(s_rgrid)::mg
    type(s_sendrecv_grid)::srgs,srg;type(s_scalar)::rhos(system%nspin),vxc(system%nspin)
    type(s_pp_info)::pp;type(s_pp_nlcc)::ppn;type(s_parallel_info)::info
    type(s_orbital)::psi;type(s_stencil)::stencil;real(8)::exc
    integer::i
    xc_calls=xc_calls+1
    if(.not.allocated(eexc_tmp))allocate(eexc_tmp,source=rhos(1)%f)
    eexc_tmp=-.45d0*rhos(1)%f**(4d0/3d0)
    vxc(1)%f=-.6d0*rhos(1)%f**(1d0/3d0)
    if(force_xc_nan)vxc(1)%f(1,1,1)=ieee_value(0d0,ieee_quiet_nan)
    exc=sum(eexc_tmp)*system%hvol
  end subroutine
end module
program dg_wpw_lda_hartree_driver
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use dg_wpw_lda_hartree
  use structures
  use hartree_sub,only:hartree_calls
  use salmon_xc,only:xc_calls,force_xc_nan
  use salmon_global,only:yn_hse
  use communication,only:forced_integer_bad,integer_sum_calls
  implicit none
  integer,parameter::np=6,nf=2
  logical::core(np,nf)
  real(8)::rho(np,nf),vxc(np,nf),exc(np,nf),weight(np),rhog(np),vh(np),excs,nvxc,eh
  integer::info,bad,i,j
  type(s_rgrid)::lg,mg;type(s_dft_system)::system;type(s_parallel_info)::pinfo
  type(s_stencil)::stencil;type(s_xc_functional)::xc;type(s_pp_info)::pp;type(s_pp_nlcc)::ppn
  type(s_orbital)::psi;type(s_sendrecv_grid)::srg,srgs;type(s_poisson)::poisson;type(s_reciprocal_grid)::fg
  type(s_scalar)::rtot,rhos(1),vh_s,vxc_s(1)
  real(8)::rho_core(np,nf,1),exc_core,nvxc_core
  core=.false.;core(1:3,1)=.true.;core(4:6,2)=.true.
  call validate_core_ownership(core,np,nf,info,bad);call req(info==0.and.bad==0,'unique ownership')
  weight=[.2d0,.3d0,.25d0,.35d0,.15d0,.4d0]
  rhog=[.4d0,.7d0,1.1d0,.8d0,.5d0,.3d0];vh=[.2d0,.5d0,.7d0,.6d0,.3d0,.1d0]
  do j=1,nf;do i=1,np
    rho(i,j)=rhog(i);vxc(i,j)=-.6d0*rhog(i)**(1d0/3d0);exc(i,j)=-.45d0*rhog(i)**(4d0/3d0)
    if(.not.core(i,j))then;rho(i,j)=1d90;vxc(i,j)=-1d90;exc(i,j)=1d90;endif
  enddo;enddo
  call integrate_core_lda_terms(rho,vxc,exc,weight,core,np,nf,excs,nvxc,info)
  call req(info==0.and.abs(excs-sum(-.45d0*rhog**(4d0/3d0)*weight))<1d-13,'core Exc')
  call req(abs(nvxc-sum(-.6d0*rhog**(4d0/3d0)*weight))<1d-13,'core nVxc')
  call hartree_energy_global(rhog,vh,weight,np,0,eh,info)
  call req(info==0.and.abs(eh-.5d0*sum(rhog*vh*weight))<1d-13,'global Hartree')
  vh(2)=ieee_value(0d0,ieee_quiet_nan);i=integer_sum_calls
  call hartree_energy_global(rhog,vh,weight,np,0,eh,info);call req(info==1.and.integer_sum_calls==i+1,'Hartree local failure handshake')
  vh(2)=.5d0;forced_integer_bad=1;i=integer_sum_calls
  call hartree_energy_global(rhog,vh,weight,np,0,eh,info);call req(info==3.and.integer_sum_calls==i+1,'Hartree peer failure handshake')
  forced_integer_bad=0
  core(3,:)=.false.;call validate_core_ownership(core,np,nf,info,bad);call req(info==2.and.bad==3,'ownership hole')
  core(3,:)=.true.;call validate_core_ownership(core,np,nf,info,bad);call req(info==2.and.bad==3,'double owner')
  core(3,2)=.false.;rho(2,1)=ieee_value(0d0,ieee_quiet_nan)
  call integrate_core_lda_terms(rho,vxc,exc,weight,core,np,nf,excs,nvxc,info);call req(info==4,'NaN rejection')
  ! Exercise the production wrapper and its fail-before-solver gates.
  core(3,2)=.false.;rho(2,1)=rhog(2);rho_core(:,:,1)=rho;mg%num=[2,3,1];system%hvol=1d0
  allocate(rtot%f(2,3,1),rhos(1)%f(2,3,1),vh_s%f(2,3,1),vxc_s(1)%f(2,3,1))
  call update_wpw_lda_hartree(lg,mg,system,pinfo,stencil,xc,pp,ppn,psi,srg,srgs,poisson,fg, &
    rtot,rhos,vh_s,vxc_s,rho_core,core,np,nf,excs,exc_core,nvxc_core,info,bad)
  call req(info==0.and.hartree_calls==1.and.xc_calls==1,'valid production update')
  call req(maxval(abs(reshape(rhos(1)%f,[np])-rhog))<1d-13,'Fortran grid order')
  call req(abs(excs-exc_core)<1d-13,'actual XC core/global equality')
  call update_wpw_owned_lda_hartree(lg,mg,system,pinfo,stencil,xc,pp,ppn,psi,srg,srgs,poisson,fg,&
    rtot,rhos,vh_s,vxc_s,reshape(rhog,[np,1]),np,excs,exc_core,nvxc_core,info)
  call req(info==0.and.hartree_calls==2.and.xc_calls==2,'owned-grid production update')
  call req(maxval(abs(reshape(rhos(1)%f,[np])-rhog))<1d-13,'owned-grid Fortran order')
  force_xc_nan=.true.
  call update_wpw_owned_lda_hartree(lg,mg,system,pinfo,stencil,xc,pp,ppn,psi,srg,srgs,poisson,fg,&
    rtot,rhos,vh_s,vxc_s,reshape(rhog,[np,1]),np,excs,exc_core,nvxc_core,info)
  call req(info/=0.and.hartree_calls==3.and.xc_calls==3,'owned-grid post-solver NaN rejected')
  force_xc_nan=.false.
  xc%use_gradient=.true.
  call update_wpw_lda_hartree(lg,mg,system,pinfo,stencil,xc,pp,ppn,psi,srg,srgs,poisson,fg, &
    rtot,rhos,vh_s,vxc_s,rho_core,core,np,nf,excs,exc_core,nvxc_core,info,bad)
  call req(info==10.and.hartree_calls==3.and.xc_calls==3,'GGA rejected before solvers');xc%use_gradient=.false.
  yn_hse='y'
  call update_wpw_lda_hartree(lg,mg,system,pinfo,stencil,xc,pp,ppn,psi,srg,srgs,poisson,fg, &
    rtot,rhos,vh_s,vxc_s,rho_core,core,np,nf,excs,exc_core,nvxc_core,info,bad)
  call req(info==10.and.hartree_calls==3.and.xc_calls==3,'HSE rejected before solvers');yn_hse='n'
  core(3,:)=.false.
  call update_wpw_lda_hartree(lg,mg,system,pinfo,stencil,xc,pp,ppn,psi,srg,srgs,poisson,fg, &
    rtot,rhos,vh_s,vxc_s,rho_core,core,np,nf,excs,exc_core,nvxc_core,info,bad)
  call req(info==2.and.hartree_calls==3.and.xc_calls==3,'ownership rejected before solvers')
  core(3,1)=.true.;rho_core(2,1,1)=ieee_value(0d0,ieee_quiet_nan)
  call update_wpw_lda_hartree(lg,mg,system,pinfo,stencil,xc,pp,ppn,psi,srg,srgs,poisson,fg, &
    rtot,rhos,vh_s,vxc_s,rho_core,core,np,nf,excs,exc_core,nvxc_core,info,bad)
  call req(info==11.and.hartree_calls==3.and.xc_calls==3,'NaN rejected before solvers')
  write(*,'(a)')'PASS DG WPW LDA Hartree Fortran driver'
contains
  subroutine req(ok,msg);logical,intent(in)::ok;character(*),intent(in)::msg
    if(.not.ok)then;write(*,*)trim(msg);error stop 1;endif
  end subroutine
end program
