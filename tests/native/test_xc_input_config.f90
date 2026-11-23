! Alternate owned functional configurations; legacy inputs must not affect evaluation.
program test_xc_input_config
 use structures,only:s_xc_config,s_xc_functional,s_pp_info
 use salmon_xc,only:init_xc_config,calc_xc,finalize_xc
 use salmon_global,only:global_xc=>xc,global_omega=>hse_omega,exx_pre_scf_active
 implicit none
 type(s_xc_config)::config
 type(s_xc_functional)::a,b
 type(s_pp_info)::pp
 real(8)::rho(2,2,2),grad(2,2,2,3),ea(2,2,2),eb(2,2,2),again(2,2,2)
 real(8)::v(2,2,2),dg(2,2,2,3)
 integer::status
 pp%flag_nlcc=.false.;rho=0.2d0;grad=0.03d0
 config%xcname='hse06';config%hse_omega=0.11d0
 call init_xc_config(a,config,status)
 if(status/=0)stop 1
 config%xcname='pbeh40'
 call init_xc_config(b,config,status)
 if(status/=0)stop 2
 if(a%config%exchange_fraction/=0.25d0.or.b%config%exchange_fraction/=0.4d0)stop 3
 call calc_xc(a,pp,rho=rho,grho=grad,eexc=ea,vxc=v,rdedd=dg)
 call calc_xc(b,pp,rho=rho,grho=grad,eexc=eb,vxc=v,rdedd=dg)
 if(maxval(abs(ea-eb))<1d-10)stop 4
 global_xc='pz';global_omega=7d0;exx_pre_scf_active=.true.
 call calc_xc(a,pp,rho=rho,grho=grad,eexc=again,vxc=v,rdedd=dg)
 if(maxval(abs(again-ea))>1d-13)stop 5
 call calc_xc(b,pp,rho=rho,grho=grad,eexc=again,vxc=v,rdedd=dg)
 if(maxval(abs(again-eb))>1d-13)stop 6
 call calc_xc(a,pp,rho=rho,grho=grad,eexc=ea,vxc=v,rdedd=dg,pre_scf_active=.true.)
 call calc_xc(b,pp,rho=rho,grho=grad,eexc=eb,vxc=v,rdedd=dg,pre_scf_active=.true.)
 if(maxval(abs(ea-eb))>1d-13)stop 7
 ! Recreate one functional with a different screening parameter.
 call finalize_xc(b)
 config%xcname='hse06';config%hse_omega=.3d0
 call init_xc_config(b,config,status)
 if(status/=0)stop 8
 call calc_xc(a,pp,rho=rho,grho=grad,eexc=ea,vxc=v,rdedd=dg)
 call calc_xc(b,pp,rho=rho,grho=grad,eexc=eb,vxc=v,rdedd=dg)
 if(maxval(abs(ea-eb))<1d-10)stop 9
 config%hse_omega=-1d0
 call init_xc_config(b,config,status)
 if(status==0)stop 10
 call calc_xc(b,pp,rho=rho,grho=grad,eexc=again,vxc=v,rdedd=dg)
 if(maxval(abs(again-eb))>1d-13)stop 11
 call finalize_xc(a);call finalize_xc(b)
 print *, 'owned XC configuration passed'
end program
