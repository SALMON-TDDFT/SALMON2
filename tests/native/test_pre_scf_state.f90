! Independent pre-SCF controllers compared with the original count/switch rules.
program probe
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
 use structures,only:s_dft_system
 use exx_native,only:exx_pre_scf_initialize,exx_pre_scf_begin,exx_pre_scf_not_ready, &
   exx_pre_scf_count_ready,exx_pre_scf,exx_set_phase
 use salmon_global,only:exx_pre_scf_active,exx_pre_scf_threshold,exx_pre_scf_steps
 implicit none
 type(s_dft_system)::a,b,disabled
 real(8)::ra(6)=[.05d0,.05d0,.2d0,.05d0,.05d0,.05d0]
 integer::i,status,ca,cb,na,nb,sw_a,sw_b
 logical::ready,active_a,active_b
 call exx_pre_scf_initialize(a,.1d0,2,status);if(status/=0)stop 1
 call exx_pre_scf_initialize(b,.2d0,3,status);if(status/=0)stop 2
 call exx_pre_scf_initialize(disabled,0d0,1,status)
 if(status/=0.or.exx_pre_scf(disabled))stop 3
 call exx_pre_scf_begin(a,0);call exx_pre_scf_begin(b,4)
 ca=0;cb=0;sw_a=0;sw_b=0;active_a=.true.;active_b=.true.
 do i=1,6
  exx_pre_scf_active=mod(i,2)==0;exx_pre_scf_threshold=99d0;exx_pre_scf_steps=99
  if(active_a)then
   na=0
   if(i==1.or..not.(ra(i)<.1d0))na=1
   if(exx_pre_scf_not_ready(a,i,ra(i))/=na)stop 4
   if(na==0)then
    ca=ca+1
   else
    ca=0
   endif
   ready=exx_pre_scf_count_ready(a,na)
   if(ready.neqv.(ca>=2))stop 5
   if(ready)then
    call exx_set_phase(a,.false.);active_a=.false.;sw_a=i
   endif
  endif
  if(active_b)then
   nb=0
   if(i==1)nb=1
   if(exx_pre_scf_not_ready(b,i+4,.1d0)/=nb)stop 6
   ! Model another DC fragment rejecting readiness at the second iteration.
   if(i==2)nb=nb+1
   if(nb==0)then
    cb=cb+1
   else
    cb=0
   endif
   ready=exx_pre_scf_count_ready(b,nb)
   if(ready.neqv.(cb>=3))stop 7
   if(ready)then
    call exx_set_phase(b,.false.);active_b=.false.;sw_b=i+4
   endif
  endif
  if(exx_pre_scf(a).neqv.active_a)stop 8
  if(exx_pre_scf(b).neqv.active_b)stop 9
 enddo
 if(sw_a/=5.or.sw_b/=9)stop 10
 call exx_pre_scf_initialize(a,ieee_value(0d0,ieee_quiet_nan),1,status)
 if(status==0.or.exx_pre_scf(a).or.a%pre_scf%threshold/=.1d0)stop 11
 call exx_pre_scf_initialize(a,.2d0,0,status)
 if(status==0.or.a%pre_scf%required_steps/=2)stop 12
 call exx_pre_scf_begin(a,10)
 if(a%pre_scf%count/=0.or.a%pre_scf%first_iteration/=11.or.exx_pre_scf(a))stop 13
 print *,'Owned pre-SCF state transitions PASS'
end program
