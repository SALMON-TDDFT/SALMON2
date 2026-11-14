program probe
use structures
use salmon_global
use current_field, only: new_field=>calc_emfields
use reference_field, only: old_field=>calc_emfields
use current_taylor, only: new_taylor=>taylor
use reference_taylor, only: old_taylor=>taylor
use, intrinsic :: ieee_arithmetic
implicit none
type(s_rt_step_config) :: config,saved
type(s_rt) :: rt,ref
type(s_rgrid) :: mg
type(s_dft_system) :: system
type(s_parallel_info) :: info
type(s_stencil) :: stencil
type(s_sendrecv_grid) :: srg
type(s_pp_grid) :: ppg
type(s_scalar) :: v(1)
type(s_orbital) :: input,output,work,rin,rout,rwork
real(8) :: invalid(5),step,eps(2),curr(3,1),nan
character(2) :: mode
integer :: status,k,j,n,m,bad
nan=ieee_value(0d0,ieee_quiet_nan)
call configure_rt_step(config,.001d0,4,'middlepoint','tr',0d0,[1d0,1d0],status)
if(status/=0)stop 1
saved=config
invalid=[nan,ieee_value(0d0,ieee_positive_inf),ieee_value(0d0,ieee_negative_inf),-1d0,0d0]
do k=1,size(invalid)
 call configure_rt_step(config,invalid(k),4,'middlepoint','tr',0d0,[1d0,1d0],status)
 if(status==0.or.config%dt/=saved%dt)stop 2
end do
call configure_rt_step(config,0d0,4,'middlepoint','tr',0d0,[1d0,1d0],status)
if(status==0)stop 3
call configure_rt_step(config,.001d0,0,'middlepoint','tr',0d0,[1d0,1d0],status)
if(status==0)stop 4
call configure_rt_step(config,.001d0,4,'bad','tr',0d0,[1d0,1d0],status)
if(status==0)stop 5
call configure_rt_step(config,.001d0,6,'hse_taylor4','tr',0d0,[1d0,1d0],status)
if(status==0)stop 13
call configure_rt_step(config,.001d0,4,'middlepoint','xx',0d0,[1d0,1d0],status)
if(status==0)stop 6
do j=1,2
 do k=1,size(invalid)
  eps=1d0;eps(j)=invalid(k)
  call configure_rt_step(config,.001d0,4,'middlepoint','tr',0d0,eps,status)
  if(status==0)stop 7
 end do
end do
do k=1,size(invalid)-1
 call configure_rt_step(config,.001d0,4,'middlepoint','tr',invalid(k),[1d0,1d0],status)
 if(status==0)stop 8
end do
if(config%dt/=saved%dt.or.config%taylor_order/=saved%taylor_order.or.config%propagator/=saved%propagator)stop 9
if(config%field_mode/=saved%field_mode.or.config%film_thickness/=saved%film_thickness.or.any(config%epsilon_em/=saved%epsilon_em))stop 14
! Actual field kernel: each rt owns a different dt/mode/material-film profile.
do k=1,2
 step=.001d0*k
 do j=1,3
  select case(j)
  case(1);mode='tr'
  case(2);mode='lo'
  case(3);mode='2d'
  end select
  call configure_rt_step(rt%step,step,4,'middlepoint',mode,2d0*k,[1d0,4d0],status)
  if(status/=0)stop 10
  allocate(rt%Ac_ext(3,0:3),rt%Ac_ind(3,0:3),rt%Ac_tot(3,0:3))
  allocate(rt%E_ext(3,0:3),rt%E_ind(3,0:3),rt%E_tot(3,0:3),rt%curr(3,0:3))
  rt%Ac_ext=.002d0;rt%Ac_ext(:,0)=0d0;rt%Ac_ind=0d0;rt%Ac_tot=rt%Ac_ext
  rt%E_ext=0d0;rt%E_ind=0d0;rt%E_tot=0d0;rt%curr=0d0
  ref=rt;curr(:,1)=[.1d0,-.02d0,.03d0]
  dt=step;trans_longi=mode;film_thickness=2d0*k;epsilon_em=[1d0,4d0]
  call old_field(1,1,curr,ref)
  dt=nan;trans_longi='xx';film_thickness=nan;epsilon_em=-1d0
  call new_field(1,1,curr,rt)
  if(any(rt%Ac_tot/=ref%Ac_tot).or.any(rt%E_tot/=ref%E_tot).or.any(rt%curr/=ref%curr))stop 11
  deallocate(rt%Ac_ext,rt%Ac_ind,rt%Ac_tot,rt%E_ext,rt%E_ind,rt%E_tot,rt%curr)
  deallocate(ref%Ac_ext,ref%Ac_ind,ref%Ac_tot,ref%E_ext,ref%E_ind,ref%E_tot,ref%curr)
 end do
end do
! Actual Taylor loops; linear hpsi fixture isolates dependence on n_hamil.
mg%is=1;mg%ie=2;system%nspin=1
info%ik_s=1;info%ik_e=1;info%io_s=1;info%io_e=1
allocate(input%zwf(2,2,2,1,1,1,1),output%zwf(2,2,2,1,1,1,1),work%zwf(2,2,2,1,1,1,1))
do k=1,2
 n=2+2*k
 call configure_rt_step(rt%step,.001d0*k,n,'middlepoint','tr',0d0,[1d0,1d0],status)
 allocate(rt%zc(n))
 do j=1,n
  rt%zc(j)=(-cmplx(0d0,1d0,8)*rt%step%dt)**j
  do m=1,j
   rt%zc(j)=rt%zc(j)/m
  end do
 end do
 input%zwf=cmplx(.1d0,.2d0,8);output%zwf=0d0;work%zwf=0d0
 rin=input;rout=output;rwork=work
 n_hamil=n
 call old_taylor(mg,system,info,stencil,srg,rin,rout,rwork,ppg,v,rt)
 n_hamil=1
 call new_taylor(mg,system,info,stencil,srg,input,output,work,ppg,v,rt)
 if(any(output%zwf/=rout%zwf))stop 12
 deallocate(rt%zc)
end do
print *, 'RT step config: scalar rejection, two dt, field tr/lo/2d and Taylor order independence PASS'
end program
