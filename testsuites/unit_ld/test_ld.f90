program test_ld
  use lorentz_drude
  use legacy_ld
  use iso_fortran_env, only: int64
  implicit none
  type(s_ld_state) :: ld, other
  type(s_legacy_ld) :: ref
  integer :: lo(3)=[0,0,0],hi(3)=[3,3,3],is(3)=[1,1,1],ie(3)=[2,2,2]
  integer :: i,j,k,p,t
  integer :: ids(1)=[1],np(0:1)=[0,3]
  real(8) :: gamma(0:1,3),omega(0:1,3),f(0:1,3),wp(0:1)
  real(8) :: a(3,1),b(3,1),c(3,1),d(3,1),dp(1),dm(1)
  real(8) :: ar(3,1),br(3,1),cr(3,1),dr(3,1),dpr(1),dmr(1)
  real(8) :: e(0:3,0:3,0:3),z(0:3,0:3,0:3),dt=0.01d0,pi
  pi=acos(-1d0)
  ld%enabled=.true.;ref%enabled=.true.
  ld%max_pole_num_ld=3;ref%max_pole_num_ld=3
  call ld_allocate(ld,lo,hi)
  call legacy_allocate(ref,lo,hi)
  gamma=0.2d0;omega=0.3d0;omega(1,1)=0d0;f=0.4d0;wp=1d0
  a=0; b=0; c=0; d=0; dp=0; dm=0
  ar=0; br=0; cr=0; dr=0; dpr=0; dmr=0
  call ld_coefficients(ids,np,gamma,omega,f,wp,dt,pi,a,b,c,d,dp,dm)
  call legacy_coefficients(ids,np,gamma,omega,f,wp,dt,pi,ar,br,cr,dr,dpr,dmr)
  if(any(a/=ar).or.any(b/=br).or.any(c/=cr).or.any(d/=dr).or.any(dp/=dpr).or.any(dm/=dmr)) &
    error stop 'coefficient mismatch'
  ld%rjx_ld=-0.1d0;ref%rjx_ld=-0.1d0
  ld%rjx_old_ld=-0.1d0;ref%rjx_old_ld=-0.1d0
  ld%rjy_ld=-0.1d0;ref%rjy_ld=-0.1d0
  ld%rjy_old_ld=-0.1d0;ref%rjy_old_ld=-0.1d0
  ld%rjz_ld=-0.1d0;ref%rjz_ld=-0.1d0
  ld%rjz_old_ld=-0.1d0;ref%rjz_old_ld=-0.1d0
  ld%px_ld=0.2d0;ref%px_ld=0.2d0
  ld%py_ld=0.2d0;ref%py_ld=0.2d0
  ld%pz_ld=0.2d0;ref%pz_ld=0.2d0
  ld%c1_jx_ld=0.8d0;ref%c1_jx_ld=0.8d0
  ld%c2_jx_ld=0.12d0;ref%c2_jx_ld=0.12d0
  ld%c3_jx_ld=0.03d0;ref%c3_jx_ld=0.03d0
  ld%c1_jy_ld=0.8d0;ref%c1_jy_ld=0.8d0
  ld%c2_jy_ld=0.12d0;ref%c2_jy_ld=0.12d0
  ld%c3_jy_ld=0.03d0;ref%c3_jy_ld=0.03d0
  ld%c1_jz_ld=0.8d0;ref%c1_jz_ld=0.8d0
  ld%c2_jz_ld=0.12d0;ref%c2_jz_ld=0.12d0
  ld%c3_jz_ld=0.03d0;ref%c3_jz_ld=0.03d0
  ld%rjx_fdtd_ld=0.07d0;ref%rjx_fdtd_ld=0.07d0
  ld%rjy_fdtd_ld=0.07d0;ref%rjy_fdtd_ld=0.07d0
  ld%rjz_fdtd_ld=0.07d0;ref%rjz_fdtd_ld=0.07d0
  ld%ex_old_ld=0d0;ref%ex_old_ld=0d0
  ld%ey_old_ld=0d0;ref%ey_old_ld=0d0
  ld%ez_old_ld=0d0;ref%ez_old_ld=0d0
  other=ld
  z=0d0;e=0d0
  do t=1,30
    call ld_store_old(ld,lo,hi,is,ie,dt,e,z,e,z,e,z)
    call legacy_store_old(ref,lo,hi,is,ie,dt,e,z,e,z,e,z)
    do k=0,3
    do j=0,3
    do i=0,3
      e(i,j,k)=sin(0.3d0*t+i+2*j+3*k)
    end do
    end do
    end do
    call ld_update(ld,lo,hi,is,ie,dt,e,z,e,z,e,z)
    call legacy_update(ref,lo,hi,is,ie,dt,e,z,e,z,e,z)
    if(any(transfer(ld%rjx_ld,[0_int64])/=transfer(ref%rjx_ld,[0_int64]))) error stop 'rjx_ld'
    if(any(transfer(ld%rjx_old_ld,[0_int64])/=transfer(ref%rjx_old_ld,[0_int64]))) error stop 'rjx_old_ld'
    if(any(transfer(ld%rjy_ld,[0_int64])/=transfer(ref%rjy_ld,[0_int64]))) error stop 'rjy_ld'
    if(any(transfer(ld%rjy_old_ld,[0_int64])/=transfer(ref%rjy_old_ld,[0_int64]))) error stop 'rjy_old_ld'
    if(any(transfer(ld%rjz_ld,[0_int64])/=transfer(ref%rjz_ld,[0_int64]))) error stop 'rjz_ld'
    if(any(transfer(ld%rjz_old_ld,[0_int64])/=transfer(ref%rjz_old_ld,[0_int64]))) error stop 'rjz_old_ld'
    if(any(transfer(ld%px_ld,[0_int64])/=transfer(ref%px_ld,[0_int64]))) error stop 'px_ld'
    if(any(transfer(ld%py_ld,[0_int64])/=transfer(ref%py_ld,[0_int64]))) error stop 'py_ld'
    if(any(transfer(ld%pz_ld,[0_int64])/=transfer(ref%pz_ld,[0_int64]))) error stop 'pz_ld'
    if(any(transfer(ld%c1_jx_ld,[0_int64])/=transfer(ref%c1_jx_ld,[0_int64]))) error stop 'c1_jx_ld'
    if(any(transfer(ld%c2_jx_ld,[0_int64])/=transfer(ref%c2_jx_ld,[0_int64]))) error stop 'c2_jx_ld'
    if(any(transfer(ld%c3_jx_ld,[0_int64])/=transfer(ref%c3_jx_ld,[0_int64]))) error stop 'c3_jx_ld'
    if(any(transfer(ld%c1_jy_ld,[0_int64])/=transfer(ref%c1_jy_ld,[0_int64]))) error stop 'c1_jy_ld'
    if(any(transfer(ld%c2_jy_ld,[0_int64])/=transfer(ref%c2_jy_ld,[0_int64]))) error stop 'c2_jy_ld'
    if(any(transfer(ld%c3_jy_ld,[0_int64])/=transfer(ref%c3_jy_ld,[0_int64]))) error stop 'c3_jy_ld'
    if(any(transfer(ld%c1_jz_ld,[0_int64])/=transfer(ref%c1_jz_ld,[0_int64]))) error stop 'c1_jz_ld'
    if(any(transfer(ld%c2_jz_ld,[0_int64])/=transfer(ref%c2_jz_ld,[0_int64]))) error stop 'c2_jz_ld'
    if(any(transfer(ld%c3_jz_ld,[0_int64])/=transfer(ref%c3_jz_ld,[0_int64]))) error stop 'c3_jz_ld'
    if(any(transfer(ld%rjx_fdtd_ld,[0_int64])/=transfer(ref%rjx_fdtd_ld,[0_int64]))) error stop 'rjx_fdtd_ld'
    if(any(transfer(ld%rjy_fdtd_ld,[0_int64])/=transfer(ref%rjy_fdtd_ld,[0_int64]))) error stop 'rjy_fdtd_ld'
    if(any(transfer(ld%rjz_fdtd_ld,[0_int64])/=transfer(ref%rjz_fdtd_ld,[0_int64]))) error stop 'rjz_fdtd_ld'
    if(any(transfer(ld%ex_old_ld,[0_int64])/=transfer(ref%ex_old_ld,[0_int64]))) error stop 'ex_old_ld'
    if(any(transfer(ld%ey_old_ld,[0_int64])/=transfer(ref%ey_old_ld,[0_int64]))) error stop 'ey_old_ld'
    if(any(transfer(ld%ez_old_ld,[0_int64])/=transfer(ref%ez_old_ld,[0_int64]))) error stop 'ez_old_ld'
    ! Interleave a second independent context with a different field.
    call ld_store_old(other,lo,hi,is,ie,dt,z,z,z,z,z,z)
    call ld_update(other,lo,hi,is,ie,dt,z,z,z,z,z,z)
  end do
  call ld_finalize(ld)
  call ld_finalize(other)
  print *, 'PASS: coefficients and all P/J/Q/history arrays bitwise match, 30 steps; interleaved contexts'
end program
