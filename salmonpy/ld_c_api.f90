! Copyright 2026 SALMON developers. Apache-2.0.
! Private C ABI for Python: one isotropic material, explicit per-pole state.
subroutine salmonpy_ld(n,nt,dt,wp,f,w,g,e,p,j,po,jo) bind(C)
  use iso_c_binding, only: c_int,c_double
  use lorentz_drude, only: s_ld_state,ld_allocate,ld_coefficients,ld_store_old,ld_update
  implicit none
  integer(c_int), value :: n,nt
  real(c_double), value :: dt,wp
  real(c_double), intent(in) :: f(n),w(n),g(n),e(3,nt)
  real(c_double), intent(inout) :: p(3,n),j(3,n)
  real(c_double), intent(out) :: po(3,nt),jo(3,nt)
  type(s_ld_state) :: ld
  real(8) :: ff(0:0,n),ww(0:0,n),gg(0:0,n),wp0(0:0)
  real(8) :: c0(n,1),c1(n,1),c2(n,1),c3(n,1),c4(1),c5(1)
  real(8) :: x(1,1,1),y(1,1,1),z(1,1,1),zero(1,1,1)
  integer :: k,q
  ld%enabled=.true.; ld%num_ld=1; ld%max_pole_num_ld=n
  call ld_allocate(ld,[1,1,1],[1,1,1])
  ff(0,:)=f; ww(0,:)=w; gg(0,:)=g; wp0(0)=wp
  c4=0; c5=0; zero=0
  call ld_coefficients([0],[n],gg,ww,ff,wp0,dt,acos(-1.d0),c0,c1,c2,c3,c4,c5)
  do q=1,n
    ld%c1_jx_ld(1,1,1,q)=c1(q,1)
    ld%c2_jx_ld(1,1,1,q)=c2(q,1)
    ld%c3_jx_ld(1,1,1,q)=c3(q,1)
    ld%px_ld(1,1,1,q)=p(1,q)
    ld%rjx_ld(1,1,1,q)=j(1,q)
    ld%c1_jy_ld(1,1,1,q)=c1(q,1)
    ld%c2_jy_ld(1,1,1,q)=c2(q,1)
    ld%c3_jy_ld(1,1,1,q)=c3(q,1)
    ld%py_ld(1,1,1,q)=p(2,q)
    ld%rjy_ld(1,1,1,q)=j(2,q)
    ld%c1_jz_ld(1,1,1,q)=c1(q,1)
    ld%c2_jz_ld(1,1,1,q)=c2(q,1)
    ld%c3_jz_ld(1,1,1,q)=c3(q,1)
    ld%pz_ld(1,1,1,q)=p(3,q)
    ld%rjz_ld(1,1,1,q)=j(3,q)
  end do
  po(:,1)=sum(p,dim=2); jo(:,1)=sum(j,dim=2)
  do k=2,nt
    x=e(1,k-1); y=e(2,k-1); z=e(3,k-1)
    call ld_store_old(ld,[1,1,1],[1,1,1],[1,1,1],[1,1,1],dt,x,zero,y,zero,z,zero)
    x=e(1,k); y=e(2,k); z=e(3,k)
    call ld_update(ld,[1,1,1],[1,1,1],[1,1,1],[1,1,1],dt,x,zero,y,zero,z,zero)
    po(1,k)=sum(ld%px_ld(1,1,1,:)); jo(1,k)=sum(ld%rjx_ld(1,1,1,:))
    po(2,k)=sum(ld%py_ld(1,1,1,:)); jo(2,k)=sum(ld%rjy_ld(1,1,1,:))
    po(3,k)=sum(ld%pz_ld(1,1,1,:)); jo(3,k)=sum(ld%rjz_ld(1,1,1,:))
  end do
  do q=1,n
    p(1,q)=ld%px_ld(1,1,1,q); j(1,q)=ld%rjx_ld(1,1,1,q)
    p(2,q)=ld%py_ld(1,1,1,q); j(2,q)=ld%rjy_ld(1,1,1,q)
    p(3,q)=ld%pz_ld(1,1,1,q); j(3,q)=ld%rjz_ld(1,1,1,q)
  end do
end subroutine salmonpy_ld
