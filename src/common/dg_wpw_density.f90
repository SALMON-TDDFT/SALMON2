module dg_wpw_density
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 public::build_wpw_density,wpw_metric_charge
contains
 subroutine build_wpw_density(psi_w,psi_p,occ,npoint,nstate,rho_direct,rho_ww,rho_wp,rho_pp,info)
  integer,intent(in)::npoint,nstate;complex(8),intent(in)::psi_w(npoint,nstate),psi_p(npoint,nstate)
  real(8),intent(in)::occ(nstate);real(8),intent(out)::rho_direct(npoint),rho_ww(npoint),rho_wp(npoint),rho_pp(npoint)
  integer,intent(out)::info;integer::i,j;complex(8)::z
  rho_direct=0;rho_ww=0;rho_wp=0;rho_pp=0;info=0
  if(any(occ<0d0).or.any(.not.ieee_is_finite(occ)))then;info=1;return;endif
  do j=1,nstate
   do i=1,npoint
    if(.not.ieee_is_finite(real(psi_w(i,j),8)).or..not.ieee_is_finite(aimag(psi_w(i,j))).or. &
       .not.ieee_is_finite(real(psi_p(i,j),8)).or..not.ieee_is_finite(aimag(psi_p(i,j))))then;info=2;return;endif
    rho_ww(i)=rho_ww(i)+occ(j)*abs(psi_w(i,j))**2
    rho_wp(i)=rho_wp(i)+2d0*occ(j)*real(conjg(psi_w(i,j))*psi_p(i,j),8)
    rho_pp(i)=rho_pp(i)+occ(j)*abs(psi_p(i,j))**2
    z=psi_w(i,j)+psi_p(i,j);rho_direct(i)=rho_direct(i)+occ(j)*abs(z)**2
   enddo
  enddo
 end subroutine
 subroutine wpw_metric_charge(c,s,occ,nbasis,nstate,charge_metric,charge_occ,orth_residual,info)
  integer,intent(in)::nbasis,nstate;complex(8),intent(in)::c(nbasis,nstate),s(nbasis,nbasis);real(8),intent(in)::occ(nstate)
  real(8),intent(out)::charge_metric,charge_occ,orth_residual;integer,intent(out)::info
  complex(8)::g(nstate,nstate);integer::i
  charge_metric=0;charge_occ=0;orth_residual=huge(1d0);info=0
  if(any(occ<0d0).or.maxval(abs(s-conjg(transpose(s))))>1d-10*max(1d0,maxval(abs(s))))then;info=1;return;endif
  g=matmul(conjg(transpose(c)),matmul(s,c));charge_occ=sum(occ)
  do i=1,nstate;charge_metric=charge_metric+occ(i)*real(g(i,i),8);g(i,i)=g(i,i)-1;enddo
  orth_residual=maxval(abs(g))
  if(.not.ieee_is_finite(charge_metric).or..not.ieee_is_finite(orth_residual))info=2
 end subroutine
end module
