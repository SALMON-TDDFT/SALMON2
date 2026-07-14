module dg_wannier_pw_scf
 use dg_wpw_fixed_operator,only:s_dg_wpw_fixed_operator,compose_dg_wpw_hamiltonian, &
  dg_wpw_fixed_dimensions,copy_dg_wpw_metric
 use dg_generalized_algebra,only:dg_generalized_eigh
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 type,public::s_dg_wpw_scf_result
  logical::converged=.false.
  integer::iterations=0
  real(8)::density_residual=huge(1d0),potential_residual=huge(1d0)
  real(8)::energy_residual=huge(1d0),q_occ_residual=huge(1d0)
  real(8)::eigen_residual=huge(1d0),orthogonality=huge(1d0)
  real(8)::charge_error=huge(1d0),fixed_point_residual=huge(1d0),total_energy=huge(1d0)
  real(8)::fixed_point_density_residual=huge(1d0),fixed_point_potential_residual=huge(1d0)
  real(8)::gap=0d0,eigenvalue_sum=0d0,e_h=0d0,n_vxc=0d0,e_xc=0d0,e_ion=0d0
 end type
 abstract interface
  subroutine dg_wpw_potential_map(c,occ,nbasis,nstate,ndensity,density,vout,eh,nvxc,exc,eion,info)
   integer,intent(in)::nbasis,nstate,ndensity
   complex(8),intent(in)::c(nbasis,nbasis);real(8),intent(in)::occ(nstate)
   real(8),intent(out)::density(ndensity),eh,nvxc,exc,eion
   complex(8),intent(out)::vout(nbasis,nbasis);integer,intent(out)::info
  end subroutine
 end interface
 public::run_dg_wpw_fixed_scf
contains
 subroutine run_dg_wpw_fixed_scf(op,occ,nstate,ndensity,v_initial,density_initial,mix,max_iter, &
   tol_density,tol_potential,tol_energy,tol_q,tol_eigen,gap_min,potential_map,c,eval,result,info)
  type(s_dg_wpw_fixed_operator),intent(in)::op
  integer,intent(in)::nstate,ndensity,max_iter
  real(8),intent(in)::occ(nstate),density_initial(ndensity),mix
  complex(8),intent(in)::v_initial(:,:)
  real(8),intent(in)::tol_density,tol_potential,tol_energy,tol_q,tol_eigen,gap_min
  procedure(dg_wpw_potential_map)::potential_map
  complex(8),intent(out)::c(:,:);real(8),intent(out)::eval(:)
  type(s_dg_wpw_scf_result),intent(out)::result;integer,intent(out)::info
  integer::n,nw,np,iter,i,discarded,map_info
  real(8)::density(ndensity),density_prev(ndensity),density_check(ndensity)
  real(8)::eh,nvxc,exc,eion,energy,energy_prev,resid,orth,cond,ne
  real(8)::eh2,nvxc2,exc2,eion2
  complex(8),allocatable::vin(:,:),vout(:,:),vcheck(:,:),h(:,:),q(:,:),qprev(:,:),ccheck(:,:),sp(:,:),s(:,:)
  real(8),allocatable::echeck(:)
  result=s_dg_wpw_scf_result();info=0;call dg_wpw_fixed_dimensions(op,nw,np,n);c=0;eval=0
  if(n<1.or.nstate<1.or.nstate>n.or.ndensity<1.or.max_iter<1.or.mix<=0d0.or.mix>1d0.or. &
     nstate>=n.or.any(occ<=0d0).or.any(occ>2d0).or.any(abs(occ-dnint(occ))>1d-12).or. &
     min(tol_density,tol_potential,tol_energy,tol_q,tol_eigen,gap_min)<=0d0.or. &
     .not.all(ieee_is_finite(occ)).or..not.ieee_is_finite(mix).or. &
     .not.all(ieee_is_finite([tol_density,tol_potential,tol_energy,tol_q,tol_eigen,gap_min])).or. &
     any(shape(v_initial)/=[n,n]).or.any(shape(c)/=[n,n]).or.size(eval)/=n.or. &
     .not.finite_complex(v_initial).or..not.all(ieee_is_finite(density_initial)))then;info=1;return;endif
  allocate(vin(n,n),vout(n,n),vcheck(n,n),h(n,n),q(n,n),qprev(n,n),ccheck(n,n),sp(n,n),s(n,n),echeck(n))
  call copy_dg_wpw_metric(op,s,info);if(info/=0)then;info=2;return;endif
  vin=v_initial;density_prev=density_initial;qprev=0;energy_prev=huge(1d0);ne=sum(occ)
  do iter=1,max_iter
   call compose_dg_wpw_hamiltonian(op,vin,h,info);if(info/=0)then;info=10+info;return;endif
   call dg_generalized_eigh(h,s,n,1d-12,eval,c,resid,orth,cond,discarded,info)
   if(info/=0)return
   result%gap=eval(nstate+1)-eval(nstate);if(.not.ieee_is_finite(result%gap).or.result%gap<gap_min)then;info=3;return;endif
   call potential_map(c,occ,n,nstate,ndensity,density,vout,eh,nvxc,exc,eion,map_info)
   if(map_info/=0)then;info=20+map_info;return;endif
   if(.not.map_output_finite(density,vout,eh,nvxc,exc,eion).or..not.hermitian(vout))then;info=24;return;endif
   call density_matrix(c,occ,n,nstate,q)
   energy=sum(occ*eval(1:nstate))-eh-nvxc+exc+eion
   result%density_residual=relative_real(density,density_prev)
   result%potential_residual=relative_complex(vout,vin)
   if(iter>1)result%energy_residual=abs(energy-energy_prev)
   result%q_occ_residual=relative_complex(q,qprev)
   result%eigen_residual=resid;result%orthogonality=orth;sp=matmul(s,q)
   result%charge_error=real(trace_matrix(sp),8)-ne
   result%iterations=iter;result%total_energy=energy;result%eigenvalue_sum=sum(occ*eval(1:nstate))
   result%e_h=eh;result%n_vxc=nvxc;result%e_xc=exc;result%e_ion=eion
   if(iter>1.and.result%density_residual<tol_density.and.result%potential_residual<tol_potential.and. &
      result%energy_residual<tol_energy.and.result%q_occ_residual<tol_q.and.resid<tol_eigen.and. &
      orth<tol_eigen.and.abs(result%charge_error)<tol_eigen)then
    result%converged=.true.;exit
   endif
   density_prev=density;qprev=q;energy_prev=energy;vin=(1d0-mix)*vin+mix*vout
  enddo
  if(.not.result%converged)then;info=30;return;endif
  call compose_dg_wpw_hamiltonian(op,vout,h,info);if(info/=0)return
  call dg_generalized_eigh(h,s,n,1d-12,echeck,ccheck,resid,orth,cond,discarded,info);if(info/=0)return
  if(echeck(nstate+1)-echeck(nstate)<gap_min)then;info=3;return;endif
  call potential_map(ccheck,occ,n,nstate,ndensity,density_check,vcheck,eh2,nvxc2,exc2,eion2,map_info)
  if(map_info/=0)then;info=20+map_info;return;endif
  if(.not.map_output_finite(density_check,vcheck,eh2,nvxc2,exc2,eion2).or..not.hermitian(vcheck))then;info=24;return;endif
  result%fixed_point_density_residual=relative_real(density_check,density)
  result%fixed_point_potential_residual=relative_complex(vcheck,vout)
  result%fixed_point_residual=max(result%fixed_point_density_residual/tol_density, &
    result%fixed_point_potential_residual/tol_potential)
  if(result%fixed_point_density_residual>=tol_density.or. &
     result%fixed_point_potential_residual>=tol_potential)then;result%converged=.false.;info=31;return;endif
  c=ccheck;eval=echeck;call density_matrix(c,occ,n,nstate,q);sp=matmul(s,q)
  result%gap=eval(nstate+1)-eval(nstate);result%eigen_residual=resid;result%orthogonality=orth
  result%charge_error=real(trace_matrix(sp),8)-ne;result%eigenvalue_sum=sum(occ*eval(1:nstate))
  result%e_h=eh2;result%n_vxc=nvxc2;result%e_xc=exc2;result%e_ion=eion2
  result%total_energy=result%eigenvalue_sum-eh2-nvxc2+exc2+eion2
 contains
  subroutine density_matrix(a,f,m,ns,p)
   integer,intent(in)::m,ns;complex(8),intent(in)::a(m,m);real(8),intent(in)::f(ns);complex(8),intent(out)::p(m,m);integer::j
   p=0;do j=1,ns;p=p+f(j)*matmul(reshape(a(:,j),[m,1]),reshape(conjg(a(:,j)),[1,m]));enddo
  end subroutine
  real(8) function relative_real(a,b);real(8),intent(in)::a(:),b(:);relative_real=sqrt(sum((a-b)**2))/max(1d-30,sqrt(sum(a**2)));end
  real(8) function relative_complex(a,b);complex(8),intent(in)::a(:,:),b(:,:);relative_complex=sqrt(sum(abs(a-b)**2))/max(1d-30,sqrt(sum(abs(a)**2)));end
  complex(8) function trace_matrix(a);complex(8),intent(in)::a(:,:);integer::j;trace_matrix=0;do j=1,min(size(a,1),size(a,2));trace_matrix=trace_matrix+a(j,j);enddo;end
  logical function finite_complex(a);complex(8),intent(in)::a(:,:);finite_complex=all(ieee_is_finite(real(a,8))).and.all(ieee_is_finite(aimag(a)));end
  logical function hermitian(a);complex(8),intent(in)::a(:,:);hermitian=maxval(abs(a-conjg(transpose(a))))<=1d-11*max(1d0,maxval(abs(a)));end
  logical function map_output_finite(d,v,a,b,e,f);real(8),intent(in)::d(:),a,b,e,f;complex(8),intent(in)::v(:,:);map_output_finite=all(ieee_is_finite(d)).and.finite_complex(v).and.all(ieee_is_finite([a,b,e,f]));end
 end subroutine
end module
