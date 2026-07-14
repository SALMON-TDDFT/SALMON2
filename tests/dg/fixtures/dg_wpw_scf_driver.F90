program dg_wpw_scf_driver
 use dg_wpw_fixed_operator
 use dg_wannier_pw_scf
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
 implicit none
 integer,parameter::nw=2,np=1,n=3,ns=2,nd=3
 type(s_dg_wpw_fixed_operator)::op
 type(s_dg_wpw_scf_result)::result
 complex(8)::sww(nw,nw),swp(nw,np),spp(np,np),hww(nw,nw),hwp(nw,np),hpp(np,np)
 complex(8)::fww(nw,nw),fwp(nw,np),v0(n,n),c(n,n),vref(n,n)
 real(8)::occ(ns),d0(nd),eval(n),dref(nd),ehref,nvxcref,excref,eionref,eref
 integer::info
 logical::emit_nan=.false.
 integer::map_calls=0,shift_call=0,converged_iterations
 sww=0;sww(1,1)=1;sww(2,2)=1;swp=0;spp=1
 hww=0;hww(1,1)=-.8d0;hww(2,2)=.15d0;hww(1,2)=.04d0;hww(2,1)=.04d0
 hwp=reshape([cmplx(.08d0,0d0,8),cmplx(-.03d0,.01d0,8)],[nw,np]);hpp=reshape([cmplx(.9d0,0d0,8)],[np,np])
 fww=0;fwp=0
 call initialize_dg_wpw_fixed_operator(op,sww,swp,spp,hww,hwp,hpp,fww,fwp,nw,np,info)
 call req(info==0,'operator');occ=[2d0,1d0];v0=0;d0=[1.5d0,1d0,.5d0]
 call run_dg_wpw_fixed_scf(op,occ,ns,nd,v0,d0,.35d0,300,1d-10,1d-12,1d-11,1d-10,1d-10, &
   1d-4,model_map,c,eval,result,info)
 call req(info==0.and.result%converged,'converged')
 call req(result%density_residual<1d-10.and.result%potential_residual<1d-10,'density potential residual')
 call req(result%energy_residual<1d-11.and.result%q_occ_residual<1d-10,'energy Q residual')
 call req(result%eigen_residual<1d-10.and.result%orthogonality<1d-10,'generalized eigensystem')
 call req(abs(result%charge_error)<1d-10.and.result%fixed_point_residual<1d0,'charge fixed point')
 converged_iterations=result%iterations
 call model_map(c,occ,n,ns,nd,dref,vref,ehref,nvxcref,excref,eionref,info)
 eref=sum(occ*eval(1:ns))-ehref-nvxcref+excref+eionref
 call req(abs(result%total_energy-eref)<1d-13.and.abs(result%eigenvalue_sum-sum(occ*eval(1:ns)))<1d-13,'independent energy')
 call req(result%gap>=1d-4,'gap')
 write(*,'(a,7es24.15)')'SCF_REF ',eval(1:3),dref(1:3),result%total_energy
 map_calls=0;shift_call=converged_iterations+1
 call run_dg_wpw_fixed_scf(op,occ,ns,nd,v0,d0,.35d0,300,1d-6,1d-12,1d-11,1d-10,1d-10, &
   1d-4,model_map,c,eval,result,info)
 call req(info==31.and.result%fixed_point_density_residual<1d-6.and. &
   result%fixed_point_potential_residual>=1d-12,'separate fixed-point tolerances')
 shift_call=0
 occ=[3d0,1d0]
 call run_dg_wpw_fixed_scf(op,occ,ns,nd,v0,d0,.35d0,20,1d-8,1d-8,1d-8,1d-8,1d-8, &
   1d-4,model_map,c,eval,result,info)
 call req(info/=0,'malformed integer occupation rejection')
 occ=[2d0,1d0]
 call run_dg_wpw_fixed_scf(op,occ,ns,nd,v0,d0,.35d0,20,1d-8,1d-8,1d-8,1d-8,1d-8, &
   10d0,model_map,c,eval,result,info)
 call req(info==3,'gap rejection')
 emit_nan=.true.
 call run_dg_wpw_fixed_scf(op,occ,ns,nd,v0,d0,.35d0,20,1d-8,1d-8,1d-8,1d-8,1d-8, &
   1d-4,model_map,c,eval,result,info)
 call req(info==24,'callback NaN rejection')
 write(*,'(a)')'PASS DG WPW fixed-point SCF'
contains
 subroutine model_map(c,occ,nbasis,nstate,ndensity,density,vout,eh,nvxc,exc,eion,info)
  integer,intent(in)::nbasis,nstate,ndensity;complex(8),intent(in)::c(nbasis,nbasis);real(8),intent(in)::occ(nstate)
  real(8),intent(out)::density(ndensity),eh,nvxc,exc,eion;complex(8),intent(out)::vout(nbasis,nbasis);integer,intent(out)::info
  complex(8)::p(nbasis,nbasis);integer::i
  p=0;do i=1,nstate;p=p+occ(i)*matmul(reshape(c(:,i),[nbasis,1]),reshape(conjg(c(:,i)),[1,nbasis]));enddo
  map_calls=map_calls+1
  density=real([(p(i,i),i=1,ndensity)],8);vout=0;do i=1,nbasis;vout(i,i)=.15d0*density(i);enddo
  eh=.5d0*sum(density*real([(vout(i,i),i=1,nbasis)],8));nvxc=0;exc=0;eion=.2d0;info=0
  if(emit_nan)density(1)=ieee_value(0d0,ieee_quiet_nan)
  if(map_calls==shift_call)vout(1,1)=vout(1,1)+1d-8
 end subroutine
 subroutine req(ok,msg);logical,intent(in)::ok;character(*),intent(in)::msg;if(.not.ok)then;write(*,*)trim(msg);error stop 1;endif;end
end program
