program wpw_density_driver
 use dg_wpw_density
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_positive_inf
 implicit none
 integer,parameter::np=4,ns=2,nb=2
 complex(8)::w(np,ns),p(np,ns),c(nb,ns),s(nb,nb),b(np,nb),psi(np,ns)
 real(8)::occ(ns),rd(np),rww(np),rwp(np),rpp(np),qmetric,qocc,orth,qgrid,weight(np)
 integer::info,i
 w=reshape([cmplx(.6d0,.1d0,8),cmplx(-.2d0,.3d0,8),cmplx(.4d0,-.1d0,8),cmplx(.1d0,.2d0,8), &
            cmplx(.2d0,-.3d0,8),cmplx(.5d0,.1d0,8),cmplx(-.1d0,.4d0,8),cmplx(.3d0,-.2d0,8)],[np,ns])
 p=reshape([cmplx(.15d0,-.05d0,8),cmplx(.08d0,.12d0,8),cmplx(-.11d0,.07d0,8),cmplx(.09d0,.03d0,8), &
            cmplx(-.07d0,.13d0,8),cmplx(.12d0,-.04d0,8),cmplx(.05d0,.09d0,8),cmplx(-.1d0,.08d0,8)],[np,ns])
 occ=[2d0,1d0]
 call build_wpw_density(w,p,occ,np,ns,rd,rww,rwp,rpp,info)
 call req(info==0.and.maxval(abs(rd-(rww+rwp+rpp)))<1d-13,'decomposition')
 call req(maxval(abs(rd-(rww+rpp)))>1d-3,'missing WP control')
 call req(maxval(abs(rd-(rww+rwp)))>1d-3,'missing PP control')
 b=reshape([cmplx(1d0,.1d0,8),cmplx(.4d0,-.2d0,8),cmplx(-.3d0,.5d0,8),cmplx(.7d0,.2d0,8), &
            cmplx(.2d0,.6d0,8),cmplx(-.5d0,.1d0,8),cmplx(.8d0,-.3d0,8),cmplx(.1d0,.4d0,8)],[np,nb])
 weight=[.2d0,.3d0,.25d0,.25d0];s=0
 do i=1,np;s=s+weight(i)*matmul(reshape(conjg(b(i,:)),[nb,1]),reshape(b(i,:),[1,nb]));enddo
 c=reshape([cmplx(.852781d0,0d0,8),cmplx(-.080701d0,.0403505d0,8), &
            cmplx(-.080701d0,-.0403505d0,8),cmplx(.970716d0,0d0,8)],[nb,ns])
 ! Orthonormalize this deterministic seed inside the fixture through the exact 2x2 Gram metric.
 call normalize_metric_columns(c,s,nb,ns)
 psi=matmul(b,c);w=.65d0*psi;p=.35d0*psi
 call build_wpw_density(w,p,occ,np,ns,rd,rww,rwp,rpp,info);qgrid=sum(weight*rd)
 call wpw_metric_charge(c,s,occ,nb,ns,qmetric,qocc,orth,info)
 call req(info==0.and.abs(qgrid-qmetric)<1d-11.and.abs(qmetric-qocc)<1d-11.and.abs(qocc-3d0)<1d-13.and.orth<1d-11,'metric charge')
 occ(1)=-1d0;call build_wpw_density(w,p,occ,np,ns,rd,rww,rwp,rpp,info);call req(info/=0,'negative occupation')
 occ=[2d0,1d0];w(1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
 call build_wpw_density(w,p,occ,np,ns,rd,rww,rwp,rpp,info);call req(info==2,'NaN amplitude')
 w(1,1)=0;occ(1)=ieee_value(0d0,ieee_positive_inf)
 call build_wpw_density(w,p,occ,np,ns,rd,rww,rwp,rpp,info);call req(info==1,'infinite occupation')
 write(*,'(a)')'PASS mixed density'
contains
 subroutine req(ok,msg);logical,intent(in)::ok;character(*),intent(in)::msg;if(.not.ok)then;write(*,*)msg;error stop 1;endif;end
 subroutine normalize_metric_columns(a,m,n,k)
  integer,intent(in)::n,k;complex(8),intent(in)::m(n,n);complex(8),intent(inout)::a(n,k);integer::i,j
  do j=1,k;do i=1,j-1;a(:,j)=a(:,j)-a(:,i)*dot_product(a(:,i),matmul(m,a(:,j)));enddo;a(:,j)=a(:,j)/sqrt(real(dot_product(a(:,j),matmul(m,a(:,j))),8));enddo
 end subroutine
end program
