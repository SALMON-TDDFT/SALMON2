program wpw_generalized_driver
 use dg_generalized_algebra
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
 implicit none
 integer,parameter::n=3
 complex(8)::s(n,n),h0(n,n),hf(n,n),h(n,n),c(n,n),x(n,n),u(n,n),c0(n),c1(n),ce(n),cnf(n),tmp(n,n),d(n,1),back(n,1)
 real(8)::eval(n),resid,orth,cond,sunit,emin,emax,nfull,neuclid
 integer::info,discarded,i,j
 complex(8)::metric_diag(4,4),metric_overlap(4,3)
 real(8)::metric_minimum,metric_maximum,minimum_ratio,retained_minimum_ratio,occupied_discarded_fraction,&
   extra_discarded_fraction,occupied_low_fraction,extra_low_fraction
 integer::effective_rank
 s=reshape([cmplx(1.4d0,0d0,8),cmplx(.12d0,-.04d0,8),cmplx(.03d0,.02d0,8), &
            cmplx(.12d0,.04d0,8),cmplx(1.1d0,0d0,8),cmplx(.08d0,-.01d0,8), &
            cmplx(.03d0,-.02d0,8),cmplx(.08d0,.01d0,8),cmplx(.9d0,0d0,8)],[n,n])
 h0=reshape([cmplx(.8d0,0d0,8),cmplx(.11d0,.03d0,8),cmplx(.02d0,-.04d0,8), &
             cmplx(.11d0,-.03d0,8),cmplx(1.3d0,0d0,8),cmplx(.07d0,.02d0,8), &
             cmplx(.02d0,.04d0,8),cmplx(.07d0,-.02d0,8),cmplx(1.8d0,0d0,8)],[n,n])
 hf=0;hf(1,2)=cmplx(-.19d0,.06d0,8);hf(2,1)=conjg(hf(1,2));hf(2,2)=.13d0;h=h0+hf
 call dg_metric_factor(s,n,1d-10,x,emin,emax,cond,discarded,info);call req(info==0.and.discarded==0,'metric')
 call dg_metric_factor(1d-12*s,n,1d-10,tmp,emin,emax,orth,discarded,info)
 call req(info==0.and.discarded==0.and.abs(orth-cond)<1d-10,'scale invariant metric cutoff')
 tmp=matmul(conjg(transpose(x)),matmul(s,x));call req(maxval(abs(tmp-identity(n)))<1d-12,'metric identity')
 call dg_generalized_eigh(h,s,n,1d-10,eval,c,resid,orth,cond,discarded,info)
 write(*,'(a,i0,3(1x,es12.4))')'generalized diagnostics ',info,resid,orth,cond
 call req(info==0.and.resid<1d-11.and.orth<1d-11,'generalized eigen')
 c0=[cmplx(.7d0,.1d0,8),cmplx(-.2d0,.3d0,8),cmplx(.1d0,-.4d0,8)]
 c0=c0/sqrt(real(dot_product(c0,matmul(s,c0)),8))
 d(:,1)=[cmplx(.2d0,.1d0,8),cmplx(-.3d0,.4d0,8),cmplx(.5d0,-.2d0,8)]
 call dg_metric_from_orth(x,d,back,n,1);call dg_metric_to_orth(x,s,back,d,n,1);call req(maxval(abs(d(:,1)-[cmplx(.2d0,.1d0,8),cmplx(-.3d0,.4d0,8),cmplx(.5d0,-.2d0,8)]))<1d-12,'transform')
 call dg_s_exp_action(h,s,n,1d-10,.17d0,c0,c1,sunit,cond,discarded,info);call req(info==0.and.sunit<1d-11,'S unitary')
 nfull=real(dot_product(c1,matmul(s,c1)),8);call req(abs(nfull-1d0)<1d-11,'S norm')
 call dg_s_exp_action(h0,s,n,1d-10,.17d0,c0,cnf,sunit,cond,discarded,info);call req(maxval(abs(cnf-c1))>1d-4,'missing face control')
 call euclidean_action(h,n,.17d0,c0,ce);neuclid=abs(real(dot_product(ce,matmul(s,ce)),8)-1d0)
 call req(neuclid>1d-5,'Euclidean negative control')
 metric_diag=0;metric_diag(1,1)=1d-12;metric_diag(2,2)=5d-10
 metric_diag(3,3)=1d-8;metric_diag(4,4)=1d0;metric_overlap=0
 metric_overlap(1,1)=3d0;metric_overlap(3,1)=4d0;metric_overlap(2,2)=1d0
 metric_overlap(3,3)=2d0;metric_overlap(4,3)=1d0
 call dg_metric_mode_residual_split(metric_diag,metric_overlap,4,3,2,1d-9,metric_minimum,metric_maximum,minimum_ratio,&
   retained_minimum_ratio,effective_rank,occupied_discarded_fraction,extra_discarded_fraction,&
   occupied_low_fraction,extra_low_fraction,info)
 call req(info==0.and.effective_rank==2.and.abs(metric_minimum-1d-12)<1d-20.and.metric_maximum==1d0.and.&
   abs(minimum_ratio-1d-12)<1d-20.and.&
   abs(retained_minimum_ratio-1d-8)<1d-16.and.&
   abs(occupied_discarded_fraction-10d0/26d0)<1d-12.and.&
   abs(extra_discarded_fraction)<1d-12.and.abs(occupied_low_fraction-16d0/26d0)<1d-12.and.&
   abs(extra_low_fraction-4d0/5d0)<1d-12,'metric mode residual split')
 metric_overlap=0
 call dg_metric_mode_residual_split(metric_diag,metric_overlap,4,3,2,1d-9,metric_minimum,metric_maximum,minimum_ratio,&
   retained_minimum_ratio,effective_rank,occupied_discarded_fraction,extra_discarded_fraction,&
   occupied_low_fraction,extra_low_fraction,info)
 call req(info==0.and.max(occupied_discarded_fraction,extra_discarded_fraction,&
   occupied_low_fraction,extra_low_fraction)==0d0,'zero metric overlap')
 metric_diag=0;do i=1,4;metric_diag(i,i)=dble(i);enddo
 call dg_metric_mode_residual_split(metric_diag,metric_overlap,4,3,2,1d-9,metric_minimum,metric_maximum,minimum_ratio,&
   retained_minimum_ratio,effective_rank,occupied_discarded_fraction,extra_discarded_fraction,&
   occupied_low_fraction,extra_low_fraction,info)
 call req(info==0.and.effective_rank==4.and.occupied_discarded_fraction==0d0,'no discarded metric mode')
 metric_overlap(1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
 call dg_metric_mode_residual_split(metric_diag,metric_overlap,4,3,2,1d-9,metric_minimum,metric_maximum,minimum_ratio,&
   retained_minimum_ratio,effective_rank,occupied_discarded_fraction,extra_discarded_fraction,&
   occupied_low_fraction,extra_low_fraction,info)
 call req(info/=0,'nonfinite metric overlap')
 metric_overlap=0;metric_overlap(4,1)=cmplx(sqrt(huge(1d0)),sqrt(huge(1d0)),8)
 call dg_metric_mode_residual_split(metric_diag,metric_overlap,4,3,2,1d-9,metric_minimum,metric_maximum,minimum_ratio,&
   retained_minimum_ratio,effective_rank,occupied_discarded_fraction,extra_discarded_fraction,&
   occupied_low_fraction,extra_low_fraction,info)
 call req(info/=0,'overflowed metric overlap weight')
 write(*,'(a,*(1x,es24.16))')'EVAL',eval
 write(*,'(a,*(1x,es24.16))')'X',((real(x(i,j),8),aimag(x(i,j)),i=1,n),j=1,n)
 write(*,'(a,*(1x,es24.16))')'C',((real(c(i,j),8),aimag(c(i,j)),i=1,n),j=1,n)
 write(*,'(a,*(1x,es24.16))')'EXP',(real(c1(i),8),aimag(c1(i)),i=1,n)
 write(*,'(a,*(1x,es24.16))')'EXPNF',(real(cnf(i),8),aimag(cnf(i)),i=1,n)
 tmp=s;tmp(3,:)=tmp(2,:);tmp(:,3)=tmp(:,2)
 call dg_metric_factor(tmp,n,1d-10,x,emin,emax,cond,discarded,info)
 call req(info/=0.and.discarded>0,'rank loss fail closed')
 tmp=h;tmp(1,2)=tmp(1,2)+cmplx(.2d0,.1d0,8);call dg_generalized_eigh(tmp,s,n,1d-10,eval,c,resid,orth,cond,discarded,info);call req(info==20,'nonhermitian fail closed')
 write(*,'(a)')'PASS generalized algebra'
contains
 subroutine req(ok,msg);logical,intent(in)::ok;character(*),intent(in)::msg
  if(.not.ok)then;write(*,'(a)')trim(msg);error stop 1;endif
 end subroutine
 subroutine euclidean_action(a,n,dt,v,w)
  integer,intent(in)::n;real(8),intent(in)::dt;complex(8),intent(in)::a(n,n),v(n);complex(8),intent(out)::w(n)
  complex(8)::q(n,n),work(64),coef(n);real(8)::e(n),rw(16);integer::i,inf
  q=a;call zheev('V','U',n,q,n,e,work,64,rw,inf);coef=matmul(conjg(transpose(q)),v)
  do i=1,n;coef(i)=coef(i)*exp(cmplx(0d0,-dt*e(i),8));enddo;w=matmul(q,coef)
 end subroutine
 function identity(m) result(a);integer,intent(in)::m;complex(8)::a(m,m);integer::k;a=0;do k=1,m;a(k,k)=1;enddo;end function
end program
