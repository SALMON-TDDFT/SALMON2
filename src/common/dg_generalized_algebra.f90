module dg_generalized_algebra
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 public::dg_metric_factor,dg_generalized_eigh,dg_reduced_generalized_eigh,dg_s_exp_action, &
  dg_metric_to_orth,dg_metric_from_orth,dg_metric_mode_residual_split
contains
 subroutine dg_metric_mode_residual_split(s,b,n,nrhs,nocc,cutoff,metric_minimum,metric_maximum,minimum_ratio,&
   retained_minimum_ratio,effective_rank,occupied_discarded_fraction,extra_discarded_fraction,&
   occupied_low_fraction,extra_low_fraction,info)
  integer,intent(in)::n,nrhs,nocc
  real(8),intent(in)::cutoff
  complex(8),intent(in)::s(n,n),b(n,nrhs)
  real(8),intent(out)::metric_minimum,metric_maximum,minimum_ratio,retained_minimum_ratio,occupied_discarded_fraction,&
    extra_discarded_fraction,occupied_low_fraction,extra_low_fraction
  integer,intent(out)::effective_rank,info
  complex(8)::v(n,n),projected(n,nrhs)
  real(8)::e(n),emax,threshold,occupied_total,extra_total
  integer::first_retained

  info=1;effective_rank=0;metric_minimum=huge(1d0);metric_maximum=huge(1d0)
  minimum_ratio=huge(1d0);retained_minimum_ratio=huge(1d0)
  occupied_discarded_fraction=huge(1d0);extra_discarded_fraction=huge(1d0)
  occupied_low_fraction=huge(1d0);extra_low_fraction=huge(1d0)
  if(n<2.or.nrhs<2.or.nocc<1.or.nocc>=nrhs.or..not.ieee_is_finite(cutoff).or.cutoff<=0d0.or.&
    .not.all(ieee_is_finite(real(s,8))).or..not.all(ieee_is_finite(aimag(s))).or.&
    .not.all(ieee_is_finite(real(b,8))).or..not.all(ieee_is_finite(aimag(b))))return
  if(maxval(abs(s-conjg(transpose(s))))>1d-10*max(1d0,maxval(abs(s))))then;info=10;return;endif
  v=.5d0*(s+conjg(transpose(s)));call heev(v,n,e,info);if(info/=0)return
  emax=e(n);if(.not.all(ieee_is_finite(e)).or.emax<=0d0.or.e(1)<-1d-12*emax)then;info=11;return;endif
  threshold=cutoff*emax;if(.not.ieee_is_finite(threshold))then;info=13;return;endif
  effective_rank=count(e>threshold)
  if(effective_rank<1)then;info=12;return;endif
  first_retained=n-effective_rank+1
  metric_minimum=e(1);metric_maximum=emax
  minimum_ratio=max(0d0,e(1)/emax);retained_minimum_ratio=e(first_retained)/emax
  projected=matmul(conjg(transpose(v)),b)
  if(.not.all(ieee_is_finite(real(projected,8))).or.&
    .not.all(ieee_is_finite(aimag(projected))))then;info=14;return;endif
  occupied_total=sum(abs(projected(:,1:nocc))**2)
  extra_total=sum(abs(projected(:,nocc+1:nrhs))**2)
  if(.not.all(ieee_is_finite([occupied_total,extra_total])))then;info=15;return;endif
  occupied_discarded_fraction=0d0;extra_discarded_fraction=0d0
  occupied_low_fraction=0d0;extra_low_fraction=0d0
  if(occupied_total>0d0)then
   occupied_discarded_fraction=sum(abs(projected(1:first_retained-1,1:nocc))**2)/occupied_total
   occupied_low_fraction=sum(abs(projected(first_retained:n,1:nocc))**2,&
     mask=spread(e(first_retained:n)<=10d0*threshold,2,nocc))/occupied_total
  endif
  if(extra_total>0d0)then
   extra_discarded_fraction=sum(abs(projected(1:first_retained-1,nocc+1:nrhs))**2)/extra_total
   extra_low_fraction=sum(abs(projected(first_retained:n,nocc+1:nrhs))**2,&
     mask=spread(e(first_retained:n)<=10d0*threshold,2,nrhs-nocc))/extra_total
  endif
  if(.not.all(ieee_is_finite([metric_minimum,metric_maximum,minimum_ratio,retained_minimum_ratio,&
    occupied_discarded_fraction,extra_discarded_fraction,occupied_low_fraction,extra_low_fraction])))return
  info=0
 end subroutine dg_metric_mode_residual_split
 subroutine heev(a,n,e,info)
  integer,intent(in)::n;complex(8),intent(inout)::a(n,n);real(8),intent(out)::e(n);integer,intent(out)::info
  complex(8),allocatable::w(:);real(8),allocatable::rw(:);complex(8)::q(1);integer::lw
  external zheev
  allocate(rw(max(1,3*n-2)));lw=-1;call zheev('V','U',n,a,n,e,q,lw,rw,info)
  if(info/=0)then;deallocate(rw);return;endif
  lw=max(1,int(real(q(1),8)));allocate(w(lw));call zheev('V','U',n,a,n,e,w,lw,rw,info);deallocate(w,rw)
 end subroutine
 subroutine dg_metric_factor(s,n,cutoff,x,emin,emax,condition,discarded,info)
  integer,intent(in)::n;real(8),intent(in)::cutoff;complex(8),intent(in)::s(n,n)
  complex(8),intent(out)::x(n,n);real(8),intent(out)::emin,emax,condition;integer,intent(out)::discarded,info
  complex(8)::v(n,n);real(8)::e(n);integer::i
  if(maxval(abs(s-conjg(transpose(s))))>1d-10*max(1d0,maxval(abs(s))))then
   x=0;emin=0;emax=0;condition=huge(1d0);discarded=n;info=10;return
  endif
  v=.5d0*(s+conjg(transpose(s)));call heev(v,n,e,info)
  if(info/=0)then;x=0;emin=0;emax=0;condition=huge(1d0);discarded=n;return;endif
  emin=e(1);emax=e(n)
  if(cutoff<=0d0.or.emax<=0d0)then;x=0;condition=huge(1d0);discarded=n;info=11;return;endif
  discarded=count(e<=cutoff*emax)
  if(emin>0)then;condition=emax/emin;else;condition=huge(1d0);endif
  if(discarded>0)then;x=0;info=100+discarded;return;endif
  do i=1,n;v(:,i)=v(:,i)/(e(i)**0.25d0);enddo;x=matmul(v,conjg(transpose(v)));info=0
 end subroutine
 subroutine dg_metric_to_orth(x,s,c,d,n,nv)
  integer,intent(in)::n,nv;complex(8),intent(in)::x(n,n),s(n,n),c(n,nv);complex(8),intent(out)::d(n,nv)
  d=matmul(conjg(transpose(x)),matmul(s,c))
 end subroutine
 subroutine dg_metric_from_orth(x,d,c,n,nv)
  integer,intent(in)::n,nv;complex(8),intent(in)::x(n,n),d(n,nv);complex(8),intent(out)::c(n,nv);c=matmul(x,d)
 end subroutine
 subroutine dg_generalized_eigh(h,s,n,cutoff,e,c,residual,orthogonality,condition,discarded,info)
  integer,intent(in)::n;real(8),intent(in)::cutoff;complex(8),intent(in)::h(n,n),s(n,n)
  real(8),intent(out)::e(n),residual,orthogonality,condition;complex(8),intent(out)::c(n,n)
  integer,intent(out)::discarded,info;complex(8)::x(n,n),ho(n,n),r(n,n),eye(n,n);real(8)::emin,emax;integer::i
  call dg_metric_factor(s,n,cutoff,x,emin,emax,condition,discarded,info)
  if(info/=0)then;c=0;e=0;residual=huge(1d0);orthogonality=huge(1d0);return;endif
  if(maxval(abs(h-conjg(transpose(h))))>1d-10*max(1d0,maxval(abs(h))))then
   c=0;e=0;residual=huge(1d0);orthogonality=huge(1d0);info=20;return
  endif
  ho=matmul(conjg(transpose(x)),matmul(h,x));ho=.5d0*(ho+conjg(transpose(ho)));call heev(ho,n,e,info)
  if(info/=0)then;c=0;residual=huge(1d0);orthogonality=huge(1d0);return;endif
  c=matmul(x,ho)
  do i=1,n;r(:,i)=matmul(h,c(:,i))-e(i)*matmul(s,c(:,i));enddo
  residual=maxval(abs(r));eye=matmul(conjg(transpose(c)),matmul(s,c));do i=1,n;eye(i,i)=eye(i,i)-1;enddo
  orthogonality=maxval(abs(eye))
 end subroutine
 subroutine dg_reduced_generalized_eigh(h,s,n,nev,cutoff,e,c,effective_rank,residual,orthogonality,info)
  integer,intent(in)::n,nev;real(8),intent(in)::cutoff;complex(8),intent(in)::h(n,n),s(n,n)
  real(8),intent(out)::e(nev),residual,orthogonality;complex(8),intent(out)::c(n,nev)
  integer,intent(out)::effective_rank,info
  complex(8),allocatable::v(:,:),x(:,:),ho(:,:),r(:,:),eye(:,:)
  real(8),allocatable::se(:),he(:)
  real(8)::emax
  integer::i,j
  c=0;e=0;effective_rank=0;residual=huge(1d0);orthogonality=huge(1d0);info=1
  if(n<1.or.nev<1.or.nev>n.or.cutoff<=0d0)return
  if(maxval(abs(s-conjg(transpose(s))))>1d-10*max(1d0,maxval(abs(s))))then;info=10;return;endif
  if(maxval(abs(h-conjg(transpose(h))))>1d-10*max(1d0,maxval(abs(h))))then;info=20;return;endif
  allocate(v(n,n),se(n));v=.5d0*(s+conjg(transpose(s)));call heev(v,n,se,info)
  if(info/=0)return
  emax=se(n)
  if(emax<=0d0)then;info=11;return;endif
  effective_rank=count(se>cutoff*emax)
  if(effective_rank<nev)then;info=100+(nev-effective_rank);return;endif
  allocate(x(n,effective_rank));j=0
  do i=1,n
   if(se(i)<=cutoff*emax)cycle
   j=j+1;x(:,j)=v(:,i)/sqrt(se(i))
  enddo
  allocate(ho(effective_rank,effective_rank),he(effective_rank))
  ho=matmul(conjg(transpose(x)),matmul(h,x));ho=.5d0*(ho+conjg(transpose(ho)))
  call heev(ho,effective_rank,he,info);if(info/=0)return
  c=matmul(x,ho(:,1:nev));e=he(1:nev)
  allocate(r(n,nev),eye(nev,nev))
  do i=1,nev;r(:,i)=matmul(h,c(:,i))-e(i)*matmul(s,c(:,i));enddo
  residual=maxval(abs(r));eye=matmul(conjg(transpose(c)),matmul(s,c))
  do i=1,nev;eye(i,i)=eye(i,i)-1;enddo
  orthogonality=maxval(abs(eye));info=0
 end subroutine
 subroutine dg_s_exp_action(h,s,n,cutoff,dt,c0,c1,sunitarity,condition,discarded,info)
  integer,intent(in)::n;real(8),intent(in)::cutoff,dt;complex(8),intent(in)::h(n,n),s(n,n),c0(n)
  complex(8),intent(out)::c1(n);real(8),intent(out)::sunitarity,condition;integer,intent(out)::discarded,info
  complex(8)::x(n,n),ho(n,n),q(n,n),d(n),u(n,n),eye(n,n);real(8)::e(n),emin,emax,cond;integer::i,disc
  call dg_metric_factor(s,n,cutoff,x,emin,emax,cond,disc,info);condition=cond;discarded=disc
  if(info/=0)then;c1=0;sunitarity=huge(1d0);return;endif
  if(maxval(abs(h-conjg(transpose(h))))>1d-10*max(1d0,maxval(abs(h))))then
   c1=0;sunitarity=huge(1d0);info=20;return
  endif
  ho=matmul(conjg(transpose(x)),matmul(h,x));ho=.5d0*(ho+conjg(transpose(ho)));call heev(ho,n,e,info)
  if(info/=0)then;c1=0;sunitarity=huge(1d0);return;endif
  q=ho;d=matmul(conjg(transpose(q)),matmul(conjg(transpose(x)),matmul(s,c0)))
  do i=1,n;d(i)=d(i)*exp(cmplx(0d0,-dt*e(i),8));enddo;c1=matmul(x,matmul(q,d))
  ! Build U robustly by applying the formula to Cartesian columns.
  u=matmul(x,matmul(q,matmul(diag_phase(e,dt,n),matmul(conjg(transpose(q)),matmul(conjg(transpose(x)),s)))))
  eye=matmul(conjg(transpose(u)),matmul(s,u))-s;sunitarity=maxval(abs(eye))
 contains
  function diag_phase(ev,t,m) result(a);integer,intent(in)::m;real(8),intent(in)::ev(m),t;complex(8)::a(m,m);integer::j
   a=0;do j=1,m;a(j,j)=exp(cmplx(0d0,-t*ev(j),8));enddo
  end function
 end subroutine
end module
