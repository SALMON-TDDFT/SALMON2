! Manual multirate coupled-timestep analysis. IEEE inquiries are scalar.
module coupled_timestep_analysis
use iso_c_binding
use,intrinsic::ieee_arithmetic,only:ieee_is_finite,ieee_value,ieee_quiet_nan,ieee_positive_inf
implicit none
contains
integer(c_int) function study_reduce(n,m,dt,weight,volume,budget,e,used,j,qe,ne,out) bind(C) result(code)
 integer(c_int),value::n,m
 real(c_double),value::dt,weight,volume
 real(c_double),intent(in)::budget(19,m*n),e(3,m*n),used(3,m*n),j(3,0:n),qe(0:n),ne(0:n)
 real(c_double),intent(out)::out(33)
 real(8)::flux,loss,h,average(3),prior(3),block(3),alpha
 integer::i,d,q,l
 code=4;out=0d0
 if(n<2.or.m<1)return
 if(.not.ieee_is_finite(dt).or..not.ieee_is_finite(weight).or..not.ieee_is_finite(volume))return
 if(dt<=0d0.or.weight<=0d0.or.volume<=0d0)return
 h=dt/m
 do q=0,n
  if(.not.ieee_is_finite(qe(q)).or..not.ieee_is_finite(ne(q)))return
  if(ne(q)<=0d0)return
  do d=1,3
   if(.not.ieee_is_finite(j(d,q)))return
  enddo
  loss=(ne(0)-ne(q))/ne(0)
  out(14)=max(out(14),loss);out(15)=max(out(15),-loss)
  out(16)=max(out(16),sqrt(dot_product(j(:,q),j(:,q))))
  if(q>3*n/4)out(20)=max(out(20),sqrt(dot_product(j(:,q),j(:,q))))
 enddo
 prior=0d0
 do i=1,m*n
  do d=1,19
   if(.not.ieee_is_finite(budget(d,i)))return
  enddo
  do d=1,3
   if(.not.ieee_is_finite(e(d,i)).or..not.ieee_is_finite(used(d,i)))return
  enddo
  do d=9,14
   out(1)=out(1)-min(budget(d,i),0d0);out(2)=out(2)+max(budget(d,i),0d0)
  enddo
  out(3)=out(3)+budget(3,i);out(4)=out(4)+budget(6,i);out(5)=out(5)+budget(8,i)
  out(17)=max(out(17),sqrt(dot_product(e(:,i),e(:,i))))
  out(27)=max(out(27),abs(budget(16,i)))
  out(29)=out(29)+h*volume*dot_product(.5d0*(prior+e(:,i)),used(:,i))
  prior=e(:,i)
 enddo
 do q=1,n
  block=0d0
  do l=1,m
   block=block+e(:,(q-1)*m+l)
  enddo
  average=block/m
  out(6)=out(6)+dt*volume*dot_product(average,.5d0*(j(:,q-1)+j(:,q)))
  out(33)=out(33)+.5d0*h*volume*dot_product(block-e(:,q*m),j(:,q)-j(:,q-1))
 enddo
 ! Electron endpoints are q*dt+h/2. Interpolate internal energy to field time n*dt.
 alpha=.5d0/m
 out(7)=weight*(alpha*qe(n-1)+(1d0-alpha)*qe(n)-qe(0))
 flux=out(2)-out(1)
 out(8)=out(3)+out(5)+flux;out(9)=out(4)+out(5)+flux
 out(10)=out(3)+out(7)+flux;out(11)=out(4)+out(7)+flux
 out(12)=out(7)-out(5);out(19)=weight*(qe(n)-qe(0));out(13)=out(19)-out(6)
 out(18)=qe(0)
 if(out(1)>0d0)out(21:24)=out(8:11)/out(1)
 if(abs(out(5))>0d0)out(25)=out(12)/abs(out(5))
 if(abs(out(6))>0d0)out(26)=out(13)/abs(out(6))
 out(28)=flux;out(30)=out(5)-out(29)
 out(31)=.5d0*h*volume*dot_product(e(:,m*n),j(:,n))
 out(32)=out(6)-out(29)-out(31)-out(33)
 do i=1,33
  if(.not.ieee_is_finite(out(i)))return
 enddo
 code=0
end function

integer(c_int) function study_compare(nc,nf,nlc,nlf,dqc,dqf,dfc,dff,jc,jf,ec,ef,out) bind(C) result(code)
 integer(c_int),value::nc,nf,nlc,nlf
 real(c_double),value::dqc,dqf,dfc,dff
 real(c_double),intent(in)::jc(3,0:nc),jf(3,0:nf),ec(3,nlc),ef(3,nlf)
 real(c_double),intent(out)::out(6)
 real(8)::a,b,t,x,alpha,dj,de,sj,sr,se,ser
 integer::i,d,k,ratio
 code=4;out=0d0;sj=0d0;sr=0d0;se=0d0;ser=0d0
 if(nc<2.or.nf<2.or.nlc<2.or.nlf<nlc)return
 if(.not.ieee_is_finite(dqc).or..not.ieee_is_finite(dqf))return
 if(.not.ieee_is_finite(dfc).or..not.ieee_is_finite(dff))return
 if(min(dqc,dqf,dfc,dff)<=0d0)return
 if(mod(nlf,nlc)/=0)return
 ratio=nlf/nlc
 if(abs(dfc-ratio*dff)>64d0*epsilon(1d0)*dfc)return
 do i=0,nc-1
  t=.5d0*dfc+i*dqc;x=(t-.5d0*dff)/dqf;k=floor(x);alpha=x-k
  if(k<0.or.k>=nf)return
  do d=1,3
   a=jc(d,i);b=(1d0-alpha)*jf(d,k)+alpha*jf(d,k+1)
   if(.not.ieee_is_finite(a).or..not.ieee_is_finite(b))return
   dj=a-b;sj=sj+dj*dj;sr=sr+b*b;out(3)=max(out(3),abs(dj))
  enddo
 enddo
 do i=1,nlc
  do d=1,3
   a=ec(d,i);b=ef(d,i*ratio)
   if(.not.ieee_is_finite(a).or..not.ieee_is_finite(b))return
   de=a-b;se=se+de*de;ser=ser+b*b;out(6)=max(out(6),abs(de))
  enddo
 enddo
 out(1)=sqrt(sj/(3*nc));out(4)=sqrt(se/(3*nlc))
 if(sr>0d0)out(2)=sqrt(sj/sr)
 if(ser>0d0)out(5)=sqrt(se/ser)
 if(sr==0d0.and.sj>0d0)return
 if(ser==0d0.and.se>0d0)return
 do i=1,6
  if(.not.ieee_is_finite(out(i)))return
 enddo
 code=0
end function

integer(c_int) function study_difference(a,b,stats) bind(C) result(code)
 real(c_double),value::a,b
 real(c_double),intent(out)::stats(2)
 integer::i
 code=4;stats=0d0
 if(.not.ieee_is_finite(a).or..not.ieee_is_finite(b))return
 code=1
 if(b==0d0)return
 stats(1)=abs(a-b);stats(2)=stats(1)/abs(b)
 code=4
 do i=1,2
  if(.not.ieee_is_finite(stats(i)))return
 enddo
 code=0
end function

integer(c_int) function study_self_test() bind(C) result(code)
 real(8)::jc(3,0:4),jf(3,0:8),ec(3,8),ef(3,16),out(33),stats(6)
 real(8)::budget(19,4),qe(0:4),ne(0:4),used(3,4)
 integer::i,status
 code=1
 status=study_difference(2d0,1d0,stats(1:2))
 if(status/=0.or.any(abs(stats(1:2)-1d0)>1d-14))return
 do i=0,4
  jc(:,i)=i*.5d0+.125d0;qe(i)=2d0*i
 enddo
 do i=0,8
  jf(:,i)=i*.5d0+.0625d0
 enddo
 ec=1d0;ef=1d0
 status=study_compare(4_c_int,4_c_int,8_c_int,16_c_int,.5d0,.5d0,.25d0,.125d0,jc,jf,ec,ef,stats)
 if(status/=0.or.maxval(abs(stats))>1d-14)return
 jf=jf+.1d0
 status=study_compare(4_c_int,4_c_int,8_c_int,16_c_int,.5d0,.5d0,.25d0,.125d0,jc,jf,ec,ef,stats)
 if(status/=0.or.stats(3)<.09d0)return
 budget=0d0;budget(8,:)=2d0;budget(9,:)=-2d0;jc=2d0;ne=32d0;used=2d0
 status=study_reduce(4_c_int,1_c_int,.5d0,1d0,1d0,budget,ec(:,1:4),used,jc,qe,ne,out)
 if(status/=0)return
 if(abs(out(5)-8d0)>1d-14.or.abs(out(7)-7d0)>1d-14.or.abs(out(8))>1d-14)return
 if(abs(out(10)+1d0)>1d-14.or.abs(out(6)-12d0)>1d-14.or.abs(out(32))>1d-14)return
 do i=0,2
  jc(:,i)=dble(i);qe(i)=4d0*i
 enddo
 used(:,1:2)=0d0;used(:,3:4)=1d0
 status=study_reduce(2_c_int,2_c_int,1d0,1d0,1d0,budget,ec(:,1:4),used,jc(:,0:2),qe(0:2),ne(0:2),out)
 if(status/=0.or.abs(out(33)-1.5d0)>1d-14.or.abs(out(32))>1d-14)return
 ne(2)=ieee_value(0d0,ieee_quiet_nan)
 status=study_reduce(2_c_int,2_c_int,1d0,1d0,1d0,budget,ec(:,1:4),used,jc(:,0:2),qe(0:2),ne(0:2),out)
 if(status/=4)return
 ne(2)=ieee_value(0d0,ieee_positive_inf)
 status=study_reduce(2_c_int,2_c_int,1d0,1d0,1d0,budget,ec(:,1:4),used,jc(:,0:2),qe(0:2),ne(0:2),out)
 if(status/=4)return
 code=0
end function
end module
