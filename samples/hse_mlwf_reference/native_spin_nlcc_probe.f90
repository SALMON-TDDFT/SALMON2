program probe
 use hybrid_semilocal, only: hybrid_density_evaluate,hse_semilocal_evaluate
 implicit none
 real(8)::r(2,2),g(2,3,2),c(2),cg(3,2),e(2),v(2,2),d(2,3,2),ep(2),em(2),vp(2,2),dp(2,3,2)
 real(8)::total(2),sig(2),ref(2),vr(2),vs(2),h
 integer::i,j,is,err
 r(1,:)=[.12d0,.3d0];r(2,:)=[.08d0,.2d0];c=[.2d0,.05d0];g=0.07d0;g(2,:,:)=.11d0;cg=-.02d0;h=1d-6
 call hybrid_density_evaluate(r,g,c,cg,e,v,d,err,.11d0)
 if(err/=0)stop 1
 do is=1,2
 do i=1,2
  r(is,i)=r(is,i)+h
  call hybrid_density_evaluate(r,g,c,cg,ep,vp,dp,err,.11d0)
  r(is,i)=r(is,i)-2*h
  call hybrid_density_evaluate(r,g,c,cg,em,vp,dp,err,.11d0)
  r(is,i)=r(is,i)+h
  if(abs((ep(i)-em(i))/(2*h)-v(is,i))>2d-8)stop 3
  do j=1,3
   g(is,j,i)=g(is,j,i)+h
   call hybrid_density_evaluate(r,g,c,cg,ep,vp,dp,err,.11d0)
   g(is,j,i)=g(is,j,i)-2*h
   call hybrid_density_evaluate(r,g,c,cg,em,vp,dp,err,.11d0)
   g(is,j,i)=g(is,j,i)+h
   if(abs((ep(i)-em(i))/(2*h)-d(is,j,i))>2d-8)stop 4
  enddo
 enddo
 enddo
 ! Equal spin densities must reproduce the unpolarized Libxc functional.
 r(1,:)=[.1d0,.2d0];r(2,:)=r(1,:);g=.03d0
 call hybrid_density_evaluate(r,g,c,cg,e,v,d,err,.11d0)
 if(err/=0)stop 5
 total=sum(r,dim=1)+c;sig=sum((g(1,:,:)+g(2,:,:)+cg)**2,dim=1)
 call hse_semilocal_evaluate(total,sig,ref,vr,vs,err,.11d0)
 if(err/=0.or.maxval(abs(e-total*ref))>1d-12)stop 6
 if(maxval(abs(v(1,:)-vr))>1d-12.or.maxval(abs(v(2,:)-vr))>1d-12)stop 7
 r=0d0;g=0d0;c=0d0;cg=0d0
 call hybrid_density_evaluate(r,g,c,cg,e,v,d,err,.11d0)
 if(err/=0.or.maxval(abs(e))>1d-14.or.maxval(abs(v))>1d-14)stop 8
 print *, 'SPIN_NLCC_VARIATIONAL_DERIVATIVES_CLOSED_SHELL_VACUUM_PASS'
end program
