program probe
 use hybrid_semilocal, only: hybrid_density_evaluate,hse_semilocal_evaluate
 implicit none
 real(8)::r(1,2),g(1,3,2),c(2),cg(3,2),e(2),v(1,2),d(1,3,2),ep(2),em(2),vp(1,2),dp(1,3,2)
 real(8)::total(2),sig(2),ref(2),vr(2),vs(2),h
 integer::i,j,err
 r(1,:)=[.12d0,.3d0];c=[.2d0,.05d0];g=0.07d0;cg=-.02d0;h=1d-6
 call hybrid_density_evaluate(r,g,c,cg,e,v,d,err,.11d0)
 if(err/=0)stop 1
 total=r(1,:)+c;sig=sum((g(1,:,:)+cg)**2,dim=1)
 call hse_semilocal_evaluate(total,sig,ref,vr,vs,err,.11d0)
 if(err/=0.or.maxval(abs(e-total*ref))>1d-13.or.maxval(abs(v(1,:)-vr))>1d-13)stop 2
 do i=1,2
  r(1,i)=r(1,i)+h
  call hybrid_density_evaluate(r,g,c,cg,ep,vp,dp,err,.11d0)
  r(1,i)=r(1,i)-2*h
  call hybrid_density_evaluate(r,g,c,cg,em,vp,dp,err,.11d0)
  r(1,i)=r(1,i)+h
  if(abs((ep(i)-em(i))/(2*h)-v(1,i))>2d-8)stop 3
  do j=1,3
   g(1,j,i)=g(1,j,i)+h
   call hybrid_density_evaluate(r,g,c,cg,ep,vp,dp,err,.11d0)
   g(1,j,i)=g(1,j,i)-2*h
   call hybrid_density_evaluate(r,g,c,cg,em,vp,dp,err,.11d0)
   g(1,j,i)=g(1,j,i)+h
   if(abs((ep(i)-em(i))/(2*h)-d(1,j,i))>2d-8)stop 4
  enddo
 enddo
 print *, 'NLCC_LIBXC_AND_VARIATIONAL_DERIVATIVES_PASS'
end program
