program probe
 use rvv10
 implicit none
 integer,parameter :: n(3)=[4,3,2],ng=24
 real(8) :: r(ng),s(ng),e(ng),v(ng),w(ng),ep(ng),vp(ng),wp(ng),em(ng),vm(ng),wm(ng)
 real(8) :: h(3)=[.7d0,.8d0,.9d0],delta,er,es,original,uniform,convergence(4),r_saved(ng),s_saved(ng)
 integer :: i,j,status
 do i=1,ng
  r(i)=.04d0+.015d0*sin(real(i,8));s(i)=.0002d0+.0001d0*cos(real(i,8))
 enddo
 call rvv10_evaluate(n,h,r,s,5.3d0,.0093d0,32,e,v,w,status)
 if(status/=0)error stop 'evaluate'
 er=0d0;es=0d0;delta=1d-7
 do i=1,ng
  original=r(i);r(i)=original+delta
  call rvv10_evaluate(n,h,r,s,5.3d0,.0093d0,32,ep,vp,wp,status)
  r(i)=original-delta
  call rvv10_evaluate(n,h,r,s,5.3d0,.0093d0,32,em,vm,wm,status)
  r(i)=original
  er=max(er,abs(sum(ep-em)/(2*delta)-v(i)))
  original=s(i);s(i)=original+delta
  call rvv10_evaluate(n,h,r,s,5.3d0,.0093d0,32,ep,vp,wp,status)
  s(i)=original-delta
  call rvv10_evaluate(n,h,r,s,5.3d0,.0093d0,32,em,vm,wm,status)
  s(i)=original
  es=max(es,abs(sum(ep-em)/(2*delta)-w(i)))
 enddo
 write(*,*)er,es
 r_saved=r;s_saved=s
 do j=1,4
  call rvv10_evaluate(n,h,r,s,5.3d0,.0093d0,8*2**j,ep,vp,wp,status)
  if(status/=0)error stop 'q mesh'
  convergence(j)=sum(ep)*product(h)
 enddo
 write(*,*)convergence
 r=reshape(cshift(reshape(r,n),1,dim=1),[ng]);s=reshape(cshift(reshape(s,n),1,dim=1),[ng])
 call rvv10_evaluate(n,h,r,s,5.3d0,.0093d0,32,ep,vp,wp,status)
 if(status/=0)error stop 'translated state'
 if(maxval(abs(ep-reshape(cshift(reshape(e,n),1,dim=1),[ng])))>1d-13)error stop 'translation energy'
 if(maxval(abs(vp-reshape(cshift(reshape(v,n),1,dim=1),[ng])))>1d-13)error stop 'translation potential'
 if(maxval(abs(wp-reshape(cshift(reshape(w,n),1,dim=1),[ng])))>1d-13)error stop 'translation sigma'
 r=r_saved;s=s_saved
 r=.04d0;s=0d0
 call rvv10_evaluate(n,h,r,s,5.3d0,.0093d0,64,e,v,w,status)
 write(*,*)maxval(abs(e/r)),maxval(abs(v)),maxval(abs(w))
 r=0d0;s=0d0
 call rvv10_evaluate(n,h,r,s,5.3d0,.0093d0,32,e,v,w,status)
 if(status/=0.or.maxval(abs(e))/=0d0)error stop 'vacuum'
 do i=1,3
  do j=0,3
   write(*,*)rvv10_kernel_fourier(.02d0,.02d0*i,.1d0*j)
  enddo
 enddo
end program
