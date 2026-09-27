program exact_pair_probe
 use iso_fortran_env,only:int64
 use, intrinsic :: ieee_arithmetic,only:ieee_is_finite
 use hse_wannier
 implicit none
 type(s_hse_wannier) :: op
 complex(8) :: target(24,12,1),action(24,12,1),dense_action(24,12,1),ref(24,12),kernel(24),metric(12,12)
 complex(8) :: density(24),potential(24),v
 real(8) :: h(3),k(3,1),angle,pi
 integer :: batch
 integer :: n(3),p(3),q(3),d(3),g,a,b,c,j,i,idx,status,trial
 integer(int64) :: expected,expected_accum
 n=[4,3,2];h=[.6d0,.8d0,.9d0];k=0d0;pi=acos(-1d0)
 call wannier_init(op,n,[1,1,1],h,k,.11d0,status)
 if(status/=0)error stop 'init'
 allocate(op%source(24,4));op%source=0d0
 op%source(1,1)=(1d0,.3d0);op%source(2,1)=(.2d0,-.1d0)
 op%source(15,2)=(.8d0,-.4d0)
 op%source(1,3)=(1d-100,2d-100)
 ! Independent inverse DFT of the Fourier multiplier, not an FFT reference.
 do g=1,24
  p=op%point(:,g);v=0d0
  do c=0,n(3)-1;do b=0,n(2)-1;do a=0,n(1)-1
   q=[a,b,c];angle=2*pi*sum(dble(p*q)/n)
   v=v+op%multiplier(a+1,b+1,c+1)*exp(cmplx(0d0,angle,8))/24d0
  enddo;enddo;enddo
  kernel(g)=v
 enddo
 do batch=1,8
 op%fft_batch_size=batch
 do trial=1,5
  target=0d0;target(:,1:4,1)=op%source
  if(trial==2)then
   target=0d0;target(1,1,1)=(.3d0,.2d0);target(15,2,1)=(.5d0,.6d0)
   target(20,3,1)=(.7d0,-.1d0);target(2,4,1)=(1d-100,-1d-100)
  endif
  if(trial==3)target=0d0
  if(trial>=4)then
   do j=1,12;do g=1,24
    target(g,j,1)=cmplx(sin(.2d0*g*j),cos(.3d0*(g+j)),8)
   enddo;enddo
   if(trial==5)target(:,6:,1)=0d0
  endif
  ref=0d0;expected=0;expected_accum=0
  do i=1,4;do j=1,12
   density=conjg(op%source(:,i))*target(:,j,1)
   if(any(density/=(0d0,0d0)))then
    expected=expected+1
    expected_accum=expected_accum+count(op%source(:,i)/=(0d0,0d0))
   endif
   potential=0d0
   do g=1,24;do idx=1,24
    d=modulo(op%point(:,g)-op%point(:,idx),n)
    a=1+d(1)+n(1)*(d(2)+n(2)*d(3))
    potential(g)=potential(g)+kernel(a)*density(idx)
   enddo;enddo
   ref(:,j)=ref(:,j)-op%source(:,i)*potential
  enddo;enddo
  op%compact_source_support=.false.
  call wannier_apply(op,target,dense_action,status)
  if(status/=0)error stop 'dense apply'
  op%compact_source_support=.true.
  call wannier_apply(op,target,action,status)
  if(op%pair_accumulation_points/=expected_accum)error stop 'compact accumulation count'
  if(op%pair_product_points/=48_int64)error stop 'compact source products not used'
  if(.not.all(ieee_is_finite(real(action))).or..not.all(ieee_is_finite(aimag(action))))error stop 'nonfinite action'
  if(maxval(abs(action-dense_action))>1d-12*max(1d0,maxval(abs(dense_action))))error stop 'dense compact mismatch'
  if(status/=0)error stop 'apply'
  if(op%fft_pairs_total/=48_int64.or.op%fft_pairs_executed/=expected)error stop 'pair count'
  if(trial==4.and.op%worker_batch==8.and.op%fft_batches_executed/=6)error stop 'full and tail batches'
  if(op%fft_batches_executed>expected)error stop 'padded empty FFTs'
  if(expected==0.and.op%fft_batches_executed/=0)error stop 'empty tile executed'
  if(trial==1.and.expected/=5_int64)error stop 'tiny nonzero pair discarded'
  if(maxval(abs(action(:,:,1)-ref))>1d-11*max(1d0,maxval(abs(ref))))error stop 'direct action'
  do j=1,12
   if(maxval(abs(ref(:,j)))>0d0)then
    if(maxval(abs(action(:,j,1)-ref(:,j)))/maxval(abs(ref(:,j)))>1d-11)error stop 'relative tiny action'
   endif
  enddo
  if(trial==1)then
   metric=matmul(conjg(transpose(target(:,:,1))),action(:,:,1))
   if(maxval(abs(metric-conjg(transpose(metric))))>1d-11)error stop 'Hermiticity'
  endif
 enddo
 enddo
 call wannier_destroy(op)
 call translated_support()
 print *, 'Exact pair FFT screening and direct convolution passed'
contains
 subroutine translated_support()
  implicit none
  type(s_hse_wannier) :: multi
  complex(8) :: targets(24,3,2),out(24,3,2),dense(24,3,2),expected_out(24,3,2)
  complex(8) :: home(48,3),result(48,3),src(48),rho(48),pot(48),kern(48),z
  real(8) :: kv(3,2),theta
  integer :: ii,jj,gg,hh,cell,ix,iy,iz,status2,rr(3),qq(3),dd(3),index2,ww
  kv=0d0;kv(1,2)=pi/(n(1)*h(1))
  call wannier_init(multi,n,[2,1,1],h,kv,.11d0,status2)
  if(status2/=0)error stop 'multi init'
  allocate(multi%source(48,2));multi%source=0d0
  multi%source(1,1)=(.7d0,.2d0);multi%source(47,1)=(.3d0,-.1d0)
  multi%source(22,2)=(1d-100,-2d-100)
  do cell=1,2;do jj=1,3;do gg=1,24
   targets(gg,jj,cell)=cmplx(sin(.2d0*gg*jj*cell),cos(.3d0*(gg+jj+cell)),8)
  enddo;enddo;enddo
  call wannier_forward(multi,targets,home)
  do gg=1,48
   rr=multi%point(:,gg);z=0d0
   do iz=0,multi%ns(3)-1;do iy=0,multi%ns(2)-1;do ix=0,multi%ns(1)-1
    qq=[ix,iy,iz];theta=2*pi*sum(dble(rr*qq)/multi%ns)
    z=z+multi%multiplier(ix+1,iy+1,iz+1)*exp(cmplx(0d0,theta,8))/48d0
   enddo;enddo;enddo
   kern(gg)=z
  enddo
  result=0d0
  do cell=1,2;do ii=1,2
   do gg=1,48
    rr=modulo(multi%point(:,gg)-multi%cell(:,cell),multi%ns)
    index2=1+rr(1)+multi%ns(1)*(rr(2)+multi%ns(2)*rr(3))
    src(gg)=multi%source(index2,ii)
   enddo
   do jj=1,3
    rho=conjg(src)*home(:,jj);pot=0d0
    do gg=1,48;do hh=1,48
     dd=modulo(multi%point(:,gg)-multi%point(:,hh),multi%ns)
     index2=1+dd(1)+multi%ns(1)*(dd(2)+multi%ns(2)*dd(3))
     pot(gg)=pot(gg)+kern(index2)*rho(hh)
    enddo;enddo
    result(:,jj)=result(:,jj)-src*pot
   enddo
  enddo;enddo
  call wannier_backward(multi,result,expected_out)
  do ww=1,3
   multi%fft_batch_size=ww;multi%compact_source_support=.false.
   call wannier_apply(multi,targets,dense,status2)
   if(status2/=0)error stop 'multi dense'
   multi%compact_source_support=.true.
   call wannier_apply(multi,targets,out,status2)
   if(status2/=0.or.multi%pair_product_points/=18_int64)error stop 'multi compact'
   if(.not.all(ieee_is_finite(real(out))).or..not.all(ieee_is_finite(aimag(out))))error stop 'multi finite'
   if(maxval(abs(out-dense))>1d-11.or.maxval(abs(out-expected_out))>1d-11)error stop 'multi convolution'
  enddo
  call wannier_destroy(multi)
 end subroutine
end program
