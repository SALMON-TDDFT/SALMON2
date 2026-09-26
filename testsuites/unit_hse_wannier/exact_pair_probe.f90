program exact_pair_probe
 use iso_fortran_env,only:int64
 use hse_wannier
 implicit none
 type(s_hse_wannier) :: op
 complex(8) :: target(24,4,1),action(24,4,1),ref(24,4),kernel(24),metric(4,4)
 complex(8) :: density(24),potential(24),v
 real(8) :: h(3),k(3,1),angle,pi
 integer :: n(3),p(3),q(3),d(3),g,a,b,c,j,i,idx,status,trial
 integer(int64) :: expected
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
 do trial=1,3
  target(:,:,1)=op%source
  if(trial==2)then
   target=0d0;target(1,1,1)=(.3d0,.2d0);target(15,2,1)=(.5d0,.6d0)
   target(20,3,1)=(.7d0,-.1d0);target(2,4,1)=(1d-100,-1d-100)
  endif
  if(trial==3)target=0d0
  ref=0d0;expected=0
  do i=1,4;do j=1,4
   density=conjg(op%source(:,i))*target(:,j,1)
   if(any(density/=(0d0,0d0)))expected=expected+1
   potential=0d0
   do g=1,24;do idx=1,24
    d=modulo(op%point(:,g)-op%point(:,idx),n)
    a=1+d(1)+n(1)*(d(2)+n(2)*d(3))
    potential(g)=potential(g)+kernel(a)*density(idx)
   enddo;enddo
   ref(:,j)=ref(:,j)-op%source(:,i)*potential
  enddo;enddo
  call wannier_apply(op,target,action,status)
  if(status/=0)error stop 'apply'
  if(op%fft_pairs_total/=16_int64.or.op%fft_pairs_executed/=expected)error stop 'pair count'
  if(trial==1.and.expected/=5_int64)error stop 'tiny nonzero pair discarded'
  if(maxval(abs(action(:,:,1)-ref))>1d-11*max(1d0,maxval(abs(ref))))error stop 'direct action'
  do j=1,4
   if(maxval(abs(ref(:,j)))>0d0)then
    if(maxval(abs(action(:,j,1)-ref(:,j)))/maxval(abs(ref(:,j)))>1d-11)error stop 'relative tiny action'
   endif
  enddo
  if(trial==1)then
   metric=matmul(conjg(transpose(target(:,:,1))),action(:,:,1))
   if(maxval(abs(metric-conjg(transpose(metric))))>1d-11)error stop 'Hermiticity'
  endif
 enddo
 call wannier_destroy(op)
 print *, 'Exact pair FFT screening and direct convolution passed'
end program
