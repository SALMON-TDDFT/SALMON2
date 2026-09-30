program probe
 use exx_local_fft
 implicit none
 integer,parameter :: n(3)=[8,7,6],ng=336
 type(s_exx_local_fft) :: plan
 real(8) :: multiplier(8,7,6),pi,angle
 complex(8) :: kernel(ng),z
 complex(8),allocatable :: density(:),potential(:),reference(:)
 integer,allocatable :: points(:,:)
 integer :: a,b,c,g,j,l,trial,status,p(3),q(3),d(3),m(3),shift(3),ns,idx
 logical :: used
 pi=acos(-1d0)
 do c=0,n(3)-1;do b=0,n(2)-1;do a=0,n(1)-1
 multiplier(a+1,b+1,c+1)=2+cos(.7d0*a)+sin(.3d0*(b+2*c))
 enddo;enddo;enddo
 ! Independent inverse DFT; include a complex (non-even) kernel deliberately.
 do g=1,ng
 p=[mod(g-1,n(1)),mod((g-1)/n(1),n(2)),(g-1)/(n(1)*n(2))];z=0
 do c=0,n(3)-1;do b=0,n(2)-1;do a=0,n(1)-1
 q=[a,b,c];angle=2*pi*sum(real(p*q,8)/n)
 z=z+multiplier(a+1,b+1,c+1)*exp(cmplx(0d0,angle,8))/ng
 enddo;enddo;enddo
 kernel(g)=z
 enddo
 call exx_local_init(plan,multiplier,status)
 if(status/=0)error stop 'init'
 do trial=1,9
 m=[2,2,2]
 if(trial==2)m=[3,2,2]
 if(trial==3)m=1
 if(trial==4)m=[2,2,2]
 if(trial==5)m=n
 if(trial==8)m=[2,6,5]
 if(trial==9)m=[2,5,4]
 shift=n-[1,1,1]
 if(trial==6)shift=0
 if(trial==7)shift=[2,3,1]
 ns=product(m)
 allocate(points(3,ns),density(ns),potential(ns),reference(ns))
 g=0
 do c=0,m(3)-1;do b=0,m(2)-1;do a=0,m(1)-1
 g=g+1;points(:,g)=modulo([a,b,c]+shift,n)
 density(g)=cmplx(sin(real(g,8)),cos(.7d0*g),8)
 enddo;enddo;enddo
 call exx_local_prepare(plan,points,used,status)
 if(status/=0)error stop 'prepare'
 if(trial>=8.and.any(plan%padded/=[3,7,6]))error stop 'mixed periodic embedding'
 if(trial==5)then
  if(used)error stop 'extended source must fall back'
 else
  if(.not.used)error stop 'local support not selected'
  if(plan%fft_points>=ng)error stop 'FFT volume not reduced'
  do l=1,3
   if(l==2)density=0
   if(l==3)then
    do j=1,ns
     density(j)=1d-100*cmplx(sin(real(j,8)),cos(.7d0*j),8)
    enddo
   endif
   call exx_local_apply(plan,density,potential,status)
   if(status/=0)error stop 'local apply'
   reference=0
   do j=1,ns;do g=1,ns
    d=modulo(points(:,j)-points(:,g),n)
    idx=1+d(1)+n(1)*(d(2)+n(2)*d(3))
    reference(j)=reference(j)+kernel(idx)*density(g)
   enddo;enddo
   if(l/=2)then
    if(maxval(abs(reference-potential))>2d-12*maxval(abs(reference)))error stop 'local/direct relative mismatch'
   else
    if(any(potential/=(0d0,0d0)))error stop 'zero density'
   endif
  enddo
 endif
 deallocate(points,density,potential,reference)
 enddo
 call exx_local_destroy(plan)
 call exx_local_destroy(plan)
 print *,'local convolution: periodic wrap, complex kernel, full fallback and cache changes passed'
end program
