program seed_probe
 use hse_wannier_gauge, only: gauge_seed,gauge_seed_gamma
 implicit none
 complex(8),allocatable :: coeff(:,:),saved(:,:),u(:,:,:),v(:,:),gram(:,:)
 real(8),allocatable :: position(:,:)
 integer :: n,ng,i,j,a,status,reference
 do a=1,4
  n=5;ng=19
  if(a==2)ng=n
  if(a==4)ng=n-1
  allocate(coeff(ng,n),saved(ng,n),u(n,n,1),v(n,n),gram(n,n),position(3,ng))
  do j=1,n;do i=1,ng
   coeff(i,j)=cmplx(sin(0.731d0*i*j),cos(0.219d0*i*(j+1)),8)
  enddo;enddo
  if(a==3)coeff(:,n)=coeff(:,1)
  saved=coeff;position=0d0
  call gauge_seed(reshape(coeff,[ng,n,1]),position,reshape([0d0,0d0,0d0],[3,1]),u,reference)
  call gauge_seed_gamma(coeff,v,status)
  if(status/=reference)error stop 'Seed status mismatch'
  if(any(coeff/=saved))error stop 'Input coefficients changed'
  if(a==3.or.a==4)then
   if(status==0)error stop 'Invalid seed accepted'
  else
   if(status/=0.or.maxval(abs(u(:,:,1)-v))>1d-12)error stop 'Seed rotation mismatch'
   gram=matmul(conjg(transpose(v)),v)
   do j=1,n;gram(j,j)=gram(j,j)-1d0;enddo
   if(maxval(abs(gram))>1d-12)error stop 'Nonunitary seed'
  endif
  deallocate(coeff,saved,u,v,gram,position)
 enddo
 allocate(coeff(3,0),v(0,0))
 call gauge_seed_gamma(coeff,v,status)
 if(status==0)error stop 'Empty seed accepted'
 deallocate(coeff,v)
 allocate(coeff(4,2),v(1,1));coeff=1d0
 call gauge_seed_gamma(coeff,v,status)
 if(status==0)error stop 'Wrong output shape accepted'
 deallocate(coeff,v)
 print *,'Gamma seed reference, unitarity, unchanged input and failure checks passed'
end program
