! Numerical diagnostics only; no orbital or Hamiltonian state is modified.
module exx_stationarity
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 public::spinor_stationarity
contains
 subroutine spinor_stationarity(psi,hpsi,f,hvol,gram_error,hermitian_error,commutator,ierr)
  complex(8),intent(in)::psi(:,:,:),hpsi(:,:,:)
  real(8),intent(in)::f(:,:),hvol
  real(8),intent(out)::gram_error,hermitian_error,commutator
  integer,intent(out)::ierr
  complex(8),allocatable::s(:,:),h(:,:),q(:,:)
  real(8)::norm2,den,local
  integer::i,j,k,n
  ierr=1;gram_error=0d0;hermitian_error=0d0;commutator=0d0
  if(any(shape(psi)/=shape(hpsi)).or.size(psi,1)<1.or.size(psi,2)<1.or.size(psi,3)<1)return
  n=size(psi,2)
  if(any(shape(f)/=[n,size(psi,3)]))return
  if(.not.ieee_is_finite(hvol))return
  if(hvol<=0d0)return
  do k=1,size(psi,3);do j=1,n
   if(.not.ieee_is_finite(f(j,k)))return
   if(f(j,k)<0d0.or.f(j,k)>1d0)return
   do i=1,size(psi,1)
    if(.not.ieee_is_finite(real(psi(i,j,k),8)).or..not.ieee_is_finite(aimag(psi(i,j,k))))return
    if(.not.ieee_is_finite(real(hpsi(i,j,k),8)).or..not.ieee_is_finite(aimag(hpsi(i,j,k))))return
   enddo
  enddo;enddo
  allocate(s(n,n),h(n,n),q(size(psi,1),n))
  do k=1,size(psi,3)
   s=matmul(conjg(transpose(psi(:,:,k))),psi(:,:,k))*hvol
   h=matmul(conjg(transpose(psi(:,:,k))),hpsi(:,:,k))*hvol
   q=hpsi(:,:,k)-matmul(psi(:,:,k),h)
   do j=1,n
    s(j,j)=s(j,j)-1d0
   enddo
   ! The commutator formula assumes an orthonormal input span and a Hermitian
   ! operator. Report those defects separately so it cannot certify either.
   norm2=0d0;den=sum(f(:,k)**2)
   do j=1,n
    norm2=norm2+2d0*f(j,k)**2*sum(abs(q(:,j))**2)*hvol
    do i=1,n
     gram_error=max(gram_error,abs(s(i,j)))
     hermitian_error=max(hermitian_error,abs(h(i,j)-conjg(h(j,i))))
     norm2=norm2+(f(i,k)-f(j,k))**2*abs(h(i,j))**2
    enddo
   enddo
   local=sqrt(norm2/max(den,tiny(1d0)))
   if(.not.ieee_is_finite(local))return
   commutator=max(commutator,local)
  enddo
  if(.not.ieee_is_finite(gram_error).or..not.ieee_is_finite(hermitian_error))return
  ierr=0
 end subroutine
end module
