! Distributed ACE action and fused projection in an orthonormal LCFO basis.
module lcfo_ace_local
 use communication,only:comm_summation
 implicit none
 private
 public :: lcfo_ace_local_action
contains
 subroutine lcfo_ace_local_action(basis,psi,hpsi,factors,grid_dv,ace_dv,comm)
  complex(8),intent(in),contiguous :: basis(:,:),psi(:,:)
  complex(8),intent(inout),contiguous :: hpsi(:,:)
  complex(8),intent(in) :: factors(:,:)
  real(8),intent(in) :: grid_dv,ace_dv
  integer,intent(in) :: comm
  complex(8),allocatable :: coeff(:,:),overlap(:,:),total(:,:),action(:,:)
  complex(8),parameter :: one=(1d0,0d0),zero=(0d0,0d0)
  integer :: ng,n,no,nf
  external :: zgemm
  ng=size(basis,1);n=size(basis,2);no=size(psi,2);nf=size(factors,2)
  if(size(psi,1)/=ng.or.any(shape(hpsi)/=shape(psi)).or.size(factors,1)/=n.or. &
     min(n,no,nf)<1.or.grid_dv<=0d0.or.ace_dv<=0d0)error stop 'LCFO local ACE: incompatible dimensions/weights'
  allocate(coeff(n,no),overlap(nf,no),total(nf,no),action(n,no))
  call zgemm('C','N',n,no,ng,one*grid_dv,basis,ng,psi,ng,zero,coeff,n)
  overlap=matmul(conjg(transpose(factors)),coeff)*ace_dv
  call comm_summation(overlap,total,size(total),comm)
  action=-matmul(factors,total)
  ! Project the non-exchange Hamiltonian and add exchange BEFORE reconstruction.
  call zgemm('C','N',n,no,ng,one*grid_dv,basis,ng,hpsi,ng,one,action,n)
  call zgemm('N','N',ng,no,n,one,basis,ng,action,n,zero,hpsi,ng)
 end subroutine
end module
