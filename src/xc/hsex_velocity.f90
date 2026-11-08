! Derivative of P(A)^dagger Sigma(A) P(A), with dP/dA = i r P.
! Occupations, k weights and volume normalization are supplied by the caller.
module hsex_velocity
  implicit none
  private
  public :: hsex_velocity_diagonal
contains
  subroutine hsex_velocity_diagonal(p,rp,sigma,derivative,velocity)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    complex(8),intent(in) :: p(:,:),rp(:,:,:),sigma(:,:),derivative(:,:,:)
    real(8),intent(out) :: velocity(:,:)
    complex(8) :: hp(size(p,1),size(p,2)),dp(size(p,1),size(p,2))
    integer :: axis,band
    interface
      subroutine zgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        implicit none
        character(1),intent(in) :: transa,transb
        integer,intent(in) :: m,n,k,lda,ldb,ldc
        complex(8),intent(in) :: alpha,beta,a(lda,*),b(ldb,*)
        complex(8),intent(inout) :: c(ldc,*)
      end subroutine
    end interface
    call zgemm('N','N',size(p,1),size(p,2),size(p,1),(1d0,0d0),sigma,size(sigma,1), &
      p,size(p,1),(0d0,0d0),hp,size(hp,1))
    do axis=1,3
      call zgemm('N','N',size(p,1),size(p,2),size(p,1),(1d0,0d0),derivative(:,:,axis),size(derivative,1), &
        p,size(p,1),(0d0,0d0),dp,size(dp,1))
      do band=1,size(p,2)
        velocity(axis,band)=real(sum(conjg(p(:,band))*dp(:,band)),8) &
          +2d0*aimag(sum(conjg(rp(:,band,axis))*hp(:,band)))
        if(.not.ieee_is_finite(velocity(axis,band)))error stop 'nonfinite HSEX projected velocity'
      enddo
    enddo
  end subroutine hsex_velocity_diagonal
end module hsex_velocity
