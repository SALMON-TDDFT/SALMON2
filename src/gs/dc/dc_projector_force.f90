! Derivative of the real core-restricted nonlocal projector expectation.
module dc_projector_force
  implicit none
  private
  public :: core_projector_force
contains
  pure subroutine core_projector_force(value,derivative,in_core,weight,force)
    complex(8),intent(in) :: value(:),derivative(:,:)
    logical,intent(in) :: in_core(:)
    real(8),intent(in) :: weight
    real(8),intent(out) :: force(3)
    complex(8) :: full,core,dfull,dcore
    integer :: axis
    full=sum(value);core=sum(value,mask=in_core)
    do axis=1,3
      dfull=sum(derivative(axis,:));dcore=sum(derivative(axis,:),mask=in_core)
      force(axis)=-weight*real(conjg(dcore)*full+conjg(core)*dfull,8)
    enddo
  end subroutine
end module
