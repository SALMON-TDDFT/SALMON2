! Cache the exact nonzero support of the fixed core basis for exchange projection.
module lcfo_projection
 implicit none
 private
 public :: s_lcfo_projection,lcfo_projection_init,lcfo_projection_apply
 type :: s_lcfo_projection
  logical :: ready=.false.
  integer :: npoints=0,ncolumns=0
  integer,allocatable :: rows(:),columns(:)
  complex(8),allocatable :: left(:,:)
 end type
contains
 subroutine lcfo_projection_init(plan,basis,scale)
  type(s_lcfo_projection),intent(inout) :: plan
  complex(8),intent(in) :: basis(:,:)
  real(8),intent(in) :: scale
  type(s_lcfo_projection) :: empty
  integer :: i
  plan=empty;plan%npoints=size(basis,1);plan%ncolumns=size(basis,2)
  plan%rows=pack([(i,i=1,plan%npoints)],any(basis/=(0d0,0d0),dim=2))
  plan%columns=pack([(i,i=1,plan%ncolumns)],any(basis/=(0d0,0d0),dim=1))
  plan%left=scale*conjg(transpose(basis(plan%rows,plan%columns)))
  plan%ready=.true.
 end subroutine
 subroutine lcfo_projection_apply(plan,action,projected)
  type(s_lcfo_projection),intent(in) :: plan
  complex(8),intent(in) :: action(:,:)
  complex(8),allocatable,intent(out) :: projected(:,:)
  complex(8),allocatable :: block(:,:),near_action(:,:)
  integer :: j
  if(.not.plan%ready)error stop 'LCFO projection: uninitialized plan'
  if(any(shape(action)/=[plan%npoints,plan%ncolumns]))error stop 'LCFO projection: incompatible action'
  allocate(projected(plan%ncolumns,plan%ncolumns));projected=0d0
  if(size(plan%rows)>0.and.size(plan%columns)>0)then
    ! Keep BLAS operands/output contiguous; scatter only after the product.
    ! Indexed-LHS MATMUL fails with GNU15/AArch64 external-BLAS optimization.
    near_action=action(plan%rows,:)
    block=matmul(plan%left,near_action)
    do j=1,size(plan%columns)
      projected(plan%columns(j),:)=block(j,:)
    enddo
  endif
  projected=.5d0*(projected+conjg(transpose(projected)))
 end subroutine
end module
