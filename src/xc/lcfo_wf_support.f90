! Fixed periodic support masks. Omit exactly zero columns before reconstruction.
module lcfo_wf_support
 implicit none
 private
 public :: s_lcfo_wf_plan,lcfo_wf_plan_init,lcfo_wf_reconstruct,lcfo_wf_total_norm
 type :: s_lcfo_wf_plan
  logical :: ready=.false.,masked=.false.
  integer :: npoints=0,ncolumns=0
  integer,allocatable :: columns(:)
  logical,allocatable :: keep(:,:)
 end type
contains
 subroutine lcfo_wf_plan_init(plan,positions,centers,length,radius,protected)
  type(s_lcfo_wf_plan),intent(inout) :: plan
  real(8),intent(in) :: positions(:,:),centers(:,:),length(3),radius
  logical,intent(in) :: protected(:)
  type(s_lcfo_wf_plan) :: empty
  logical,allocatable :: full_mask(:,:)
  real(8) :: delta(3)
  integer :: j,g,n
  plan=empty;n=size(centers,2)
  if(size(positions,1)/=3.or.size(centers,1)/=3.or.size(protected)/=n) &
    error stop 'LCFO WF support: incompatible geometry'
  plan%npoints=size(positions,2);plan%ncolumns=n
  if(radius<=0d0.or.all(protected))then
   plan%columns=[(j,j=1,n)]
  else
   allocate(full_mask(plan%npoints,n))
   do j=1,n;do g=1,plan%npoints
    delta=modulo(positions(:,g)-centers(:,j)+.5d0*length,length)-.5d0*length
    full_mask(g,j)=protected(j).or.sum(delta**2)<=radius**2
   enddo;enddo
   plan%columns=pack([(j,j=1,n)],any(full_mask,dim=1))
   plan%keep=full_mask(:,plan%columns)
   plan%masked=size(plan%columns)/=n.or..not.all(plan%keep)
  endif
  plan%ready=.true.
 end subroutine
 subroutine lcfo_wf_reconstruct(plan,basis,frame,wf)
  type(s_lcfo_wf_plan),intent(in) :: plan
  complex(8),intent(in) :: basis(:,:),frame(:,:)
  complex(8),allocatable,intent(out) :: wf(:,:)
  if(.not.plan%ready.or.size(basis,1)/=plan%npoints.or.size(frame,2)/=plan%ncolumns) &
    error stop 'LCFO WF support: incompatible reconstruction'
  if(.not.plan%masked)then
   wf=matmul(basis,frame)
  else
   allocate(wf(plan%npoints,size(plan%columns)))
   if(size(plan%columns)==0)return
   wf=matmul(basis,frame(:,plan%columns))
   where(.not.plan%keep)wf=0d0
  endif
 end subroutine
 real(8) function lcfo_wf_total_norm(gram,frame) result(norm)
  ! Exact identity for any basis, including a nonorthogonal one: ||B F||_F^2.
  complex(8),intent(in) :: gram(:,:),frame(:,:)
  norm=real(sum(conjg(frame)*matmul(gram,frame)),8)
 end function
end module
