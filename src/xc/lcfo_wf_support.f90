! Fixed periodic support masks. Omit exactly zero columns before reconstruction.
module lcfo_wf_support
 use iso_fortran_env, only: int64
 implicit none
 private
 public :: s_lcfo_wf_plan,lcfo_wf_plan_init,lcfo_wf_reconstruct,lcfo_wf_total_norm
 type :: s_lcfo_wf_plan
  logical :: ready=.false.,masked=.false.
  integer :: npoints=0,ncolumns=0
  integer,allocatable :: columns(:)
  logical,allocatable :: keep(:,:)
 end type
 public :: s_lcfo_wf_kernel,lcfo_wf_kernel_init,lcfo_wf_kernel_apply
 type :: s_wf_block
  integer,allocatable :: rows(:),columns(:),basis_columns(:)
  complex(8),allocatable :: basis(:,:)
 end type
 type :: s_lcfo_wf_kernel
  integer :: npoints=0,ncolumns=0,nbasis=0
  integer(int64) :: products=0
  logical :: ready=.false.
  type(s_wf_block),allocatable :: blocks(:)
 end type
contains
 subroutine lcfo_wf_kernel_init(kernel,plan,basis)
  type(s_lcfo_wf_kernel),intent(out) :: kernel
  type(s_lcfo_wf_plan),intent(in) :: plan
  complex(8),intent(in) :: basis(:,:)
  logical,allocatable :: keep(:,:),nonzero(:,:)
  integer,allocatable :: groups(:),representative(:)
  integer(int64),allocatable :: hashes(:)
  integer(int64) :: hash
  integer :: ng,nc,nb,g,j,k,n,match,r
  if(.not.plan%ready.or.size(basis,1)/=plan%npoints)error stop 'LCFO WF kernel: invalid basis'
  ng=size(basis,1);nb=size(basis,2);nc=size(plan%columns)
  kernel%npoints=ng;kernel%ncolumns=nc;kernel%nbasis=nb
  allocate(keep(nc,ng),nonzero(nb,ng),groups(ng),representative(ng),hashes(ng))
  keep=.true.
  if(plan%masked)keep=transpose(plan%keep)
  nonzero=transpose(basis/=(0d0,0d0))
  groups=0;n=0
  ! Rows with identical WF masks AND exact basis support form disjoint blocks.
  ! Each basis row is stored at most once; no per-WF duplication of the basis.
  do g=1,ng
   if(.not.any(keep(:,g)).or..not.any(nonzero(:,g)))cycle
   hash=0_int64
   do j=1,nc
    hash=ieor(ishftc(hash,1),int(merge(j,0,keep(j,g)),int64))
   enddo
   do j=1,nb
    hash=ieor(ishftc(hash,1),int(merge(j,0,nonzero(j,g)),int64))
   enddo
   match=0
   do k=1,n
    if(hashes(k)/=hash)cycle
    r=representative(k)
    if(.not.all(keep(:,r).eqv.keep(:,g)))cycle
    if(.not.all(nonzero(:,r).eqv.nonzero(:,g)))cycle
    match=k;exit
   enddo
   if(match==0)then
    n=n+1;match=n;representative(n)=g;hashes(n)=hash
   endif
   groups(g)=match
  enddo
  allocate(kernel%blocks(n))
  do k=1,n
   r=representative(k)
   associate(b=>kernel%blocks(k))
    b%rows=pack([(g,g=1,ng)],groups==k)
    b%columns=pack([(j,j=1,nc)],keep(:,r))
    b%basis_columns=pack([(j,j=1,nb)],nonzero(:,r))
    b%basis=basis(b%rows,b%basis_columns)
    kernel%products=kernel%products+int(size(b%rows),int64)*size(b%columns)*size(b%basis_columns)
   end associate
  enddo
  kernel%ready=.true.
 end subroutine
 subroutine lcfo_wf_kernel_apply(kernel,frame,wf)
  type(s_lcfo_wf_kernel),intent(in) :: kernel
  complex(8),intent(in) :: frame(:,:)
  complex(8),allocatable,intent(out) :: wf(:,:)
  complex(8),allocatable :: local_frame(:,:),result(:,:)
  integer :: k,j,i
  if(.not.kernel%ready.or.size(frame,1)/=kernel%nbasis.or.size(frame,2)/=kernel%ncolumns) &
   error stop 'LCFO WF kernel: incompatible frame'
  allocate(wf(kernel%npoints,kernel%ncolumns));wf=0d0
  do k=1,size(kernel%blocks)
   associate(b=>kernel%blocks(k))
    local_frame=frame(b%basis_columns,b%columns)
    ! Keep the BLAS result contiguous (indexed MATMUL LHS miscompiled on GNU15/AArch64).
    result=matmul(b%basis,local_frame)
    do j=1,size(b%columns);do i=1,size(b%rows)
     wf(b%rows(i),b%columns(j))=result(i,j)
    enddo;enddo
   end associate
  enddo
 end subroutine
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
 subroutine lcfo_wf_reconstruct(plan,basis,frame,wf,compact)
  type(s_lcfo_wf_plan),intent(in) :: plan
  complex(8),intent(in) :: basis(:,:),frame(:,:)
  complex(8),allocatable,intent(out) :: wf(:,:)
  logical,intent(in),optional :: compact
  logical :: selected
  integer :: nc
  selected=.false.
  if(present(compact))selected=compact
  nc=plan%ncolumns
  if(selected)nc=size(plan%columns)
  if(.not.plan%ready.or.size(basis,1)/=plan%npoints.or.size(frame,2)/=nc) &
    error stop 'LCFO WF support: incompatible reconstruction'
  if(.not.plan%masked)then
   wf=matmul(basis,frame)
  else
   allocate(wf(plan%npoints,size(plan%columns)))
   if(size(plan%columns)==0)return
   if(selected)then
    wf=matmul(basis,frame)
   else
    wf=matmul(basis,frame(:,plan%columns))
   endif
   where(.not.plan%keep)wf=0d0
  endif
 end subroutine
 real(8) function lcfo_wf_total_norm(gram,frame) result(norm)
  ! Exact identity for any basis, including a nonorthogonal one: ||B F||_F^2.
  complex(8),intent(in) :: gram(:,:),frame(:,:)
  norm=real(sum(conjg(frame)*matmul(gram,frame)),8)
 end function
end module
