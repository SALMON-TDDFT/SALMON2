! Exact restriction of a discrete periodic convolution to compact grid support.
! Padding removes local circular aliasing; the kernel remains the GLOBAL one.
! A plan belongs to one immutable multiplier and is not concurrently reentrant.
module exx_local_fft
 use iso_c_binding
 use iso_fortran_env, only: int64
 implicit none
 private
 include 'fftw3.f03'
 public :: exx_local_prepare_compact,smooth_size
 public :: s_exx_local_fft,exx_local_init,exx_local_prepare,exx_local_apply,exx_local_destroy
 type s_exx_local_fft
  integer :: n(3)=0,box(3)=0,padded(3)=0,fft_points=0
  logical :: ready=.false.
  complex(c_double_complex),allocatable :: kernel(:,:,:),filter(:,:,:),work(:,:,:)
  integer,allocatable :: indices(:)
  type(c_ptr) :: forward=c_null_ptr,backward=c_null_ptr
 end type
contains
 subroutine clear_box(plan)
  type(s_exx_local_fft),intent(inout) :: plan
  if(c_associated(plan%forward))call fftw_destroy_plan(plan%forward)
  if(c_associated(plan%backward))call fftw_destroy_plan(plan%backward)
  plan%forward=c_null_ptr;plan%backward=c_null_ptr
  if(allocated(plan%filter))deallocate(plan%filter,plan%work)
  if(allocated(plan%indices))deallocate(plan%indices)
  plan%box=0;plan%padded=0;plan%fft_points=0;plan%ready=.false.
 end subroutine

 subroutine exx_local_destroy(plan)
  type(s_exx_local_fft),intent(inout) :: plan
  call clear_box(plan)
  if(allocated(plan%kernel))deallocate(plan%kernel)
  plan%n=0
 end subroutine

 subroutine exx_local_init(plan,multiplier,status)
  type(s_exx_local_fft),intent(inout) :: plan
  real(8),intent(in) :: multiplier(:,:,:)
  integer,intent(out) :: status
  type(c_ptr) :: inverse
  call exx_local_destroy(plan)
  status=1;plan%n=shape(multiplier)
  if(any(plan%n<1))return
  allocate(plan%kernel(plan%n(1),plan%n(2),plan%n(3)))
  plan%kernel=cmplx(multiplier,0d0,c_double_complex)
  inverse=fftw_plan_dft_3d(plan%n(3),plan%n(2),plan%n(1),plan%kernel,plan%kernel,FFTW_BACKWARD,FFTW_ESTIMATE)
  if(.not.c_associated(inverse))then
   call exx_local_destroy(plan);return
  endif
  call fftw_execute_dft(inverse,plan%kernel,plan%kernel)
  call fftw_destroy_plan(inverse)
  plan%kernel=plan%kernel/real(product(int(plan%n,int64)),8)
  status=0
 end subroutine

 integer function smooth_size(lower) result(value)
  integer,intent(in) :: lower
  integer :: remainder,factor
  value=lower
  do
   remainder=value
   do factor=2,5
    if(factor==4)cycle
    do while(mod(remainder,factor)==0)
     remainder=remainder/factor
    enddo
   enddo
   if(remainder==1)return
   value=value+1
  enddo
 end function

 subroutine exx_local_prepare(plan,points,used,status)
  type(s_exx_local_fft),intent(inout) :: plan
  integer,intent(in) :: points(:,:) ! unique zero-based global grid coordinates
  logical,intent(out) :: used
  integer,intent(out) :: status
  integer :: axis,j,ns,gap,best_gap,start,origin(3),box(3),padded(3),a,b,c,p(3),q(3)
  integer,allocatable :: occupied(:)
  logical,allocatable :: seen(:)
  status=1;used=.false.;plan%ready=.false.
  if(.not.allocated(plan%kernel).or.size(points,1)/=3)return
  ns=size(points,2)
  if(ns==0)then
   status=0;return
  endif
  do axis=1,3
   if(any(points(axis,:)<0).or.any(points(axis,:)>=plan%n(axis)))return
   allocate(seen(plan%n(axis)));seen=.false.
   seen(points(axis,:)+1)=.true.
   occupied=pack([(j-1,j=1,plan%n(axis))],seen)
   best_gap=-1;start=1
   do j=1,size(occupied)
    if(j<size(occupied))then
     gap=occupied(j+1)-occupied(j)
    else
     gap=occupied(1)+plan%n(axis)-occupied(j)
    endif
    if(gap>best_gap)then
     best_gap=gap;start=mod(j,size(occupied))+1
    endif
   enddo
   origin(axis)=occupied(start)
   box(axis)=plan%n(axis)-best_gap+1
   padded(axis)=smooth_size(2*box(axis)-1)
   deallocate(seen,occupied)
  enddo
  status=0
  if(product(int(padded,int64))>=product(int(plan%n,int64)))return
  if(any(box/=plan%box))then
   call clear_box(plan)
   plan%box=box;plan%padded=padded;plan%fft_points=product(padded)
   allocate(plan%filter(padded(1),padded(2),padded(3)),plan%work(padded(1),padded(2),padded(3)))
   plan%forward=fftw_plan_dft_3d(padded(3),padded(2),padded(1),plan%work,plan%work,FFTW_FORWARD,FFTW_ESTIMATE)
   plan%backward=fftw_plan_dft_3d(padded(3),padded(2),padded(1),plan%work,plan%work,FFTW_BACKWARD,FFTW_ESTIMATE)
   if(.not.c_associated(plan%forward).or..not.c_associated(plan%backward))then
    call clear_box(plan);status=1;return
   endif
   plan%work=0d0
   do c=1-box(3),box(3)-1;do b=1-box(2),box(2)-1;do a=1-box(1),box(1)-1
    p=modulo([a,b,c],padded)+1;q=modulo([a,b,c],plan%n)+1
    plan%work(p(1),p(2),p(3))=plan%kernel(q(1),q(2),q(3))
   enddo;enddo;enddo
   call fftw_execute_dft(plan%forward,plan%work,plan%work)
   plan%filter=plan%work
  endif
  if(allocated(plan%indices))deallocate(plan%indices)
  allocate(plan%indices(ns),seen(plan%fft_points));seen=.false.
  do j=1,ns
   p=modulo(points(:,j)-origin,plan%n)
   if(any(p>=box))then
    status=1;return
   endif
   a=1+p(1)+padded(1)*(p(2)+padded(2)*p(3))
   if(seen(a))then
    status=1;return
   endif
   seen(a)=.true.;plan%indices(j)=a
  enddo
  used=.true.;plan%ready=.true.
 end subroutine

 ! Kernel entries are indexed by displacement + box, over -(box-1):box-1.
 ! This interface never stores a full global kernel; it accepts a compact tile.
 subroutine exx_local_prepare_compact(plan,n,box,kernel,status)
  type(s_exx_local_fft),intent(inout) :: plan
  integer,intent(in) :: n(3),box(3)
  complex(8),intent(in) :: kernel(:,:,:)
  integer,intent(out) :: status
  integer :: padded(3),a,b,c,p(3),j
  status=1
  if(any(n<1).or.any(box<1).or.any(box>n))return
  if(any(shape(kernel)/=2*box-1))return
  padded=[smooth_size(2*box(1)-1),smooth_size(2*box(2)-1),smooth_size(2*box(3)-1)]
  if(product(int(padded,int64))>=product(int(n,int64)))return
  if(any(box/=plan%box).or..not.c_associated(plan%forward))then
   call clear_box(plan)
   plan%box=box;plan%padded=padded;plan%fft_points=product(padded)
   allocate(plan%filter(padded(1),padded(2),padded(3)),plan%work(padded(1),padded(2),padded(3)))
   plan%forward=fftw_plan_dft_3d(padded(3),padded(2),padded(1),plan%work,plan%work,FFTW_FORWARD,FFTW_ESTIMATE)
   plan%backward=fftw_plan_dft_3d(padded(3),padded(2),padded(1),plan%work,plan%work,FFTW_BACKWARD,FFTW_ESTIMATE)
   if(.not.c_associated(plan%forward).or..not.c_associated(plan%backward))then
    call clear_box(plan);return
   endif
  endif
  plan%n=n;plan%work=0d0
  do c=1-box(3),box(3)-1;do b=1-box(2),box(2)-1;do a=1-box(1),box(1)-1
   p=modulo([a,b,c],padded)+1
   plan%work(p(1),p(2),p(3))=kernel(a+box(1),b+box(2),c+box(3))
  enddo;enddo;enddo
  call fftw_execute_dft(plan%forward,plan%work,plan%work);plan%filter=plan%work
  if(allocated(plan%indices))deallocate(plan%indices)
  allocate(plan%indices(product(box)));j=0
  do c=0,box(3)-1;do b=0,box(2)-1;do a=0,box(1)-1
   j=j+1;plan%indices(j)=1+a+padded(1)*(b+padded(2)*c)
  enddo;enddo;enddo
  plan%ready=.true.;status=0
 end subroutine

 subroutine exx_local_apply(plan,density,potential,status)
  type(s_exx_local_fft),intent(inout),target :: plan
  complex(8),intent(in) :: density(:)
  complex(8),intent(out) :: potential(:)
  integer,intent(out) :: status
  complex(c_double_complex),pointer :: flat(:)
  status=1;potential=0d0
  if(.not.plan%ready)return
  if(size(density)/=size(plan%indices).or.size(potential)/=size(density))return
  flat(1:plan%fft_points)=>plan%work
  flat=0d0;flat(plan%indices)=density
  call fftw_execute_dft(plan%forward,plan%work,plan%work)
  plan%work=plan%work*plan%filter
  call fftw_execute_dft(plan%backward,plan%work,plan%work)
  potential=flat(plan%indices)/plan%fft_points
  status=0
 end subroutine
end module exx_local_fft
