module ptcn_rectangular_fft
 use iso_c_binding
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 include 'fftw3.f03'
 public::ptcn_fft_plan,ptcn_fft_create,ptcn_fft_apply,ptcn_fft_destroy
 type ptcn_fft_plan
  integer::grid(3)=0,ng=0,no=0
  type(c_ptr)::forward=c_null_ptr,backward=c_null_ptr
  complex(c_double_complex),allocatable::buffer(:,:)
 end type
contains
 subroutine ptcn_fft_create(p,grid,no,err)
  type(ptcn_fft_plan),intent(inout)::p
  integer,intent(in)::grid(3),no
  integer,intent(out)::err
  integer::dims(3),s
  integer(c_int64_t)::points
  err=1
  call ptcn_fft_destroy(p)
  if(minval(grid)<1.or.maxval(grid)>1000000.or.no<0)return
  points=int(grid(1),c_int64_t)*grid(2)*grid(3)
  if(points>huge(1))return
  if(points*int(no,c_int64_t)>huge(1))return
  p%grid=grid;p%ng=int(points);p%no=no
  allocate(p%buffer(p%ng,no),stat=s);if(s/=0)return
  if(no==0)then
   err=0;return
  endif
  ! FFTW C dimension order is reversed relative to Fortran x-fast storage.
  dims=grid(3:1:-1)
  p%forward=fftw_plan_many_dft(3,dims,no,p%buffer,dims,1,p%ng,p%buffer,dims,1,p%ng, &
    FFTW_FORWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
  p%backward=fftw_plan_many_dft(3,dims,no,p%buffer,dims,1,p%ng,p%buffer,dims,1,p%ng, &
    FFTW_BACKWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
  if(.not.c_associated(p%forward).or..not.c_associated(p%backward))then
   call ptcn_fft_destroy(p);return
  endif
  err=0
 end subroutine
 subroutine ptcn_fft_apply(p,u,denominator,v,err)
  type(ptcn_fft_plan),intent(inout)::p
  complex(8),intent(in)::u(:,:),denominator(:)
  complex(8),intent(out)::v(:,:)
  integer,intent(out)::err
  integer::i,j
  err=1
  if(.not.allocated(p%buffer))return
  if(any(shape(u)/=[p%ng,p%no]).or.any(shape(v)/=shape(u)).or.size(denominator)/=p%ng)return
  do i=1,p%ng
   if(.not.ieee_is_finite(real(denominator(i),8)).or..not.ieee_is_finite(aimag(denominator(i))))return
   if(abs(denominator(i))<=tiny(1d0))return
  enddo
  do j=1,p%no;do i=1,p%ng
   if(.not.ieee_is_finite(real(u(i,j),8)).or..not.ieee_is_finite(aimag(u(i,j))))return
  enddo;enddo
  if(p%no==0)then
   err=0;return
  endif
  p%buffer=u
  call fftw_execute_dft(p%forward,p%buffer,p%buffer)
  do j=1,p%no;do i=1,p%ng
   p%buffer(i,j)=p%buffer(i,j)/denominator(i)
  enddo;enddo
  call fftw_execute_dft(p%backward,p%buffer,p%buffer)
  p%buffer=p%buffer/real(p%ng,8)
  do j=1,p%no;do i=1,p%ng
   if(.not.ieee_is_finite(real(p%buffer(i,j),8)).or..not.ieee_is_finite(aimag(p%buffer(i,j))))return
  enddo;enddo
  v=p%buffer;err=0
 end subroutine
 subroutine ptcn_fft_destroy(p)
  type(ptcn_fft_plan),intent(inout)::p
  if(c_associated(p%forward))call fftw_destroy_plan(p%forward)
  if(c_associated(p%backward))call fftw_destroy_plan(p%backward)
  p%forward=c_null_ptr;p%backward=c_null_ptr
  if(allocated(p%buffer))deallocate(p%buffer)
  p%grid=0;p%ng=0;p%no=0
 end subroutine
end module
