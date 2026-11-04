! Atomic units. Columns of step are the Cartesian real-space grid increments.
! Recenter each orbital in fractional BvK grid coordinates before zero padding.
module kli_self_hartree
  use iso_c_binding
  use ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  include 'fftw3.f03'
  public :: isolated_box,prepare_box,self_hartree,release_box
  type isolated_box
    integer :: n(3)=0, m(3)=0
    real(8) :: dv=0d0,rc=0d0
    complex(c_double_complex),allocatable :: buf(:,:,:)
    real(8),allocatable :: kernel(:,:,:)
    type(c_ptr) :: fw=c_null_ptr,bw=c_null_ptr
  end type
contains
  pure function cross(a,b) result(c)
    real(8),intent(in)::a(3),b(3)
    real(8)::c(3)
    c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
  end function
  subroutine release_box(box)
    type(isolated_box),intent(inout)::box
    if(c_associated(box%fw))call fftw_destroy_plan(box%fw)
    if(c_associated(box%bw))call fftw_destroy_plan(box%bw)
    box%fw=c_null_ptr;box%bw=c_null_ptr
    if(allocated(box%buf))deallocate(box%buf)
    if(allocated(box%kernel))deallocate(box%kernel)
    box%n=0
  end subroutine
  subroutine prepare_box(box,n,step,status,scale)
    type(isolated_box),intent(inout)::box
    integer,intent(in)::n(3)
    real(8),intent(in)::step(3,3)
    real(8),intent(in),optional::scale
    integer,intent(out)::status
    integer::i,j,x,y,z,s(3),m(3)
    real(8)::dual(3,3),det,pi,q(3),q2,rc
    status=1
    call release_box(box)
    if(any(n<1))return
    do i=1,3
      do j=1,3
        if(.not.ieee_is_finite(step(i,j)))return
      enddo
      dual(:,i)=cross(step(:,mod(i,3)+1),step(:,mod(i+1,3)+1))
    enddo
    det=dot_product(step(:,1),dual(:,1))
    if(abs(det)<tiny(det))return
    pi=acos(-1d0);dual=dual/det
    ! Radius exceeds every pair separation in the unwrapped support box.
    rc=0d0
    do x=-1,1,2;do y=-1,1,2;do z=-1,1,2
      rc=max(rc,norm2(matmul(step,real((n-1)*[x,y,z],8))))
    enddo;enddo;enddo
    if(rc<=0d0)return
    ! Each padded face separation exceeds twice the cutoff, excluding images.
    do i=1,3
      m(i)=max(n(i),ceiling(2*rc*norm2(dual(:,i)))+1)
      if(present(scale))m(i)=ceiling(scale*m(i))
      dual(:,i)=dual(:,i)*(2*pi/m(i))
    enddo
    box%n=n;box%m=m;box%dv=abs(det);box%rc=rc
    allocate(box%buf(0:m(1)-1,0:m(2)-1,0:m(3)-1),box%kernel(0:m(1)-1,0:m(2)-1,0:m(3)-1))
    do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
      s=[x,y,z];where(2*s>=m)s=s-m
      q=matmul(dual,real(s,8));q2=sum(q*q)
      if(q2==0d0)then
        box%kernel(x,y,z)=2*pi*rc**2
      else
        box%kernel(x,y,z)=8*pi*sin(.5d0*rc*sqrt(q2))**2/q2
      endif
    enddo;enddo;enddo
    box%fw=fftw_plan_dft_3d(m(3),m(2),m(1),box%buf,box%buf,FFTW_FORWARD,FFTW_ESTIMATE)
    box%bw=fftw_plan_dft_3d(m(3),m(2),m(1),box%buf,box%buf,FFTW_BACKWARD,FFTW_ESTIMATE)
    if(.not.c_associated(box%fw).or..not.c_associated(box%bw))then
      call release_box(box)
      return
    endif
    status=0
  end subroutine
  subroutine self_hartree(box,rho,v,eh,boundary,status,frame_shift,track)
    type(isolated_box),intent(inout)::box
    real(8),intent(in)::rho(0:,0:,0:)
    real(8),intent(out)::v(0:,0:,0:),eh,boundary
    integer,intent(out)::status
    integer,optional,intent(inout)::frame_shift(3)
    logical,optional,intent(in)::track
    integer::x,y,z,n(3),s(3),shift(3),i
    complex(8)::moment(3)
    real(8)::twopi,chg
    status=1;n=box%n
    if(any(shape(rho)/=n).or.any(shape(v)/=n).or.any(n==0))return
    twopi=2*acos(-1d0);moment=0d0
    do z=0,n(3)-1;do y=0,n(2)-1;do x=0,n(1)-1
      if(.not.ieee_is_finite(rho(x,y,z)).or.rho(x,y,z)<0d0)return
      s=[x,y,z]
      do i=1,3
        moment(i)=moment(i)+rho(x,y,z)*exp(cmplx(0d0,twopi*s(i)/n(i),8))
      enddo
    enddo;enddo;enddo
    chg=sum(rho)*box%dv
    if(abs(chg-1d0)>1d-7)return
    shift=n/2-nint(atan2(aimag(moment),real(moment,8))*n/twopi)
    ! Keep the chosen image cut fixed throughout a variational minimization.
    if(present(frame_shift))then
      if(any(frame_shift==huge(1)))frame_shift=shift
      if(present(track))then
        if(track)then
          ! A one-grid-point hysteresis prevents cut chatter at half-grid centers.
          where(abs(modulo(shift-frame_shift+n/2,n)-n/2)>1)frame_shift=shift
        endif
      endif
      shift=frame_shift
    endif
    box%buf=0d0;boundary=0d0
    do z=0,n(3)-1;do y=0,n(2)-1;do x=0,n(1)-1
      s=modulo([x,y,z]+shift,n)
      box%buf(s(1),s(2),s(3))=rho(x,y,z)
      if(any(s==0).or.any(s==n-1))boundary=boundary+rho(x,y,z)*box%dv
    enddo;enddo;enddo
    call fftw_execute_dft(box%fw,box%buf,box%buf)
    box%buf=box%buf*box%kernel
    call fftw_execute_dft(box%bw,box%buf,box%buf)
    do z=0,n(3)-1;do y=0,n(2)-1;do x=0,n(1)-1
      s=modulo([x,y,z]+shift,n)
      v(x,y,z)=real(box%buf(s(1),s(2),s(3)),8)/product(box%m)
    enddo;enddo;enddo
    eh=.5d0*sum(v*rho)*box%dv
    if(.not.ieee_is_finite(eh).or.eh<=0d0)return
    status=0
  end subroutine
end module
