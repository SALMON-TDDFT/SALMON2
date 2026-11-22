program test_screen_forms
  use hsex_screen_math, only: screen_forms,screen_finite
  implicit none
  integer :: n(3),boxes,trial,i,j,case_index
  integer(8) :: started,finished,rate
  complex(8),allocatable :: actual(:,:),expected(:,:)
  real(8) :: q(3,64),error,fast(3),slow(3),tmp
  integer,parameter :: choices(5)=[1,3,4,6,12]
  character(32) :: argument
  do j=1,64
    do i=1,3
      q(i,j)=3d0*sin(dble(7*i+3*j))
    enddo
  enddo
  q(:,1)=0d0;q(:,2)=[6d0,-12d0,0d0]
  do case_index=1,size(choices)
    boxes=choices(case_index);n=[12,24,12]
    allocate(actual(boxes**3,64),expected(boxes**3,64))
    call reference_forms(n,boxes,q,expected)
    call screen_forms(n,boxes,q,actual)
    if(.not.screen_finite(actual).or..not.screen_finite(expected))error stop 'nonfinite box forms'
    error=maxval(abs(actual-expected))
    if(error>2d-14)error stop 'cached box forms differ from the direct finite-grid sums'
    print *, 'SCREEN_FORMS_MAX_ERROR',boxes,error
    deallocate(actual,expected)
  enddo
  call get_command_argument(1,argument)
  if(argument=='--benchmark')then
    boxes=6;n=[12,12,12]
    allocate(actual(boxes**3,64),expected(boxes**3,64))
    do trial=1,3
      call system_clock(started,rate)
      do i=1,100
        call reference_forms(n,boxes,q,expected)
      enddo
      call system_clock(finished)
      slow(trial)=dble(finished-started)/rate
      call system_clock(started)
      do i=1,100
        call screen_forms(n,boxes,q,actual)
      enddo
      call system_clock(finished)
      fast(trial)=dble(finished-started)/rate
      if(.not.screen_finite(actual).or..not.screen_finite(expected))error stop 'nonfinite benchmark forms'
      if(maxval(abs(actual-expected))>2d-14)error stop 'benchmark forms mismatch'
    enddo
    do i=1,2
      do j=i+1,3
        if(slow(j)<slow(i))then
          tmp=slow(i);slow(i)=slow(j);slow(j)=tmp
        endif
        if(fast(j)<fast(i))then
          tmp=fast(i);fast(i)=fast(j);fast(j)=tmp
        endif
      enddo
    enddo
    print *, 'SCREEN_FORMS_MEDIAN reference/cached/speedup',slow(2),fast(2),slow(2)/fast(2)
    if(slow(2)<2d0*fast(2))error stop 'box-form reuse did not reach the local performance target'
  endif
  print *, 'SCREEN_FORMS_PASS'
contains
  subroutine reference_forms(n,boxes,q,forms)
    implicit none
    integer,intent(in) :: n(3),boxes
    real(8),intent(in) :: q(:,:)
    complex(8),intent(out) :: forms(boxes**3,size(q,2))
    integer :: i,g,a,b,j,m
    complex(8) :: factor
    real(8) :: angle
    forms=1d0
    do g=1,size(q,2)
      do i=0,boxes**3-1
        do a=1,3
          b=mod(i/boxes**(a-1),boxes);m=n(a)/boxes
          factor=0d0
          do j=b*m,(b+1)*m-1
            angle=-2d0*acos(-1d0)*j*q(a,g)/n(a)
            factor=factor+cmplx(cos(angle),sin(angle),8)
          enddo
          forms(i+1,g)=forms(i+1,g)*factor/m
        enddo
      enddo
    enddo
  end subroutine reference_forms
end program test_screen_forms
