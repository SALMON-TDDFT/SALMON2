program test_response_symmetry
  use hsex_response_symmetry, only: s_response_stars,response_q_stars
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  implicit none
  integer,parameter :: ng=64,nk=8,nv=2,nq=3
  type(s_response_stars) :: stars
  complex(8) :: occupied(ng,nv,nk)
  real(8) :: cell(3,3),rec(3,3),frac(3,nk),levels(nv,nk),q(3,nq),atoms(3,1),r(3)
  integer :: g(3,1,nq),counts(nq),sources(nq),i,j,code,species(1),sites(3,nk)
  logical :: broken_h
  integer :: invalid_h
  complex(8) :: saved
  cell=0d0;rec=0d0;broken_h=.false.;invalid_h=0
  do i=1,3
    cell(i,i)=8d0;rec(i,i)=2d0*acos(-1d0)/8d0
  enddo
  do j=0,nk-1
    sites(:,j+1)=[mod(j,2),mod(j/2,2),j/4]
    frac(:,j+1)=-.25d0+.5d0*sites(:,j+1)
  enddo
  do i=0,ng-1
    r=dble([mod(i,4),mod(i/4,4),i/16])/4d0
    occupied(i+1,1,:)=1d0/8d0
    occupied(i+1,2,:)=exp(cmplx(0d0,2d0*acos(-1d0)*sum(r),8))/8d0
  enddo
  ! Equal occupied energies permit arbitrary complex unitary orbital gauges.
  do j=2,nk
    occupied(:,1,j)=(occupied(:,1,j)+(0d0,1d0)*occupied(:,2,j))/sqrt(2d0)
    occupied(:,2,j)=(0d0,1d0)*occupied(:,1,1)/sqrt(2d0)+occupied(:,2,1)/sqrt(2d0)
  enddo
  levels=-1d0;atoms=0d0;species=1;g=0;counts=1;q=0d0
  q(1,1)=.5d0;q(2,2)=.5d0;q(3,3)=.5d0;sources=[2,3,5]
  call response_q_stars([4,4,4],[2,2,2],[2,2,2],cell,rec,atoms,species,frac,levels,occupied, &
    1,q,g,counts,sources,apply,stars,code)
  if(code/=0.or.stars%nsym/=6.or.stars%nrep/=1)error stop 'axis permutation symmetry not found'
  if(any(stars%representative/=1))error stop 'q star representative'
  print *,'RESPONSE_SYMMETRY_GAUGE_STAR_PASS',stars%nsym,stars%nrep
  atoms(:,1)=[.13d0,.21d0,.37d0]
  call response_q_stars([4,4,4],[2,2,2],[2,2,2],cell,rec,atoms,species,frac,levels,occupied, &
    1,q,g,counts,sources,apply,stars,code)
  if(code/=0.or.stars%nsym/=1.or.stars%nrep/=nq)error stop 'broken atom symmetry retained'
  print *,'RESPONSE_SYMMETRY_BROKEN_ATOMS_PASS'
  atoms=0d0;broken_h=.true.
  call response_q_stars([4,4,4],[2,2,2],[2,2,2],cell,rec,atoms,species,frac,levels,occupied, &
    1,q,g,counts,sources,apply,stars,code)
  if(code/=0.or.stars%nsym/=2.or.stars%nrep/=2)error stop 'broken H symmetry retained'
  print *,'RESPONSE_SYMMETRY_BROKEN_H_PASS'
  broken_h=.false.;levels(1,3)=-.8d0
  call response_q_stars([4,4,4],[2,2,2],[2,2,2],cell,rec,atoms,species,frac,levels,occupied, &
    1,q,g,counts,sources,apply,stars,code)
  if(code/=0.or.stars%nsym/=6.or.stars%nrep/=2)error stop 'split source energies retained'
  print *,'RESPONSE_SYMMETRY_SPLIT_ENERGY_PASS'
  levels=-1d0
  do i=0,ng-1
    occupied(i+1,2,3)=exp(cmplx(0d0,2d0*acos(-1d0)*mod(i,4)/4d0,8))/8d0
  enddo
  call response_q_stars([4,4,4],[2,2,2],[2,2,2],cell,rec,atoms,species,frac,levels,occupied, &
    1,q,g,counts,sources,apply,stars,code)
  if(code/=0.or.stars%nsym/=6.or.stars%nrep/=2)error stop 'different source projector retained'
  print *,'RESPONSE_SYMMETRY_SOURCE_SPAN_PASS'
  occupied(:,:,3)=occupied(:,:,2)
  ! Equivalent reciprocal representations exercise nonzero target folding
  ! and q_delta, with a corresponding nonzero G relabeling.
  do j=1,nk
    frac(1,j)=frac(1,j)+1d0
    do i=0,ng-1
      occupied(i+1,:,j)=occupied(i+1,:,j)*exp(cmplx(0d0,-2d0*acos(-1d0)*mod(i,4)/4d0,8))
    enddo
  enddo
  q(2,2)=-.5d0;g(:,1,2)=[0,1,0]
  call response_q_stars([4,4,4],[2,2,2],[2,2,2],cell,rec,atoms,species,frac,levels,occupied, &
    1,q,g,counts,sources,apply,stars,code)
  if(code/=0.or.stars%nsym/=6.or.stars%nrep/=1)error stop 'target fold or reciprocal G shift failed'
  print *,'RESPONSE_SYMMETRY_TARGET_FOLD_G_SHIFT_PASS'
  saved=occupied(ng,nv,nk)
  do j=1,2
    if(j==1)occupied(ng,nv,nk)=ieee_value(0d0,ieee_quiet_nan)
    if(j==2)occupied(ng,nv,nk)=ieee_value(0d0,ieee_positive_inf)
    call response_q_stars([4,4,4],[2,2,2],[2,2,2],cell,rec,atoms,species,frac,levels,occupied, &
      1,q,g,counts,sources,apply,stars,code)
    if(code==0)error stop 'nonfinite symmetry occupied state accepted'
  enddo
  occupied(ng,nv,nk)=saved
  do invalid_h=1,2
    call response_q_stars([4,4,4],[2,2,2],[2,2,2],cell,rec,atoms,species,frac,levels,occupied, &
      1,q,g,counts,sources,apply,stars,code)
    if(code==0)error stop 'nonfinite symmetry H action accepted'
  enddo
  print *,'RESPONSE_SYMMETRY_IEEE_REJECTION_PASS'
contains
  subroutine apply(x,y,status)
    implicit none
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: status
    integer :: point
    y=2d0*x;status=0
    if(invalid_h==1)y(size(x))=ieee_value(0d0,ieee_quiet_nan)
    if(invalid_h==2)y(size(x))=ieee_value(0d0,ieee_positive_inf)
    if(broken_h)then
      do point=1,size(x)
        y(point)=y(point)+.1d0*mod(point-1,4)*x(point)
      enddo
    endif
  end subroutine
end program
