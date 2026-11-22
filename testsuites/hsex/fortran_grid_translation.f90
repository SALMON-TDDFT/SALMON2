program test_grid_translation
  use hsex_screen_math, only: screen_shift
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  implicit none
  complex(8),allocatable :: psi(:),v(:),p(:),vp(:),product(:),shifted_product(:)
  real(8) :: shift(3),x,error,energy_change,pi
  integer :: n(3),ng,nx,g,status,site(3)
  complex(8) :: phase
  pi=acos(-1d0)
  do nx=8,16,8
    n=[nx,2,2];ng=product_size(n)
    allocate(psi(ng),v(ng),p(ng),vp(ng),product(ng),shifted_product(ng))
    do g=0,ng-1
      x=dble(mod(g,nx))/nx
      psi(g+1)=cos(6d0*pi*x);v(g+1)=cos(4d0*pi*x)
    enddo
    psi=psi/sqrt(sum(abs(psi)**2));product=v*psi
    shift=[1d0/nx,0d0,0d0]
    call check(error,energy_change)
    if(error>1d-12.or.abs(energy_change)>1d-12)error stop 'integer translation covariance'
    shift=[.5d0/nx,0d0,0d0]
    call check(error,energy_change)
    if(abs(sum(abs(p)**2)-1d0)>1d-12)error stop 'translation is not unitary'
    if(nx==8)then
      if(error<.1d0.or.abs(energy_change)<.1d0)error stop 'collocation pinning not detected'
    else
      if(error>1d-12.or.abs(energy_change)>1d-12)error stop 'refined band-limited covariance'
    endif
    print '(a,i4,2es18.9)','GRID_TRANSLATION Nx/defect/pinning_energy ',nx,error,energy_change
    shift(1)=ieee_value(0d0,ieee_quiet_nan)
    call screen_shift(n,shift,psi,p,status)
    if(status==0)error stop 'NaN translation accepted'
    shift=0d0;psi(1)=cmplx(ieee_value(0d0,ieee_positive_inf),0d0,8)
    call screen_shift(n,shift,psi,p,status)
    if(status==0)error stop 'infinite translation input accepted'
    deallocate(psi,v,p,vp,product,shifted_product)
  enddo
  n=[5,6,7];ng=product_size(n);allocate(psi(ng),p(ng))
  do g=0,ng-1
    site=[mod(g,n(1)),mod(g/n(1),n(2)),g/(n(1)*n(2))]
    psi(g+1)=exp(cmplx(0d0,2d0*pi*sum(dble(site)*[1d0,-2d0,3d0]/n),8))
  enddo
  shift=[.13d0,.27d0,-.19d0]
  phase=exp(cmplx(0d0,-2d0*pi*sum(shift*[1d0,-2d0,3d0]),8))
  call screen_shift(n,shift,psi,p,status)
  if(status/=0.or.maxval(abs(p-phase*psi))>1d-12)error stop '3D non-power-two translation'
  call screen_shift(n,shift,psi,p,status,real_field=.true.)
  if(status==0)error stop 'complex real field accepted'
  deallocate(psi,p)
  n=[4,1,1];allocate(psi(4),p(4));psi=[1d0,-1d0,1d0,-1d0]
  shift=[.125d0,0d0,0d0]
  call screen_shift(n,shift,psi,p,status,real_field=.true.)
  if(status/=0.or.maxval(abs(p))>1d-12)error stop 'real Nyquist cosine convention'
  call screen_shift(n,shift,psi,p,status)
  if(status/=0.or.abs(sum(abs(p)**2)-4d0)>1d-12)error stop 'complex Nyquist norm'
  print *, 'GRID_TRANSLATION_COVARIANCE_PINNING_REFINEMENT_IEEE_PASS'
contains
  integer function product_size(n)
    implicit none
    integer,intent(in) :: n(3)
    product_size=n(1)*n(2)*n(3)
  end function
  subroutine check(error,energy)
    implicit none
    real(8),intent(out) :: error,energy
    call screen_shift(n,shift,psi,p,status)
    if(status/=0)error stop 'occupied translation rejected'
    call screen_shift(n,shift,v,vp,status,real_field=.true.)
    if(status/=0)error stop 'real potential translation rejected'
    call screen_shift(n,shift,product,shifted_product,status)
    if(status/=0)error stop 'product translation rejected'
    error=sqrt(sum(abs(vp*p-shifted_product)**2))
    energy=real(dot_product(p,vp*p)-dot_product(psi,product),8)
  end subroutine
end program
