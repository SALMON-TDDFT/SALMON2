! Discrete XC variation with prescribed valence/core density; no GS data.
program test_nlcc_xc_variation
  use salmon_xc, only: init_xc, finalize_xc, calc_xc
  use structures, only: s_xc_functional, s_pp_info
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer,parameter :: n=8
  real(8),parameter :: dx=0.7d0,dv=dx**3,eps=1d-4,pi=acos(-1d0)
  type(s_xc_functional) :: xc
  type(s_pp_info) :: pp
  real(8) :: rho(n,n,n),core(n,n,n),direction(n,n,n),grad(n,n,n,3)
  real(8) :: e(n,n,n),v(n,n,n),rd(n,n,n,3),vp(n,n,n),tmp(n,n,n)
  real(8) :: eplus,eminus,fd,inner,err
  character(16) :: names(3)=[character(16)::'pz','libxc_pz','libxc_pbe']
  integer :: i,j,k,case_,nc,failures
  failures=0
  do nc=0,1
    pp%flag_nlcc=nc==1
    do k=1,n
    do j=1,n
    do i=1,n
      rho(i,j,k)=.3d0+.025d0*cos(2*pi*i/n)+.018d0*sin(2*pi*j/n)
      core(i,j,k)=nc*(.07d0+.015d0*cos(2*pi*(j+k)/n))
      direction(i,j,k)=.02d0+.007d0*cos(2*pi*(i+j)/n)
    end do
    end do
    end do
    do case_=1,3
      call init_xc(xc,'unpolarized',0d0,trim(names(case_)),'none','none')
      call evaluate(rho,e,v,rd)
      vp=v
      do i=1,3
        vp=vp+(cshift(rd(:,:,:,i),1,i)-cshift(rd(:,:,:,i),-1,i))/(2*dx)
      end do
      inner=sum(vp*direction)*dv
      call evaluate(rho+eps*direction,e,v,rd)
      eplus=sum(e)*dv
      call evaluate(rho-eps*direction,e,v,rd)
      eminus=sum(e)*dv
      fd=(eplus-eminus)/(2*eps)
      err=abs(fd-inner)
      write(*,'(a,1x,a,1x,i1,3(1x,es24.16))') 'XC_VARIATION',trim(names(case_)),nc,fd,inner,err
      if(.not.ieee_is_finite(fd) .or. .not.ieee_is_finite(inner)) failures=failures+1
      if(err>1d-9*max(1d0,abs(inner))) failures=failures+1
      call finalize_xc(xc)
    end do
  end do
  if(failures/=0) stop 1
  write(*,*) 'NLCC_XC_VARIATION_PASS'
contains
  subroutine evaluate(density,energy,potential,gradient_term)
    real(8),intent(in) :: density(n,n,n)
    real(8),intent(out) :: energy(n,n,n),potential(n,n,n),gradient_term(n,n,n,3)
    integer :: a
    tmp=density+core
    do a=1,3
      grad(:,:,:,a)=(cshift(tmp,1,a)-cshift(tmp,-1,a))/(2*dx)
    end do
    call calc_xc(xc,pp,rho=density,rho_nlcc=core,grho=grad,eexc=energy,vxc=potential,rdedd=gradient_term)
  end subroutine
end program
