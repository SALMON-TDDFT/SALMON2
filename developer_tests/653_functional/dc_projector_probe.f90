program probe
  use dc_projector_force, only: core_projector_force
  implicit none
  complex(8) :: v(27),dv(3,27),vp(27),vm(27),ignored(3,27)
  real(8) :: r(3),rp(3),rm(3),f(3),fd,ep,em,h,w
  logical :: core(27)
  integer :: channel,axis,point,mode
  r=[.21d0,-.17d0,.08d0];h=1d-5
  do mode=1,3
    do point=1,27
      core(point)=(modulo(point,3)/=0)
    enddo
    if(mode==2)core=.true.
    if(mode==3)core=.false.
    do channel=0,2
      w=(-1d0)**channel*.37d0
      call sample(r,channel,v,dv)
      call core_projector_force(v,dv,core,w,f)
      do axis=1,3
        rp=r;rm=r;rp(axis)=rp(axis)+h;rm(axis)=rm(axis)-h
        call sample(rp,channel,vp,ignored);call sample(rm,channel,vm,ignored)
        ep=w*real(conjg(sum(vp,mask=core))*sum(vp),8)
        em=w*real(conjg(sum(vm,mask=core))*sum(vm),8)
        fd=-(ep-em)/(2*h)
        if(abs(fd-f(axis))>1d-7)error stop 'truncated projector derivative mismatch'
      enddo
    enddo
  enddo
contains
  subroutine sample(center,l,value,derivative)
    real(8),intent(in) :: center(3)
    integer,intent(in) :: l
    complex(8),intent(out) :: value(27),derivative(3,27)
    real(8) :: x(3),radial,harmonic,dh(3),k(3)
    complex(8) :: psi,phase
    integer :: ix,iy,iz,i
    k=[.13d0,.29d0,-.11d0];i=0
    do iz=-1,1;do iy=-1,1;do ix=-1,1
      i=i+1;x=[real(ix,8),real(iy,8),real(iz,8)]-center
      radial=exp(-sum(x*x)*.4d0);dh=0d0
      select case(l)
      case(0);harmonic=1d0
      case(1);harmonic=x(1);dh(1)=1d0
      case(2);harmonic=x(1)*x(2);dh=[x(2),x(1),0d0]
      end select
      psi=cmplx(sin(real(i,8)),cos(.7d0*i),8)
      phase=exp(cmplx(0d0,sum(k*x),8))
      value(i)=radial*harmonic*phase*psi
      derivative(:,i)=(radial*(.8d0*x*harmonic-dh)-cmplx(0d0,1d0,8)*k*radial*harmonic)*phase*psi
    enddo;enddo;enddo
  end subroutine
end program
