! Native Poisson FFTE layout: x-complete pencils, distributed over y and z.
! x-only grid decomposition replicates pencils; use y/z decomposition to scale.
module rvv10_distributed
  use communication, only: comm_summation,comm_get_max
  use rvv10, only: rvv10_evaluate,rvv10_kernel_fourier
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: rvv10_evaluate_distributed
contains
  subroutine rvv10_evaluate_distributed(n,lo,m,dims,coords,comm,comm_r,h,rho,sigma,b,c,nq, &
      energy,vrho,vsigma,used,status)
    integer,intent(in) :: n(3),lo(3),m(3),dims(3),coords(3),comm(3),comm_r,nq
    real(8),intent(in) :: h(3),rho(:),sigma(:),b,c
    real(8),intent(out) :: energy(:),vrho(:),vsigma(:)
    logical,intent(out) :: used
    integer,intent(out) :: status
    integer :: bad,i,j,t,tile(3),nt
    real(8),allocatable :: part(:,:,:),whole(:,:,:),r(:),s(:),e(:),v(:),w(:)
    used=.false.;status=0;bad=0
    ! FFTE fixed workspace and transpose divisibility constraints. Fall back
    ! collectively before entering any axis-communicator operation.
    if(any(n<2).or.any(n>4096).or.any(dims<1).or.any(dims>65536))bad=1
    if(bad==0)then
      if(any(modulo(n,dims)/=0))bad=1
      if(modulo(n(1),dims(2))/=0.or.modulo(n(2),dims(3))/=0)bad=1
      if(any(m/=n/dims).or.any(lo/=coords*m+1))bad=1
      if(any(coords<0).or.any(coords>=dims))bad=1
      do i=1,3
        t=n(i)
        do j=2,5
          if(j==4)cycle
          do while(modulo(t,j)==0)
            t=t/j
          enddo
        enddo
        if(t/=1)bad=1
      enddo
    endif
    call comm_get_max(bad,comm_r)
    if(bad/=0)return
    used=.true.;status=1;bad=0
    if(size(rho)/=product(m).or.size(sigma)/=product(m).or.size(energy)/=product(m).or. &
      size(vrho)/=product(m).or.size(vsigma)/=product(m))bad=1
    if(any(rho<0d0).or.any(sigma<0d0).or..not.all(ieee_is_finite(rho)).or. &
       .not.all(ieee_is_finite(sigma)))bad=1
    if(any(h<=0d0).or..not.all(ieee_is_finite(h)).or.nq<8.or.nq>128)bad=1
    if(.not.ieee_is_finite(b).or..not.ieee_is_finite(c).or.b<=0d0.or.c<0d0)bad=1
    call comm_get_max(bad,comm_r)
    if(bad/=0)return
    tile=[n(1),m(2),m(3)];nt=product(tile)
    allocate(part(tile(1),tile(2),tile(3)),whole(tile(1),tile(2),tile(3)))
    allocate(r(nt),s(nt),e(nt),v(nt),w(nt))
    part=0d0;part(lo(1):lo(1)+m(1)-1,:,:)=reshape(rho,m)
    call comm_summation(part,whole,nt,comm(1));r=reshape(whole,[nt])
    part=0d0;part(lo(1):lo(1)+m(1)-1,:,:)=reshape(sigma,m)
    call comm_summation(part,whole,nt,comm(1));s=reshape(whole,[nt])
    call rvv10_evaluate(tile,h,r,s,b,c,nq,e,v,w,status,convolve)
    call comm_get_max(status,comm_r)
    if(status/=0)return
    whole=reshape(e,tile);energy=reshape(whole(lo(1):lo(1)+m(1)-1,:,:),[product(m)])
    whole=reshape(v,tile);vrho=reshape(whole(lo(1):lo(1)+m(1)-1,:,:),[product(m)])
    whole=reshape(w,tile);vsigma=reshape(whole(lo(1):lo(1)+m(1)-1,:,:),[product(m)])
  contains
    subroutine convolve(theta,u,mesh,ierr)
      complex(8),intent(in) :: theta(:,:)
      complex(8),intent(out) :: u(:,:)
      real(8),intent(in) :: mesh(:)
      integer,intent(out) :: ierr
      complex(8),allocatable :: transformed(:,:),a(:),bb(:)
      integer :: q,d,x,y,z,p(3),gidx
      real(8) :: g,phi,pi
      allocate(transformed(nt,nq),a(nt),bb(nt));pi=acos(-1d0)
      do q=1,nq
        a=theta(:,q)
        call pzfft3dv_rvv10(a,bb,n(1),n(2),n(3),dims(2),dims(3),-1,comm(2),comm(3))
        transformed(:,q)=bb
      enddo
      u=0d0
!$omp parallel do collapse(3) private(x,y,z,p,gidx,g,q,d,phi)
      do z=0,tile(3)-1;do y=0,tile(2)-1;do x=0,tile(1)-1
        p=[x,y+lo(2)-1,z+lo(3)-1];where(p>=(n+1)/2)p=p-n
        gidx=1+x+tile(1)*(y+tile(2)*z);g=sqrt(sum((2*pi*p/(n*h))**2))
        do q=1,nq;do d=1,q
          phi=rvv10_kernel_fourier(mesh(q),mesh(d),g)
          u(gidx,q)=u(gidx,q)+phi*transformed(gidx,d)
          if(q/=d)u(gidx,d)=u(gidx,d)+phi*transformed(gidx,q)
        enddo;enddo
      enddo;enddo;enddo
!$omp end parallel do
      do q=1,nq
        a=u(:,q)
        ! Native FFTE inverse already divides by the GLOBAL number of points.
        call pzfft3dv_rvv10(a,bb,n(1),n(2),n(3),dims(2),dims(3),1,comm(2),comm(3))
        u(:,q)=bb
      enddo
      ierr=0
    end subroutine
  end subroutine
end module
