! Full-support Gamma exchange on x-complete y/z pencils.
! All band matrices are replicated, but orbital grid rows and FFT work are local.
! Collective contract: n/h/dims/band counts/radius/omega/maxiter and call order agree;
! coords and local grid rows vary. Communicators follow spatial coordinate order.
module hse_spatial
  use communication, only: comm_summation,comm_get_max
  use fftw_pencils, only: pencil_transform
  use hse_wannier_gauge, only: gauge_transport,gauge_minimize
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: spatial_exx_state,spatial_exx_refresh,spatial_exx_apply
  type spatial_exx_state
    integer :: updates=0,iterations=0,localization_status=1
    real(8) :: spread=0d0,gradient=0d0,min_singular=0d0
    complex(8),allocatable :: gauge(:,:,:),previous(:,:,:),source(:,:)
  end type
contains
  subroutine spatial_exx_refresh(op,n,h,dims,coords,comm,comm_r,psi,maxiter,tolerance,status,occupation)
    type(spatial_exx_state),intent(inout) :: op
    integer,intent(in) :: n(3),dims(2),coords(2),comm(2),comm_r,maxiter
    real(8),intent(in) :: h(3),tolerance
    real(8),intent(in),optional :: occupation(:,:)
    complex(8),intent(in) :: psi(:,:,:)
    integer,intent(out) :: status
    complex(8),allocatable :: raw(:,:,:,:),shifted(:,:),phase(:)
    real(8) :: b(3,6),weights(6),delta,pi,position(3)
    integer :: no,ng,m(3),lo(3),x,y,z,g,j,axis,neighbors(6,1),bad
    no=size(psi,2);ng=size(psi,1);status=1;bad=0
    if(no<1.or.size(psi,3)/=1.or.any(n<1).or.any(dims<1))bad=1
    if(any(h<=0d0).or..not.all(ieee_is_finite(h)))bad=1
    if(.not.all(ieee_is_finite(real(psi))).or..not.all(ieee_is_finite(aimag(psi))))bad=1
    if(maxiter<0.or.tolerance<=0d0.or..not.ieee_is_finite(tolerance))bad=1
    if(present(occupation))then
      if(any(shape(occupation)/=[no,1]))bad=1
      if(any(occupation<0d0).or.any(occupation>2d0).or..not.all(ieee_is_finite(occupation)))bad=1
    endif
    call comm_get_max(bad,comm_r)
    if(bad/=0)return
    m=[n(1),n(2)/dims(1),n(3)/dims(2)];lo=[0,coords(1)*m(2),coords(2)*m(3)]
    if(ng/=product(m).or.any(coords<0).or.any(coords>=dims))bad=1
    if(modulo(n(1),dims(1))/=0.or.modulo(n(2),dims(1))/=0.or. &
       modulo(n(2),dims(2))/=0.or.modulo(n(3),dims(2))/=0)bad=1
    if(allocated(op%gauge))then
      if(any(shape(op%gauge)/=[no,no,1]).or.any(shape(op%previous)/=shape(psi)))bad=1
    endif
    call comm_get_max(bad,comm_r)
    if(bad/=0)return
    if(.not.allocated(op%gauge))then
      allocate(op%gauge(no,no,1));op%gauge=0d0
      do j=1,no
        op%gauge(j,j,1)=1d0
      enddo
    else
      call gauge_transport(psi,op%previous,product(h),op%gauge,op%min_singular,status,sum_grid)
      if(status/=0)then
        op%gauge=0d0
        do j=1,no
          op%gauge(j,j,1)=1d0
        enddo
      endif
    endif
    op%iterations=0;op%localization_status=2;op%spread=-1d0;op%gradient=-1d0
    if(maxiter>0)then
      allocate(raw(no,no,6,1),shifted(ng,no),phase(ng))
      pi=acos(-1d0);b=0d0;neighbors=1
      do axis=1,3
        delta=2*pi/(n(axis)*h(axis));b(axis,axis)=delta;b(axis,axis+3)=-delta
        weights(axis)=1d0/(2*delta**2);weights(axis+3)=weights(axis)
        g=0
        do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
          g=g+1;position=([x,y,z]+lo)*h
          phase(g)=exp(cmplx(0d0,-position(axis)*delta,8))
        enddo;enddo;enddo
        do j=1,no
          shifted(:,j)=psi(:,j,1)*phase
        enddo
        raw(:,:,axis,1)=matmul(conjg(transpose(psi(:,:,1))),shifted)*product(h)
        call sum_grid(raw(:,:,axis,1))
        raw(:,:,axis+3,1)=conjg(transpose(raw(:,:,axis,1)))
      enddo
      call gauge_minimize(op%gauge,raw,neighbors,b,weights,maxiter,tolerance,op%spread,op%gradient, &
        op%iterations,op%localization_status)
    endif
    ! A non-converged unitary gauge preserves the full-support exchange operator.
    if(.not.all(ieee_is_finite(real(op%gauge))).or..not.all(ieee_is_finite(aimag(op%gauge))))bad=1
    call comm_get_max(bad,comm_r)
    if(bad/=0)return
    op%source=matmul(psi(:,:,1),op%gauge(:,:,1))
    op%previous=reshape(op%source,[ng,no,1])
    if(present(occupation))then
      if(.not.allocated(shifted))allocate(shifted(ng,no))
      do j=1,no
        shifted(:,j)=psi(:,j,1)*sqrt(occupation(j,1)/2d0)
      enddo
      op%source=matmul(shifted,op%gauge(:,:,1))
    endif
    op%updates=op%updates+1;status=0
  contains
    subroutine sum_grid(a)
      complex(8),intent(inout) :: a(:,:)
      complex(8) :: total(size(a,1),size(a,2))
      call comm_summation(a,total,size(a),comm_r)
      a=total
    end subroutine
  end subroutine

  subroutine spatial_exx_apply(op,n,h,dims,coords,comm,comm_r,radius_input,target,action,status,omega)
    type(spatial_exx_state),intent(in) :: op
    integer,intent(in) :: n(3),dims(2),coords(2),comm(2),comm_r
    real(8),intent(in) :: h(3),radius_input
    real(8),intent(in),optional :: omega
    complex(8),intent(in) :: target(:,:,:)
    complex(8),intent(out) :: action(:,:,:)
    integer,intent(out) :: status
    complex(8),allocatable :: density(:,:),spectrum(:,:),potential(:,:)
    real(8),allocatable :: multiplier(:)
    real(8) :: radius,pi,q(3),q2,screening
    integer :: ng,nt,m(3),lo(3),x,y,z,g,p(3),i,j,first,nb,bad
    status=1;action=0d0;bad=0
    screening=0d0
    if(present(omega))screening=omega
    if(.not.ieee_is_finite(screening).or.screening<0d0)bad=1
    if(.not.allocated(op%source))bad=1
    if(any(n<1).or.any(dims<1).or.any(h<=0d0).or..not.all(ieee_is_finite(h)))bad=1
    if(screening==0d0)then
      if(.not.ieee_is_finite(radius_input).or.radius_input<0d0)bad=1
    endif
    call comm_get_max(bad,comm_r)
    if(bad/=0)return
    m=[n(1),n(2)/dims(1),n(3)/dims(2)];lo=[0,coords(1)*m(2),coords(2)*m(3)]
    ng=product(m);nt=size(target,2)
    if(size(target,1)/=ng.or.size(target,3)/=1.or.nt<1.or.any(shape(action)/=shape(target)))bad=1
    if(size(op%source,1)/=ng)bad=1
    if(.not.all(ieee_is_finite(real(target))).or..not.all(ieee_is_finite(aimag(target))))bad=1
    radius=.5d0*minval(n*h)
    if(radius_input>0d0)radius=radius_input
    if(screening==0d0.and.radius>.5d0*minval(n*h)*(1d0+1d-12))bad=1
    call comm_get_max(bad,comm_r)
    if(bad/=0)return
    allocate(multiplier(ng));pi=acos(-1d0);g=0
    ! X-pencil spectral output has the same local ordering as real-space input.
    do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
      g=g+1;p=[x,y,z]+lo
      where(p>=(n+1)/2)p=p-n
      q=2*pi*p/(n*h);q2=sum(q*q)
      if(screening>0d0)then
        if(q2<1d-24)then
          multiplier(g)=pi/screening**2
        else
          multiplier(g)=4*pi*(1d0-exp(-q2/(4*screening**2)))/q2
        endif
      else if(q2<1d-24)then
        multiplier(g)=2*pi*radius**2
      else
        multiplier(g)=8*pi*sin(.5d0*sqrt(q2)*radius)**2/q2
      endif
    enddo;enddo;enddo
    allocate(density(ng,min(4,nt)),spectrum(ng,min(4,nt)),potential(ng,min(4,nt)))
    do i=1,size(op%source,2)
      do first=1,nt,4
        nb=min(4,nt-first+1)
        do j=1,nb
          density(:,j)=conjg(op%source(:,i))*target(:,first+j-1,1)
        enddo
        call pencil_transform(n,dims,coords,comm,density(:,:nb),spectrum(:,:nb),-1,status)
        if(status/=0)return
        do j=1,nb
          spectrum(:,j)=spectrum(:,j)*multiplier
        enddo
        ! Inverse pencil_transform already includes 1/product(n).
        call pencil_transform(n,dims,coords,comm,spectrum(:,:nb),potential(:,:nb),1,status)
        if(status/=0)return
        do j=1,nb
          action(:,first+j-1,1)=action(:,first+j-1,1)-op%source(:,i)*potential(:,j)
        enddo
      enddo
    enddo
    status=0
  end subroutine
end module
