! Full-support Gamma exchange on x-complete y/z pencils.
! Refresh and apply accept orbital-local columns through optional comm_o.
! Grid rows and FFT work are spatially local.
! Collective contract: n/h/dims/radius/omega/maxiter and call order agree;
! band counts agree within spatial groups, and may differ across orbital groups.
! coords and local grid rows vary. Communicators follow spatial coordinate order.
module hse_spatial
  use communication, only: comm_summation,comm_get_max,comm_bcast,comm_get_groupinfo
  use exx_orbitals, only: orbital_layout,orbital_check,orbital_overlap,orbital_rotate
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
  subroutine spatial_exx_refresh(op,n,h,dims,coords,comm,comm_r,psi,maxiter,tolerance,status,occupation,comm_o)
    type(spatial_exx_state),intent(inout) :: op
    integer,intent(in) :: n(3),dims(2),coords(2),comm(2),comm_r,maxiter
    integer,intent(in),optional :: comm_o
    real(8),intent(in) :: h(3),tolerance
    real(8),intent(in),optional :: occupation(:,:)
    complex(8),intent(in) :: psi(:,:,:)
    integer,intent(out) :: status
    complex(8),allocatable :: raw(:,:,:,:),shifted(:,:),phase(:)
    real(8) :: b(3,6),weights(6),delta,pi,position(3)
    integer :: no,ng,m(3),lo(3),x,y,z,g,j,axis,neighbors(6,1),bad
    if(present(comm_o))then
      call refresh_orbitals(op,n,h,dims,coords,comm_r,comm_o,psi,maxiter,tolerance,status,occupation)
      return
    endif
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

  subroutine refresh_orbitals(op,n,h,dims,coords,comm_r,comm_o,psi,maxiter,tolerance,status,occupation)
    type(spatial_exx_state),intent(inout) :: op
    integer,intent(in) :: n(3),dims(2),coords(2),comm_r,comm_o,maxiter
    real(8),intent(in) :: h(3),tolerance
    complex(8),intent(in) :: psi(:,:,:)
    real(8),intent(in),optional :: occupation(:,:)
    integer,intent(out) :: status
    integer,allocatable :: counts(:)
    complex(8),allocatable :: raw(:,:,:,:),phase(:),overlap(:,:),left(:,:),right(:,:),work(:)
    real(8),allocatable :: singular(:),rwork(:),occupation_weights(:)
    real(8) :: b(3,6),weights(6),delta,pi,position(3)
    integer :: no,ng,nlocal,first,m(3),lo(3),axis,g,x,y,z,j,bad,neighbors(6,1),initialized,total_initialized
    ng=size(psi,1);nlocal=size(psi,2);status=1;bad=0
    if(size(psi,3)/=1.or.any(n<1).or.any(dims<1))bad=1
    if(any(h<=0d0).or..not.all(ieee_is_finite(h)))bad=1
    if(.not.all(ieee_is_finite(real(psi))).or..not.all(ieee_is_finite(aimag(psi))))bad=1
    if(maxiter<0.or.tolerance<=0d0.or..not.ieee_is_finite(tolerance))bad=1
    if(present(occupation))then
      if(any(shape(occupation)/=[nlocal,1]))bad=1
      if(any(occupation<0d0).or.any(occupation>2d0).or..not.all(ieee_is_finite(occupation)))bad=1
    endif
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    m=[n(1),n(2)/dims(1),n(3)/dims(2)];lo=[0,coords(1)*m(2),coords(2)*m(3)]
    if(ng/=product(m).or.any(coords<0).or.any(coords>=dims))bad=1
    if(modulo(n(1),dims(1))/=0.or.modulo(n(2),dims(1))/=0.or. &
       modulo(n(2),dims(2))/=0.or.modulo(n(3),dims(2))/=0)bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    call orbital_layout(nlocal,comm_r,comm_o,counts,first,bad)
    if(bad/=0.or.sum(counts)<1)return
    no=sum(counts)
    initialized=0
    if(allocated(op%gauge))initialized=1
    total_initialized=initialized
    call orbital_check(total_initialized,comm_r,comm_o)
    if(initialized/=total_initialized)bad=1
    if(allocated(op%gauge))then
      if(any(shape(op%gauge)/=[no,no,1]).or..not.allocated(op%previous))then
        bad=1
      else
        if(any(shape(op%previous)/=shape(psi)))bad=1
        if(.not.all(ieee_is_finite(real(op%previous))).or..not.all(ieee_is_finite(aimag(op%previous))))bad=1
      endif
    endif
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    if(.not.allocated(op%gauge))then
      allocate(op%gauge(no,no,1));op%gauge=0d0
      do j=1,no
        op%gauge(j,j,1)=1d0
      enddo
    else
      allocate(overlap(no,no),left(no,no),right(no,no),work(8*no),singular(no),rwork(5*no))
      call orbital_overlap(psi(:,:,1),op%previous(:,:,1),product(h),comm_r,comm_o,counts,first,overlap)
      call zgesvd('A','A',no,no,overlap,no,singular,left,no,right,no,work,size(work),rwork,bad)
      if(bad/=0)bad=1
      call orbital_check(bad,comm_r,comm_o)
      if(bad==0)then
        op%min_singular=minval(singular)
        if(.not.all(ieee_is_finite(singular)).or.op%min_singular<1d-8)bad=1
      endif
      call orbital_check(bad,comm_r,comm_o)
      if(bad==0)then
        op%gauge(:,:,1)=matmul(left,right)
      else
        op%gauge=0d0
        do j=1,no
          op%gauge(j,j,1)=1d0
        enddo
        bad=0
      endif
    endif
    op%iterations=0;op%localization_status=2;op%spread=-1d0;op%gradient=-1d0
    if(maxiter>0)then
      allocate(raw(no,no,6,1),phase(ng))
      pi=acos(-1d0);b=0d0;neighbors=1
      do axis=1,3
        delta=2*pi/(n(axis)*h(axis));b(axis,axis)=delta;b(axis,axis+3)=-delta
        weights(axis)=1d0/(2*delta**2);weights(axis+3)=weights(axis)
        g=0
        do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
          g=g+1;position=([x,y,z]+lo)*h
          phase(g)=exp(cmplx(0d0,-position(axis)*delta,8))
        enddo;enddo;enddo
        call orbital_overlap(psi(:,:,1),psi(:,:,1),product(h),comm_r,comm_o,counts,first,raw(:,:,axis,1),phase)
        raw(:,:,axis+3,1)=conjg(transpose(raw(:,:,axis,1)))
      enddo
      call gauge_minimize(op%gauge,raw,neighbors,b,weights,maxiter,tolerance,op%spread,op%gradient, &
        op%iterations,op%localization_status)
    endif
    if(.not.all(ieee_is_finite(real(op%gauge))).or..not.all(ieee_is_finite(aimag(op%gauge))))bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    if(allocated(op%source))deallocate(op%source)
    if(allocated(op%previous))deallocate(op%previous)
    allocate(op%source(ng,nlocal),op%previous(ng,nlocal,1),occupation_weights(nlocal))
    call orbital_rotate(psi(:,:,1),op%gauge(:,:,1),comm_o,counts,first,op%previous(:,:,1))
    occupation_weights=1d0
    if(present(occupation))occupation_weights=sqrt(occupation(:,1)/2d0)
    call orbital_rotate(psi(:,:,1),op%gauge(:,:,1),comm_o,counts,first,op%source,occupation_weights)
    op%updates=op%updates+1;status=0
  end subroutine

  ! Optional comm_o joins matching grid pencils across orbital groups. Source and
  ! target columns may have different (including zero) local counts. Counts must
  ! agree within comm_r, and the communicators form a spatial/orbital product.
  subroutine spatial_exx_apply(op,n,h,dims,coords,comm,comm_r,radius_input,target,action,status,omega,comm_o)
    type(spatial_exx_state),intent(in) :: op
    integer,intent(in) :: n(3),dims(2),coords(2),comm(2),comm_r
    integer,intent(in),optional :: comm_o
    real(8),intent(in) :: h(3),radius_input
    real(8),intent(in),optional :: omega
    complex(8),intent(in) :: target(:,:,:)
    complex(8),intent(out) :: action(:,:,:)
    integer,intent(out) :: status
    complex(8),allocatable :: density(:,:),spectrum(:,:),source_column(:)
    real(8),allocatable :: multiplier(:)
    real(8) :: radius,pi,q(3),q2,screening
    integer :: ng,nt,m(3),lo(3),x,y,z,g,p(3),i,j,first,nb,bad,owner,orb_rank,orb_size,count,counts_max,nt_max
    status=1;action=0d0;bad=0
    screening=0d0
    if(present(omega))screening=omega
    if(.not.ieee_is_finite(screening).or.screening<0d0)bad=1
    if(.not.allocated(op%source))bad=1
    if(any(n<1).or.any(dims<1).or.any(h<=0d0).or..not.all(ieee_is_finite(h)))bad=1
    if(screening==0d0)then
      if(.not.ieee_is_finite(radius_input).or.radius_input<0d0)bad=1
    endif
    call collective_bad()
    if(bad/=0)return
    m=[n(1),n(2)/dims(1),n(3)/dims(2)];lo=[0,coords(1)*m(2),coords(2)*m(3)]
    ng=product(m);nt=size(target,2)
    if(size(target,1)/=ng.or.size(target,3)/=1.or.nt<0.or.any(shape(action)/=shape(target)))bad=1
    if(size(op%source,1)/=ng)bad=1
    if(.not.all(ieee_is_finite(real(target))).or..not.all(ieee_is_finite(aimag(target))))bad=1
    radius=.5d0*minval(n*h)
    if(radius_input>0d0)radius=radius_input
    if(screening==0d0.and.radius>.5d0*minval(n*h)*(1d0+1d-12))bad=1
    call collective_bad()
    if(bad/=0)return
    orb_rank=0;orb_size=1
    if(present(comm_o))call comm_get_groupinfo(comm_o,orb_rank,orb_size)
    ! FFT peers must use identical batch sizes and source broadcast counts.
    nt_max=nt
    call comm_get_max(nt_max,comm_r)
    if(nt_max/=nt)bad=1
    counts_max=size(op%source,2)
    call comm_get_max(counts_max,comm_r)
    if(counts_max/=size(op%source,2))bad=1
    if(.not.all(ieee_is_finite(real(op%source))).or..not.all(ieee_is_finite(aimag(op%source))))bad=1
    call collective_bad()
    if(bad/=0)return
    if(present(comm_o))allocate(source_column(ng))
    allocate(multiplier(ng));pi=acos(-1d0);g=0
    ! Keep the forward FFT in Z pencils: local storage order is (z,x,y).
    m=[n(1)/dims(1),n(2)/dims(2),n(3)]
    lo=[coords(1)*m(1),coords(2)*m(2),0]
    do y=0,m(2)-1;do x=0,m(1)-1;do z=0,m(3)-1
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
    allocate(density(ng,min(4,nt)),spectrum(ng,min(4,nt)))
    ! Grid rows match across comm_o; stream one source column from its owner.
    ! Empty source/target partitions participate in all orbital collectives.
    do owner=0,orb_size-1
      count=size(op%source,2)
      if(present(comm_o))call comm_bcast(count,comm_o,owner)
      do i=1,count
        if(present(comm_o))then
          if(orb_rank==owner)source_column=op%source(:,i)
          call comm_bcast(source_column,comm_o,owner)
        endif
        do first=1,nt,4
          nb=min(4,nt-first+1)
          do j=1,nb
            if(present(comm_o))then
              density(:,j)=conjg(source_column)*target(:,first+j-1,1)
            else
              density(:,j)=conjg(op%source(:,i))*target(:,first+j-1,1)
            endif
          enddo
          call pencil_transform(n,dims,coords,comm,density(:,:nb),spectrum(:,:nb),-1,status,spectral_z=.true.)
          bad=status
          if(present(comm_o))call comm_get_max(bad,comm_r)
          if(bad/=0)exit
          do j=1,nb
            spectrum(:,j)=spectrum(:,j)*multiplier
          enddo
          ! Inverse pencil_transform already includes 1/product(n).
          call pencil_transform(n,dims,coords,comm,spectrum(:,:nb),density(:,:nb),1,status,spectral_z=.true.)
          bad=status
          if(present(comm_o))call comm_get_max(bad,comm_r)
          if(bad/=0)exit
          do j=1,nb
            if(present(comm_o))then
              action(:,first+j-1,1)=action(:,first+j-1,1)-source_column*density(:,j)
            else
              action(:,first+j-1,1)=action(:,first+j-1,1)-op%source(:,i)*density(:,j)
            endif
          enddo
        enddo
        if(present(comm_o))call collective_bad()
        if(bad/=0)return
      enddo
    enddo
    status=0
  contains
    subroutine collective_bad()
      call comm_get_max(bad,comm_r)
      if(present(comm_o))call comm_get_max(bad,comm_o)
      ! Preserve a nonzero return on every peer, even after a successful local FFT.
      if(bad/=0)status=1
    end subroutine
  end subroutine
end module
