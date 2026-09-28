! Gamma exchange on x-complete y/z pencils, with optional compact source action.
! Refresh and apply accept orbital-local columns through optional comm_o.
! Grid rows and FFT work are spatially local.
! Collective contract: n/h/dims/radius/omega/maxiter and call order agree;
! band counts agree within spatial groups, and may differ across orbital groups.
! coords and local grid rows vary. Communicators follow spatial coordinate order.
module hse_spatial
  use iso_fortran_env, only: int64
  use exx_pair_candidates, only: exx_pair_catalog,pair_catalog_build,pair_catalog_query,pair_source_box
  use exx_spatial_local, only: s_exx_spatial_local,spatial_local_init,spatial_local_apply,spatial_local_destroy
  use communication, only: comm_summation,comm_get_max,comm_bcast,comm_get_groupinfo
  use exx_orbitals, only: orbital_layout,orbital_check,orbital_overlap,orbital_rotate
  use fftw_pencils, only: pencil_transform
  use hse_wannier_gauge, only: gauge_transport,gauge_minimize_gamma_inplace
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: spatial_exx_state,spatial_exx_refresh,spatial_exx_apply,spatial_exx_canonical_source
  type spatial_exx_state
    integer :: updates=0,iterations=0,localization_status=1,last_localization_status=1
    logical :: compact=.false.,seed_localized=.false.,seed_needed=.true.,retained_gauge=.false.
    logical :: retain_accepted_gauge=.false.
    integer :: screen_mode=0 ! 0 off, 1 diagnose, 2 omit
    real(8) :: screen_tolerance=0d0,screen_bound=0d0,screen_cpu_seconds=0d0
    integer(int64) :: screen_candidates=0,screen_skipped=0
    integer(int64) :: pair_products=0,pair_catalog_entries=0
    integer(int64) :: local_pairs=0,local_points=0,global_pairs=0
    real(8) :: spread=0d0,gradient=0d0,min_singular=0d0
    complex(8),allocatable :: gauge(:,:,:),previous(:,:,:),source(:,:)
  end type
contains
  subroutine spatial_exx_canonical_source(op,psi,occupation,comm_r,status,comm_o)
    type(spatial_exx_state),intent(inout) :: op
    complex(8),intent(in) :: psi(:,:,:)
    real(8),intent(in) :: occupation(:,:)
    integer,intent(in) :: comm_r
    integer,intent(in),optional :: comm_o
    integer,intent(out) :: status
    integer :: bad,j
    bad=0;status=1
    if(size(psi,3)/=1.or.any(shape(occupation)/=[size(psi,2),1]))bad=1
    if(.not.all(ieee_is_finite(real(psi))).or..not.all(ieee_is_finite(aimag(psi))))bad=1
    if(any(occupation<0d0).or.any(occupation>2d0).or..not.all(ieee_is_finite(occupation)))bad=1
    call comm_get_max(bad,comm_r)
    if(present(comm_o))call comm_get_max(bad,comm_o)
    if(bad/=0)return
    ! No gauge, overlap transport, or extra previous-state grid is needed.
    if(allocated(op%gauge))deallocate(op%gauge)
    if(allocated(op%previous))deallocate(op%previous)
    op%source=psi(:,:,1)
    do j=1,size(psi,2)
      op%source(:,j)=op%source(:,j)*sqrt(occupation(j,1)/2d0)
    enddo
    op%iterations=0;op%localization_status=2;op%last_localization_status=2
    op%spread=-1d0;op%gradient=-1d0;op%min_singular=0d0
    op%compact=.false.;op%retained_gauge=.false.;op%screen_mode=0
    op%updates=op%updates+1;status=0
  end subroutine spatial_exx_canonical_source

  subroutine spatial_exx_refresh(op,n,h,dims,coords,comm,comm_r,psi,maxiter,tolerance,status,occupation,comm_o)
    type(spatial_exx_state),intent(inout) :: op
    integer,intent(in) :: n(3),dims(2),coords(2),comm(2),comm_r,maxiter
    integer,intent(in),optional :: comm_o
    real(8),intent(in) :: h(3),tolerance
    real(8),intent(in),optional :: occupation(:,:)
    complex(8),intent(in) :: psi(:,:,:)
    integer,intent(out) :: status
    complex(8),allocatable :: raw(:,:,:,:),shifted(:,:),phase(:),transported(:,:,:)
    logical :: accepted
    real(8) :: b(3,6),weights(6),delta,pi,position(3)
    integer :: no,ng,m(3),lo(3),x,y,z,g,j,axis,bad
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
        op%last_localization_status=1
        op%seed_needed=.true.
        op%gauge=0d0
        do j=1,no
          op%gauge(j,j,1)=1d0
        enddo
      endif
    endif
    accepted=op%retain_accepted_gauge.and.op%seed_localized.and.op%last_localization_status==0
    op%retained_gauge=.false.
    if(accepted)transported=op%gauge
    op%iterations=0;op%localization_status=2;op%spread=-1d0;op%gradient=-1d0
    if(maxiter>0)then
      allocate(raw(no,no,6,1),shifted(ng,no),phase(ng))
      pi=acos(-1d0);b=0d0
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
      if(op%seed_localized.and.op%seed_needed)then
        call projected_position_seed(raw,op%gauge(:,:,1),bad)
        call comm_get_max(bad,comm_r)
        if(bad/=0)return
        op%seed_needed=.false.
      endif
      ! This backend constructs Gamma +/- links; consume them with the shared Jacobi solver.
      call gauge_minimize_gamma_inplace(op%gauge,raw,b,weights,maxiter,tolerance,op%spread,op%gradient, &
        op%iterations,op%localization_status)
      if(op%localization_status/=0.and.accepted)then
        ! Keep the accepted gauge transported into the current occupied space.
        ! The failed minimization remains visible in localization_status.
        op%gauge=transported;op%retained_gauge=.true.
      else
        op%last_localization_status=op%localization_status
      endif
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
    complex(8),allocatable :: raw(:,:,:,:),phase(:),overlap(:,:),left(:,:),right(:,:),work(:),transported(:,:,:)
    logical :: accepted
    real(8),allocatable :: singular(:),rwork(:),occupation_weights(:)
    real(8) :: b(3,6),weights(6),delta,pi,position(3)
    integer :: no,ng,nlocal,first,m(3),lo(3),axis,g,x,y,z,j,bad,initialized,total_initialized
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
        op%last_localization_status=1
        op%seed_needed=.true.
        op%gauge=0d0
        do j=1,no
          op%gauge(j,j,1)=1d0
        enddo
        bad=0
      endif
    endif
    accepted=op%retain_accepted_gauge.and.op%seed_localized.and.op%last_localization_status==0
    op%retained_gauge=.false.
    if(accepted)transported=op%gauge
    op%iterations=0;op%localization_status=2;op%spread=-1d0;op%gradient=-1d0
    if(maxiter>0)then
      allocate(raw(no,no,6,1),phase(ng))
      pi=acos(-1d0);b=0d0
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
      if(op%seed_localized.and.op%seed_needed)then
        call projected_position_seed(raw,op%gauge(:,:,1),bad)
        call orbital_check(bad,comm_r,comm_o)
        if(bad/=0)return
        op%seed_needed=.false.
      endif
      ! This backend constructs Gamma +/- links; consume them with the shared Jacobi solver.
      call gauge_minimize_gamma_inplace(op%gauge,raw,b,weights,maxiter,tolerance,op%spread,op%gradient, &
        op%iterations,op%localization_status)
      if(op%localization_status/=0.and.accepted)then
        ! Keep the accepted gauge transported into the current occupied space.
        ! The failed minimization remains visible in localization_status.
        op%gauge=transported;op%retained_gauge=.true.
      else
        op%last_localization_status=op%localization_status
      endif
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

  ! A symmetry-adapted occupied basis can be a stationary saddle of the spread
  ! functional. Diagonalize a fixed Hermitian projected periodic-position
  ! combination before the first minimization to select localized directions.
  ! The links are already reduced No x No matrices; no grid/WF gather is needed.
  subroutine projected_position_seed(raw,gauge,status)
    complex(8),intent(in) :: raw(:,:,:,:)
    complex(8),intent(out) :: gauge(:,:)
    integer,intent(out) :: status
    complex(8),allocatable :: projected(:,:),work(:)
    real(8),allocatable :: eigenvalues(:),rwork(:)
    real(8) :: cosine_weight(3),sine_weight(3)
    complex(8) :: coefficient,phase
    integer :: no,axis,j,pivot
    no=size(gauge,1);status=1
    allocate(projected(no,no),eigenvalues(no),work(max(1,2*no)),rwork(max(1,3*no-2)))
    cosine_weight=sqrt([2d0,3d0,5d0]);sine_weight=sqrt([7d0,11d0,13d0])
    projected=0d0
    do axis=1,3
      coefficient=cmplx(cosine_weight(axis),-sine_weight(axis),8)/2d0
      projected=projected+coefficient*raw(:,:,axis,1)+conjg(coefficient)*raw(:,:,axis+3,1)
    enddo
    ! Remove only floating-point anti-Hermitian roundoff before LAPACK.
    projected=(projected+conjg(transpose(projected)))/2d0
    call zheev('V','U',no,projected,no,eigenvalues,work,size(work),rwork,status)
    if(status/=0)then
      status=1;return
    endif
    if(.not.all(ieee_is_finite(real(projected))).or..not.all(ieee_is_finite(aimag(projected))))then
      status=1;return
    endif
    do j=1,no
      pivot=maxloc(abs(projected(:,j)),dim=1)
      phase=projected(pivot,j)/abs(projected(pivot,j))
      projected(:,j)=projected(:,j)*conjg(phase)
    enddo
    gauge=projected
  end subroutine projected_position_seed

  ! Optional comm_o joins matching grid pencils across orbital groups. Source and
  ! target columns may have different (including zero) local counts. Counts must
  ! agree within comm_r, and the communicators form a spatial/orbital product.
  subroutine spatial_exx_apply(op,n,h,dims,coords,comm,comm_r,radius_input,target,action,status,omega,comm_o)
    type(spatial_exx_state),intent(inout) :: op
    integer,intent(in) :: n(3),dims(2),coords(2),comm(2),comm_r
    integer,intent(in),optional :: comm_o
    real(8),intent(in) :: h(3),radius_input
    real(8),intent(in),optional :: omega
    complex(8),intent(in) :: target(:,:,:)
    complex(8),intent(out) :: action(:,:,:)
    integer,intent(out) :: status
    complex(8),allocatable :: density(:,:),spectrum(:,:),source_column(:)
    real(8),allocatable :: multiplier(:)
    type(exx_pair_catalog) :: catalogue
    type(s_exx_spatial_local) :: compact_plan
    complex(8),allocatable :: compact_action(:,:)
    logical :: compact_used
    integer,allocatable :: selected(:),broad_kept(:)
    real(8),allocatable :: omitted(:),pair_norms(:,:),pair_totals(:,:),source_norms(:),global_norms(:)
    real(8) :: kernel_local(2),kernel_sum(2),lambda,kzero,krms,budget,normq,qmax,qnorm_local
    real(8) :: candidate_bound,summary_local(3),summary_total(3),cpu_start,cpu_end,bound_scale
    integer :: source_total,target_total,nselected,ncandidate,k,broad_sources,box_lower(3),box_upper(3)
    integer :: mesh_local(3),mesh_lo(3)
    real(8) :: envelope_max,envelope_norm,envelope_factor,threshold_floor,threshold_pair,factor
    real(8) :: pair_stats(2),global_pair_stats(2)
    integer(int64) :: compact_pairs,compact_points
    real(8) :: radius,pi,q(3),q2,screening
    integer :: ng,nt,m(3),lo(3),x,y,z,g,p(3),i,j,first,nb,bad,owner,orb_rank,orb_size,count,counts_max,nt_max
    status=1;action=0d0;bad=0
    op%screen_candidates=0;op%screen_skipped=0;op%screen_bound=0d0;op%screen_cpu_seconds=0d0
    op%pair_products=0;op%pair_catalog_entries=0
    screening=0d0
    if(present(omega))screening=omega
    if(.not.ieee_is_finite(screening).or.screening<0d0)bad=1
    if(.not.allocated(op%source))bad=1
    if(op%screen_mode<0.or.op%screen_mode>2)bad=1
    if(.not.ieee_is_finite(op%screen_tolerance).or.op%screen_tolerance<0d0)bad=1
    if(any(n<1).or.any(dims<1).or.any(h<=0d0).or..not.all(ieee_is_finite(h)))bad=1
    if(screening==0d0)then
      if(.not.ieee_is_finite(radius_input).or.radius_input<0d0)bad=1
    endif
    call collective_bad()
    if(bad/=0)return
    m=[n(1),n(2)/dims(1),n(3)/dims(2)];lo=[0,coords(1)*m(2),coords(2)*m(3)]
    ng=product(m);nt=size(target,2);mesh_local=m;mesh_lo=lo
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
    allocate(source_column(ng))
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
    allocate(selected(nt),omitted(nt),broad_kept(nt));omitted=0d0;broad_kept=0;broad_sources=0
    if(op%screen_mode/=0)then
      call cpu_time(cpu_start)
      ! Norms of the actual discrete convolution, including the G=0 mode.
      kernel_local=[sum(abs(multiplier)),sum(abs(multiplier)**2)]
      call comm_summation(kernel_local,kernel_sum,2,comm_r)
      lambda=maxval(abs(multiplier));call max_scalar(lambda,comm_r)
      kzero=kernel_sum(1)/real(product(int(n,int64)),8)
      krms=sqrt(kernel_sum(2)/real(product(int(n,int64)),8))
      source_total=size(op%source,2);target_total=nt
      if(present(comm_o))then
        call comm_summation(size(op%source,2),source_total,comm_o)
        call comm_summation(nt,target_total,comm_o)
      endif
      budget=op%screen_tolerance/(real(max(1,source_total),8)*sqrt(real(max(1,target_total),8)))
      allocate(pair_norms(2,nt),pair_totals(2,nt))
      ! A common floor retains every block that any source query might need.
      ! Upper bounds on max|q| and ||q||_2 avoid one collective per source.
      envelope_max=0d0
      if(size(op%source)>0)envelope_max=maxval(abs(op%source))
      call max_scalar(envelope_max,comm_r)
      if(present(comm_o))call max_scalar(envelope_max,comm_o)
      allocate(source_norms(size(op%source,2)),global_norms(size(op%source,2)))
      source_norms=0d0
      if(envelope_max>0d0)source_norms=sum((abs(op%source)/envelope_max)**2,dim=1)
      if(size(source_norms)>0)call comm_summation(source_norms,global_norms,size(source_norms),comm_r)
      envelope_norm=0d0
      if(size(global_norms)>0)envelope_norm=envelope_max*sqrt(product(h)*maxval(global_norms))
      if(present(comm_o))call max_scalar(envelope_norm,comm_o)
      envelope_factor=envelope_max*lambda*envelope_norm*(1d0+512d0*epsilon(1d0))
      threshold_floor=0d0
      if(ieee_is_finite(envelope_factor).and.envelope_factor>tiny(1d0)) &
        threshold_floor=(budget/envelope_factor)*(1d0-512d0*epsilon(1d0))
      if(.not.ieee_is_finite(threshold_floor))threshold_floor=0d0
      call pair_catalog_build(catalogue,n,mesh_lo,mesh_local,comm_r,target(:,:,1),threshold_floor,status)
      call collective_bad_status()
      if(bad/=0)return
      op%pair_catalog_entries=catalogue%entries
      call cpu_time(cpu_end);op%screen_cpu_seconds=cpu_end-cpu_start
    endif
    op%local_pairs=0;op%local_points=0;op%global_pairs=0
    if(op%compact)then
      call spatial_local_init(compact_plan,n,dims,coords,comm,multiplier,status)
      call collective_bad_status()
      if(status/=0)return
    endif
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
        else
          source_column=op%source(:,i)
        endif
        if(op%screen_mode/=0)then
          call cpu_time(cpu_start)
          qmax=maxval(abs(source_column));call max_scalar(qmax,comm_r)
          qnorm_local=0d0
          if(qmax>0d0)qnorm_local=sum((abs(source_column)/qmax)**2)
          call comm_summation(qnorm_local,normq,comm_r)
          normq=qmax*sqrt(product(h)*normq)
          factor=qmax*lambda*normq*(1d0+256d0*epsilon(1d0))
          threshold_pair=threshold_floor
          if(ieee_is_finite(factor).and.factor>tiny(1d0)) &
            threshold_pair=max(threshold_floor,(budget/factor)*(1d0-512d0*epsilon(1d0)))
          if(.not.ieee_is_finite(threshold_pair))threshold_pair=threshold_floor
          call pair_source_box(n,mesh_lo,mesh_local,comm_r,source_column,box_lower,box_upper,status)
          if(status==0)call pair_catalog_query(catalogue,box_lower,box_upper,threshold_pair,selected,ncandidate,status)
          call collective_bad_status()
          if(bad/=0)return
          ! All absent pairs have action norm <= budget. Count their complement
          ! without an Nsource-by-Ntarget skip table or per-source dense update.
          broad_sources=broad_sources+1
          do k=1,ncandidate
            j=selected(k);broad_kept(j)=broad_kept(j)+1
            pair_norms(1,k)=sum(abs(conjg(source_column)*target(:,j,1))**2)
            pair_norms(2,k)=sum(abs(conjg(source_column)*target(:,j,1)))
          enddo
          op%pair_products=op%pair_products+int(ncandidate,int64)
          op%screen_candidates=op%screen_candidates+int(nt-ncandidate,int64)
          if(ncandidate>0)call comm_summation(pair_norms(:,:ncandidate),pair_totals(:,:ncandidate),2*ncandidate,comm_r)
          nselected=0
          do k=1,ncandidate
            j=selected(k)
            candidate_bound=min(qmax*lambda*sqrt(product(h)*pair_totals(1,k)), &
              normq*kzero*pair_totals(2,k),normq*krms*sqrt(pair_totals(1,k)))
            ! Do not turn a subnormal squared norm into a zero error certificate.
            if(pair_totals(1,k)<tiny(1d0).and.pair_totals(2,k)>0d0) &
              candidate_bound=normq*kzero*pair_totals(2,k)
            candidate_bound=candidate_bound*(1d0+128d0*epsilon(1d0))
            if(ieee_is_finite(candidate_bound))then
              if(candidate_bound<=budget.and.(budget>0d0.or.pair_totals(2,k)==0d0))then
                op%screen_candidates=op%screen_candidates+1_int64
                omitted(j)=omitted(j)+candidate_bound
                if(op%screen_mode==2)cycle
              endif
            endif
            nselected=nselected+1;selected(nselected)=j
          enddo
          call cpu_time(cpu_end);op%screen_cpu_seconds=op%screen_cpu_seconds+cpu_end-cpu_start
        endif
        if(op%screen_mode/=2)then
          nselected=nt
          do j=1,nt
            selected(j)=j
          enddo
        endif
        op%screen_skipped=op%screen_skipped+int(nt-nselected,int64)
        if(op%compact)then
          allocate(compact_action(ng,nselected))
          call spatial_local_apply(compact_plan,comm_r,source_column,target(:,selected(:nselected),1),compact_action, &
            compact_used,status,compact_pairs,compact_points)
          call collective_bad_status()
          if(status/=0)then
            call spatial_local_destroy(compact_plan)
            return
          endif
          if(compact_used)then
            do j=1,nselected
              action(:,selected(j),1)=action(:,selected(j),1)+compact_action(:,j)
            enddo
            deallocate(compact_action)
            op%local_pairs=op%local_pairs+compact_pairs
            op%local_points=op%local_points+compact_points
            if(present(comm_o))call collective_bad()
            if(bad/=0)return
            cycle
          endif
          deallocate(compact_action)
        endif
        op%global_pairs=op%global_pairs+nselected
        do first=1,nselected,4
          nb=min(4,nselected-first+1)
          do j=1,nb
            if(present(comm_o))then
              density(:,j)=conjg(source_column)*target(:,selected(first+j-1),1)
            else
              density(:,j)=conjg(op%source(:,i))*target(:,selected(first+j-1),1)
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
              action(:,selected(first+j-1),1)=action(:,selected(first+j-1),1)-source_column*density(:,j)
            else
              action(:,selected(first+j-1),1)=action(:,selected(first+j-1),1)-op%source(:,i)*density(:,j)
            endif
          enddo
        enddo
        if(present(comm_o))call collective_bad()
        if(bad/=0)return
      enddo
    enddo
    call spatial_local_destroy(compact_plan)
    if(op%screen_mode/=0)then
      omitted=omitted+budget*real(broad_sources-broad_kept,8)
      pair_stats=real([op%pair_products,op%pair_catalog_entries],8);global_pair_stats=pair_stats
      if(present(comm_o))call comm_summation(pair_stats,global_pair_stats,2,comm_o)
      op%pair_products=int(global_pair_stats(1),int64);op%pair_catalog_entries=int(global_pair_stats(2),int64)
      bound_scale=0d0
      if(nt>0)bound_scale=maxval(omitted)
      if(present(comm_o))call max_scalar(bound_scale,comm_o)
      summary_local=[0d0,real(op%screen_candidates,8),real(op%screen_skipped,8)]
      if(bound_scale>0d0)summary_local(1)=sum((omitted/bound_scale)**2)
      summary_total=summary_local
      if(present(comm_o))call comm_summation(summary_local,summary_total,3,comm_o)
      op%screen_bound=bound_scale*sqrt(summary_total(1))
      op%screen_candidates=int(summary_total(2),int64);op%screen_skipped=int(summary_total(3),int64)
      call max_scalar(op%screen_cpu_seconds,comm_r)
      if(present(comm_o))call max_scalar(op%screen_cpu_seconds,comm_o)
    endif
    status=0
  contains
    subroutine max_scalar(value,group)
      real(8),intent(inout) :: value
      integer,intent(in) :: group
      real(8) :: result(1)
      call comm_get_max([value],result,1,group)
      value=result(1)
    end subroutine
    subroutine collective_bad_status()
      bad=status
      call collective_bad()
    end subroutine
    subroutine collective_bad()
      call comm_get_max(bad,comm_r)
      if(present(comm_o))call comm_get_max(bad,comm_o)
      ! Preserve a nonzero return on every peer, even after a successful local FFT.
      if(bad/=0)status=1
    end subroutine
  end subroutine
end module
