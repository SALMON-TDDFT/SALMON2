! Gamma exchange on native Cartesian domains, with optional compact source action.
! Refresh and apply accept orbital-local columns through optional comm_o.
! Grid rows and FFT work are spatially local.
! Collective contract: n/h/dims/radius/omega/maxiter and call order agree;
! band counts agree within spatial groups, and may differ across orbital groups.
! coords and local grid rows vary. Communicators follow spatial coordinate order.
module exx_spatial
  use exx_sparse_orbitals, only: s_sparse_orbitals,sparse_valid,sparse_column,sparse_norms,sparse_clear
  use iso_fortran_env, only: int64
  use exx_pair_candidates, only: exx_pair_catalog,pair_catalog_build,pair_catalog_query,pair_source_box
  use exx_batch_backend, only: local_backend_factory
  use exx_spatial_local, only: s_exx_spatial_local,spatial_local_init,spatial_local_apply,spatial_local_destroy, &
    local_batch_action,s_exx_sr,sr_prepare,sr_apply,sr_destroy,sr_mask_kernel
  use communication, only: comm_summation,comm_get_max,comm_bcast,comm_get_groupinfo
  use exx_orbitals, only: orbital_layout,orbital_check,orbital_overlap,orbital_rotate
  use exx_distributed_metric, only: distributed_metric_available
  use exx_distributed_gauge, only: s_exx_gauge,gauge_tiles_clear,gauge_tiles_refresh,gauge_tiles_rotate
  use fftw_blocks, only: mesh_transform,block_layout
  use exx_wannier_gauge, only: gauge_transport,gauge_minimize_gamma_inplace
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: spatial_exx_state,spatial_exx_refresh,spatial_exx_apply,spatial_exx_canonical_source
  type spatial_exx_state
    integer :: updates=0,iterations=0,localization_status=1,last_localization_status=1
    logical :: compact=.false.,seed_localized=.false.,seed_needed=.true.,retained_gauge=.false.
    logical :: retain_accepted_gauge=.false.
    procedure(local_batch_action),pointer,nopass :: local_batch=>null()
    procedure(local_backend_factory),pointer,nopass :: create_local_backend=>null()
    integer :: local_batch_size=8
    real(8) :: sr_tolerance=0d0
    integer :: screen_mode=0 ! 0 off, 1 diagnose, 2 omit
    real(8) :: screen_tolerance=0d0,screen_bound=0d0,screen_cpu_seconds=0d0
    integer(int64) :: screen_candidates=0,screen_skipped=0
    integer(int64) :: pair_products=0,pair_catalog_entries=0,pair_product_points=0
    integer(int64) :: local_pairs=0,local_points=0,global_pairs=0
    real(8) :: spread=0d0,gradient=0d0,min_singular=0d0
    type(s_exx_gauge) :: gauge_tiles
    type(s_sparse_orbitals) :: sparse_source
    complex(8),allocatable :: gauge(:,:,:),previous(:,:,:),source(:,:)
  end type
  ! Rank-specific loops keep IEEE inquiries scalar on Fujitsu compilers.
  private :: salmon_all_finite,finite_real_1d,finite_real_2d,finite_real_3d
  interface salmon_all_finite
    module procedure finite_real_1d,finite_real_2d,finite_real_3d
  end interface
contains
  subroutine spatial_exx_canonical_source(op,psi,occupation,comm_r,status,comm_o)
    implicit none
    type(spatial_exx_state),intent(inout) :: op
    complex(8),intent(in) :: psi(:,:,:)
    real(8),intent(in) :: occupation(:,:)
    integer,intent(in) :: comm_r
    integer,intent(in),optional :: comm_o
    integer,intent(out) :: status
    integer :: bad,j
    bad=0;status=1
    if(size(psi,3)/=1.or.any(shape(occupation)/=[size(psi,2),1]))bad=1
    if(.not.salmon_all_finite(real(psi)).or..not.salmon_all_finite(aimag(psi)))bad=1
    if(any(occupation<0d0).or.any(occupation>2d0).or..not.salmon_all_finite(occupation))bad=1
    call comm_get_max(bad,comm_r)
    if(present(comm_o))call comm_get_max(bad,comm_o)
    if(bad/=0)return
    ! No gauge, overlap transport, or extra previous-state grid is needed.
    call gauge_tiles_clear(op%gauge_tiles)
    if(allocated(op%gauge))deallocate(op%gauge)
    if(allocated(op%previous))deallocate(op%previous)
    call sparse_clear(op%sparse_source)
    op%source=psi(:,:,1)
    do j=1,size(psi,2)
      op%source(:,j)=op%source(:,j)*sqrt(occupation(j,1)/2d0)
    enddo
    op%iterations=0;op%localization_status=2;op%last_localization_status=2
    op%spread=-1d0;op%gradient=-1d0;op%min_singular=0d0
    op%compact=.false.;op%retained_gauge=.false.;op%screen_mode=0
    op%updates=op%updates+1;status=0
  end subroutine spatial_exx_canonical_source

  subroutine spatial_exx_refresh(op,n,h,dims,coords,comm,comm_r,psi,maxiter,tolerance,status,occupation,comm_o,comm_matrix)
    implicit none
    integer,intent(in),optional :: comm_matrix
    type(spatial_exx_state),intent(inout) :: op
    integer,intent(in) :: n(3),dims(:),coords(:),comm(:),comm_r,maxiter
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
      call refresh_orbitals(op,n,h,dims,coords,comm_r,comm_o,psi,maxiter,tolerance,status,occupation,comm_matrix)
      return
    endif
    no=size(psi,2);ng=size(psi,1);status=1;bad=0
    if(allocated(op%gauge_tiles%matrix))bad=1
    if(no<1.or.size(psi,3)/=1.or.any(n<1).or.any(dims<1))bad=1
    if(any(h<=0d0).or..not.salmon_all_finite(h))bad=1
    if(.not.salmon_all_finite(real(psi)).or..not.salmon_all_finite(aimag(psi)))bad=1
    if(maxiter<0.or.tolerance<=0d0.or..not.ieee_is_finite(tolerance))bad=1
    if(present(occupation))then
      if(any(shape(occupation)/=[no,1]))bad=1
      if(any(occupation<0d0).or.any(occupation>2d0).or..not.salmon_all_finite(occupation))bad=1
    endif
    call comm_get_max(bad,comm_r)
    if(bad/=0)return
    call block_layout(n,dims,coords,m,lo,bad)
    if(ng/=product(m).or.any(coords<0).or.any(coords>=dims))bad=1
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
    if(.not.salmon_all_finite(real(op%gauge)).or..not.salmon_all_finite(aimag(op%gauge)))bad=1
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
      implicit none
      complex(8),intent(inout) :: a(:,:)
      complex(8) :: total(size(a,1),size(a,2))
      call comm_summation(a,total,size(a),comm_r)
      a=total
    end subroutine
  end subroutine

  subroutine refresh_orbitals(op,n,h,dims,coords,comm_r,comm_o,psi,maxiter,tolerance,status,occupation,comm_matrix)
    implicit none
    integer,intent(in),optional :: comm_matrix
    type(spatial_exx_state),intent(inout) :: op
    integer,intent(in) :: n(3),dims(:),coords(:),comm_r,comm_o,maxiter
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
    if(any(h<=0d0).or..not.salmon_all_finite(h))bad=1
    if(.not.salmon_all_finite(real(psi)).or..not.salmon_all_finite(aimag(psi)))bad=1
    if(maxiter<0.or.tolerance<=0d0.or..not.ieee_is_finite(tolerance))bad=1
    if(present(occupation))then
      if(any(shape(occupation)/=[nlocal,1]))bad=1
      if(any(occupation<0d0).or.any(occupation>2d0).or..not.salmon_all_finite(occupation))bad=1
    endif
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    call block_layout(n,dims,coords,m,lo,bad)
    if(ng/=product(m).or.any(coords<0).or.any(coords>=dims))bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    call orbital_layout(nlocal,comm_r,comm_o,counts,first,bad)
    if(bad/=0.or.sum(counts)<1)return
    no=sum(counts)
    if(present(comm_matrix))then
      if(distributed_metric_available(comm_matrix))then
        call refresh_tiles(op,n,h,m,lo,comm_r,comm_o,comm_matrix,psi,maxiter,tolerance,status,occupation)
        return
      endif
    endif
    if(allocated(op%gauge_tiles%matrix))bad=1
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
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
        if(.not.salmon_all_finite(real(op%previous)).or..not.salmon_all_finite(aimag(op%previous)))bad=1
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
        if(.not.salmon_all_finite(singular).or.op%min_singular<1d-8)bad=1
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
    if(.not.salmon_all_finite(real(op%gauge)).or..not.salmon_all_finite(aimag(op%gauge)))bad=1
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

  subroutine refresh_tiles(op,n,h,m,lo,comm_r,comm_o,comm_matrix,psi,maxiter,tolerance,status,occupation)
    implicit none
    type(spatial_exx_state),intent(inout) :: op
    integer,intent(in) :: n(3),m(3),lo(3),comm_r,comm_o,comm_matrix,maxiter
    complex(8),intent(in) :: psi(:,:,:)
    real(8),intent(in) :: h(3),tolerance
    real(8),intent(in),optional :: occupation(:,:)
    integer,intent(out) :: status
    complex(8),allocatable :: phase(:,:)
    real(8),allocatable :: occupation_weights(:)
    real(8) :: b(3,6),weights(6),delta,pi,position(3)
    integer :: ng,nlocal,axis,g,x,y,z,bad
    ng=size(psi,1);nlocal=size(psi,2)
    bad=0
    ! A state must not silently change matrix distribution after initialization.
    if(allocated(op%gauge))bad=1
    call orbital_check(bad,comm_r,comm_o)
    status=1
    if(bad/=0)return
    allocate(phase(ng,3));phase=0d0
    pi=acos(-1d0);b=0d0
    do axis=1,3
      delta=2*pi/(n(axis)*h(axis));b(axis,axis)=delta;b(axis,axis+3)=-delta
      weights(axis)=1d0/(2*delta**2);weights(axis+3)=weights(axis)
      if(maxiter==0)cycle
      g=0
      do z=0,m(3)-1
        do y=0,m(2)-1
          do x=0,m(1)-1
            g=g+1;position=([x,y,z]+lo)*h
            phase(g,axis)=exp(cmplx(0d0,-position(axis)*delta,8))
          enddo
        enddo
      enddo
    enddo
    call gauge_tiles_refresh(op%gauge_tiles,psi(:,:,1),op%previous,product(h),phase,b,weights,comm_r,comm_o, &
      comm_matrix,maxiter,tolerance,op%seed_localized,op%seed_needed,op%retain_accepted_gauge, &
      op%last_localization_status,op%retained_gauge,op%min_singular,op%spread,op%gradient,op%iterations, &
      op%localization_status,status)
    if(status/=0)return
    if(allocated(op%source))deallocate(op%source)
    if(allocated(op%previous))deallocate(op%previous)
    allocate(op%source(ng,nlocal),op%previous(ng,nlocal,1),occupation_weights(nlocal))
    call gauge_tiles_rotate(op%gauge_tiles,psi(:,:,1),comm_r,comm_o,op%previous(:,:,1),status)
    if(status/=0)return
    occupation_weights=1d0
    if(present(occupation))occupation_weights=sqrt(occupation(:,1)/2d0)
    call gauge_tiles_rotate(op%gauge_tiles,psi(:,:,1),comm_r,comm_o,op%source,status,weights=occupation_weights)
    if(status==0)op%updates=op%updates+1
  end subroutine

  ! A symmetry-adapted occupied basis can be a stationary saddle of the spread
  ! functional. Diagonalize a fixed Hermitian projected periodic-position
  ! combination before the first minimization to select localized directions.
  ! The links are already reduced No x No matrices; no grid/WF gather is needed.
  subroutine projected_position_seed(raw,gauge,status)
    implicit none
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
    if(.not.salmon_all_finite(real(projected)).or..not.salmon_all_finite(aimag(projected)))then
      status=1;return
    endif
    do j=1,no
      pivot=maxloc(abs(projected(:,j)),dim=1)
      phase=projected(pivot,j)/abs(projected(pivot,j))
      projected(:,j)=projected(:,j)*conjg(phase)
    enddo
    gauge=projected
  end subroutine projected_position_seed

  ! Optional comm_o joins matching grid blocks across orbital groups. Source and
  ! target columns may have different (including zero) local counts. Counts must
  ! agree within comm_r, and the communicators form a spatial/orbital product.
  subroutine spatial_exx_apply(op,n,h,dims,coords,comm,comm_r,radius_input,target,action,status,omega,comm_o, &
                               screen_target_count)
!$  use omp_lib, only: omp_get_max_threads
    implicit none
    type(spatial_exx_state),intent(inout) :: op
    integer,intent(in) :: n(3),dims(:),coords(:),comm(:),comm_r
    integer,intent(in),optional :: comm_o,screen_target_count
    real(8),intent(in) :: h(3),radius_input
    real(8),intent(in),optional :: omega
    complex(8),intent(in) :: target(:,:,:)
    complex(8),intent(out) :: action(:,:,:)
    integer,intent(out) :: status
    complex(8),allocatable :: density(:,:),spectrum(:,:),source_column(:)
    real(8),allocatable :: multiplier(:)
    type(exx_pair_catalog) :: catalogue
    type(s_exx_spatial_local) :: compact_plan
    type(s_exx_sr) :: sr
    complex(8),allocatable :: compact_action(:,:)
    logical :: compact_used,packed_source
    integer :: no_source,nnz,k0,point_workers
    integer,allocatable :: wire_rows(:)
    complex(8),allocatable :: wire_values(:)
    integer,allocatable :: selected(:),broad_kept(:),active_rows(:)
    real(8),allocatable :: omitted(:),pair_norms(:,:),pair_totals(:,:),source_norms(:),global_norms(:),pair_values(:)
    real(8) :: kernel_local(2),kernel_sum(2),lambda,kzero,krms,budget,normq,qmax,qnorm_local
    real(8) :: candidate_bound,summary_local(3),summary_total(3),cpu_start,cpu_end,bound_scale
    integer :: source_total,target_total,nselected,ncandidate,k,broad_sources,box_lower(3),box_upper(3)
    integer :: mesh_local(3),mesh_lo(3),nactive,stride(3)
    real(8) :: envelope_max,envelope_norm,envelope_factor,threshold_floor,threshold_pair,factor
    real(8) :: pair_stats(3),global_pair_stats(3),point_total
    integer(int64) :: compact_pairs,compact_points
    real(8) :: radius,pi,q(3),q2,screening,sr_cost,fft_volume
    integer(int64) :: sr_pairs,wf_pairs
    integer :: ng,nt,m(3),lo(3),x,y,z,g,p(3),i,j,first,nb,bad,owner,orb_rank,orb_size,count,counts_max,nt_max
    status=1;action=0d0;bad=0
    op%screen_candidates=0;op%screen_skipped=0;op%screen_bound=0d0;op%screen_cpu_seconds=0d0
    op%pair_products=0;op%pair_catalog_entries=0;op%pair_product_points=0
    screening=0d0
    if(present(omega))screening=omega
    if(.not.ieee_is_finite(screening).or.screening<0d0)bad=1
    packed_source=.not.allocated(op%source).and.allocated(op%sparse_source%offset)
    no_source=0
    if(allocated(op%source))then
      no_source=size(op%source,2)
    else if(packed_source)then
      no_source=op%sparse_source%no
    else
      bad=1
    endif
    if(op%screen_mode<0.or.op%screen_mode>2)bad=1
    if(.not.ieee_is_finite(op%screen_tolerance).or.op%screen_tolerance<0d0)bad=1
    if(.not.ieee_is_finite(op%sr_tolerance).or.op%sr_tolerance<0d0.or.op%sr_tolerance>=1d0)bad=1
    if(op%sr_tolerance>0d0)then
      if(.not.op%compact.or.screening<=0d0.or.associated(op%create_local_backend))bad=1
    endif
    if(any(n<1).or.any(dims<1).or.any(h<=0d0).or..not.salmon_all_finite(h))bad=1
    if(screening==0d0)then
      if(.not.ieee_is_finite(radius_input).or.radius_input<0d0)bad=1
    endif
    call collective_bad()
    if(bad/=0)return
    call block_layout(n,dims,coords,m,lo,bad)
    ng=product(m);nt=size(target,2);mesh_local=m;mesh_lo=lo
    if(size(target,1)/=ng.or.size(target,3)/=1.or.nt<0.or.any(shape(action)/=shape(target)))bad=1
    if(packed_source)then
      if(.not.sparse_valid(op%sparse_source,ng,no_source))bad=1
    else
      if(size(op%source,1)/=ng)bad=1
    endif
    if(.not.salmon_all_finite(real(target)).or..not.salmon_all_finite(aimag(target)))bad=1
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
    counts_max=no_source
    call comm_get_max(counts_max,comm_r)
    if(counts_max/=no_source)bad=1
    if(.not.packed_source)then
      if(.not.salmon_all_finite(real(op%source)).or..not.salmon_all_finite(aimag(op%source)))bad=1
    endif
    call collective_bad()
    if(bad/=0)return
    allocate(source_column(ng))
    allocate(multiplier(ng));pi=acos(-1d0);g=0
    ! Every FFT preserves Cartesian xyz ownership, including two-axis callers.
    stride=[1,m(1),m(1)*m(2)]
    if(screening>0d0)then
!$omp parallel do collapse(2) default(none) schedule(static) &
!$omp private(y,x,z,g,p,q,q2) shared(m,lo,n,h,pi,screening,radius,multiplier,stride)
      do y=0,m(2)-1;do x=0,m(1)-1;do z=0,m(3)-1
        g=1+x*stride(1)+y*stride(2)+z*stride(3)
        p=[x,y,z]+lo
        where(p>=(n+1)/2)p=p-n
        q=2*pi*p/(n*h);q2=sum(q*q)
        if(q2<1d-24)then
          multiplier(g)=pi/screening**2
        else
          multiplier(g)=4*pi*(1d0-exp(-q2/(4*screening**2)))/q2
        endif
      enddo;enddo;enddo
!$omp end parallel do
    else
!$omp parallel do collapse(2) default(none) schedule(static) &
!$omp private(y,x,z,g,p,q,q2) shared(m,lo,n,h,pi,screening,radius,multiplier,stride)
      do y=0,m(2)-1;do x=0,m(1)-1;do z=0,m(3)-1
        g=1+x*stride(1)+y*stride(2)+z*stride(3)
        p=[x,y,z]+lo
        where(p>=(n+1)/2)p=p-n
        q=2*pi*p/(n*h);q2=sum(q*q)
        if(q2<1d-24)then
          multiplier(g)=2*pi*radius**2
        else
          multiplier(g)=8*pi*sin(.5d0*sqrt(q2)*radius)**2/q2
        endif
      enddo;enddo;enddo
!$omp end parallel do
    endif
    allocate(selected(nt),omitted(nt),broad_kept(nt));omitted=0d0;broad_kept=0;broad_sources=0
    if(op%screen_mode/=0)then
      call cpu_time(cpu_start)
      ! Norms of the actual discrete convolution, including the G=0 mode.
      kernel_local=[sum(abs(multiplier)),sum(abs(multiplier)**2)]
      call comm_summation(kernel_local,kernel_sum,2,comm_r)
      lambda=maxval(abs(multiplier));call max_scalar(lambda,comm_r)
      kzero=kernel_sum(1)/real(product(int(n,int64)),8)
      krms=sqrt(kernel_sum(2)/real(product(int(n,int64)),8))
      source_total=no_source;target_total=nt
      if(present(comm_o))then
        call comm_summation(no_source,source_total,comm_o)
        call comm_summation(nt,target_total,comm_o)
      endif
      ! Blocked callers retain the full-target error budget and catalogue floor.
      if(present(screen_target_count))then
        if(screen_target_count<target_total)bad=1
        target_total=screen_target_count
      endif
      call collective_bad()
      if(bad/=0)return
      budget=op%screen_tolerance/(real(max(1,source_total),8)*sqrt(real(max(1,target_total),8)))
      allocate(pair_norms(2,nt),pair_totals(2,nt))
      ! A common floor retains every block that any source query might need.
      ! Upper bounds on max|q| and ||q||_2 avoid one collective per source.
      envelope_max=0d0
      if(packed_source)then
        if(size(op%sparse_source%value)>0)envelope_max=maxval(abs(op%sparse_source%value))
      else
        if(size(op%source)>0)envelope_max=maxval(abs(op%source))
      endif
      call max_scalar(envelope_max,comm_r)
      if(present(comm_o))call max_scalar(envelope_max,comm_o)
      allocate(source_norms(no_source),global_norms(no_source))
      source_norms=0d0
      if(packed_source)then
        call sparse_norms(op%sparse_source,envelope_max,source_norms)
      else
        if(envelope_max>0d0)source_norms=sum((abs(op%source)/envelope_max)**2,dim=1)
      endif
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
      compact_plan%batch_action=>op%local_batch
      compact_plan%batch_size=op%local_batch_size
      if(associated(op%create_local_backend))call op%create_local_backend(compact_plan%backend)
      call collective_bad_status()
      if(status/=0)return
    endif
    sr_cost=huge(1d0);sr_pairs=0;wf_pairs=0
    if(op%sr_tolerance>0d0)then
      call sr_prepare(sr,compact_plan,h,screening,op%sr_tolerance,comm_r,status)
      if(status/=0)return
      if(sr%ready)then
        call sr_mask_kernel(sr,compact_plan,h)
        fft_volume=real(product(int(sr%length,int64)),8)
        sr_cost=real(sr%peers,8)*fft_volume*log(max(2d0,fft_volume))
      endif
      if(sr%rank==0)write(*,'(a,l2,3es18.9,3i8)')'EXX_SR active/radius/tail L1/L2/FFT shape: ', &
        sr%ready,sr%radius,sr%tail_l1,sr%tail_l2,sr%length
    endif
    point_workers=1
!$  point_workers=omp_get_max_threads()
    allocate(density(ng,min(4,nt)),spectrum(ng,min(4,nt)))
    ! Grid rows match across comm_o; stream one source column from its owner.
    ! Empty source/target partitions participate in all orbital collectives.
    do owner=0,orb_size-1
      count=no_source
      if(present(comm_o))call comm_bcast(count,comm_o,owner)
      do i=1,count
        ! Send support entries only; materialize one spatial column for FFTs.
        nnz=-1
        if(orb_rank==owner.and.packed_source)nnz=op%sparse_source%offset(i+1)-op%sparse_source%offset(i)
        if(present(comm_o))call comm_bcast(nnz,comm_o,owner)
        if(nnz>=0)then
          allocate(wire_rows(nnz),wire_values(nnz))
          if(orb_rank==owner)then
            k0=op%sparse_source%offset(i)
            wire_rows=op%sparse_source%row(k0:k0+nnz-1)
            wire_values=op%sparse_source%value(k0:k0+nnz-1)
          endif
          if(present(comm_o))then
            call comm_bcast(wire_rows,comm_o,owner)
            call comm_bcast(wire_values,comm_o,owner)
          endif
          source_column=0d0
          source_column(wire_rows)=wire_values
          deallocate(wire_rows,wire_values)
        else
          if(orb_rank==owner)source_column=op%source(:,i)
          if(present(comm_o))call comm_bcast(source_column,comm_o,owner)
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
          if(bad/=0)then
            call sr_destroy(sr)
            call spatial_local_destroy(compact_plan)
            return
          endif
          ! All absent pairs have action norm <= budget. Count their complement
          ! without an Nsource-by-Ntarget skip table or per-source dense update.
          broad_sources=broad_sources+1
          nactive=0
          do g=1,ng
            if(source_column(g)/=(0d0,0d0))nactive=nactive+1
          enddo
          allocate(active_rows(nactive),pair_values(nactive));k=0
          do g=1,ng
            if(source_column(g)==(0d0,0d0))cycle
            k=k+1;active_rows(k)=g
          enddo
          do k=1,ncandidate
            j=selected(k);broad_kept(j)=broad_kept(j)+1
            pair_values=abs(source_column(active_rows))*abs(target(active_rows,j,1))
            pair_norms(1,k)=sum(pair_values**2)
            pair_norms(2,k)=sum(pair_values)
          enddo
          deallocate(active_rows,pair_values)
          op%pair_product_points=op%pair_product_points+int(ncandidate,int64)*int(nactive,int64)
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
            compact_used,status,compact_pairs,compact_points,fft_cost_limit=sr_cost)
          call collective_bad_status()
          if(status/=0)then
            call sr_destroy(sr)
            call spatial_local_destroy(compact_plan)
            return
          endif
          if(compact_used)then
            do j=1,nselected
              action(:,selected(j),1)=action(:,selected(j),1)+compact_action(:,j)
            enddo
            deallocate(compact_action)
            wf_pairs=wf_pairs+compact_pairs
            op%local_pairs=op%local_pairs+compact_pairs
            op%local_points=op%local_points+compact_points
            if(present(comm_o))call collective_bad()
            if(bad/=0)then
              call sr_destroy(sr)
              call spatial_local_destroy(compact_plan)
              return
            endif
            cycle
          endif
          deallocate(compact_action)
        endif
        if(sr%ready)then
          allocate(compact_action(ng,nselected))
          call sr_apply(sr,comm_r,source_column,target(:,selected(:nselected),1),compact_action)
          do j=1,nselected
            action(:,selected(j),1)=action(:,selected(j),1)+compact_action(:,j)
          enddo
          deallocate(compact_action)
          sr_pairs=sr_pairs+nselected
          op%local_pairs=op%local_pairs+nselected
          op%local_points=op%local_points+int(nselected,int64)*int(sr%peers,int64)*product(int(sr%length,int64))
          cycle
        endif
        op%global_pairs=op%global_pairs+nselected
        do first=1,nselected,4
          nb=min(size(density,2),nselected-first+1)
          ! Pad the final batch to reuse a bounded set of FFT plans/work arrays.
          if(nb<size(density,2))density(:,nb+1:)=0d0
!$omp parallel do default(none) schedule(static) private(j) &
!$omp shared(nb,ng,density,source_column,target,selected,first) &
!$omp num_threads(min(nb,point_workers)) if(ng*nb>=262144.and.point_workers>1)
          do j=1,nb
            density(:,j)=conjg(source_column)*target(:,selected(first+j-1),1)
          enddo
!$omp end parallel do
          call mesh_transform(n,dims,coords,comm,density,spectrum,-1,status)
          bad=status
          if(present(comm_o))call comm_get_max(bad,comm_r)
          if(bad/=0)exit
!$omp parallel do default(none) schedule(static) private(j) &
!$omp shared(nb,ng,spectrum,multiplier) &
!$omp num_threads(min(nb,point_workers)) if(ng*nb>=262144.and.point_workers>1)
          do j=1,nb
            spectrum(:,j)=spectrum(:,j)*multiplier
          enddo
!$omp end parallel do
          ! Inverse mesh_transform already includes 1/product(n).
          call mesh_transform(n,dims,coords,comm,spectrum,density,1,status)
          bad=status
          if(present(comm_o))call comm_get_max(bad,comm_r)
          if(bad/=0)exit
!$omp parallel do default(none) schedule(static) private(j) &
!$omp shared(nb,ng,action,selected,first,source_column,density) &
!$omp num_threads(min(nb,point_workers)) if(ng*nb>=262144.and.point_workers>1)
          do j=1,nb
            action(:,selected(first+j-1),1)=action(:,selected(first+j-1),1)-source_column*density(:,j)
          enddo
!$omp end parallel do
        enddo
        if(present(comm_o))call collective_bad()
        if(bad/=0)then
          call sr_destroy(sr)
          call spatial_local_destroy(compact_plan)
          return
        endif
      enddo
    enddo
    if(op%sr_tolerance>0d0.and.sr%rank==0) &
      write(*,'(a,2i18)')'EXX_SR WF-local/neighborhood pairs: ',wf_pairs,sr_pairs
    call sr_destroy(sr)
    call spatial_local_destroy(compact_plan,status)
    call collective_bad_status()
    if(status/=0)return
    if(op%screen_mode/=0)then
      omitted=omitted+budget*real(broad_sources-broad_kept,8)
      call comm_summation(real(op%pair_product_points,8),point_total,comm_r)
      pair_stats=[real(op%pair_products,8),real(op%pair_catalog_entries,8),point_total];global_pair_stats=pair_stats
      if(present(comm_o))call comm_summation(pair_stats,global_pair_stats,3,comm_o)
      op%pair_products=int(global_pair_stats(1),int64);op%pair_catalog_entries=int(global_pair_stats(2),int64)
      op%pair_product_points=int(global_pair_stats(3),int64)
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
      implicit none
      real(8),intent(inout) :: value
      integer,intent(in) :: group
      real(8) :: result(1)
      call comm_get_max([value],result,1,group)
      value=result(1)
    end subroutine
    subroutine collective_bad_status()
      implicit none
      ! Backends use signed error codes; MPI_MAX must see a Boolean failure.
      bad=merge(1,0,status/=0)
      call collective_bad()
    end subroutine
    subroutine collective_bad()
      implicit none
      call comm_get_max(bad,comm_r)
      if(present(comm_o))call comm_get_max(bad,comm_o)
      ! Preserve a nonzero return on every peer, even after a successful local FFT.
      if(bad/=0)status=1
    end subroutine
  end subroutine


  pure logical function finite_real_1d(values) result(finite)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in) :: values(:)
    real(8) :: value
    integer :: i
    finite=.false.
    do i=1,size(values,1)
      value=values(i)
      if(.not.ieee_is_finite(value))return
    enddo
    finite=.true.
  end function

  pure logical function finite_real_2d(values) result(finite)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in) :: values(:,:)
    real(8) :: value
    integer :: i,j
    finite=.false.
    do j=1,size(values,2)
      do i=1,size(values,1)
        value=values(i,j)
        if(.not.ieee_is_finite(value))return
      enddo
    enddo
    finite=.true.
  end function

  pure logical function finite_real_3d(values) result(finite)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in) :: values(:,:,:)
    real(8) :: value
    integer :: i,j,k
    finite=.false.
    do k=1,size(values,3)
      do j=1,size(values,2)
        do i=1,size(values,1)
          value=values(i,j,k)
          if(.not.ieee_is_finite(value))return
        enddo
      enddo
    enddo
    finite=.true.
  end function
end module
