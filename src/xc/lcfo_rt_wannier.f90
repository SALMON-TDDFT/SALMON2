! Initial distributed MV gauge, then polar transport of U against prior WFs.
! Positive spherical radius is an explicit source-mask approximation, not a
! variational energy functional. Global density/Hartree are never masked.
module lcfo_rt_wannier
  use iso_fortran_env, only: int32
  use lcfo_rt_basis
  use communication, only: comm_summation,comm_bcast
  use hse_wannier_gauge, only: gauge_seed,gauge_minimize_gamma
  use lcfo_dist_rows, only: s_lcfo_halo,lcfo_gather_root,lcfo_halo_get
  use lcfo_dist_dense, only: lcfo_distributed_polar
  use lcfo_wf_support, only: s_lcfo_wf_plan,lcfo_wf_plan_init,lcfo_wf_reconstruct,lcfo_wf_total_norm
  use salmon_global, only: hse_mlwf_maxiter,hse_mlwf_tolerance
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: lcfo_mlwf_enabled,lcfo_mlwf_configure,lcfo_mlwf_source,lcfo_mlwf_stage,lcfo_mlwf_accept_cached,lcfo_mlwf_track
  logical,save :: lcfo_mlwf_enabled=.false.
  real(8),save :: radius=0d0
  complex(8),allocatable,save :: rotation(:,:)
  real(8),allocatable,save :: centers(:,:)
  logical,allocatable,save :: protected(:)
  integer,save :: uses=0
  complex(8),allocatable,save :: previous_wf(:,:),current_frame(:,:),step_reference(:,:)
  logical,save :: step_active=.false.
  type(s_lcfo_wf_plan),save :: fragment_support,core_support
  complex(8),allocatable,save :: core_gram(:,:)
contains
  subroutine lcfo_mlwf_configure()
    character(64) :: value
    integer :: status,ios,bad
    bad=0
    if(lcfo_rank==0)then
      call get_environment_variable('SALMON_LCFO_RT_MLWF',value,status=status)
      lcfo_mlwf_enabled=status==0.and.trim(value)=='1'
      call get_environment_variable('SALMON_LCFO_RT_RADIUS',value,status=status)
      if(status==0.and.len_trim(value)>0)then
        read(value,*,iostat=ios)radius
        if(ios/=0)bad=1
      endif
      if(.not.ieee_is_finite(radius).or.radius<0d0)bad=1
      if(radius>0d0.and..not.lcfo_mlwf_enabled)bad=1
    endif
    call comm_bcast(bad,lcfo_comm,0)
    if(bad/=0)error stop 'LCFO MLWF: finite nonnegative radius requires SALMON_LCFO_RT_MLWF=1'
    call comm_bcast(lcfo_mlwf_enabled,lcfo_comm,0)
    call comm_bcast(radius,lcfo_comm,0)
    if(lcfo_mlwf_enabled.and.lcfo_rank==0)then
      write(*,*) 'LCFO MLWF: initial U with polar temporal transport; global periodic 3D radius (0=full) =',radius
      if(radius>0d0.and.radius<.5d0*sqrt(sum((lcfo_grid*lcfo_h)**2))) &
        write(*,*) 'LCFO MLWF source-mask approximation: no renormalization; trace energy is diagnostic'
    endif
  end subroutine

  subroutine initialize_rotation(coeff)
    complex(8),intent(in) :: coeff(:,:)
    complex(8),allocatable :: grid(:,:),shifted(:,:),raw_local(:,:,:,:),raw(:,:,:,:),u(:,:,:),overlap(:,:)
    complex(8),allocatable :: moment_local(:,:),moment(:,:),full_coeff(:,:)
    real(8),allocatable :: position(:,:),norm_local(:),norms(:),seed_position(:,:)
    real(8) :: b(3,6),weights(6),pi,length(3),wf_spread,gradient,delta,unitary_error
    integer :: no,ng,g,x,y,z,a,j,status,iterations,neighbors(6,1),seed_status,iu
    no=size(coeff,2);ng=product(lcfo_core);pi=acos(-1d0);length=lcfo_grid*lcfo_h
    grid=matmul(lcfo_basis,coeff)
    allocate(position(3,ng));g=0
    do z=0,lcfo_core(3)-1;do y=0,lcfo_core(2)-1;do x=0,lcfo_core(1)-1
      g=g+1;position(:,g)=(lcfo_origins(:,lcfo_rank+1)+[x,y,z])*lcfo_h
    enddo;enddo;enddo
    allocate(raw_local(no,no,6,1),raw(no,no,6,1),shifted(ng,no),u(no,no,1))
    b=0d0;neighbors=1
    do a=1,3
      delta=2*pi/length(a);b(a,a)=delta;b(a,a+3)=-delta
      weights(a)=1d0/(2*delta**2);weights(a+3)=weights(a)
      do j=1,no
        shifted(:,j)=grid(:,j)*exp(cmplx(0d0,-delta*position(a,:),8))
      enddo
      raw_local(:,:,a,1)=matmul(conjg(transpose(grid)),shifted)*lcfo_dv
      raw_local(:,:,a+3,1)=conjg(transpose(raw_local(:,:,a,1)))
    enddo
    call comm_summation(raw_local,raw,size(raw),lcfo_comm)
    call lcfo_gather_root(coeff,lcfo_counts,lcfo_comm,full_coeff)
    if(lcfo_rank==0)then
      ! LCFO functions have disjoint compact cores. Pivoted coefficient rows
      ! supply localized trial functions without gathering the global grid.
      allocate(seed_position(3,size(full_coeff,1)));seed_position=0d0
      call gauge_seed(reshape(full_coeff,[size(full_coeff,1),no,1]),seed_position,reshape([0d0,0d0,0d0],[3,1]), &
                      u,seed_status)
      if(seed_status/=0)then
        u=0d0
        do j=1,no;u(j,j,1)=1d0;enddo
      endif
      open(newunit=iu,file='lcfo_mlwf_links.bin',access='stream',form='unformatted',status='replace')
      write(iu)int([16909060,1,no],int32),b,weights,u,raw
      close(iu)
      call gauge_minimize_gamma(u,raw,b,weights,hse_mlwf_maxiter,hse_mlwf_tolerance, &
                          wf_spread,gradient,iterations,status)
      write(*,'(a,3i7,2es17.8)') 'LCFO MLWF initial iterations/status/seed/spread/gradient:', &
        iterations,status,seed_status,wf_spread,gradient
    endif
    call comm_bcast(status,lcfo_comm,0)
    if(status/=0.and.radius>0d0.and.radius<.5d0*sqrt(sum(length**2))) &
      error stop 'LCFO MLWF: initial localization unconverged; support comparison requires converged U'
    call comm_bcast(u,lcfo_comm,0)
    rotation=u(:,:,1)
    overlap=matmul(conjg(transpose(rotation)),rotation)
    do j=1,no;overlap(j,j)=overlap(j,j)-1d0;enddo
    unitary_error=maxval(abs(overlap))
    if(.not.all(ieee_is_finite(real(rotation))).or..not.all(ieee_is_finite(aimag(rotation))).or. &
       unitary_error>1d-10)error stop 'LCFO MLWF: invalid initial U'
    grid=matmul(grid,rotation)
    allocate(moment_local(3,no),moment(3,no),norm_local(no),norms(no),centers(3,no),protected(no))
    do j=1,no
      norm_local(j)=sum(abs(grid(:,j))**2)*lcfo_dv
      do a=1,3
        moment_local(a,j)=sum(abs(grid(:,j))**2*exp(cmplx(0d0,2*pi*position(a,:)/length(a),8)))*lcfo_dv
      enddo
    enddo
    call comm_summation(norm_local,norms,no,lcfo_comm)
    call comm_summation(moment_local,moment,size(moment),lcfo_comm)
    if(any(norms<=0d0))error stop 'LCFO MLWF: empty occupied orbital'
    do a=1,3
      centers(a,:)=modulo(atan2(aimag(moment(a,:)),real(moment(a,:)))*length(a)/(2*pi),length(a))
    enddo
    ! A poorly defined center on any axis makes a spherical cut unreliable.
    protected=any(abs(moment)/spread(norms,1,3)<.1d0,dim=1)
    previous_wf=matmul(coeff,rotation);current_frame=previous_wf
    if(lcfo_rank==0)then
      open(newunit=iu,file='lcfo_mlwf_initial.bin',access='stream',form='unformatted',status='replace')
      write(iu)int([16909060,2,no,size(full_coeff,1),lcfo_grid],int32),lcfo_h,full_coeff,rotation,centers,norms
      write(iu)int(merge(1,0,protected),int32),int([iterations,status],int32),wf_spread,gradient
      close(iu)
    endif
    if(lcfo_rank==0)write(*,'(a,es12.4,a,i6,a,es12.4)')'LCFO MLWF U error ',unitary_error, &
      ' protected factors ',count(protected),' minimum xyz center reliability ',minval(abs(moment)/spread(norms,1,3))
  end subroutine

  subroutine lcfo_mlwf_track(coeff)
    complex(8),intent(in) :: coeff(:,:)
    complex(8),allocatable :: new_u(:,:,:)
    real(8) :: minimum_overlap
    integer :: no,status
    if(.not.allocated(rotation))then
      call initialize_rotation(coeff)
      return
    endif
    no=size(coeff,2)
    if(allocated(current_frame))then
      allocate(new_u(no,no,1))
      if(step_active)then
        call lcfo_distributed_polar(coeff,step_reference,lcfo_comm,new_u(:,:,1),minimum_overlap,status)
      else
        call lcfo_distributed_polar(coeff,previous_wf,lcfo_comm,new_u(:,:,1),minimum_overlap,status)
      endif
      if(status/=0)error stop 'LCFO MLWF: occupied subspace overlap lost; relocalization required'
      rotation=new_u(:,:,1)
      if(lcfo_rank==0)write(*,'(a,es14.6)')'LCFO MLWF transport minimum overlap ',minimum_overlap
    endif
    current_frame=matmul(coeff,rotation)
    previous_wf=current_frame
  end subroutine

  subroutine lcfo_mlwf_source(coeff,basis,selected,source,halo)
    complex(8),intent(in) :: coeff(:,:),basis(:,:)
    integer,intent(in) :: selected(:)
    complex(8),allocatable,intent(out) :: source(:,:,:)
    type(s_lcfo_halo),intent(in),optional :: halo
    complex(8),allocatable :: core_wf(:,:),near_frame(:,:),fragment_wf(:,:)
    real(8) :: length(3),local_loss(2),loss(2),mask_radius
    real(8),allocatable :: positions(:,:)
    integer :: ns(3),p(3),point(3),g,x,y,z,no,active_sources,active_sum,ncore
    logical :: cut
    call lcfo_mlwf_track(coeff)
    no=size(coeff,2);length=lcfo_grid*lcfo_h
    cut=radius>0d0.and.radius<.5d0*sqrt(sum(length**2))
    if(.not.fragment_support%ready)then
      ! Centers and support geometry are fixed during polar transport.
      ns=lcfo_core+2*lcfo_buffer
      allocate(positions(3,size(basis,1)));positions=0d0;mask_radius=0d0
      if(cut)then
        if(size(basis,1)/=product(ns))error stop 'LCFO WF support: fragment grid mismatch'
        mask_radius=radius;g=0
        do z=0,ns(3)-1;do y=0,ns(2)-1;do x=0,ns(1)-1
          g=g+1;p=[x,y,z]
          where(p>=lcfo_core+lcfo_buffer)p=p-ns
          point=modulo(lcfo_origins(:,lcfo_rank+1)+p,lcfo_grid)
          positions(:,g)=point*lcfo_h
        enddo;enddo;enddo
      endif
      call lcfo_wf_plan_init(fragment_support,positions,centers,length,mask_radius,protected)
      deallocate(positions);ncore=no
      if(cut)then
        allocate(positions(3,product(lcfo_core)));g=0
        do z=0,lcfo_core(3)-1;do y=0,lcfo_core(2)-1;do x=0,lcfo_core(1)-1
          g=g+1;positions(:,g)=(lcfo_origins(:,lcfo_rank+1)+[x,y,z])*lcfo_h
        enddo;enddo;enddo
        call lcfo_wf_plan_init(core_support,positions,centers,length,radius,protected)
        core_gram=matmul(conjg(transpose(lcfo_basis)),lcfo_basis)
        ncore=size(core_support%columns)
      endif
      write(*,'(a,4i8)')'LCFO WF reconstruction rank/fragment/core/total columns:', &
        lcfo_rank,size(fragment_support%columns),ncore,no
    endif
    if(present(halo))then
      call lcfo_halo_get(halo,current_frame,near_frame)
    else
      if(size(lcfo_counts)/=1)error stop 'LCFO MLWF: distributed source requires halo plan'
      near_frame=current_frame(selected,:)
    endif
    call lcfo_wf_reconstruct(fragment_support,basis,near_frame,fragment_wf)
    allocate(source(size(basis,1),size(fragment_wf,2),1));source(:,:,1)=fragment_wf
    local_loss=0d0
    if(cut)then
      local_loss(2)=lcfo_wf_total_norm(core_gram,current_frame)
      if(core_support%masked)then
        call lcfo_wf_reconstruct(core_support,lcfo_basis,current_frame,core_wf)
        local_loss(1)=max(0d0,local_loss(2)-sum(abs(core_wf)**2))
      endif
    endif
    call comm_summation(local_loss,loss,2,lcfo_comm)
    active_sources=count(any(abs(source(:,:,1))>0d0,dim=1))
    call comm_summation(active_sources,active_sum,lcfo_comm)
    uses=uses+1
    if(lcfo_rank==0)write(*,'(a,2i9)') 'LCFO MLWF active/possible fragment sources: ', &
      active_sum,no*size(lcfo_counts)
    if(lcfo_rank==0)write(*,'(a,i8,a,es14.6)')'LCFO MLWF reuse ',uses, &
      ' discarded global source norm fraction ',loss(1)/max(loss(2),tiny(1d0))
  end subroutine
  subroutine lcfo_mlwf_stage(stage)
    integer,intent(in) :: stage
    if(.not.lcfo_mlwf_enabled)return
    select case(stage)
    case(0)
      step_reference=previous_wf;step_active=.true.
    case(1)
      ! Predictor and corrected orbitals must use the same accepted reference.
    case(2)
      previous_wf=step_reference;step_active=.false.
      if(allocated(step_reference))deallocate(step_reference)
    case default
      error stop 'LCFO MLWF: invalid Taylor stage'
    end select
  end subroutine

  subroutine lcfo_mlwf_accept_cached()
    if(.not.lcfo_mlwf_enabled)return
    ! A cache hit after predictor rollback still accepts the frame corresponding
    ! to those exact coefficients, without repeating FFT exchange or localization.
    if(allocated(current_frame))previous_wf=current_frame
  end subroutine
end module
