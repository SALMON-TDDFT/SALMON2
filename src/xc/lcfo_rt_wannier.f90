! Initial distributed MV gauge, then polar transport of U against prior WFs.
! Positive spherical radius is an explicit source-mask approximation, not a
! variational energy functional. Global density/Hartree are never masked.
module lcfo_rt_wannier
  use iso_fortran_env, only: int32,int64
  use lcfo_mlwf_links, only: lcfo_initial_links
  use lcfo_rt_basis, only: lcfo_direct_wf,lcfo_basis,lcfo_counts,lcfo_origins,lcfo_grid, &
    lcfo_core,lcfo_buffer,lcfo_rank,lcfo_comm,lcfo_h, &
    lcfo_dv
  use communication, only: comm_summation,comm_bcast
  use exx_wannier_gauge, only: gauge_minimize_gamma_inplace
  use lcfo_seed, only: lcfo_seed_gamma
  use lcfo_dist_rows, only: s_lcfo_halo
  use lcfo_dist_rows, only: s_lcfo_column_halo,lcfo_column_halo_init,lcfo_column_halo_get
  use lcfo_dist_dense, only: lcfo_distributed_polar
  use lcfo_wf_support, only: s_lcfo_wf_plan,lcfo_wf_plan_init,lcfo_wf_reconstruct,lcfo_wf_total_norm
  use lcfo_wf_support, only: s_lcfo_wf_kernel,lcfo_wf_kernel_init,lcfo_wf_kernel_apply
  use lcfo_wf_support, only: lcfo_wf_sphere_norm
  use salmon_global, only: exx_mlwf_maxiter,exx_mlwf_tolerance,hse_lcfo_wf_radius, &
    yn_hse_wannier,hse_lcfo_u_interval,yn_hse_lcfo_seed_distributed
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: lcfo_mlwf_enabled,lcfo_mlwf_configure,lcfo_mlwf_source,lcfo_mlwf_stage,lcfo_mlwf_accept_cached,lcfo_mlwf_track
  public :: lcfo_mlwf_rebase
  logical,save :: lcfo_mlwf_enabled=.false.
  real(8),save :: radius=0d0
  complex(8),allocatable,save :: rotation(:,:)
  real(8),allocatable,save :: centers(:,:)
  logical,allocatable,save :: protected(:)
  integer,save :: uses=0,u_interval=1,physical_step=0
  complex(8),allocatable,save :: step_rotation(:,:),frame_rotation(:,:)
  complex(8),allocatable,save :: transport_anchor(:,:),step_anchor(:,:)
  logical,save :: frame_transported=.false.
  complex(8),allocatable,save :: previous_wf(:,:),current_frame(:,:),step_reference(:,:)
  logical,save :: step_active=.false.
  type(s_lcfo_wf_plan),save :: fragment_support,core_support
  type(s_lcfo_column_halo),save :: source_halo
  type(s_lcfo_wf_kernel),save :: fragment_kernel,core_kernel
  complex(8),allocatable,save :: core_gram(:,:)
contains
  subroutine lcfo_mlwf_configure()
    implicit none
    integer :: bad
    bad=0
    if(lcfo_rank==0)then
      lcfo_mlwf_enabled=yn_hse_wannier=='y'
      radius=hse_lcfo_wf_radius
      u_interval=hse_lcfo_u_interval
      if(u_interval<1.or.(u_interval>1.and..not.lcfo_mlwf_enabled))bad=1
      if(.not.ieee_is_finite(radius).or.radius<0d0)bad=1
      if(radius>0d0.and..not.lcfo_mlwf_enabled)bad=1
      if(lcfo_direct_wf.and..not.lcfo_mlwf_enabled)bad=1
    endif
    call comm_bcast(bad,lcfo_comm,0)
    if(bad/=0)error stop 'LCFO MLWF: invalid radius or U interval (requires MLWF, finite radius >=0, integer interval >=1)'
    call comm_bcast(lcfo_mlwf_enabled,lcfo_comm,0)
    call comm_bcast(radius,lcfo_comm,0)
    call comm_bcast(u_interval,lcfo_comm,0)
    if(lcfo_mlwf_enabled.and.lcfo_rank==0)then
      write(*,*) 'LCFO MLWF: initial U with polar temporal transport; global periodic 3D radius (0=full) =',radius
      if(radius>0d0.and.radius<.5d0*sqrt(sum((lcfo_grid*lcfo_h)**2))) &
        write(*,*) 'LCFO MLWF source-mask approximation: no renormalization; trace energy is diagnostic'
    endif
  end subroutine

  subroutine initialize_rotation(coeff)
    implicit none
    complex(8),intent(in) :: coeff(:,:)
    complex(8),allocatable :: grid(:,:),raw(:,:,:,:),u(:,:,:),overlap(:,:)
    complex(8),allocatable :: moment_local(:,:),moment(:,:)
    real(8),allocatable :: position(:,:),norm_local(:),norms(:)
    real(8) :: b(3,6),weights(6),pi,length(3),wf_spread,gradient,delta,unitary_error
    integer :: no,ng,g,x,y,z,a,j,status,iterations,seed_status,iu
    integer(int64) :: link_scratch
    no=size(coeff,2);ng=product(lcfo_core);pi=acos(-1d0);length=lcfo_grid*lcfo_h
    grid=matmul(lcfo_basis,coeff)
    allocate(position(3,ng));g=0
    do z=0,lcfo_core(3)-1;do y=0,lcfo_core(2)-1;do x=0,lcfo_core(1)-1
      g=g+1;position(:,g)=(lcfo_origins(:,lcfo_rank+1)+[x,y,z])*lcfo_h
    enddo;enddo;enddo
    allocate(u(no,no,1))
    b=0d0
    do a=1,3
      delta=2*pi/length(a);b(a,a)=delta;b(a,a+3)=-delta
      weights(a)=1d0/(2*delta**2);weights(a+3)=weights(a)
    enddo
    ! Stream the unchanged coefficient snapshot while gathering QR layout.
    ! No root full_coeff copy survives alongside the QR matrix.
    iu=0
    if(lcfo_rank==0)then
      open(newunit=iu,file='lcfo_mlwf_initial.bin',access='stream',form='unformatted',status='replace')
      write(iu)int([16909060,2,no,sum(lcfo_counts),lcfo_grid],int32),lcfo_h
    endif
    call lcfo_seed_gamma(coeff,lcfo_counts,lcfo_comm,u(:,:,1),seed_status,snapshot_unit=iu, &
      distributed=yn_hse_lcfo_seed_distributed=='y')
    if(lcfo_rank==0)then
      close(iu)
      if(seed_status/=0)then
        u=0d0
        do j=1,no;u(j,j,1)=1d0;enddo
      endif
    endif
    call lcfo_initial_links(grid,position,length,lcfo_dv,lcfo_comm,raw,scratch_elements=link_scratch)
    if(lcfo_rank==0)write(*,'(a,2i18)')'LCFO initial links root/scratch complex elements: ' , &
      size(raw,kind=int64),link_scratch
    if(lcfo_rank==0)then
      open(newunit=iu,file='lcfo_mlwf_links.bin',access='stream',form='unformatted',status='replace')
      write(iu)int([16909060,1,no],int32),b,weights,u,raw
      close(iu)
      call gauge_minimize_gamma_inplace(u,raw,b,weights,exx_mlwf_maxiter,exx_mlwf_tolerance, &
                          wf_spread,gradient,iterations,status)
      write(*,'(a,3i7,2es17.8)') 'LCFO MLWF initial iterations/status/seed/spread/gradient:', &
        iterations,status,seed_status,wf_spread,gradient
    endif
    deallocate(raw)
    call comm_bcast(status,lcfo_comm,0)
    if(status/=0.and.radius>0d0.and.radius<.5d0*sqrt(sum(length**2))) &
      error stop 'LCFO MLWF: initial localization unconverged; support comparison requires converged U'
    call comm_bcast(u,lcfo_comm,0)
    rotation=u(:,:,1)
    deallocate(u)
    overlap=matmul(conjg(transpose(rotation)),rotation)
    do j=1,no;overlap(j,j)=overlap(j,j)-1d0;enddo
    unitary_error=maxval(abs(overlap))
    if(.not.all(ieee_is_finite(real(rotation))).or..not.all(ieee_is_finite(aimag(rotation))).or. &
       unitary_error>1d-10)error stop 'LCFO MLWF: invalid initial U'
    deallocate(overlap)
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
    call report_radius_coverage(grid,position,length,norms)
    previous_wf=matmul(coeff,rotation);current_frame=previous_wf;frame_rotation=rotation
    transport_anchor=current_frame;frame_transported=.true.
    if(lcfo_rank==0)then
      open(newunit=iu,file='lcfo_mlwf_initial.bin',access='stream',form='unformatted',status='old',position='append')
      write(iu)rotation,centers,norms
      write(iu)int(merge(1,0,protected),int32),int([iterations,status],int32),wf_spread,gradient
      close(iu)
    endif
    if(lcfo_rank==0)write(*,'(a,es12.4,a,i6,a,es12.4)')'LCFO MLWF U error ',unitary_error, &
      ' protected factors ',count(protected),' minimum xyz center reliability ',minval(abs(moment)/spread(norms,1,3))
  end subroutine

  subroutine report_radius_coverage(grid,position,length,norms)
    implicit none
    complex(8),intent(in) :: grid(:,:)
    real(8),intent(in) :: position(:,:),length(3),norms(:)
    real(8),allocatable :: local_norm(:),inside(:),fraction(:)
    integer :: n,j,iu,worst(1),below
    n=size(norms)
    if(radius==0d0)then
      ! Full support already has a globally reduced norm; no grid scan or MPI needed.
      if(lcfo_rank/=0)return
      allocate(inside(n));inside=norms
    else
      allocate(local_norm(n),inside(n))
      call lcfo_wf_sphere_norm(grid,position,centers,length,lcfo_dv,radius,local_norm)
      call comm_summation(local_norm,inside,n,lcfo_comm)
      if(lcfo_rank/=0)return
    endif
    allocate(fraction(n));fraction=inside/norms
    open(newunit=iu,file='lcfo_mlwf_radius.dat',status='replace')
    write(iu,'(a,es24.16)')'# Initial WF geometric sphere coverage; radius_bohr (0=full): ',radius
    write(iu,'(a)')'# wf  total_norm  sphere_norm  sphere_fraction  protected_uncut'
    do j=1,n
      write(iu,'(i10,3es25.16,i4)')j,norms(j),inside(j),fraction(j),merge(1,0,protected(j))
    enddo
    close(iu)
    worst=minloc(fraction);below=count(fraction<.999d0)
    write(*,'(a,es16.8,a,i8)')'LCFO MLWF initial minimum sphere norm fraction: ', &
      fraction(worst(1)),' WF ',worst(1)
    if(below>0)then
      write(*,'(a,i8,a,i8)')'WARNING LCFO MLWF radius: sphere norm below 99.9% for ',below,' of ',n
      write(*,'(a,es16.8,a,i8)')'WARNING LCFO MLWF radius unchanged (bohr): ',radius, &
        '; protected WFs remain uncut: ',count(protected)
    endif
  end subroutine

  subroutine lcfo_mlwf_track(coeff)
    implicit none
    complex(8),intent(in) :: coeff(:,:)
    complex(8),allocatable :: new_u(:,:,:)
    real(8) :: minimum_overlap
    integer :: no,status
    if(.not.allocated(rotation))then
      call initialize_rotation(coeff)
      return
    endif
    no=size(coeff,2);frame_transported=.false.
    if(allocated(current_frame).and.(physical_step<=1.or.mod(physical_step-1,u_interval)==0))then
      allocate(new_u(no,no,1))
      if(u_interval>1)then
        ! Do not turn held-U dephasing into permanent motion of the gauge anchor.
        if(step_active)then
          call lcfo_distributed_polar(coeff,step_anchor,lcfo_comm,new_u(:,:,1),minimum_overlap,status)
        else
          call lcfo_distributed_polar(coeff,transport_anchor,lcfo_comm,new_u(:,:,1),minimum_overlap,status)
        endif
      else if(step_active)then
        call lcfo_distributed_polar(coeff,step_reference,lcfo_comm,new_u(:,:,1),minimum_overlap,status)
      else
        call lcfo_distributed_polar(coeff,previous_wf,lcfo_comm,new_u(:,:,1),minimum_overlap,status)
      endif
      if(status/=0)error stop 'LCFO MLWF: occupied subspace overlap lost; relocalization required'
      rotation=new_u(:,:,1);frame_transported=.true.
      if(lcfo_rank==0)then
        write(*,'(a,es14.6)')'LCFO MLWF transport minimum overlap ',minimum_overlap
        write(*,'(a,2i8)')'LCFO MLWF U refreshed step/interval:',physical_step,u_interval
      endif
    endif
    if(physical_step>1.and.mod(physical_step-1,u_interval)/=0.and.lcfo_rank==0) &
      write(*,'(a,2i8)')'LCFO MLWF U held step/interval:',physical_step,u_interval
    current_frame=matmul(coeff,rotation)
    frame_rotation=rotation
    if(frame_transported)transport_anchor=current_frame
    previous_wf=current_frame
  end subroutine

  subroutine lcfo_mlwf_source(coeff,basis,selected,source,halo)
    implicit none
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
      if(.not.source_halo%ready)then
        call lcfo_column_halo_init(source_halo,halo,fragment_support%columns,no)
        write(*,'(a,3i12)')'LCFO WF halo rank/selected/full values:', &
          lcfo_rank,sum(source_halo%recv_count),halo%nselected*no
      endif
      call lcfo_column_halo_get(source_halo,halo,current_frame,near_frame)
    else
      if(size(lcfo_counts)/=1)error stop 'LCFO MLWF: distributed source requires halo plan'
      near_frame=current_frame(selected,fragment_support%columns)
    endif
    if(.not.fragment_kernel%ready)then
      call lcfo_wf_kernel_init(fragment_kernel,fragment_support,basis)
      write(*,'(a,i6,2i14)')'LCFO WF local kernel rank/blocks/products:', &
        lcfo_rank,size(fragment_kernel%blocks),fragment_kernel%products
    endif
    call lcfo_wf_kernel_apply(fragment_kernel,near_frame,fragment_wf)
    allocate(source(size(basis,1),size(fragment_wf,2),1));source(:,:,1)=fragment_wf
    local_loss=0d0
    if(cut)then
      local_loss(2)=lcfo_wf_total_norm(core_gram,current_frame)
      if(core_support%masked)then
        if(.not.core_kernel%ready)call lcfo_wf_kernel_init(core_kernel,core_support,lcfo_basis)
        call lcfo_wf_kernel_apply(core_kernel,current_frame(:,core_support%columns),core_wf)
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
  subroutine lcfo_mlwf_rebase(coeff)
    implicit none
    ! Change the native occupied-orbital coordinates, not the physical WF frame.
    complex(8),intent(inout) :: coeff(:,:)
    integer :: j
    if(.not.lcfo_mlwf_enabled.or..not.allocated(rotation).or.step_active) &
      error stop 'LCFO direct WF: rebase requires an accepted MLWF frame'
    coeff=matmul(coeff,rotation)
    current_frame=coeff;previous_wf=coeff
    if(frame_transported)transport_anchor=coeff
    rotation=0d0
    do j=1,size(rotation,1);rotation(j,j)=1d0;enddo
    frame_rotation=rotation
  end subroutine

  subroutine lcfo_mlwf_stage(stage)
    implicit none
    integer,intent(in) :: stage
    if(.not.lcfo_mlwf_enabled)return
    select case(stage)
    case(0)
      physical_step=physical_step+1
      step_reference=previous_wf;step_rotation=rotation;step_anchor=transport_anchor;step_active=.true.
    case(1)
      ! Predictor and corrected orbitals must use the same accepted reference.
    case(2)
      previous_wf=step_reference;rotation=step_rotation;transport_anchor=step_anchor;step_active=.false.
      if(allocated(step_anchor))deallocate(step_anchor)
      if(allocated(step_rotation))deallocate(step_rotation)
      if(allocated(step_reference))deallocate(step_reference)
    case default
      error stop 'LCFO MLWF: invalid Taylor stage'
    end select
  end subroutine

  subroutine lcfo_mlwf_accept_cached()
    implicit none
    if(.not.lcfo_mlwf_enabled)return
    ! A cache hit after predictor rollback still accepts the frame corresponding
    ! to those exact coefficients, without repeating FFT exchange or localization.
    if(allocated(current_frame))then
      previous_wf=current_frame
      rotation=frame_rotation
      if(frame_transported)transport_anchor=current_frame
    endif
  end subroutine
end module
