#include "config.h"
! SALMON adapter: legacy full-grid/orbital layout distributes k points.
! Gamma HSE SCF and DC-initialized hybrid mesh RT support spatial y/z FFTW pencils.
! Its source and ACE factors retain only local grid rows; overlaps are reduced.
module exx_native
  use exx_functional, only: exchange_fraction
  use iso_fortran_env, only: int64
  use lcfo_rt_basis, only: lcfo_rt_active
  use exx_lcfo_rt, only: lcfo_exx_refresh,lcfo_exx_add_action,lcfo_exx_stage
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use structures
  use plusU_global, only: PLUS_U_ON
  use hse_exchange
  use exx_ace
  use exx_orbitals, only: orbital_ace_build,orbital_ace_apply,orbital_layout,orbital_rotate,orbital_hermitian_action
  use exx_adaptive_support, only: adaptive_source_mask
  use exx_spatial
  use exx_wannier
  use exx_symmetry
  use sym_sub, only: use_symmetry,SymMatA,SymMatB
  use communication, only: comm_summation,comm_alltoall,comm_get_max
  use salmon_global, only: xc,yn_periodic,yn_spinorbit,yn_jm,yn_dc,yn_md,yn_symmetrized_stencil,propagator,num_kgrid,hse_omega, &
    pbeh_coulomb_radius,theory,yn_conventional_from_dcdft,num_rgrid, &
    yn_hse_wannier,exx_mlwf_interval,exx_mlwf_maxiter,exx_mlwf_tolerance,exx_mlwf_radius,exx_mlwf_norm_fraction,exx_local_fft, &
    yn_exx_dc_mlwf,exx_pre_scf_active,exx_ace_support,exx_pair_screening,exx_pair_tolerance,hse_block_rows, &
    yn_hse_profile,hse_fft_layout,yn_hse_eigen_diagnostic,yn_hse_solver_diagnostic,yn_hse_wannier_snapshot
  implicit none
  private
  public :: exx_export_snapshot,exx_eigen_diagnostic_enabled,exx_export_eigen_pair
  public :: exx_check_localization
  public :: exx_enabled,exx_refresh,exx_add_action,exx_exchange_energy,exx_freeze
  public :: exx_pack,exx_unpack,exx_timings,exx_walltime
  public :: exx_taylor_stage,exx_core_exchange,exx_force_full_action
  type(s_exx_symmetry_map),save :: symmetry_map
  type(hse_kernel),save :: kernel
  type(spatial_exx_state),target,save :: spatial
  type(s_exx_wannier),save :: wannier
  real(8),allocatable,save :: cached_occupation(:,:)
  complex(8),allocatable,save :: cached_action(:,:,:)
  logical,save :: finite_support_localized=.true.,radius_warning_reported=.false.
  logical,save,public :: exx_adaptive_ready=.false.,exx_support_changed=.false.
  logical,save :: adaptive_active=.false.,cached_adaptive_ready=.false.
  logical,save :: exx_force_full_action=.false.
  type(s_exx_ace),save :: ace
  type(s_exx_ace),save :: initial_ace,midpoint_ace
  complex(8),allocatable,save :: full_source(:,:,:),initial_source(:,:,:),midpoint_source(:,:,:)
  logical,save :: taylor_active=.false.,taylor_midpoint=.false.
  complex(8),allocatable,save :: cached_source(:,:,:),target_work(:,:,:),action_work(:,:,:)
  real(8),save :: exx_exchange_energy=0d0
  real(8),save :: exx_timings(4)=0d0 ! full EXX, ACE build, ACE apply, EXX collectives
  logical,save :: exx_freeze=.false.,reported_team=.false.,timing_enabled=.false.
contains
  subroutine exx_check_localization()
    implicit none
    if(dc_canonical())return
    ! A density criterion cannot certify a gauge-dependent truncated operator.
    if(exx_mlwf_norm_fraction>0d0.and.exx_mlwf_radius==0d0.and..not.adaptive_active) &
      error stop 'Adaptive EXX: localized support not established; SCF result rejected'
    if(exx_mlwf_radius>0d0.and..not.finite_support_localized) &
      error stop 'EXX MLWF finite support: localization not converged; SCF result rejected'
  end subroutine exx_check_localization

  logical function exx_eigen_diagnostic_enabled(info,solver) result(enabled)
    use communication, only: comm_bcast
    implicit none
    type(s_parallel_info),intent(in) :: info
    logical,optional,intent(in) :: solver
    integer :: flag
    enabled=.false.
    if(.not.exx_enabled())return
    flag=0
    if(info%id_rko==0)then
      if(yn_hse_eigen_diagnostic=='y')flag=1
      if(present(solver))then
        if(solver)then
          flag=0
          if(yn_hse_solver_diagnostic=='y')flag=1
        endif
      endif
    endif
    call comm_bcast(flag,info%icomm_rko,0)
    enabled=flag==1
  end function

  subroutine exx_export_eigen_pair(system,mg,info,psi,hpsi,tag)
    use iso_fortran_env, only: int32
    use salmon_global, only: base_directory
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: psi,hpsi
    character(*),optional,intent(in) :: tag
    character(:),allocatable :: filename
    complex(8),allocatable :: p(:,:,:),hp(:,:,:)
    integer :: iu,status,closed,total,ng,no
    ! Diagnostic format deliberately supports only a full Gamma fragment.
    if(system%nk/=1.or.info%isize_r/=1.or.info%isize_o/=1) &
      error stop 'HSE eigen diagnostic requires Gamma and full grid/orbitals'
    ng=product(mg%num);no=system%no
    allocate(p(ng,no,1),hp(ng,no,1))
    call exx_pack(psi,mg,info,p);call exx_pack(hpsi,mg,info,hp)
    filename=trim(base_directory)//'hse_eigen_pair.bin'
    if(present(tag))filename=trim(base_directory)//'hse_eigen_'//tag//'.bin'
    open(newunit=iu,file=filename, &
      status='replace',access='stream',form='unformatted',iostat=status)
    if(status==0)then
      write(iu,iostat=status)int([16909060,1,ng,no,1],int32)
      if(status==0)write(iu,iostat=status)system%hvol,system%rocc(:,:,1),p,hp
      close(iu,iostat=closed)
      if(status==0)status=closed
    endif
    call comm_summation(abs(status),total,info%icomm_rko)
    if(total/=0)error stop 'HSE eigen diagnostic write failed'
  end subroutine

  subroutine exx_export_snapshot(system,mg,info,psi,iteration,residual,converged)
    use salmon_global, only: base_directory
    use communication, only: comm_bcast
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: psi
    integer,intent(in) :: iteration
    real(8),intent(in) :: residual
    logical,intent(in) :: converged
    integer :: status,enabled
    if(.not.exx_enabled().or..not.use_wannier_exchange())return
    enabled=0
    if(info%id_k==0)then
      if(yn_hse_wannier_snapshot=='y')enabled=1
    endif
    call comm_bcast(enabled,info%icomm_k,0)
    if(enabled==0)return
    call exx_refresh(system,mg,info,psi)
    status=0
    if(info%id_k==0)call wannier_snapshot(wannier,wannier%source_occupation,hse_omega,exx_exchange_energy, &
      residual,iteration,converged,trim(base_directory)//'hse_wannier_snapshot.bin',status)
    call comm_bcast(status,info%icomm_k,0)
    if(status/=0)error stop 'HSE Wannier snapshot: write failed'
  end subroutine
  subroutine exx_taylor_stage(stage,system,mg,info,psi)
    implicit none
    integer,intent(in) :: stage
    type(s_dft_system),intent(in),optional :: system
    type(s_rgrid),intent(in),optional :: mg
    type(s_parallel_info),intent(in),optional :: info
    type(s_orbital),intent(in),optional :: psi
    integer :: ierr,no
    if(lcfo_rt_active)then
      call lcfo_exx_stage(stage,system,mg,info,psi)
      return
    endif
    select case(stage)
    case(0)
      taylor_active=.true.;taylor_midpoint=.false.
      if(propagator=='hse_taylor4_full')then
        initial_source=full_source
      else
        initial_ace=ace
      endif
    case(1)
      if(propagator=='hse_taylor4_full')then
        no=size(full_source,2)
        if(allocated(midpoint_source))deallocate(midpoint_source)
        allocate(midpoint_source(size(full_source,1),2*no,size(full_source,3)))
        midpoint_source(:,:no,:)=initial_source/sqrt(2d0)
        midpoint_source(:,no+1:,:)=full_source/sqrt(2d0)
      endif
      ! Apply the two endpoint operators with half weights; no doubled factors.
      taylor_midpoint=.true.
    case(2)
      taylor_active=.false.;taylor_midpoint=.false.
      call exx_ace_clear(initial_ace)
      call exx_ace_clear(midpoint_ace)
      if(allocated(initial_source))deallocate(initial_source)
      if(allocated(midpoint_source))deallocate(midpoint_source)
    case default
      error stop 'Invalid HSE Taylor stage'
    end select
  end subroutine

  real(8) function exx_walltime()
    implicit none
    integer(int64) :: count,rate
    call system_clock(count,rate)
    exx_walltime=real(count,8)/real(rate,8)
  end function

  subroutine warn_fixed_radius(max_loss)
    implicit none
    real(8),intent(in) :: max_loss
    real(8) :: target
    if(exx_mlwf_radius<=0d0.or.radius_warning_reported)return
    target=merge(exx_mlwf_norm_fraction,.999d0,exx_mlwf_norm_fraction>0d0)
    if(1d0-max_loss>=target-64d0*epsilon(1d0))return
    write(*,'(a,3es18.9)') 'WARNING EXX fixed radius retains less than target; R(bohr)/min retained/target: ', &
      exx_mlwf_radius,1d0-max_loss,target
    radius_warning_reported=.true.
  end subroutine warn_fixed_radius

  logical function dc_canonical()
    implicit none
    dc_canonical=yn_dc=='y'.and.yn_exx_dc_mlwf=='n'
  end function

  logical function exx_enabled()
    use exx_functional, only: is_hybrid
    implicit none
    exx_enabled=(is_hybrid(xc)).and..not.exx_pre_scf_active
  end function

  subroutine exx_pack(psi,mg,info,a)
    implicit none
    type(s_orbital),intent(in) :: psi
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    complex(8),intent(out) :: a(:,:,:)
    integer :: ik,io,is(3),ie(3)
    is=mg%is;ie=mg%ie
    do ik=info%ik_s,info%ik_e;do io=info%io_s,info%io_e
      a(:,io-info%io_s+1,ik-info%ik_s+1)=reshape( &
        psi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),1,io,ik,1),[product(mg%num)])
    enddo;enddo
  end subroutine

  subroutine exx_unpack(a,psi,mg,info)
    implicit none
    complex(8),intent(in) :: a(:,:,:)
    type(s_orbital),intent(inout) :: psi
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    integer :: ik,io,is(3),ie(3)
    is=mg%is;ie=mg%ie
    do ik=info%ik_s,info%ik_e;do io=info%io_s,info%io_e
      psi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),1,io,ik,1)= &
        reshape(a(:,io-info%io_s+1,ik-info%ik_s+1),mg%num)
    enddo;enddo
    psi%update_zwf_overlap=.false.
  end subroutine

  subroutine exx_refresh(system,mg,info,psi)
    use exx_functional, only: is_global_hybrid
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: psi
    complex(8),allocatable :: w(:,:,:),local(:,:,:)
    real(8) :: ex,offdiag(3,3),tick,communication_before
    integer :: ierr,total_error,ng,nk,no,n,mesh,j,first_full,count_full
    if(.not.exx_enabled().or.exx_freeze)return
    if(yn_periodic/='y'.or.system%nspin/=1.or..not.allocated(psi%zwf)) &
      error stop 'HSE06: periodic complex unpolarized orbitals required'
    if(yn_md=='y'.and.theory/='dft_md')then
      if((.not.is_global_hybrid(xc)).or. &
         (theory/='tddft_response'.and.theory/='tddft_pulse').or.yn_conventional_from_dcdft/='y'.or.lcfo_rt_active) &
        error stop 'Hybrid: unsupported real-time ionic extension'
    endif
    if(yn_spinorbit/='n'.or.yn_jm/='n'.or.yn_symmetrized_stencil=='y') &
      error stop 'HSE06: unsupported Hamiltonian/ionic extension'
    if(PLUS_U_ON)error stop 'HSE06: DFT+U combination unsupported'
    if(allocated(system%Ac_micro%v))error stop 'HSE06: microscopic vector potential unsupported'
    if(lcfo_rt_active)then
      call lcfo_exx_refresh(system,mg,info,psi,exx_exchange_energy)
      return
    endif
    if((info%isize_r>1.or.info%isize_o>1.or.(exx_mlwf_norm_fraction>0d0.and.exx_mlwf_radius==0d0)).and. &
       ((theory=='dft').or. &
        (yn_dc=='n'.and.(yn_conventional_from_dcdft=='y'.or.is_global_hybrid(xc)).and. &
         (theory=='tddft_response'.or.theory=='tddft_pulse'))))then
      call refresh_spatial(system,mg,info,psi)
      return
    endif
    if(info%isize_r/=1.or.info%isize_o/=1.or.info%numm/=1) &
      error stop 'HSE06: initial native support requires k-only MPI distribution'

    if(use_wannier_exchange())then
      call refresh_wannier(system,mg,info,psi)
      return
    endif
    ng=product(mg%num);nk=system%nk;no=system%no;n=mg%num(1);mesh=nint(real(nk,8)**(1d0/3d0))
    if(use_symmetry)mesh=num_kgrid(1)
    if(any(mg%num/=n).or.(.not.use_symmetry.and.mesh**3/=nk).or. &
      maxval(abs(system%hgs-system%hgs(1)))>1d-12) &
      error stop 'HSE06: cubic grid and full cubic k mesh required'
    offdiag=system%primitive_a
    do j=1,3;offdiag(j,j)=0;enddo
    if(maxval(abs(offdiag))>1d-12)error stop 'HSE06: orthogonal cell required'
    if((.not.use_symmetry.and.maxval(abs(system%wtk-1d0/nk))>1d-12).or. &
      maxval(abs(system%rocc-2d0))>1d-12) &
      error stop 'HSE06: uniform k weights and fully occupied spin pairs required'
    if(info%io_s/=1.or.info%io_e/=no.or.info%numk<1)error stop 'HSE06: unsupported orbital layout'
    if(kernel%n==0)then
      ! Keep representatives persistent; expand only rank-local stars for EXX.
      if(use_symmetry)then
        if(any(num_kgrid/=mesh))error stop 'HSE symmetry: cubic full mesh required'
        call symmetry_init(symmetry_map,mg%num,system%hgs,system%vec_k(:,:nk),system%wtk(:nk), &
          SymMatA,SymMatB,mesh,ierr)
        call comm_summation(ierr,total_error,info%icomm_rko)
        if(total_error/=0)error stop 'HSE symmetry: invalid grid, stars or weights'
        call symmetry_validate_atoms(symmetry_map,system%Rion,system%kion,ierr)
        call comm_summation(ierr,total_error,info%icomm_rko)
        if(total_error/=0)error stop 'HSE symmetry: operation does not preserve atoms'
        first_full=symmetry_map%first(info%ik_s)
        count_full=symmetry_map%first(info%ik_e+1)-first_full
        call hse_kernel_init(kernel,n,mesh,system%hgs(1),symmetry_map%full_k,hse_omega, &
          max(1,min(16,64/info%isize_k)),ierr,first_full,count_full, &
          block_rows=hse_block_rows,profile=yn_hse_profile=='y',fft_layout=hse_fft_layout)
      else
        call hse_kernel_init(kernel,n,mesh,system%hgs(1),system%vec_k,hse_omega, &
          max(1,min(16,64/info%isize_k)),ierr,info%ik_s,info%numk, &
          block_rows=hse_block_rows,profile=yn_hse_profile=='y',fft_layout=hse_fft_layout)
      endif
      call comm_summation(ierr,total_error,info%icomm_rko)
      if(total_error/=0)error stop 'HSE06: kernel initialization failed'
      timing_enabled=kernel%profile.or.propagator=='hse_ptcn'
    endif
    allocate(local(ng,no,info%numk))
    call exx_pack(psi,mg,info,local)
    ierr=1
    if(allocated(cached_source))then
      if(all(shape(cached_source)==shape(local)))then
        if(all(cached_source==local))ierr=0
      endif
    endif
    call comm_summation(ierr,total_error,info%icomm_k)
    if(total_error==0)return
    allocate(w(ng,no,info%numk))
    if(propagator=='hse_taylor4_full')full_source=local
    if(timing_enabled)then
      tick=exx_walltime();communication_before=exx_timings(4)
    endif
    call apply_distributed(local,local,w,info,ierr)
    if(timing_enabled)exx_timings(1)=exx_timings(1)+exx_walltime()-tick-(exx_timings(4)-communication_before)
    call comm_summation(ierr,total_error,info%icomm_k)
    if(total_error/=0)error stop 'HSE06: distributed exchange action failed'
    if(kernel%profile.and.info%id_k==0)write(*,'(a,i0,a,6es14.5)') &
      'HSE_PROFILE rank0 block=',kernel%block,' density comm fft kernel packing action=',kernel%seconds
    if(.not.reported_team.and.info%id_k==0)write(*,'(a,i0)')'HSE_OPENMP threads=',kernel%threads_used
    if(.not.reported_team.and.info%id_k==0)write(*,'(a,l1)')'HSE_FFT_CONTIGUOUS=',kernel%contiguous_fft
    if(.not.reported_team.and.info%id_k==0.and.kernel%auto_fft) &
      write(*,'(a,2es14.5)')'HSE_FFT_AUTO trial strided contiguous seconds=',kernel%fft_trial_seconds
    reported_team=.true.
    if(timing_enabled)tick=exx_walltime()
    call exx_ace_build(ace,local,w,system%hvol,ierr)
    if(timing_enabled)exx_timings(2)=exx_timings(2)+exx_walltime()-tick
    call comm_summation(ierr,total_error,info%icomm_k)
    if(total_error/=0)error stop 'HSE06: ACE metric failed'
    cached_source=local
    ex=0d0
    do j=1,info%numk
      ex=ex+exchange_fraction()*real(sum(conjg(local(:,:,j))*w(:,:,j)),8)*system%hvol*system%wtk(info%ik_s+j-1)
    enddo
    call comm_summation(ex,exx_exchange_energy,info%icomm_k)
  end subroutine

  subroutine apply_distributed(source,target,action,info,ierr)
    implicit none
    complex(8),intent(in) :: source(:,:,:),target(:,:,:)
    complex(8),intent(out) :: action(:,:,:)
    type(s_parallel_info),intent(in) :: info
    integer,intent(out) :: ierr
    integer :: local_layout(2*info%isize_k),layout(2*info%isize_k),np,j,index,total_error
    complex(8),allocatable :: expanded_target(:,:,:),expanded_action(:,:,:),rotated(:,:)
    np=info%isize_k;local_layout=0
    local_layout(info%id_k+1)=info%ik_s
    local_layout(np+info%id_k+1)=info%numk
    if(use_symmetry)then
      local_layout(info%id_k+1)=symmetry_map%first(info%ik_s)
      local_layout(np+info%id_k+1)=symmetry_map%first(info%ik_e+1)-symmetry_map%first(info%ik_s)
    endif
    call comm_summation(local_layout,layout,size(layout),info%icomm_k)
    if(use_symmetry)then
      ierr=0
      if(.not.all(ieee_is_finite(real(source))).or..not.all(ieee_is_finite(aimag(source))))ierr=1
      call comm_summation(ierr,total_error,info%icomm_k)
      if(total_error/=0)then
        ierr=1;return
      endif
      allocate(expanded_target(size(target,1),size(target,2),layout(np+info%id_k+1)))
      allocate(rotated(size(source,1),size(source,2)))
      do j=1,size(expanded_target,3)
        index=layout(info%id_k+1)+j-1
        call symmetry_transform(symmetry_map,index,1,target(:,:,symmetry_map%owner(index)-info%ik_s+1), &
          expanded_target(:,:,j),ierr)
        if(ierr/=0)error stop 'HSE symmetry target transformation failed'
      enddo
      allocate(expanded_action(size(expanded_target,1),size(expanded_target,2),size(expanded_target,3)))
      ! Source argument is only a valid layout placeholder when the density callback is supplied.
      call hse_kernel_apply_distributed(kernel,expanded_target,expanded_target,expanded_action,layout(:np), &
        layout(np+1:),info%id_k,transpose_tiles,ierr,fill_density)
      if(ierr==0)then
        do j=1,info%numk
          index=symmetry_map%first(info%ik_s+j-1)-symmetry_map%first(info%ik_s)+1
          action(:,:,j)=expanded_action(:,:,index)
        enddo
      endif
      return
    endif
    call hse_kernel_apply_distributed(kernel,source,target,action,layout(:np),layout(np+1:), &
      info%id_k,transpose_tiles,ierr)
  contains
    subroutine fill_density(j,lo,rows,density)
      implicit none
      integer,intent(in) :: j,lo,rows
      complex(8),intent(out) :: density(:,:)
      integer :: full_index,rep,op,g,stat,ng,no
      complex(8) :: coefficient
      external :: zgemm
      ng=size(source,1);no=size(source,2)
      full_index=layout(info%id_k+1)+j-1
      rep=symmetry_map%owner(full_index)-info%ik_s+1
      coefficient=cmplx(1d0/symmetry_map%multiplicity(full_index),0d0,8)
      density=0d0
      ! Average little-group projectors one orbital block at a time, not copies of all orbitals.
      do op=1,symmetry_map%multiplicity(full_index)
        call symmetry_transform(symmetry_map,full_index,op,source(:,:,rep),rotated,stat)
        if(stat/=0)error stop 'HSE symmetry source transformation failed'
        do g=1,no
          rotated(:,g)=rotated(:,g)*kernel%phase(:,j)
        enddo
        call zgemm('N','C',rows,ng,no,coefficient,rotated(lo,1),ng,rotated(1,1),ng, &
          (1d0,0d0),density(1,1),size(density,1))
      enddo
    end subroutine
    subroutine transpose_tiles(send,recv,count)
      implicit none
      complex(8),intent(in) :: send(:)
      complex(8),intent(out) :: recv(:)
      integer,intent(in) :: count
      real(8) :: start
      if(timing_enabled)start=exx_walltime()
      call comm_alltoall(send,recv,info%icomm_k,count)
      if(timing_enabled)exx_timings(4)=exx_timings(4)+exx_walltime()-start
    end subroutine
  end subroutine

  subroutine exx_add_action(psi,hpsi,system,mg,info,lcfo_coeff,lcfo_action)
    implicit none
    type(s_orbital),intent(in) :: psi
    type(s_orbital),intent(inout) :: hpsi
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    complex(8),intent(in),optional :: lcfo_coeff(:,:)
    complex(8),intent(out),optional :: lcfo_action(:,:)
    integer :: ierr,ng,total_error,info_error
    integer,allocatable :: orbital_comm
    real(8) :: tick,communication_before,action_scale
    if(.not.exx_enabled())return
    if(lcfo_rt_active)then
      call lcfo_exx_add_action(psi,hpsi,system,mg,info,lcfo_coeff,lcfo_action)
      return
    endif
    if(.not.exx_ace_ready(ace))error stop 'HSE06: occupied exchange source is not initialized'
    if(info%isize_o>1)orbital_comm=info%icomm_o
    ng=product(mg%num)
    if(allocated(target_work))then
      if(any(shape(target_work)/=[ng,info%numo,info%numk]))deallocate(target_work,action_work)
    endif
    if(.not.allocated(target_work))allocate(target_work(ng,info%numo,info%numk), &
      action_work(ng,info%numo,info%numk))
    call exx_pack(psi,mg,info,target_work)
    action_scale=exchange_fraction()
    if(timing_enabled)tick=exx_walltime()
    if(use_wannier_exchange().and.exx_force_full_action)then
      if(info%isize_r>1.or.info%isize_o>1.or.(exx_mlwf_norm_fraction>0d0.and.exx_mlwf_radius==0d0))then
        call spatial_exx_apply(spatial,num_rgrid,system%hgs,[info%isize_y,info%isize_z], &
          [info%id_y,info%id_z],[info%icomm_y,info%icomm_z],info%icomm_r, &
          pbeh_coulomb_radius,target_work,action_work,info_error,omega=merge(hse_omega,0d0,xc=='hse06'),comm_o=orbital_comm)
        if(info_error/=0)error stop 'Spatial EXX: full action failed'
      else
        call apply_wannier_collective(target_work,action_work,info)
      endif
      ierr=0
    else if(taylor_active.and.propagator=='hse_taylor4_full')then
      communication_before=exx_timings(4)
      if(taylor_midpoint)then
        call apply_distributed(midpoint_source,target_work,action_work,info,ierr)
      else
        call apply_distributed(initial_source,target_work,action_work,info,ierr)
      endif
      call comm_summation(ierr,total_error,info%icomm_k)
      if(total_error/=0)error stop 'HSE Taylor full target action failed'
      if(timing_enabled)exx_timings(1)=exx_timings(1)+exx_walltime()-tick-(exx_timings(4)-communication_before)
    else
      if(taylor_active.and.taylor_midpoint)then
        call apply_endpoint(initial_ace)
        if(ierr/=0)error stop 'HSE06: initial endpoint ACE application failed'
        call add_mesh_action(.5d0*exchange_fraction())
        action_scale=.5d0*exchange_fraction()
      endif
      call apply_endpoint(ace)
      if(timing_enabled)exx_timings(3)=exx_timings(3)+exx_walltime()-tick
    endif
    if(ierr/=0)error stop 'HSE06: ACE application failed'
    call add_mesh_action(action_scale)
  contains
    subroutine apply_endpoint(state)
      implicit none
      type(s_exx_ace),intent(in) :: state
      if(info%isize_o>1.or.state%packed)then
        call orbital_ace_apply(state,target_work,action_work,info%icomm_r,info%icomm_o,ierr)
      else if(info%isize_r>1)then
        call exx_ace_apply(state,target_work,action_work,ierr,sum_spatial)
      else
        call exx_ace_apply(state,target_work,action_work,ierr)
      endif
    end subroutine
    subroutine add_mesh_action(weight)
      implicit none
      real(8),intent(in) :: weight
      integer :: ix,iy,iz,io,ik,g
      do ik=info%ik_s,info%ik_e;do io=info%io_s,info%io_e
        g=0
        do iz=mg%is(3),mg%ie(3);do iy=mg%is(2),mg%ie(2);do ix=mg%is(1),mg%ie(1)
          g=g+1
          hpsi%zwf(ix,iy,iz,1,io,ik,1)=hpsi%zwf(ix,iy,iz,1,io,ik,1) &
            +weight*action_work(g,io-info%io_s+1,ik-info%ik_s+1)
        enddo;enddo;enddo
      enddo;enddo
      hpsi%update_zwf_overlap=.false.
    end subroutine
    subroutine sum_spatial(a)
      implicit none
      complex(8),intent(inout) :: a(:,:)
      complex(8) :: total(size(a,1),size(a,2))
      call comm_summation(a,total,size(a),info%icomm_r)
      a=total
    end subroutine
  end subroutine
  logical function use_wannier_exchange()
    implicit none
    use_wannier_exchange=yn_dc=='y'.or.yn_hse_wannier=='y'
  end function

  subroutine refresh_spatial(system,mg,info,psi)
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: psi
    complex(8),allocatable :: local(:,:,:),w(:,:,:)
    real(8) :: ex,offdiag(3,3)
    integer :: status,total,changed,maxiter,j,io
    integer,allocatable :: orbital_comm
    real(8),allocatable :: radii(:),loss(:)
    logical,allocatable :: protected(:)
    integer(int64) :: fft_work(3)
    logical :: was_active,screen_fallback,radius_covers_cell,support_accepted
    integer :: requested_screen_mode
    real(8) :: correction_norm,accepted_bound
    integer :: adaptive_bad
    real(8) :: mask_diagnostic(2),mask_maximum(2)
    if((exx_mlwf_norm_fraction>0d0.or.exx_mlwf_radius>0d0).and.theory/='dft')exx_adaptive_ready=.true.
    if(info%isize_x/=1.or.info%isize_k/=1.or.info%numm/=1.or.system%nk/=1) &
      error stop 'Spatial EXX: Gamma y/z pencils required'
    if(info%isize_o>system%no)error stop 'Spatial EXX: each orbital group must own at least one state'
    if(any(num_kgrid/=1).or.maxval(abs(system%vec_k))>1d-12.or.use_symmetry) &
      error stop 'Spatial EXX: unshifted Gamma required'
    if(theory/='dft'.and.maxval(abs(system%rocc-2d0))>1d-12) &
      error stop 'Spatial EXX: occupied spin pairs required'
    if(theory/='dft'.and.propagator/='hse_taylor4')error stop 'Spatial EXX: Taylor4 ACE required'
    if(any(mg%num/=num_rgrid/[1,info%isize_y,info%isize_z]).or. &
       any(mg%is/=[1,info%id_y*mg%num(2)+1,info%id_z*mg%num(3)+1])) &
      error stop 'Spatial EXX: mesh pencil layout mismatch'
    offdiag=system%primitive_a
    do j=1,3
      offdiag(j,j)=0d0
    enddo
    if(maxval(abs(offdiag))>1d-12)error stop 'Spatial EXX: orthogonal cell required'
    if(info%isize_o>1)orbital_comm=info%icomm_o
    allocate(local(product(mg%num),info%numo,1))
    call exx_pack(psi,mg,info,local)
    changed=1
    if(allocated(cached_source).and.allocated(cached_occupation))then
      if(all(shape(cached_source)==shape(local)).and.all(shape(cached_occupation)==shape(system%rocc(:,:,1))))then
        if(all(cached_source==local).and.all(cached_occupation==system%rocc(:,:,1)))changed=0
      endif
    endif
    if((exx_mlwf_norm_fraction>0d0.or.exx_mlwf_radius>0d0).and.(exx_adaptive_ready.neqv.cached_adaptive_ready))changed=1
    call comm_summation(changed,total,info%icomm_ro)
    if(total==0)return
    allocate(w(product(mg%num),info%numo,1))
    maxiter=0
    if(mod(spatial%updates,exx_mlwf_interval)==0)maxiter=exx_mlwf_maxiter
    if((exx_mlwf_norm_fraction>0d0.or.exx_mlwf_radius>0d0).and.exx_adaptive_ready.and.spatial%last_localization_status/=0) &
      maxiter=exx_mlwf_maxiter
    spatial%seed_localized=exx_mlwf_norm_fraction>0d0.or.exx_mlwf_radius>0d0
    spatial%retain_accepted_gauge=theory/='dft'
    if(dc_canonical())then
      maxiter=0
      call spatial_exx_canonical_source(spatial,local,system%rocc(info%io_s:info%io_e,:,1), &
        info%icomm_r,status,comm_o=orbital_comm)
      if(info%id_ro==0.and.spatial%updates==1)write(*,'(a)') 'EXX_DC canonical full-fragment source (spatial)'
    else
    call spatial_exx_refresh(spatial,num_rgrid,system%hgs,[info%isize_y,info%isize_z], &
      [info%id_y,info%id_z],[info%icomm_y,info%icomm_z],info%icomm_r,local, &
      maxiter,exx_mlwf_tolerance,status,occupation=system%rocc(info%io_s:info%io_e,:,1),comm_o=orbital_comm)
    endif
    if(status/=0)error stop 'Spatial EXX: source refresh failed'
    if((exx_mlwf_norm_fraction>0d0.or.exx_mlwf_radius>0d0).and.exx_mlwf_norm_fraction<1d0.and.theory/='dft')then
      if(spatial%last_localization_status/=0) &
        error stop 'Adaptive RT requires an accepted transported MLWF gauge; refine initial localization'
    endif
    if(spatial%retained_gauge.and.info%id_ro==0)write(*,'(a,i10)') &
      'EXX_SPATIAL retained accepted transported gauge at refresh: ',spatial%updates
    was_active=adaptive_active
    ! A sphere this large cannot truncate any point under the periodic metric;
    ! as for fraction=1, its exact action does not require a converged gauge.
    radius_covers_cell=exx_mlwf_radius>0d0.and. &
      exx_mlwf_radius>=sqrt(sum((.5d0*num_rgrid*system%hgs)**2))
    adaptive_active=.not.dc_canonical().and.(exx_mlwf_norm_fraction>0d0.or.exx_mlwf_radius>0d0).and.exx_adaptive_ready.and. &
      (spatial%last_localization_status==0.or.radius_covers_cell.or. &
       (exx_mlwf_norm_fraction==1d0.and.exx_mlwf_radius==0d0))
    if(adaptive_active)then
      allocate(radii(info%numo),loss(info%numo),protected(info%numo))
      call adaptive_source_mask(num_rgrid,system%hgs,mg%is-1,mg%num,info%icomm_r,spatial%source, &
        merge(exx_mlwf_norm_fraction,.999d0,exx_mlwf_norm_fraction>0d0),radii,loss,protected,status, &
        fixed_radius=exx_mlwf_radius)
      call comm_summation(status,adaptive_bad,info%icomm_ro)
      if(adaptive_bad/=0)error stop 'Adaptive EXX: source mask failed'
      mask_diagnostic=[maxval(loss),maxval(radii)]
      call comm_get_max(mask_diagnostic,mask_maximum,2,info%icomm_ro)
      if(info%id_ro==0)then
        if(exx_mlwf_radius>0d0)then
          write(*,'(a,3es18.9)') 'EXX_FIXED radius/min retained norm/max loss: ', &
            exx_mlwf_radius,1d0-mask_maximum(1),mask_maximum(1)
          call warn_fixed_radius(mask_maximum(1))
        else
          write(*,'(a,3es18.9)')'EXX_ADAPTIVE fraction/max radius/max norm loss: ', &
            exx_mlwf_norm_fraction,mask_maximum(2),mask_maximum(1)
        endif
      endif
    endif
    if(exx_mlwf_radius>0d0)finite_support_localized=adaptive_active
    if(adaptive_active.neqv.was_active)exx_support_changed=.true.
    cached_adaptive_ready=exx_adaptive_ready
    spatial%compact=adaptive_active.and.exx_local_fft=='auto'
    requested_screen_mode=0
    if(exx_pair_screening=='diagnose')requested_screen_mode=1
    if(exx_pair_screening=='on')requested_screen_mode=2
    if(dc_canonical())requested_screen_mode=0
    spatial%screen_mode=requested_screen_mode;spatial%screen_tolerance=exx_pair_tolerance/2d0
    fft_work=0_int64
    support_accepted=.false.
    if(exx_ace_support=='source')then
      if(.not.adaptive_active)error stop 'Source-support ACE requires active MLWF support'
      call build_source_support_ace(support_accepted)
      ! Fixed-ion source ACE never uses the full-action DC route. Its next
      ! refresh regenerates source from the mesh; only transport previous persists.
      if(support_accepted)deallocate(spatial%source)
      if(info%id_ro==0)write(*,'(a,l1,a)')'EXX_SUPPORT_ACE accepted: ',support_accepted, &
        ' (failure falls back to occupied-vector ACE)'
    endif
    if(.not.support_accepted)then
      spatial%screen_mode=requested_screen_mode;spatial%screen_tolerance=exx_pair_tolerance/2d0
      call apply_exchange_action()
      if(status/=0)error stop 'Spatial EXX: exchange action failed'
      if(requested_screen_mode/=0.and.info%id_ro==0)write(*,'(a,i2,2i18,2es18.9)') &
        'EXX_PAIR mode/candidates/skipped/action bound/max rank CPU seconds: ',requested_screen_mode, &
        spatial%screen_candidates,spatial%screen_skipped,spatial%screen_bound,spatial%screen_cpu_seconds
      if(requested_screen_mode/=0.and.info%id_ro==0)write(*,'(a,2i18)') &
        'EXX_PAIR generated grid products/catalogue entries: ',spatial%pair_products,spatial%pair_catalog_entries
      if(requested_screen_mode/=0.and.info%id_ro==0)write(*,'(a,i18)') &
        'EXX_PAIR evaluated product points: ',spatial%pair_product_points
      correction_norm=0d0;accepted_bound=0d0
      if(requested_screen_mode==2.and.spatial%screen_skipped>0)then
        call orbital_hermitian_action(local,w,system%hvol,info%icomm_r,info%icomm_o, &
          exx_pair_tolerance-spatial%screen_bound,correction_norm,status)
        accepted_bound=spatial%screen_bound+correction_norm
      endif
      if(status==0)call build_exchange_ace()
      screen_fallback=status/=0.and.requested_screen_mode==2.and.spatial%screen_skipped>0
      if(screen_fallback)then
        ! Pair-dependent omissions need not define a Hermitian input metric.
        ! Never relax the existing Hermitian/positive ACE validation to accept them.
        spatial%screen_mode=0;accepted_bound=0d0
        call apply_exchange_action()
        if(status/=0)error stop 'Spatial EXX: unscreened fallback action failed'
        call build_exchange_ace()
      endif
      if(status/=0)error stop 'Spatial EXX: ACE build failed'
      if(requested_screen_mode/=0.and.info%id_ro==0)write(*,'(a,l1,a,es18.9)') &
        'EXX_PAIR unscreened ACE fallback: ',screen_fallback,' accepted action bound: ', &
        accepted_bound
    endif ! occupied-vector ACE, including support-ACE fallback
    if(adaptive_active.and.info%id_ro==0)write(*,'(a,3i18)') &
      'EXX_ADAPTIVE local/global pairs/local FFT points (orbital group 0): ', &
      fft_work
    cached_occupation=system%rocc(:,:,1)
    if(spatial%updates==1)write(*,'(a,4i10)')'EXX_ORBITALS rank/local/global/grid: ', &
      info%id_ro,info%numo,system%no,product(mg%num)
    ex=0d0
    do j=1,info%numo
      io=info%io_s+j-1
      ex=ex+.5d0*exchange_fraction()*system%hvol*system%rocc(io,1,1)*system%wtk(1) &
        *real(sum(conjg(local(:,j,1))*w(:,j,1)),8)
    enddo
    call comm_summation(ex,exx_exchange_energy,info%icomm_ro)
    call move_alloc(local,cached_source)
    if(yn_dc=='y')then
      call move_alloc(w,cached_action)
    else if(allocated(cached_action))then
      deallocate(cached_action)
    endif
    if(.not.dc_canonical().and.info%id_ro==0.and.(spatial%updates==1.or.maxiter>0)) &
      write(*,'(a,3i8,3es16.7)')'EXX_SPATIAL refresh/iterations/status/spread/gradient/overlap: ', &
      spatial%updates,spatial%iterations,spatial%localization_status,spatial%spread,spatial%gradient,spatial%min_singular
  contains
    subroutine build_source_support_ace(accepted)
      implicit none
      logical,intent(out) :: accepted
      complex(8),pointer :: training(:,:,:)
      accepted=.false.
      ! Alias the existing finite source; do not allocate another full-grid copy.
      training(1:size(spatial%source,1),1:size(spatial%source,2),1:1)=>spatial%source
      spatial%screen_mode=2;spatial%screen_tolerance=0d0
      call spatial_exx_apply(spatial,num_rgrid,system%hgs,[info%isize_y,info%isize_z], &
        [info%id_y,info%id_z],[info%icomm_y,info%icomm_z],info%icomm_r, &
        pbeh_coulomb_radius,training,w,status,omega=merge(hse_omega,0d0,xc=='hse06'),comm_o=orbital_comm)
      if(status/=0)error stop 'Source-support ACE: exchange action failed'
      call record_fft_work()
      if(info%id_ro==0)write(*,'(a,4i18)')'EXX_SUPPORT_ACE products/skipped/catalogue/product points: ', &
        spatial%pair_products,spatial%screen_skipped,spatial%pair_catalog_entries,spatial%pair_product_points
      ! The ACE metric is -S^H K_S S. S need not be orthonormal; the existing
      ! Hermitian/positive metric checks and conditioning threshold still apply.
      call orbital_ace_build(ace,training,w,system%hvol,info%icomm_r,info%icomm_o,status,packed=.true.)
      call comm_summation(status,adaptive_bad,info%icomm_ro)
      if(adaptive_bad/=0)return
      call orbital_ace_apply(ace,local,w,info%icomm_r,info%icomm_o,status)
      if(status/=0)error stop 'Source-support ACE: occupied mesh action failed'
      accepted=.true.
    end subroutine
    subroutine apply_exchange_action()
      implicit none
      complex(8),allocatable :: localized_action(:,:,:),adjoint(:,:)
      integer,allocatable :: counts(:)
      integer :: first
      if(requested_screen_mode==0)then
        call spatial_exx_apply(spatial,num_rgrid,system%hgs,[info%isize_y,info%isize_z], &
          [info%id_y,info%id_z],[info%icomm_y,info%icomm_z],info%icomm_r, &
          pbeh_coulomb_radius,local,w,status,omega=merge(hse_omega,0d0,xc=='hse06'),comm_o=orbital_comm)
      else
        allocate(localized_action(size(w,1),size(w,2),1))
        call spatial_exx_apply(spatial,num_rgrid,system%hgs,[info%isize_y,info%isize_z], &
          [info%id_y,info%id_z],[info%icomm_y,info%icomm_z],info%icomm_r, &
          pbeh_coulomb_radius,spatial%previous,localized_action,status, &
          omega=merge(hse_omega,0d0,xc=='hse06'),comm_o=orbital_comm)
        if(status/=0)return
        adjoint=conjg(transpose(spatial%gauge(:,:,1)))
        if(info%isize_o>1)then
          call orbital_layout(info%numo,info%icomm_r,info%icomm_o,counts,first,status)
          if(status/=0)return
          call orbital_rotate(localized_action(:,:,1),adjoint,info%icomm_o,counts,first,w(:,:,1))
        else
          w(:,:,1)=matmul(localized_action(:,:,1),adjoint)
        endif
      endif
      call record_fft_work()
    end subroutine
    subroutine record_fft_work()
      ! Count rejected ACE attempts too: each exchange call resets its counters.
      implicit none
      fft_work=fft_work+[spatial%local_pairs,spatial%global_pairs,spatial%local_points]
    end subroutine
    subroutine build_exchange_ace()
      implicit none
      if(info%isize_o>1)then
        call orbital_ace_build(ace,local,w,system%hvol,info%icomm_r,info%icomm_o,status)
      else
        call exx_ace_build(ace,local,w,system%hvol,status,sum_spatial)
      endif
      call comm_summation(status,adaptive_bad,info%icomm_ro)
      status=adaptive_bad
    end subroutine
    subroutine sum_spatial(a)
      implicit none
      complex(8),intent(inout) :: a(:,:)
      complex(8) :: total(size(a,1),size(a,2))
      call comm_summation(a,total,size(a),info%icomm_r)
      a=total
    end subroutine
  end subroutine

  subroutine refresh_wannier(system,mg,info,psi)
    use communication, only: comm_bcast
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: psi
    complex(8),allocatable :: local(:,:,:),send(:,:,:),allpsi(:,:,:),allw(:,:,:),w(:,:,:)
    real(8) :: ex,offdiag(3,3)
    integer :: ng,no,nk,ik,j,status,changed,total_changed,maxiter
    ng=product(mg%num);no=system%no;nk=system%nk
    if(use_symmetry.or.nk/=product(num_kgrid))error stop 'HSE Wannier: full k mesh required'
    if(maxval(abs(system%wtk-1d0/nk))>1d-12)error stop 'HSE Wannier: uniform k weights required'
    offdiag=system%primitive_a
    do j=1,3
      offdiag(j,j)=0d0
    enddo
    if(maxval(abs(offdiag))>1d-12)error stop 'HSE Wannier: orthogonal cell required'
    if(info%io_s/=1.or.info%io_e/=no.or.info%numk<1)error stop 'HSE Wannier: invalid orbital layout'
    allocate(local(ng,no,info%numk))
    call exx_pack(psi,mg,info,local)
    changed=1
    if(allocated(cached_source).and.allocated(cached_occupation))then
      if(all(shape(cached_source)==shape(local)).and.all(shape(cached_occupation)==shape(system%rocc(:,:,1))))then
        if(all(cached_source==local).and.all(cached_occupation==system%rocc(:,:,1)))changed=0
      endif
    endif
    call comm_summation(changed,total_changed,info%icomm_k)
    if(total_changed==0)return
    allocate(send(ng,no,nk),allpsi(ng,no,nk),allw(ng,no,nk),w(ng,no,info%numk))
    send=0d0;send(:,:,info%ik_s:info%ik_e)=local
    call comm_summation(send,allpsi,size(send),info%icomm_k)
    status=0
    if(info%id_k==0)then
      if(wannier%ng==0)then
        if(xc=='hse06')then
          call wannier_init(wannier,mg%num,num_kgrid,system%hgs,system%vec_k,hse_omega,status)
        else
          call wannier_init(wannier,mg%num,num_kgrid,system%hgs,system%vec_k,0d0,status,pbeh_coulomb_radius)
          write(*,*)'PBEh40 Coulomb radius (bohr): ', &
            merge(pbeh_coulomb_radius,.5d0*minval(mg%num*num_kgrid*system%hgs),pbeh_coulomb_radius>0d0)
        endif
      endif
      if(status==0)then
        maxiter=0
        if(mod(wannier%updates,exx_mlwf_interval)==0)maxiter=exx_mlwf_maxiter
        if(dc_canonical())maxiter=0
        call wannier_refresh_source(wannier,allpsi,system%rocc(:,:,1),maxiter,exx_mlwf_tolerance,status, &
          localize=.not.dc_canonical())
        if(dc_canonical().and.wannier%updates==1)write(*,'(a)') 'EXX_DC canonical full-fragment source (k mesh)'
      endif
      if(status==0)then
        wannier%use_local_fft=exx_local_fft=='auto'
        call wannier_truncate_source(wannier,exx_mlwf_radius,status)
        call warn_fixed_radius(wannier%max_discarded_norm_fraction)
        finite_support_localized=wannier%discarded_norm_fraction==0d0.or.wannier%last_localization_status==0
      endif
      if(status==0)call wannier_apply(wannier,allpsi,allw,status)
      if(.not.dc_canonical().and.status==0.and.(wannier%updates==1.or.maxiter>0))then
        write(*,'(a,3i7,3es16.7)')'EXX_WANNIER refresh/iterations/status/spread/gradient/overlap: ', &
        wannier%updates,wannier%localization_iterations,wannier%localization_status, &
        wannier%spread,wannier%gradient,wannier%min_singular
        if(wannier%use_local_fft)then
          write(*,'(a,4i18)')'EXX_FFT local/global pairs/actual/full pair grid points: ', &
            wannier%local_fft_pairs_executed,wannier%fft_pairs_executed-wannier%local_fft_pairs_executed, &
            wannier%fft_pair_grid_points,wannier%fft_pairs_executed*int(wannier%ngs,int64)
        endif
        if(exx_mlwf_radius>0d0)then
          write(*,'(a,es16.7,i8,2es16.7)')'EXX_MLWF radius (bohr)/protected/total loss/max factor loss: ', &
            exx_mlwf_radius,wannier%protected_sources,wannier%discarded_norm_fraction,wannier%max_discarded_norm_fraction
          if(wannier%localization_status/=0) &
            write(*,'(a)')'EXX_MLWF: localization not converged; finite-radius source approximation is gauge dependent.'
        else if(wannier%localization_status/=0)then
          write(*,'(a)')'EXX_WANNIER: localization not converged; retaining full-support exact exchange.'
        endif
      endif
    endif
    call comm_bcast(status,info%icomm_k,0)
    if(status/=0)error stop 'HSE Wannier: collective exchange refresh failed'
    call comm_bcast(finite_support_localized,info%icomm_k,0)
    call comm_bcast(allw,info%icomm_k,0)
    w=allw(:,:,info%ik_s:info%ik_e)
    call exx_ace_build(ace,local,w,system%hvol,status)
    call comm_summation(status,total_changed,info%icomm_k)
    if(total_changed/=0)error stop 'HSE Wannier: ACE construction metric failed'
    cached_source=local;cached_occupation=system%rocc(:,:,1)
    if(yn_dc=='y')then
      cached_action=w
    else if(allocated(cached_action))then
      deallocate(cached_action)
    endif
    ex=0d0
    do ik=1,info%numk
      do j=1,no
        ex=ex+.5d0*exchange_fraction()*system%rocc(j,info%ik_s+ik-1,1)*system%wtk(info%ik_s+ik-1)*system%hvol &
          *real(sum(conjg(local(:,j,ik))*w(:,j,ik)),8)
      enddo
    enddo
    call comm_summation(ex,exx_exchange_energy,info%icomm_k)
  end subroutine

  subroutine apply_wannier_collective(target,action,info)
    use communication, only: comm_bcast
    implicit none
    complex(8),intent(in) :: target(:,:,:)
    complex(8),intent(out) :: action(:,:,:)
    type(s_parallel_info),intent(in) :: info
    complex(8),allocatable :: send(:,:,:),alltarget(:,:,:),allaction(:,:,:)
    integer :: ng,nt,nk,status
    ng=size(target,1);nt=size(target,2);nk=product(num_kgrid)
    allocate(send(ng,nt,nk),alltarget(ng,nt,nk),allaction(ng,nt,nk))
    send=0d0;send(:,:,info%ik_s:info%ik_e)=target
    call comm_summation(send,alltarget,size(send),info%icomm_k)
    status=0
    if(info%id_k==0)call wannier_apply(wannier,alltarget,allaction,status)
    call comm_bcast(status,info%icomm_k,0)
    if(status/=0)error stop 'HSE Wannier: full trial action failed'
    call comm_bcast(allaction,info%icomm_k,0)
    action=allaction(:,:,info%ik_s:info%ik_e)
  end subroutine

  subroutine exx_core_exchange(system,mg,info,psi,core,energy)
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: psi
    integer,intent(in) :: core(3)
    real(8),intent(out) :: energy
    real(8) :: local
    integer :: ix,iy,iz,io,ik,g
    if(.not.allocated(cached_action))error stop 'DC HSE: refreshed exchange action missing'
    local=0d0
    do ik=info%ik_s,info%ik_e
      do io=info%io_s,info%io_e
        do iz=mg%is(3),min(mg%ie(3),core(3))
          do iy=mg%is(2),min(mg%ie(2),core(2))
            do ix=mg%is(1),min(mg%ie(1),core(1))
              g=1+(ix-mg%is(1))+mg%num(1)*((iy-mg%is(2))+mg%num(2)*(iz-mg%is(3)))
              local=local+.5d0*exchange_fraction()*system%rocc(io,ik,1)*system%wtk(ik)*system%hvol &
                *real(conjg(psi%zwf(ix,iy,iz,1,io,ik,1))*cached_action(g,io-info%io_s+1,ik-info%ik_s+1),8)
            enddo
          enddo
        enddo
      enddo
    enddo
    call comm_summation(local,energy,info%icomm_rko)
  end subroutine
end module
