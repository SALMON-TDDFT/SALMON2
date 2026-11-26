#include "config.h"
! Optional SALMONpy process-isolated endpoint-A protocol; ordinary runs are unchanged.
module tddft_worker
  use filesystem,only:resolve_path
  implicit none
  private
  public :: worker_requested, worker_configure, worker_initial_A, worker_next_A, worker_reply
  public :: worker_performance_requested,worker_spin_supported
  public :: worker_solid_requested,worker_solid_load,worker_near_requested,worker_near_load
  public :: worker_E_requested, worker_reply_E,worker_norm_requested,worker_electron_number,worker_state_reply
  public :: worker_read_line,worker_channel_context,worker_channel_begin,worker_channel_end,worker_host_requested
  type worker_channel_context
    integer::transport=0,controller=0,decomposition(5)=1
  end type
  ! Borrowed only during one blocking host call; the call frame owns the context.
  type(worker_channel_context),pointer::host_channel=>null()
  ! Embedded mode is set only by the explicit host-MPI ABI, once per process.
  logical,public,save :: worker_embedded=.false.
  real(8),public,save :: embedded_dt
  integer,public,save :: embedded_steps,embedded_decomposition(5)
  real(8), save :: worker_initial_A(3)=0d0
  real(8),save :: worker_electron_number=0d0
contains
  logical function worker_host_requested()
    worker_host_requested=associated(host_channel)
  end function
  subroutine worker_channel_begin(owner)
    type(worker_channel_context),target,intent(inout)::owner
    if(associated(host_channel))error stop 'molecular host channel already active'
    host_channel=>owner
  end subroutine
  subroutine worker_channel_end()
    nullify(host_channel)
  end subroutine
  subroutine channel_read(line,ios)
#ifdef USE_MPI
    use mpi,only:MPI_Recv,MPI_CHARACTER,MPI_STATUS_SIZE,MPI_SUCCESS
#endif
    character(*),intent(inout)::line
    integer,intent(out)::ios
#ifdef USE_MPI
    integer::packet_status(MPI_STATUS_SIZE),ierr
    call MPI_Recv(line,len(line),MPI_CHARACTER,host_channel%controller,4101, &
                  host_channel%transport,packet_status,ierr)
    ios=0
    if(ierr/=MPI_SUCCESS)ios=1
#else
    ios=1
#endif
  end subroutine
  subroutine channel_observation(step,values)
#ifdef USE_MPI
    use mpi,only:MPI_Send,MPI_DOUBLE_PRECISION
#endif
    integer,intent(in)::step
    real(8),intent(in)::values(9)
#ifdef USE_MPI
    real(8)::packet(10)
    integer::ierr
    packet(1)=real(step,8);packet(2:10)=values
    call MPI_Send(packet,10,MPI_DOUBLE_PRECISION,host_channel%controller,4102,host_channel%transport,ierr)
#endif
  end subroutine
  subroutine channel_status(operation,status,step,updates)
#ifdef USE_MPI
    use mpi,only:MPI_Send,MPI_INTEGER
#endif
    integer,intent(in)::operation,status,step,updates
#ifdef USE_MPI
    integer::packet(4),ierr
    packet=[operation,status,step,updates]
    call MPI_Send(packet,4,MPI_INTEGER,host_channel%controller,4103,host_channel%transport,ierr)
#endif
  end subroutine

  ! Control input belongs to the caller-selected material communicator.
  ! Only its root owns the stream; no MPI_COMM_WORLD operation is introduced.
  subroutine worker_read_line(line,ios,stream_unit)
    use iso_fortran_env,only:input_unit
    use parallelization,only:nproc_id_global,nproc_group_global
    use communication,only:comm_is_root,comm_bcast
    character(*),intent(out)::line
    integer,intent(out)::ios
    integer,optional,intent(in)::stream_unit
    integer::unit
    line='';ios=0;unit=input_unit
    if(present(stream_unit))unit=stream_unit
    if(comm_is_root(nproc_id_global))then
      if(associated(host_channel))then
        call channel_read(line,ios)
      else
        read(unit,'(a)',iostat=ios)line
      endif
    endif
    call comm_bcast(ios,nproc_group_global)
    if(ios==0)then
      call comm_bcast(line,nproc_group_global)
    else
      line=''
    endif
  end subroutine

  logical function worker_performance_requested()
    character(16)::value
    integer::status
    ! One explicitly owned child process; never select embedded-host behavior.
    call get_environment_variable('SALMONPY_TDDFT_EXECUTION_MODE',value,status=status)
    worker_performance_requested=.not.worker_embedded.and..not.associated(host_channel).and. &
      status==0.and.trim(value)=='performance'
  end function
  logical function worker_norm_requested()
    character(8) :: value
    integer :: status
    call get_environment_variable("SALMONPY_TDDFT_NORM",value,status=status)
    worker_norm_requested=status==0.and.trim(value)=="1"
  end function
  logical function worker_requested()
    character(8) :: value
    integer :: status
    call get_environment_variable('SALMONPY_TDDFT_WORKER',value,status=status)
    worker_requested=associated(host_channel).or.worker_embedded.or.(status==0.and.( &
      trim(value)=='1'.or.trim(value)=='E3'.or.trim(value)=='E3N'.or. &
      trim(value)=='E3Y'.or.trim(value)=='E3M'))
  end function

  logical function worker_E_requested()
    character(8) :: value
    integer :: status
    call get_environment_variable('SALMONPY_TDDFT_WORKER',value,status=status)
    worker_E_requested=associated(host_channel).or.(.not.worker_embedded.and.status==0.and. &
      (trim(value)=='E3'.or.trim(value)=='E3N'.or.trim(value)=='E3Y'.or.trim(value)=='E3M'))
  end function

  logical function worker_nearfield_requested()
    character(8) :: value
    integer :: status
    call get_environment_variable('SALMONPY_TDDFT_WORKER',value,status=status)
    worker_nearfield_requested=associated(host_channel).or.(.not.worker_embedded.and.status==0.and. &
      (trim(value)=='E3N'.or.trim(value)=='E3Y'.or.trim(value)=='E3M'))
  end function

  logical function worker_solid_requested()
    character(8)::value
    integer::status
    call get_environment_variable('SALMONPY_TDDFT_WORKER',value,status=status)
    worker_solid_requested=.not.worker_embedded.and.status==0.and.trim(value)=='E3Y'
  end function

  logical function worker_near_requested()
    character(8)::value
    integer::status
    call get_environment_variable('SALMONPY_TDDFT_WORKER',value,status=status)
    worker_near_requested=associated(host_channel).or. &
      (.not.worker_embedded.and.status==0.and.trim(value)=='E3M')
  end function

  subroutine worker_near_load(rt)
    use structures,only:s_rt
    use salmonpy_solid_harmonics,only:solid_lmax
    use,intrinsic::ieee_arithmetic,only:ieee_is_finite
    type(s_rt),intent(inout)::rt
    character(64)::magic,units
    integer::u,ios,lmax,frames,n,k,d,p,allocation_status
    real(8)::geometry(7),step
    if(.not.worker_near_requested())return
    open(newunit=u,file='salmonpy-near-eh.dat',status='old',action='read',iostat=ios)
    if(ios/=0)error stop 'SALMONpy missing salmonpy-near-eh.dat'
    read(u,'(a)',iostat=ios)magic
    if(ios/=0)error stop 'SALMONpy invalid E/H field header'
    if(magic=='SALMONPY_NEAR_GRID_V1')then
      call worker_grid_load(rt,u)
      close(u,iostat=ios)
      if(ios/=0)error stop 'SALMONpy direct-grid input close failed'
      return
    endif
    if(magic/='SALMONPY_NEAR_EH_V3')error stop 'SALMONpy invalid E/H field header'
    read(u,*,iostat=ios)lmax,frames
    if(ios/=0)error stop 'SALMONpy invalid E/H dimensions'
    if(lmax<0.or.lmax>solid_lmax.or.int(frames,8)<int(rt%max_steps,8)+2_8) &
      error stop 'SALMONpy E/H field history is insufficient'
    read(u,*,iostat=ios)geometry,step
    if(ios/=0)error stop 'SALMONpy invalid E/H geometry'
    do k=1,7
      if(.not.ieee_is_finite(geometry(k)))error stop 'SALMONpy nonfinite E/H geometry'
    enddo
    if(.not.ieee_is_finite(step))error stop 'SALMONpy nonfinite E/H timestep'
    if(geometry(4)<=0d0.or.abs(step-rt%time_step)>1d-12*rt%time_step)error stop 'SALMONpy E/H timestep differs'
    read(u,'(a)',iostat=ios)units
    if(ios/=0.or.units/='E:au H:au bohr atomic_time vacuum')error stop 'SALMONpy E/H molecular vacuum units required'
    p=(lmax+1)*(lmax+2)*(lmax+3)/6
    allocate(rt%near_Ecoeff(3,p,0:rt%max_steps+1),rt%near_Hcoeff(3,p,0:rt%max_steps+1),stat=allocation_status)
    if(allocation_status/=0)error stop 'SALMONpy E/H allocation failed'
    rt%near_magnetic=.false.
    do n=0,rt%max_steps+1
      read(u,*,iostat=ios)rt%near_Ecoeff(:,:,n),rt%near_Hcoeff(:,:,n)
      if(ios/=0)error stop 'SALMONpy E/H coefficients are incomplete'
      do k=1,p
      do d=1,3
        if(.not.ieee_is_finite(rt%near_Ecoeff(d,k,n)))error stop 'SALMONpy nonfinite E coefficient'
        if(.not.ieee_is_finite(rt%near_Hcoeff(d,k,n)))error stop 'SALMONpy nonfinite H coefficient'
        if(rt%near_Hcoeff(d,k,n)/=0d0)rt%near_magnetic=.true.
      enddo
      enddo
    enddo
    close(u,iostat=ios)
    if(ios/=0)error stop 'SALMONpy E/H close failed'
    rt%near_geometry=geometry;rt%near_lmax=lmax
  end subroutine

  subroutine worker_grid_load(rt,u)
    use structures,only:s_rt
    use salmonpy_grid_field,only:grid_validate
    use,intrinsic::ieee_arithmetic,only:ieee_is_finite
    type(s_rt),intent(inout)::rt
    integer,intent(in)::u
    integer::ios,frames,k,n,allocation_status,status
    integer(8)::elements
    real(8)::step
    character(64)::units
    read(u,*,iostat=ios)rt%near_grid_counts,frames
    if(ios/=0)error stop 'SALMONpy invalid direct-grid dimensions'
    if(any(rt%near_grid_counts<2).or.int(frames,8)<int(rt%max_steps,8)+2_8) &
      error stop 'SALMONpy direct-grid history is insufficient'
    elements=6_8*(int(rt%max_steps,8)+2_8)
    do k=1,3
      if(elements>huge(elements)/rt%near_grid_counts(k))error stop 'SALMONpy direct-grid size overflow'
      elements=elements*rt%near_grid_counts(k)
    enddo
    read(u,*,iostat=ios)rt%near_geometry(1:3),rt%near_grid_center,step
    if(ios/=0)error stop 'SALMONpy invalid direct-grid molecular positions'
    do k=1,3
      if(.not.ieee_is_finite(rt%near_geometry(k)).or..not.ieee_is_finite(rt%near_grid_center(k))) &
        error stop 'SALMONpy nonfinite direct-grid molecular positions'
    enddo
    if(.not.ieee_is_finite(step))error stop 'SALMONpy nonfinite direct-grid timestep'
    if(abs(step-rt%time_step)>1d-12*rt%time_step)error stop 'SALMONpy direct-grid timestep differs'
    read(u,*,iostat=ios)rt%near_grid_origins
    if(ios/=0)error stop 'SALMONpy invalid E/H component origins'
    read(u,*,iostat=ios)rt%near_grid_spacing
    if(ios/=0)error stop 'SALMONpy invalid direct-grid spacing'
    read(u,'(a)',iostat=ios)units
    if(ios/=0.or.units/='E:au H:au bohr atomic_time vacuum')error stop 'SALMONpy direct-grid vacuum units required'
    allocate(rt%near_grid_values(3,rt%near_grid_counts(1),rt%near_grid_counts(2),rt%near_grid_counts(3), &
      2,0:rt%max_steps+1),stat=allocation_status)
    if(allocation_status/=0)error stop 'SALMONpy direct-grid allocation failed'
    rt%near_magnetic=.false.
    do n=0,rt%max_steps+1
      read(u,*,iostat=ios)rt%near_grid_values(:,:,:,:,:,n)
      if(ios/=0)error stop 'SALMONpy direct-grid field history is incomplete'
      call grid_validate(rt%near_grid_counts,rt%near_grid_origins,rt%near_grid_spacing, &
        rt%near_grid_values(:,:,:,:,:,n),status)
      if(status/=0)error stop 'SALMONpy invalid direct-grid fields or metadata'
      if(any(rt%near_grid_values(:,:,:,:,2,n)/=0d0))rt%near_magnetic=.true.
    enddo
  end subroutine

  subroutine worker_solid_load(rt)
    use structures,only:s_rt
    use salmonpy_solid_harmonics,only:solid_lmax
    use,intrinsic::ieee_arithmetic,only:ieee_is_finite
    type(s_rt),intent(inout)::rt
    character(32)::magic
    integer::u,ios,lmax,frames,n,k,p,allocation_status
    real(8)::center(3),radius,step
    if(.not.worker_solid_requested())return
    open(newunit=u,file='salmonpy-ylm.dat',status='old',action='read',iostat=ios)
    if(ios/=0)error stop 'SALMONpy missing salmonpy-ylm.dat'
    read(u,'(a)',iostat=ios)magic
    if(ios/=0.or.magic/='SALMONPY_SOLID_FIELD_V1')error stop 'SALMONpy invalid spatial field header'
    read(u,*,iostat=ios)lmax,frames
    if(ios/=0)error stop 'SALMONpy invalid spatial field dimensions'
    if(lmax<1.or.lmax>solid_lmax.or.frames<rt%max_steps+2)error stop 'SALMONpy spatial field history is insufficient'
    read(u,*,iostat=ios)center,radius,step
    if(ios/=0)error stop 'SALMONpy invalid spatial field geometry'
    do k=1,3
      if(.not.ieee_is_finite(center(k)))error stop 'SALMONpy nonfinite spatial field center'
    enddo
    if(.not.ieee_is_finite(radius).or..not.ieee_is_finite(step))error stop 'SALMONpy nonfinite spatial field scale'
    if(radius<=0d0.or.abs(step-rt%time_step)>1d-12*rt%time_step)error stop 'SALMONpy spatial field timestep differs'
    p=lmax*(lmax+2)
    allocate(rt%solid_coeff(p,0:rt%max_steps+1),stat=allocation_status)
    if(allocation_status/=0)error stop 'SALMONpy spatial field allocation failed'
    do n=0,rt%max_steps+1
      read(u,*,iostat=ios)rt%solid_coeff(:,n)
      if(ios/=0)error stop 'SALMONpy spatial field coefficients are incomplete'
      do k=1,p
        if(.not.ieee_is_finite(rt%solid_coeff(k,n)))error stop 'SALMONpy nonfinite spatial field coefficient'
      enddo
    enddo
    close(u,iostat=ios)
    if(ios/=0)error stop 'SALMONpy spatial field close failed'
    rt%solid_center=center;rt%solid_radius=radius;rt%solid_lmax=lmax
  end subroutine

  logical function worker_spin_supported()
    use salmon_global,only:spin,yn_spinorbit,xc
    worker_spin_supported=spin=='unpolarized'.and.yn_spinorbit=='n'
    if(worker_embedded.and.spin=='noncollinear'.and.yn_spinorbit=='y') &
      worker_spin_supported=xc=='pz'.or.xc=='libxc_pz'.or.xc=='libxc_pbe'
  end function

  subroutine worker_configure()
    use salmon_global
    use parallelization, only: nproc_size_global
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    character(32) :: command
    character(131072)::initial_line
    integer :: ios,i,j
    logical :: plus_u_file
    if(.not.worker_requested()) return
    if(.not.worker_embedded.and..not.associated(host_channel).and.nproc_size_global/=1) &
      error stop 'SALMONpy process worker requires one MPI rank'
    if(worker_E_requested()) then
      if(theory/='tddft_pulse'.or.yn_periodic/='n'.or.yn_md/='n') &
        error stop 'SALMONpy E worker requires isolated fixed-ion tddft_pulse'
      if(xc/='pz') error stop 'SALMONpy E worker requires PZ'
      if(worker_nearfield_requested())then
        if(num_dipole_source<0.or.num_dipole_source>2) error stop 'SALMONpy dipole source count must be 0..2'
        if(num_dipole_source>0)then
          if(.not.ieee_is_finite(rad_dipole_source)) error stop 'SALMONpy nonfinite dipole source radius'
          if(rad_dipole_source<=0d0) error stop 'SALMONpy positive dipole source radius required'
          do i=1,num_dipole_source
            do j=1,3
              if(.not.ieee_is_finite(vec_dipole_source(j,i))) error stop 'SALMONpy nonfinite dipole source vector'
              if(.not.ieee_is_finite(cood_dipole_source(j,i))) error stop 'SALMONpy nonfinite dipole source position'
            enddo
          enddo
          call check_source_pulse(ae_shape1,tw1,omega1,phi_CEP1)
          if(ae_shape2/='none')then
            call check_source_pulse(ae_shape2,tw2,omega2,phi_CEP2)
          else
            if(.not.ieee_is_finite(tw2)) error stop 'SALMONpy nonfinite second pulse width'
            if(tw2/=0d0) error stop 'SALMONpy absent second source pulse requires tw2=0'
          endif
          if(.not.ieee_is_finite(t1_t2)) error stop 'SALMONpy nonfinite source pulse delay'
        endif
      else
        if(num_dipole_source/=0) error stop 'SALMONpy dipole sources require explicit near-field mode'
      endif
    else
    if(theory/='tddft_pulse'.or.yn_periodic/='y'.or.yn_md/='n') &
      error stop 'SALMONpy TDDFT requires periodic fixed-ion tddft_pulse'
    end if
    if(propagator/='middlepoint'.or.trans_longi/='tr') &
      error stop 'SALMONpy TDDFT requires middlepoint and transverse supplied total A'
    if(yn_restart/='n'.or.yn_dc/='n'.or.yn_conventional_from_dcdft/='n') &
      error stop 'SALMONpy TDDFT requires a fresh conventional GS initial state'
    if(.not.worker_spin_supported().or.yn_jm/='n') &
      error stop 'SALMONpy TDDFT requires unpolarized scalar pseudopotentials'
    if(xname/='none'.or.cname/='none') error stop 'SALMONpy: separate exchange/correlation unsupported'
    select case(xc)
    case('pz','tbmbj')
    case('libxc_pz','libxc_pbe')
#ifndef USE_LIBXC
      error stop 'SALMONpy TDDFT: Libxc functional requires USE_LIBXC=ON'
#endif
    case default
      error stop 'SALMONpy TDDFT supports pz, tbmbj, libxc_pz and libxc_pbe only'
    end select
    if(projection_option/='no'.or.dm_unfold_option/='no') &
      error stop 'SALMONpy TDDFT worker does not support projection/unfolding'
    if(write_rt_wfn_k/='n') &
      error stop 'SALMONpy TDDFT worker does not support RT wavefunction checkpoints'
    inquire(file=resolve_path('hubbard_param.dat'),exist=plus_u_file)
    if(plus_u_file) error stop 'SALMONpy TDDFT worker does not support +U inputs'
    if(len_trim(sysname)==0.or.index(sysname,'/')/=0.or.index(sysname,achar(92))/=0.or. &
       index(sysname,':')/=0.or.index(sysname,'..')/=0) error stop 'SALMONpy: sysname must be a filename stem'
    call relative_file(file_kw)
#ifdef USE_OPENACC
    error stop 'SALMONpy TDDFT worker currently supports CPU execution'
#endif
    call relative_file(file_atom_coor)
    call relative_file(file_atom_red_coor)
    do i=1,nelem
      call relative_file(file_pseudo(i))
    end do
    if(worker_embedded) then
      dt=embedded_dt; nt=embedded_steps
      nproc_k=embedded_decomposition(1); nproc_ob=embedded_decomposition(2)
      nproc_rgrid=embedded_decomposition(3:5)
      ios=0; command='INIT'
    else
      call worker_read_line(initial_line,ios)
      if(ios==0)read(initial_line,*,iostat=ios)command,dt,nt,worker_initial_A
    end if
    if(worker_E_requested()) then
      if(worker_performance_requested())then
        if(ios/=0.or.trim(command)/='INIT_E_PERFORMANCE')error stop 'SALMONpy: expected INIT_E_PERFORMANCE dt capacity E'
      else
        if(ios/=0.or.trim(command)/='INIT_E') error stop 'SALMONpy: expected INIT_E dt capacity E'
      endif
    else if(worker_performance_requested())then
      if(ios/=0.or.trim(command)/='INIT_PERFORMANCE')error stop 'SALMONpy: expected INIT_PERFORMANCE dt capacity A'
    else
      if(ios/=0.or.trim(command)/='INIT') error stop 'SALMONpy: expected INIT dt capacity A'
    end if
    if(.not.ieee_is_finite(dt)) error stop 'SALMONpy: invalid dt'
    if(dt<=0d0.or.nt<1) error stop 'SALMONpy: invalid dt/capacity'
    do i=1,3
      if(.not.ieee_is_finite(worker_initial_A(i))) error stop 'SALMONpy: invalid initial A'
    end do
    if(associated(host_channel))then
      nproc_k=host_channel%decomposition(1);nproc_ob=host_channel%decomposition(2)
      nproc_rgrid=host_channel%decomposition(3:5)
    endif
    ! No external field or restart-output path can own the evolving state.
    if(worker_near_requested())yn_symmetrized_stencil='y'
    base_directory='./'
    directory_read_data='./data_for_restart/'
    if(.not.worker_nearfield_requested().or.num_dipole_source==0)then
      ae_shape1='none'; ae_shape2='none'
    endif
    checkpoint_interval=-1; time_shutdown=-1d0; write_rt_wfn_k='n'
  end subroutine

  subroutine check_source_pulse(shape,width,omega,phase)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    character(*),intent(in) :: shape
    real(8),intent(in) :: width,omega,phase
    if(shape/='Acos2'.and.shape/='Ecos2') error stop 'SALMONpy source pulse requires Acos2 or Ecos2'
    if(.not.ieee_is_finite(width)) error stop 'SALMONpy nonfinite source pulse width'
    if(.not.ieee_is_finite(omega)) error stop 'SALMONpy nonfinite source frequency'
    if(.not.ieee_is_finite(phase)) error stop 'SALMONpy nonfinite source phase'
    if(width<=0d0.or.omega<=0d0) error stop 'SALMONpy positive source pulse width and frequency required'
  end subroutine

  subroutine relative_file(path)
    character(*), intent(in) :: path
    if(len_trim(path)==0) return
    if(path(1:1)=='/'.or.index(path,'..')/=0.or.index(path,':')/=0.or.index(path,achar(92))/=0) &
      error stop 'SALMONpy: input assets must use relative bundle paths'
  end subroutine

  logical function worker_next_A(A,operation,rt,next_step)
    use structures,only:s_rt
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8), intent(out) :: A(3)
    integer,optional,intent(out)::operation
    type(s_rt),optional,intent(inout)::rt
    integer,optional,intent(in)::next_step
    real(8),allocatable::ce(:,:,:),ch(:,:,:)
    integer::p,n,k,d,side,allocation_status
    character(131072) :: line
    character(32) :: command
    integer :: ios,i
    worker_next_A=.false.
    if(present(operation))operation=0
    call worker_read_line(line,ios)
    if(ios<0) return
    if(ios/=0) error stop 'SALMONpy: command read failed'
    read(line,*,iostat=ios) command
    if(ios/=0) error stop 'SALMONpy: missing command'
    if(trim(command)=='CLOSE') return
    if(worker_E_requested().and.present(operation))then
      select case(trim(command))
      case('CAPTURE_E');operation=1
      case('RESTORE_E');operation=2
      case('RELEASE_E');operation=3
      case('RESUME_E');operation=4
      case('PROBE_E');operation=5
      case('SAVE_E');operation=6
      case('LOAD_E');operation=7
      case('DENSITY_E');operation=8
      case('GRID_E');operation=9
      case('ORBITAL_E');operation=10
      case('ORBITAL_AUDIT_E');operation=11
      case('OBSERVE_E');operation=12
      case('REPARTITION_E');operation=13
      end select
      if(operation/=0)then
        A=0d0;worker_next_A=.true.;return
      endif
    endif
    if(worker_near_requested().and.trim(command)=='STEP_EH')then
      if(present(rt))then
        if(allocated(rt%near_grid_values))error stop 'SALMONpy direct-grid histories use STEP_E'
      endif
      if(.not.present(rt).or..not.present(next_step))error stop 'SALMONpy spatial step requires owning state'
      p=(rt%near_lmax+1)*(rt%near_lmax+2)*(rt%near_lmax+3)/6
      allocate(ce(3,p,2),ch(3,p,2),stat=allocation_status)
      if(allocation_status/=0)error stop 'SALMONpy spatial command allocation failed'
      read(line,*,iostat=ios)command,A,n,ce(:,:,1),ce(:,:,2),ch(:,:,1),ch(:,:,2)
      if(ios/=0.or.n/=next_step)error stop 'SALMONpy spatial coefficient step differs'
      if(n<1.or.n>rt%max_steps)error stop 'SALMONpy spatial coefficient step outside capacity'
      do side=1,2
      do k=1,p
      do d=1,3
        if(.not.ieee_is_finite(ce(d,k,side)).or..not.ieee_is_finite(ch(d,k,side))) &
          error stop 'SALMONpy nonfinite spatial coefficient'
      enddo
      enddo
      enddo
      do d=1,3
        if(.not.ieee_is_finite(A(d)))error stop 'SALMONpy invalid additional uniform E'
      enddo
      rt%near_Ecoeff(:,:,n:n+1)=ce;rt%near_Hcoeff(:,:,n:n+1)=ch
      ! A live stream may introduce H after an initially electric-only table.
      rt%near_magnetic=.true.
      worker_next_A=.true.;return
    endif
    read(line,*,iostat=ios) command,A
    if(worker_E_requested()) then
      if(ios/=0.or.trim(command)/='STEP_E') error stop 'SALMONpy: expected STEP_E E'
    else
      if(ios/=0.or.trim(command)/='STEP') error stop 'SALMONpy: expected STEP A'
    end if
    if(.not.worker_performance_requested())then
    do i=1,3
      if(.not.ieee_is_finite(A(i))) error stop 'SALMONpy: invalid A'
    end do
    endif
    worker_next_A=.true.
  end function

  subroutine worker_state_reply(operation,status,step,updates)
    use parallelization,only:nproc_id_global
    use communication,only:comm_is_root,comm_get_max
    use parallelization,only:nproc_group_global
    use iso_fortran_env,only:output_unit
    integer,intent(in)::operation,status,step,updates
    integer::agreed
    agreed=status
    if(associated(host_channel))call comm_get_max(agreed,nproc_group_global)
    if(.not.comm_is_root(nproc_id_global))return
    if(associated(host_channel))then
      call channel_status(operation,agreed,step,updates);return
    endif
    write(output_unit,'(a,4(1x,i0))')'SALMONPY_E_STATE_V1',operation,status,step,updates
    flush(output_unit)
  end subroutine

  subroutine worker_reply(step,J)
    use parallelization,only:nproc_id_global
    use communication,only:comm_is_root
    use iso_fortran_env, only: output_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer,intent(in) :: step
    real(8),intent(in) :: J(3)
    integer :: i
    character(32)::prefix
    prefix="SALMONPY_TDDFT_V1"
    if(worker_performance_requested())prefix="SALMONPY_TDDFT_PERFORMANCE_V1"
    if(.not.worker_performance_requested())then
    do i=1,3
      if(.not.ieee_is_finite(J(i)))then
        write(output_unit,'(a)') 'SALMONPY_NUMERICAL_FAILURE nonfinite current'
        flush(output_unit)
        error stop 'SALMONpy: nonfinite current'
      endif
    end do
    endif
    if(.not.comm_is_root(nproc_id_global))return
    if(worker_norm_requested())then
      write(output_unit,'(a,1x,i0,4(1x,es26.17e3))') trim(prefix),step,J,worker_electron_number
    else
      write(output_unit,'(a,1x,i0,3(1x,es26.17e3))') trim(prefix),step,J
    endif
    flush(output_unit)
  end subroutine
  subroutine worker_reply_E(step,delta_dipole,total_dipole,electron_number,internal_energy,energy_fresh)
    use parallelization,only:nproc_id_global
    use communication,only:comm_is_root
    use iso_fortran_env, only: output_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    integer,intent(in) :: step
    real(8),intent(in) :: delta_dipole(3),total_dipole(3),electron_number,internal_energy
    logical,intent(in)::energy_fresh
    integer :: i
    character(40)::prefix
    prefix="SALMONPY_TDDFT_E_V3"
    if(worker_performance_requested())prefix="SALMONPY_TDDFT_E_PERFORMANCE_V3"
    if(.not.worker_performance_requested())then
    do i=1,3
      if(.not.ieee_is_finite(delta_dipole(i)).or..not.ieee_is_finite(total_dipole(i))) &
        error stop 'SALMONpy: nonfinite dipole'
    end do
    if(.not.ieee_is_finite(electron_number))error stop 'SALMONpy: nonfinite electron number'
    if(.not.ieee_is_finite(internal_energy))error stop 'SALMONpy: nonfinite internal energy'
    endif
    if(.not.comm_is_root(nproc_id_global))return
    if(associated(host_channel))then
      call channel_observation(step,[delta_dipole,total_dipole,electron_number,internal_energy, &
        merge(1d0,0d0,energy_fresh)])
      return
    endif
    write(output_unit,'(a,1x,i0,9(1x,es26.17e3))') trim(prefix), &
      step,delta_dipole,total_dipole,electron_number,internal_energy,merge(1d0,0d0,energy_fresh)
    flush(output_unit)
  end subroutine
end module
