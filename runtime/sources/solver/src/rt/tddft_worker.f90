#include "config.h"
! Optional SALMONpy process-isolated endpoint-A protocol; ordinary runs are unchanged.
module tddft_worker
  use filesystem,only:resolve_path
  implicit none
  private
  public :: worker_requested, worker_configure, worker_initial_A, worker_next_A, worker_reply
  public :: worker_E_requested, worker_reply_E,worker_norm_requested,worker_electron_number
  ! Embedded mode is set only by the explicit host-MPI ABI, once per process.
  logical,public,save :: worker_embedded=.false.
  real(8),public,save :: embedded_dt
  integer,public,save :: embedded_steps,embedded_decomposition(5)
  real(8), save :: worker_initial_A(3)=0d0
  real(8),save :: worker_electron_number=0d0
contains
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
    worker_requested=worker_embedded.or.(status==0.and.(trim(value)=='1'.or.trim(value)=='E3'))
  end function

  logical function worker_E_requested()
    character(8) :: value
    integer :: status
    call get_environment_variable('SALMONPY_TDDFT_WORKER',value,status=status)
    worker_E_requested=.not.worker_embedded.and.status==0.and.trim(value)=='E3'
  end function

  subroutine worker_configure()
    use salmon_global
    use parallelization, only: nproc_size_global
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    character(32) :: command
    integer :: ios,i
    logical :: plus_u_file
    if(.not.worker_requested()) return
    if(.not.worker_embedded.and.nproc_size_global/=1) error stop 'SALMONpy TDDFT requires one MPI rank'
    if(worker_E_requested()) then
      if(theory/='tddft_pulse'.or.yn_periodic/='n'.or.yn_md/='n') &
        error stop 'SALMONpy E worker requires isolated fixed-ion tddft_pulse'
      if(xc/='pz'.or.num_dipole_source/=0) error stop 'SALMONpy E worker requires PZ without dipole sources'
    else
    if(theory/='tddft_pulse'.or.yn_periodic/='y'.or.yn_md/='n') &
      error stop 'SALMONpy TDDFT requires periodic fixed-ion tddft_pulse'
    end if
    if(propagator/='middlepoint'.or.trans_longi/='tr') &
      error stop 'SALMONpy TDDFT requires middlepoint and transverse supplied total A'
    if(yn_restart/='n'.or.yn_dc/='n'.or.yn_conventional_from_dcdft/='n') &
      error stop 'SALMONpy TDDFT requires a fresh conventional GS initial state'
    if(yn_spinorbit/='n'.or.yn_jm/='n'.or.spin/='unpolarized') &
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
      read(*,*,iostat=ios) command,dt,nt,worker_initial_A
    end if
    if(worker_E_requested()) then
      if(ios/=0.or.trim(command)/='INIT_E') error stop 'SALMONpy: expected INIT_E dt capacity E'
    else
      if(ios/=0.or.trim(command)/='INIT') error stop 'SALMONpy: expected INIT dt capacity A'
    end if
    if(.not.ieee_is_finite(dt)) error stop 'SALMONpy: invalid dt'
    if(dt<=0d0.or.nt<1) error stop 'SALMONpy: invalid dt/capacity'
    do i=1,3
      if(.not.ieee_is_finite(worker_initial_A(i))) error stop 'SALMONpy: invalid initial A'
    end do
    ! No external field or restart-output path can own the evolving state.
    base_directory='./'
    directory_read_data='./data_for_restart/'
    ae_shape1='none'; ae_shape2='none'
    checkpoint_interval=-1; time_shutdown=-1d0; write_rt_wfn_k='n'
  end subroutine

  subroutine relative_file(path)
    character(*), intent(in) :: path
    if(len_trim(path)==0) return
    if(path(1:1)=='/'.or.index(path,'..')/=0.or.index(path,':')/=0.or.index(path,achar(92))/=0) &
      error stop 'SALMONpy: input assets must use relative bundle paths'
  end subroutine

  logical function worker_next_A(A)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8), intent(out) :: A(3)
    character(512) :: line
    character(32) :: command
    integer :: ios,i
    worker_next_A=.false.
    read(*,'(a)',iostat=ios) line
    if(ios<0) return
    if(ios/=0) error stop 'SALMONpy: command read failed'
    read(line,*,iostat=ios) command
    if(ios/=0) error stop 'SALMONpy: missing command'
    if(trim(command)=='CLOSE') return
    read(line,*,iostat=ios) command,A
    if(worker_E_requested()) then
      if(ios/=0.or.trim(command)/='STEP_E') error stop 'SALMONpy: expected STEP_E E'
    else
      if(ios/=0.or.trim(command)/='STEP') error stop 'SALMONpy: expected STEP A'
    end if
    do i=1,3
      if(.not.ieee_is_finite(A(i))) error stop 'SALMONpy: invalid A'
    end do
    worker_next_A=.true.
  end function

  subroutine worker_reply(step,J)
    use iso_fortran_env, only: output_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer,intent(in) :: step
    real(8),intent(in) :: J(3)
    integer :: i
    do i=1,3
      if(.not.ieee_is_finite(J(i)))then
        write(output_unit,'(a)') 'SALMONPY_NUMERICAL_FAILURE nonfinite current'
        flush(output_unit)
        error stop 'SALMONpy: nonfinite current'
      endif
    end do
    if(worker_norm_requested())then
      write(output_unit,'(a,1x,i0,4(1x,es26.17e3))') 'SALMONPY_TDDFT_V1',step,J,worker_electron_number
    else
      write(output_unit,'(a,1x,i0,3(1x,es26.17e3))') 'SALMONPY_TDDFT_V1',step,J
    endif
    flush(output_unit)
  end subroutine
  subroutine worker_reply_E(step,delta_dipole,total_dipole,electron_number,internal_energy,energy_fresh)
    use iso_fortran_env, only: output_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    integer,intent(in) :: step
    real(8),intent(in) :: delta_dipole(3),total_dipole(3),electron_number,internal_energy
    logical,intent(in)::energy_fresh
    integer :: i
    do i=1,3
      if(.not.ieee_is_finite(delta_dipole(i)).or..not.ieee_is_finite(total_dipole(i))) &
        error stop 'SALMONpy: nonfinite dipole'
    end do
    if(.not.ieee_is_finite(electron_number))error stop 'SALMONpy: nonfinite electron number'
    if(.not.ieee_is_finite(internal_energy))error stop 'SALMONpy: nonfinite internal energy'
    write(output_unit,'(a,1x,i0,9(1x,es26.17e3))') 'SALMONPY_TDDFT_E_V3', &
      step,delta_dipole,total_dipole,electron_number,internal_energy,merge(1d0,0d0,energy_fresh)
    flush(output_unit)
  end subroutine
end module
