#include "config.h"
  program main
  use salmon_global
  use parallelization
  use communication, only: comm_is_root
  use inputoutput
  use math_constants
  use timer
  use plusU_global, only: read_Hubbard_parameters
  use sym_sub, only: read_sw_symmetry
  implicit none


  call set_math_constants
  call setup_parallel
  if (nproc_id_global == 0) then
    call print_software_version
  endif

  call read_input
#if defined(USE_HSE) && defined(USE_MPI)
  call start_exchange_helpers
#endif
  call set_basic_flag(theory)

  call read_Hubbard_parameters  !should move into read_input subroutine!
  call read_sw_symmetry( yn_symmetry )

  call timer_initialize

  quiet = .false. ! enable to output message into stdout

  !GCEED: (main)
  if(nproc_id_global==0) write(*,*)"  theory= ", trim(theory)
  select case(theory)
  case('dft','dft_band')              ; call main_dft
  case('dft_md')                      ; call main_dft_md
  case('tddft_response','tddft_pulse'); call main_tddft
  case('dft_k_expand')                ; call main_dft_k_expand !convert DFT/k-points data to supercell/gammma DFT
  case('single_scale_maxwell_tddft'  ); call main_tddft
  case('multi_scale_maxwell_tddft'   ); call main_ms
  case('maxwell')                     ; call main_maxwell
  case('sbe', 'maxwell_sbe')          ; call main_ssbe(nproc_group_global)
 !case('maxwell_sbe')                 ; call main_maxwell_sbe
 !case('ttm')                         ; call main_ttm
 !case('maxwell_ttm')                 ; call main_maxwell_ttm
  case default ; stop 'invalid theory @ main'
  end select

  if (yn_out_perflog == 'y') call write_perflog

  if (nproc_id_global == 0) print '(A)',"end SALMON"

#if defined(USE_HSE) && defined(USE_MPI)
  call stop_exchange_helpers
#endif
  call end_parallel
contains

#if defined(USE_HSE) && defined(USE_MPI)
  subroutine start_exchange_helpers
    use iso_fortran_env, only: output_unit
    use communication, only: comm_select_application
    use exx_worker_service, only: service_enter
    use mpi, only: MPI_COMM_WORLD,MPI_INTEGER,MPI_MAX,MPI_Allreduce,MPI_Abort,MPI_Barrier
    implicit none
    integer :: bad,total_bad,ierr
    logical :: supported
    ! Ordinary orbital/spatial layouts retain their original decomposition.
    if(nproc_k<1.or.nproc_k>=nproc_size_global.or.nproc_ob/=1.or.any(nproc_rgrid/=1))return
    supported=theory=='dft'.and.xc=='hse06'.and.yn_spinorbit=='y'.and.spin=='noncollinear'.and. &
      yn_dc=='n'.and.yn_md=='n'.and.yn_opt=='n'.and.yn_periodic=='y'.and.exx_factor_mode=='off'.and. &
      propagator/='hse_taylor4_full'.and.yn_exx_spinor_pair_fft=='y'.and.index(yn_symmetry,'y')/=0.and. &
      all(hse_q_reduction==1).and.exx_kpoint_backend=='cpu'.and.temperature<0d0.and. &
      yn_jm=='n'.and.yn_symmetrized_stencil/='y'.and.yn_hse_wannier_snapshot/='y'.and. &
      exx_pair_screening=='off'.and.dm_unfold_option=='no'.and. &
      exx_mlwf_radius==0d0.and.exx_mlwf_norm_fraction==0d0.and.hse_sr_radius==0d0.and.hse_sr_tolerance==0d0.and. &
      al_vec1(2)==0d0.and.al_vec1(3)==0d0.and.al_vec2(1)==0d0.and. &
      al_vec2(3)==0d0.and.al_vec3(1)==0d0.and.al_vec3(2)==0d0
    supported=supported.and.modulo(nproc_size_global,exx_worker_rprime_groups)==0
    bad=merge(0,1,supported)
    call MPI_Allreduce(bad,total_bad,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(total_bad/=0)then
      if(nproc_id_global==0)print *, &
        'EXX_SERVICE_UNSUPPORTED: helpers require static full-support CPU SOC HSE06 pair exchange with symmetry'
      flush(output_unit)
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      call MPI_Abort(MPI_COMM_WORLD,86,ierr)
      return
    endif
    call comm_select_application(nproc_k)
    call service_enter(nproc_k,exx_worker_tile_divisor,exx_worker_rprime_groups)
    ! Only SCF owners return; refresh their rank context and clock offset.
    call setup_parallel
  end subroutine start_exchange_helpers

  subroutine stop_exchange_helpers
    use exx_worker_service, only: service_stop
    implicit none
    call service_stop
  end subroutine stop_exchange_helpers
#endif


  subroutine set_basic_flag(theory)
    implicit none
    character(32)  :: theory

    if(theory=="maxwell ") return

    select case(theory)
    case('dft','dft2tddft')
       calc_mode='GS'
    case('dft_md')
       yn_md='y'
       calc_mode='GS'
    case('dft_band')
       calc_mode='DFT_BAND'
    case('tddft_response','tddft_pulse')
       calc_mode='RT'
    case('multi_scale_maxwell_tddft')
       calc_mode='RT'
    case('single_scale_maxwell_tddft')
       calc_mode='RT'
    case('dft_tddft')
       calc_mode='GS_RT'  !legacy-- this is not supported officially now
       write(*,*) "theory=dft_tddft is not supported officially !!"
    end select

  end subroutine set_basic_flag

  subroutine print_software_version
    use salmon_xc, only: print_xc_info
    implicit none
    include 'versionf.h'
    print '(A)',         '##############################################################################'
    print '(A)',         '# SALMON: Scalable Ab-initio Light-Matter simulator for Optics and Nanoscience'
    print '(A)',         '#'
    print '(A,I1,".",I1,".",I1)', &
    &                    '#                             Version ', SALMON_VER_MAJOR, SALMON_VER_MINOR, SALMON_VER_MICRO
    print '(A)',         '##############################################################################'
    
    call print_xc_info()    
  end subroutine

  subroutine write_perflog
    use misc_routines, only: gen_logfilename
    use filesystem, only: get_filehandle
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use iso_fortran_env, only: output_unit
    use perflog
    implicit none
    integer :: fh, imode
    logical :: is_opened

    fh = -1

    select case(format_perflog)
    case default
      fh    = output_unit
      imode = write_mode_readable
    case('text')
      if (comm_is_root(nproc_id_global)) then
        fh = get_filehandle()
        open(fh, file=gen_logfilename('perflog','txt'))
      end if
      imode = write_mode_readable
    case('csv')
      if (comm_is_root(nproc_id_global)) then
        fh = get_filehandle()
        open(fh, file=gen_logfilename('perflog','csv'))
      end if
      imode = write_mode_csv
    end select

    call write_performance(fh,imode)

    if (fh /= output_unit) then
      inquire(fh, opened=is_opened)
      if (is_opened) close(fh)
    end if
  end subroutine write_perflog
end program main
