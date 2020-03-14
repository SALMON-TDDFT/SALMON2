program main
  use salmon_global
  use parallelization
  use communication, only: comm_is_root
  use inputoutput
  use math_constants
  use timer
  use plusU_global, only: read_Hubbard_parameters
  implicit none


  call set_math_constants
  call setup_parallel
  if (nproc_id_global == 0) then
    call print_software_version
  endif

  call read_input
  call set_basic_flag(theory)

  call read_Hubbard_parameters  !should move into read_input subroutine!

  call timer_initialize

  !ARTED: (legacy: only in the case of iperiodic=3 + domain parallel=y)
  select case(yn_domain_parallel)  
  case('n')
     select case(iperiodic)
     case(3) 
        if( theory/='maxwell' ) then
          call arted
          stop
        end if
     end select
  end select

  !GCEED: (main)
  select case(theory)
  case('dft','dft_band')              ; call main_dft
  case('dft_md')                      ; call main_dft_md
  case('tddft_response','tddft_pulse'); call main_tddft
  case('dft2tddft')                   ; call main_dft2tddft ! DFT data redistributor to use TDDFT
  case('dft_k_expand')                ; call main_dft_k_expand !convert DFT/k-points data to supercell/gammma DFT
  case('single_scale_maxwell_tddft'  ); call main_tddft
  case('multi_scale_maxwell_tddft'   ); call arted      !temporally
  case('multiscale_experiment' )      ; call main_ms    ! experimental
  case('maxwell')                     ; call main_maxwell
 !case('sbe')                         ; call main_sbe
 !case('maxwell_sbe')                 ; call main_maxwell_sbe
 !case('ttm')                         ; call main_ttm
 !case('maxwell_ttm')                 ; call main_maxwell_ttm
  case default ; stop 'invalid theory'
  end select


  call write_perflog_csv
  
  if (nproc_id_global == 0) print '(A)',"end SALMON"

  call end_parallel
contains

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
       use_adiabatic_md='y'
    case('tddft_response','tddft_pulse')
       calc_mode='RT'
    case('multi_scale_maxwell_tddft')
       calc_mode='RT'
       use_ms_maxwell='y'
    case('single_scale_maxwell_tddft')
       calc_mode='RT'
       use_singlescale='y'
    case('dft_tddft')
       calc_mode='GS_RT'  !legacy-- this is not supported officially now
       write(*,*) "theory=dft_tddft is not supported officially !!"
    end select

    select case(yn_md)
    case('y') ; use_ehrenfest_md='y'
    end select

    select case(yn_opt)
    case('y') ; use_geometry_opt='y'
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
    if (GIT_FOUND) then 
      print '(A)',       '#'
      print '(A,A,A,A)', '#   [Git revision] ', GIT_COMMIT_HASH, ' in ', GIT_BRANCH
    endif
    print '(A)',         '##############################################################################'
    
    call print_xc_info()    
  end subroutine

  subroutine write_perflog_csv
    use perflog
    use misc_routines, only: gen_logfilename
    use filesystem, only: get_filehandle
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use iso_fortran_env, only: output_unit
    implicit none
    integer :: fh

    if (comm_is_root(nproc_id_global)) then
      fh = get_filehandle()
      open(fh, file=gen_logfilename('perflog','csv'))
    end if

    call write_performance(fh,write_mode_csv)

    if (comm_is_root(nproc_id_global)) then
      close(fh)
    end if
  end subroutine
end program main
