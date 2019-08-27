program main !test
  use salmon_global
  use salmon_parallel
  use inputoutput
  use math_constants
  use timer
  implicit none

  call set_math_constants

  call setup_parallel
  if (nproc_id_global == 0) then
    call print_software_version
  endif

  call read_input

  call timer_initialize

  select case(theory)
  case('TDDFT')
    select case(iperiodic)
    case(0)
      call gceed
    case(3)
      select case(yn_domain_parallel)
      case('y')
        call gceed
      case('n')
        call arted
      case default
        stop 'invalid yn_domain_parallel'
      end select
    case default
      stop 'invalid iperiodic'
    end select
  case('Maxwell')
    call classic_em
  case default
    stop 'invalid theory'
  end select

  call write_perflog_csv

  call end_parallel
contains
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
    use salmon_file, only: get_filehandle
    use salmon_parallel, only: nproc_id_global
    use salmon_communication, only: comm_is_root
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
