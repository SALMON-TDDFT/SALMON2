program main
  use salmon_global
  use salmon_parallel
  use salmon_communication, only: comm_is_root
  use inputoutput
  use math_constants
  use timer
  implicit none
  character(32)  :: theory_org


  call set_math_constants
  call setup_parallel
  if (nproc_id_global == 0) then
    call print_software_version
  endif

  call read_input

  !convert old keyword of "theory" to new keyword (if it is old)
  theory_org = theory
  if(theory=="TDDFT ") then
     call convert_theory_to_new_keyword(theory_org,theory)
     if(comm_is_root(nproc_id_global)) then
        write(*,'(a)') "# theory keyword was converted to a new one:"
        write(*,'(2a)')"# theory=",trim(theory)
     endif
  else if(theory.ne."Maxwell ") then
     call set_old_input_keyword(theory)  !remove later...
  endif

  call timer_initialize

  !ARTED: (legacy: only in the case of iperiodic=3 + domain parallel=y)
  select case(yn_domain_parallel)  
  case('n')
     select case(iperiodic)
     case(3) 
        if( theory/='Maxwell' ) then
          call arted
          stop
        end if
     end select
  end select

  !GCEED: (main)
  select case(theory)
  case('DFT')                         ; call main_dft
  case('DFT_MD')                      ; call arted   !temporally
  case('TDDFT_response','TDDFT_pulse'); call main_tddft
  case('Single_scale_Maxwell_TDDFT'  ); call main_tddft
  case('Multi_scale_Maxwell_TDDFT'   ); call arted   !temporally
  case('Maxwell')                     ; call classic_em
 !case('SBE')                         ; call main_sbe
 !case('Maxwell_SBE')                 ; call main_maxwell_sbe
 !case('TTM')                         ; call main_ttm
 !case('Maxwell_TTM')                 ; call main_maxwell_ttm
  case default ; stop 'invalid theory'
  end select


!===== Theory branch of OLDER VERSION =====
!  select case(theory)
!  case('TDDFT')
!     select case(iperiodic)
!     case(0);  call gceed
!     case(3)
!        select case(yn_domain_parallel)
!        case('y')   ; call gceed
!        case('n')   ; call arted
!        case default; stop 'invalid yn_domain_parallel'
!        end select
!     case default; stop 'invalid iperiodic'
!     end select
!  case('Maxwell'); call classic_em
!  case default   ; stop 'invalid theory'
!  end select

  call write_perflog_csv

  call end_parallel
contains

  subroutine convert_theory_to_new_keyword(theory_org,theory)
    implicit none
    character(32)  :: theory_org, theory

    if(theory_org=="TDDFT   ") then

       select case(calc_mode)
       case('GS')
          select case(use_adiabatic_md)
          case('y') ; theory="DFT_MD"
          case('n') ; theory="DFT"
          end select
       case('RT')
          select case(use_ms_maxwell)
          case('y') ; theory="Multi_scale_Maxwell_TDDFT"
          case('n') 
             select case(use_singlescale)
             case('y') ; theory="Single_scale_Maxwell_TDDFT"
             case('n') 
                select case(ae_shape1)
                case('impulse') ; theory="TDDFT_response"
                case default    ; theory="TDDFT_pulse"
                end select
             end select
          end select
       case('GS_RT') 
          theory="DFT_TDDFT" !legacy-- this is not supported officially now
          write(*,*) "calc_mode=GS_RT is not supported now!!!"
       end select

       select case(use_ehrenfest_md)
       case('y') ; yn_md='y'
       end select

       select case(use_geometry_opt)
       case('y') ; yn_opt='y'
       end select


    endif
  end subroutine convert_theory_to_new_keyword

  subroutine set_old_input_keyword(theory)
    character(32)  :: theory

    if(theory.ne."TDDFT " .and. theory.ne."Maxwell ") then

       select case(theory)
       case('DFT')
          calc_mode='GS'
       case('DFT_MD')
          yn_md='y'
          calc_mode='GS'
          use_adiabatic_md='y'
       case('TDDFT_response','TDDFT_pulse')
          calc_mode='RT'
       case('Multi_scale_Maxwell_TDDFT')
          calc_mode='RT'
          use_ms_maxwell='y'
       case('Single_scale_Maxwell_TDDFT')
          calc_mode='RT'
          use_singlescale='y'
       case('DFT_TDDFT')
          calc_mode='GS_RT'  !legacy-- this is not supported officially now
          write(*,*) "theory=DFT_TDDFT is not supported officially !!"
       end select

       select case(yn_md)
       case('y') ; use_ehrenfest_md='y'
       end select

       select case(yn_opt)
       case('y') ; use_geometry_opt='y'
       end select

    endif
  end subroutine set_old_input_keyword

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
