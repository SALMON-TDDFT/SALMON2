subroutine main_ssbe(icomm)
    ! use mpi
    use omp_lib
    use communication
    use multiscale_ssbe
    use realtime_ssbe
    use salmon_global
    implicit none
    integer, intent(in) :: icomm

    select case(trim(theory))
    case ("sbe")
        call main_realtime_ssbe(icomm)
    case ("maxwell_sbe")
        call main_multiscale_ssbe(icomm)
    end select

    call comm_sync_all(icomm)

    return
end subroutine main_ssbe 
