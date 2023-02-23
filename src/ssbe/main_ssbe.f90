subroutine main_ssbe(icomm)
    ! use mpi
    use omp_lib
    use communication
    use multiscale_ssbe
    use realtime_ssbe
    use salmon_global
    implicit none
    integer, intent(in) :: icomm

    select case (trim(theory))
    case ("sbe")
        call realtime_main_ssbe(icomm)
    case ("maxwell_sbe")
        call multiscale_main_ssbe(icomm)
    end select

    stop "Bye!"
end subroutine main_ssbe 
