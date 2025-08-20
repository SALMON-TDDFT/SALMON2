subroutine main_ssbe(icomm)
    ! use mpi
    use omp_lib
    use communication
    use multiscale_vg_ssbe
    use realtime_vg_ssbe
    use salmon_global
    implicit none
    integer, intent(in) :: icomm

    select case(trim(theory))
    case ("vg_sbe")
        call main_realtime_vg_ssbe(icomm)
    case ("maxwell_vg_sbe")
        call main_multiscale_vg_ssbe(icomm)
    end select

    call comm_sync_all(icomm)

    return
end subroutine main_ssbe 
