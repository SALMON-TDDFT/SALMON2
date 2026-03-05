module util_ssbe
    implicit none
contains
! Almost equally split integer n over array of size m...
subroutine split_num(n, m, ntbl)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(out) :: ntbl(m)
    integer :: ntmp, itmp
    ntmp = n / m
    itmp = mod(n, m)
    ntbl(1:itmp) = ntmp + 1
    ntbl(itmp+1:m) = ntmp
    return
end subroutine split_num

! Almost equally split range of imin:imax to m sub-regions...
subroutine split_range(imin, imax, m, itbl_min, itbl_max)
    implicit none
    integer, intent(in) :: imin
    integer, intent(in) :: imax
    integer, intent(in) :: m
    integer, intent(out) :: itbl_min(m)
    integer, intent(out) :: itbl_max(m)
    integer :: ntbl(m), itmp, i
    call split_num(imax-imin+1, m, ntbl)
    itmp = imin
    do i = 1, m
        itbl_min(i) = itmp
        itmp = itmp + ntbl(i)
        itbl_max(i) = itmp - 1
    end do
    return
end subroutine split_range
end module util_ssbe

