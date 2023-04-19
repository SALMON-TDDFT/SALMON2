module shaper_ssbe
    implicit none

contains

subroutine generate_macropoint(fh)
    use salmon_global
    implicit none
    integer, intent(in) :: fh

    integer :: ishape, m_id, ix, iy, iz
    integer :: imedia_grid(1:nx_m, 1:ny_m, 1:nz_m)
    integer :: ibuff(4, nx_m*ny_m*nz_m), icount, npoint

    imedia_grid(:, :, :) = 0

    do ishape = 1, n_s
        select case (trim(media_type(id_s(ishape))))
        case ("multiscale", "ms")
            m_id = -1
        case ("vacuum")
            m_id = 0
        case ("constant media")
            m_id = id_s(ishape)
        case default
            stop "Unsupported media"
        end select

        select case(trim(typ_s(ishape)))
        case('ellipsoid')
            call fill_ellipsoid(ori_s(ishape, 1), ori_s(ishape, 2), ori_s(ishape, 3), &
                & inf_s(ishape, 1), inf_s(ishape, 2), inf_s(ishape, 3), inf_s(ishape, 4), inf_s(ishape, 5), m_id)
        case('half-ellipsoid')
            call fill_half_ellipsoid(ori_s(ishape, 1), ori_s(ishape, 2), ori_s(ishape, 3), &
                & inf_s(ishape, 1), inf_s(ishape, 2), inf_s(ishape, 3), inf_s(ishape, 4), inf_s(ishape, 5), m_id)
        case('elliptic-cylinder')
            call fill_elliptic_cylinder(ori_s(ishape, 1), ori_s(ishape, 2), ori_s(ishape, 3), &
                & inf_s(ishape, 1), inf_s(ishape, 2), inf_s(ishape, 3), inf_s(ishape, 4), inf_s(ishape, 5), m_id)
        case('triangular-cylinder')
            call fill_triangular_cylinder(ori_s(ishape, 1), ori_s(ishape, 2), ori_s(ishape, 3), &
                & inf_s(ishape, 1), inf_s(ishape, 2), inf_s(ishape, 3), inf_s(ishape, 4), inf_s(ishape, 5), m_id)
        case('rectangular-cylinder')
            call fill_rectangular_cylinder(ori_s(ishape, 1), ori_s(ishape, 2), ori_s(ishape, 3), &
                & inf_s(ishape, 1), inf_s(ishape, 2), inf_s(ishape, 3), inf_s(ishape, 4), inf_s(ishape, 5), m_id)
        case('elliptic-cone')
            call fill_elliptic_cone(ori_s(ishape, 1), ori_s(ishape, 2), ori_s(ishape, 3), &
                & inf_s(ishape, 1), inf_s(ishape, 2), inf_s(ishape, 3), inf_s(ishape, 4), inf_s(ishape, 5), m_id)
        case('triangular-cone')
            call fill_triangular_cone(ori_s(ishape, 1), ori_s(ishape, 2), ori_s(ishape, 3), &
                & inf_s(ishape, 1), inf_s(ishape, 2), inf_s(ishape, 3), inf_s(ishape, 4), inf_s(ishape, 5), m_id)
        case('rectangular-cone')
            call fill_rectangular_cone(ori_s(ishape, 1), ori_s(ishape, 2), ori_s(ishape, 3), &
                & inf_s(ishape, 1), inf_s(ishape, 2), inf_s(ishape, 3), inf_s(ishape, 4), inf_s(ishape, 5), m_id)
        case('elliptic-ring')
            call fill_elliptic_ring(ori_s(ishape, 1), ori_s(ishape, 2), ori_s(ishape, 3), &
                & inf_s(ishape, 1), inf_s(ishape, 2), inf_s(ishape, 3), inf_s(ishape, 4), inf_s(ishape, 5), m_id)
        case default
            stop "Unsupported shape"
        end select
    end do

    icount = 0
    do iz = 1, nz_m
        do iy = 1, ny_m
            do ix = 1, nx_m
                if (imedia_grid(ix, iy, iz) .ne. 0) then
                    icount = icount + 1
                    ibuff(1, icount) = ix
                    ibuff(2, icount) = iy
                    ibuff(3, icount) = iz
                    ibuff(4, icount) = imedia_grid(ix, iy, iz)
                end if
            end do
        end do
    end do

    npoint = icount
    write(fh, "(i9)") npoint
    do icount = 1, npoint
        write(fh, "(5i9)") icount, ibuff(1:4, icount)
    end do

contains


    subroutine fill_ellipsoid(x0, y0, z0, fx, fy, fz, a, b, mm)
        implicit none
        real(8), intent(in) :: x0, y0, z0, fx, fy, fz, a, b
        integer, intent(in) :: mm
        integer :: jx, jy, jz
        real(8) :: x, y, z, tmp
        logical :: flag
        do jz = 1, nz_m
            do jy = 1, ny_m
                do jx = 1, nx_m
                    x = hx_m * jx - x0
                    y = hy_m * jy - y0
                    z = hz_m * jz - z0
                    tmp = (x / (fx*0.5d0))**2 + (y / (fy*0.5d0))**2 + (z / (fz*0.5d0))**2
                    flag = (tmp <= 1.0)
                    if (flag) imedia_grid(jx, jy, jz) = mm
                end do
            end do
        end do
        return
    end subroutine fill_ellipsoid


    subroutine fill_half_ellipsoid(x0, y0, z0, fx, fy, fz, a, b, mm)
        implicit none
        real(8), intent(in) :: x0, y0, z0, fx, fy, fz, a, b
        integer, intent(in) :: mm
        integer :: jx, jy, jz
        real(8) :: x, y, z, tmp
        logical :: flag
        do jz = 1, nz_m
            do jy = 1, ny_m
                do jx = 1, nx_m
                    x = hx_m * jx - x0
                    y = hy_m * jy - y0
                    z = hz_m * jz - z0
                    tmp = (x / (fx*0.5d0))**2 + (y / (fy*0.5d0))**2 + (z / (fz))**2
                    flag = ((tmp <= 1.0) .and. (z >= 0))
                    if (flag) imedia_grid(jx, jy, jz) = mm
                end do
            end do
        end do
        return
    end subroutine fill_half_ellipsoid


    subroutine fill_elliptic_cylinder(x0, y0, z0, fx, fy, fz, a, b, mm)
        implicit none
        real(8), intent(in) :: x0, y0, z0, fx, fy, fz, a, b
        integer, intent(in) :: mm
        integer :: jx, jy, jz
        real(8) :: x, y, z, tmp
        logical :: flag
        do jz = 1, nz_m
            do jy = 1, ny_m
                do jx = 1, nx_m
                    x = hx_m * jx - x0
                    y = hy_m * jy - y0
                    z = hz_m * jz - z0
                    tmp = (x / (fx*0.5d0))**2 + (y / (fy*0.5d0))**2
                    flag = ((tmp <= 1.0) .and. (z >= -fz*0.5d0) .and. (z <= fz*0.5d0))
                    if (flag) imedia_grid(jx, jy, jz) = mm
                end do
            end do
        end do
        return
    end subroutine fill_elliptic_cylinder


    subroutine fill_triangular_cylinder(x0, y0, z0, fx, fy, fz, a, b, mm)
        implicit none
        real(8), intent(in) :: x0, y0, z0, fx, fy, fz, a, b
        integer, intent(in) :: mm
        integer :: jx, jy, jz
        real(8) :: x, y, z, tmp
        logical :: flag
        do jz = 1, nz_m
            do jy = 1, ny_m
                do jx = 1, nx_m
                    x = hx_m * jx - x0
                    y = hy_m * jy - y0
                    z = hz_m * jz - z0
                    flag = ((x >= -fx*0.5d0) .and. (x <= fx*0.5d0) .and. (y >= -fy / 3.0) &
                        .and. (y <= (fy / (fx*0.5d0)*x + fy*2.0 / 3.0)) &
                        .and. (y <= (fy / (fx*0.5d0)*x + fy*2.0 / 3.0)) &
                        .and. (z >= -fz*0.5d0) .and. (z <= fz*0.5d0))
                    if (flag) imedia_grid(jx, jy, jz) = mm
                end do
            end do
        end do
        return
    end subroutine fill_triangular_cylinder


    subroutine fill_rectangular_cylinder(x0, y0, z0, fx, fy, fz, a, b, mm)
        implicit none
        real(8), intent(in) :: x0, y0, z0, fx, fy, fz, a, b
        integer, intent(in) :: mm
        integer :: jx, jy, jz
        real(8) :: x, y, z, tmp
        logical :: flag
        do jz = 1, nz_m
            do jy = 1, ny_m
                do jx = 1, nx_m
                    x = hx_m * jx - x0
                    y = hy_m * jy - y0
                    z = hz_m * jz - z0
                    flag = ((x >= -fx*0.5d0) .and. (x <= fx*0.5d0) &
                        .and. (y >= -fy*0.5d0) .and. (y <= fy*0.5d0) &
                        .and. (z >= -fz*0.5d0) .and. (z <= fz*0.5d0))
                    if (flag) imedia_grid(jx, jy, jz) = mm
                end do
            end do
        end do
        return
    end subroutine fill_rectangular_cylinder


    subroutine fill_elliptic_cone(x0, y0, z0, fx, fy, fz, a, b, mm)
        implicit none
        real(8), intent(in) :: x0, y0, z0, fx, fy, fz, a, b
        integer, intent(in) :: mm
        integer :: jx, jy, jz
        real(8) :: x, y, z, tmp
        logical :: flag
        do jz = 1, nz_m
            do jy = 1, ny_m
                do jx = 1, nx_m
                    x = hx_m * jx - x0
                    y = hy_m * jy - y0
                    z = hz_m * jz - z0
                    tmp = (x / (fx*0.5d0*(fz-z) / fz))**2 + (y / (fy*0.5d0*(fz-z) / fz))**2
                    flag = ((tmp <= 1.0) .and. (z >= 0) .and. (z <= fz))
                    if (flag) imedia_grid(jx, jy, jz) = mm
                end do
            end do
        end do
        return
    end subroutine fill_elliptic_cone


    subroutine fill_triangular_cone(x0, y0, z0, fx, fy, fz, a, b, mm)
        implicit none
        real(8), intent(in) :: x0, y0, z0, fx, fy, fz, a, b
        integer, intent(in) :: mm
        integer :: jx, jy, jz
        real(8) :: x, y, z, tmp
        logical :: flag
        do jz = 1, nz_m
            do jy = 1, ny_m
                do jx = 1, nx_m
                    x = hx_m * jx - x0
                    y = hy_m * jy - y0
                    z = hz_m * jz - z0
                    flag = ((x >= -fx*0.5d0*(fz-z) / fz) &
                        .and. (x <= fx*0.5d0*(fz-z) / fz) &
                        .and. (y >= -fy / 3.0*(fz-z) / fz) &
                        .and. (y <= (fy / (fx*0.5d0)*x + fy*2.0 / 3.0*(fz-z) / fz)) &
                        .and. (y <= (-fy / (fx*0.5d0)*x + fy*2.0 / 3.0*(fz-z) / fz)) &
                        .and. (z >= 0) .and. (z <= fz))
                    if (flag) imedia_grid(jx, jy, jz) = mm
                end do
            end do
        end do
        return
    end subroutine fill_triangular_cone


    subroutine fill_rectangular_cone(x0, y0, z0, fx, fy, fz, a, b, mm)
        implicit none
        real(8), intent(in) :: x0, y0, z0, fx, fy, fz, a, b
        integer, intent(in) :: mm
        integer :: jx, jy, jz
        real(8) :: x, y, z, tmp
        logical :: flag
        do jz = 1, nz_m
            do jy = 1, ny_m
                do jx = 1, nx_m
                    x = hx_m * jx - x0
                    y = hy_m * jy - y0
                    z = hz_m * jz - z0
                    flag = ((x >= -fx*0.5d0*(fz-z) / fz) &
                        .and. (x <= fx*0.5d0*(fz-z) / fz) &
                        .and. (y >= -fy*0.5d0*(fz-z) / fz) &
                        .and. (y <= fy*0.5d0*(fz-z) / fz) &
                        .and. (z >= 0) .and. (z <= fz))
                    if (flag) imedia_grid(jx, jy, jz) = mm
                end do
            end do
        end do
        return
    end subroutine fill_rectangular_cone


    subroutine fill_elliptic_ring(x0, y0, z0, fx, fy, fz, a, b, mm)
        implicit none
        real(8), intent(in) :: x0, y0, z0, fx, fy, fz, a, b
        integer, intent(in) :: mm
        integer :: jx, jy, jz
        real(8) :: x, y, z, tmp
        logical :: flag
        do jz = 1, nz_m
            do jy = 1, ny_m
                do jx = 1, nx_m
                    x = hx_m * jx - x0
                    y = hy_m * jy - y0
                    z = hz_m * jz - z0
                    tmp = (x / (fx*0.5d0))**2 + (y / (fy*0.5d0))**2
                    if((tmp <= 1.0) .and. (z >= -fz*0.5d0) .and. (z <= fz*0.5d0)) then
                        tmp = (x / (a*0.5d0))**2.0 + (y / (b*0.5d0))**2.0
                        flag = (tmp >= 1.0)
                    else
                        flag = .false.
                    end if
                    if (flag) imedia_grid(jx, jy, jz) = mm
                end do
            end do
        end do
        return
    end subroutine fill_elliptic_ring

end subroutine generate_macropoint
end module shaper_ssbe




