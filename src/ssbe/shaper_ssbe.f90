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
        end select

        select case(trim(typ_s(ishape)))
        case('ellipsoid')            
            call fill_ellipsoid(ori_s(ishape, 1), ori_s(ishape, 2), ori_s(ishape, 3), &
                & inf_s(ishape, 1)*0.5, inf_s(ishape, 2)*0.5, inf_s(ishape, 3)*0.5, &
                &  .false., .false., .false., m_id)
        case('half-ellipsoid')
            call fill_ellipsoid(ori_s(ishape, 1), ori_s(ishape, 2), ori_s(ishape, 3), &
                & inf_s(ishape, 1)*0.5, inf_s(ishape, 2)*0.5, inf_s(ishape, 3), &
                &  .false., .false., .true., m_id)
        case('elliptic-cylinder')
            call fill_elliptic_cylinder(ori_s(ishape, 1), ori_s(ishape, 2), ori_s(ishape, 3), &
                & inf_s(ishape, 1)*0.5, inf_s(ishape, 2)*0.5, inf_s(ishape, 3)*0.5, m_id)
        case('rectangular-cylinder') 
            call fill_rectangular_cylinder(ori_s(ishape, 1), ori_s(ishape, 2), ori_s(ishape, 3), &
                & inf_s(ishape, 1)*0.5, inf_s(ishape, 2)*0.5, inf_s(ishape, 3)*0.5, m_id)
        ! case('triangular-cylinder')  
        ! case('elliptic-cone')        
        ! case('triangular-cone')      
        ! case('rectangular-cone')     
        ! case('elliptic-ring')
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
        write(fh, "(4i9)") ibuff(1:4, icount)
    end do

contains

    subroutine fill_ellipsoid(x0, y0, z0, rx, ry, rz, xh, yh, zh, mm)
        implicit none
        real(8), intent(in) :: x0, y0, z0, rx, ry, rz
        logical, intent(in) :: xh, yh, zh
        integer, intent(in) :: mm
        integer :: jx, jx_min, jx_max
        integer :: jy, jy_min, jy_max
        integer :: jz, jz_min, jz_max
        real(8) :: x, y, z, t
        jx_min = max(1, floor((x0 - rx) / hx_m))
        jy_min = max(1, floor((y0 - ry) / hy_m))
        jz_min = max(1, floor((z0 - rz) / hz_m))
        jx_max = min(nx_m, ceiling((x0 + rx) / hx_m))
        jy_max = min(ny_m, ceiling((y0 + ry) / hy_m))
        jz_max = min(nz_m, ceiling((z0 + rz) / hz_m))
        do jz = jz_min, jz_max
            do jy = jy_min, jy_max
                do jx = jx_min, jx_max
                    x = hx_m * jx
                    y = hy_m * jy
                    z = hz_m * jz
                    if ((zh == .false.) .or. (z > z0)) then
                        t = ((x - x0) / rx) ** 2 &
                            & + ((y - y0) / ry) ** 2 &
                            & + ((z - z0) / rz) ** 2
                        if (t <= 1.0d0) imedia_grid(jx, jy, jz) = mm
                    end if
                end do
            end do
        end do
        return
    end subroutine fill_ellipsoid

    subroutine fill_elliptic_cylinder(x0, y0, z0, rx, ry, rz, mm)
        implicit none
        real(8), intent(in) :: x0, y0, z0, rx, ry, rz
        integer, intent(in) :: mm
        integer :: jx, jx_min, jx_max
        integer :: jy, jy_min, jy_max
        integer :: jz, jz_min, jz_max
        real(8) :: x, y, z, t
        jx_min = max(1, floor((x0 - rx) / hx_m))
        jy_min = max(1, floor((y0 - ry) / hy_m))
        jz_min = max(1, floor((z0 - rz) / hz_m))
        jx_max = min(nx_m, ceiling((x0 + rx) / hx_m))
        jy_max = min(ny_m, ceiling((y0 + ry) / hy_m))
        jz_max = min(nz_m, ceiling((z0 + rz) / hz_m))
        do jz = jz_min, jz_max
            do jy = jy_min, jy_max
                do jx = jx_min, jx_max
                    x = hx_m * jx
                    y = hy_m * jy
                    z = hz_m * jz
                    if ((z < (z0 - rz)) .or. (z > (z0 + rz))) then
                        t = ((x - x0) / rx) ** 2 + ((y - y0) / ry) ** 2 
                        if (t <= 1.0d0) imedia_grid(jx, jy, jz) = mm
                    end if
                end do
            end do
        end do
        return
    end subroutine fill_elliptic_cylinder
                    
    subroutine fill_rectangular_cylinder(x0, y0, z0, rx, ry, rz, mm)
        implicit none
        real(8), intent(in) :: x0, y0, z0, rx, ry, rz
        integer, intent(in) :: mm
        integer :: jx, jx_min, jx_max
        integer :: jy, jy_min, jy_max
        integer :: jz, jz_min, jz_max
        real(8) :: x, y, z, t
        jx_min = max(1, floor((x0 - rx) / hx_m))
        jy_min = max(1, floor((y0 - ry) / hy_m))
        jz_min = max(1, floor((z0 - rz) / hz_m))
        jx_max = min(nx_m, ceiling((x0 + rx) / hx_m))
        jy_max = min(ny_m, ceiling((y0 + ry) / hy_m))
        jz_max = min(nz_m, ceiling((z0 + rz) / hz_m))
        do jz = jz_min, jz_max
            do jy = jy_min, jy_max
                do jx = jx_min, jx_max
                    x = hx_m * jx
                    y = hy_m * jy
                    z = hz_m * jz
                    if ((abs(x-x0)<=rx) .and. (abs(y-y0)<=ry) .and. (abs(z-z0)<=rz)) then
                        imedia_grid(jx, jy, jz) = mm
                    end if
                end do
            end do
        end do
        return
    end subroutine fill_rectangular_cylinder

end subroutine generate_macropoint
end module shaper_ssbe




