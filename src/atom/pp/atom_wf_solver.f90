module solver
    use salmon_math, only: pi
    implicit none

contains

subroutine normalize(nr, dr, u)
    implicit none
    integer, intent(in) :: nr
    real(8), intent(in) :: dr
    real(8), intent(inout) :: u(nr)
    real(8) :: tmp_s, tmp_c
    tmp_s = 4.0 * pi * sum(u(:) ** 2) * dr
    tmp_c = 1.0 / sqrt(tmp_s)
    if (sum(u(1:nr)) < 0) then
        tmp_c = -tmp_c
    end if
    u = tmp_c * u
end subroutine normalize

subroutine calc_eigen_real_sym(n, m, hmat, vec, eig)
    implicit none
    integer, intent(in)  :: n
    integer, intent(in)  :: m
    real(8), intent(in)  :: hmat(n, n)
    real(8), intent(out) :: vec(n, m)
    real(8), intent(out) :: eig(m)
    real(8) :: tmp(n, n)
    real(8) :: work(8*n)
    integer :: iwork(5*n), ifail(n)
    integer :: info, m_out
    real(8), parameter :: vl = 0.0d0, vu = 0.0d0, abstol = -1.0d0
    tmp(:, :) = hmat(:, :)
    call DSYEVX('V', 'I', 'U', &
                n, tmp, n, &
                vl, vu, 1, m, abstol, &
                m_out, eig, vec, n, &
                work, size(work), iwork, ifail, info)
    if ((info .ne. 0) .or. (m_out .ne. m)) then
        stop "Error: DSYEVX in calc_eigen_real_sym failed!"
    end if
    return
end subroutine calc_eigen_real_sym


subroutine calc_eigen_real_sym_tri(n, m, alpha, beta, vec, eig)
    implicit none
    integer, intent(in)  :: n
    integer, intent(in)  :: m
    real(8), intent(in)  :: alpha(n), beta(n-1)
    real(8), intent(out) :: vec(n, m)
    real(8), intent(out) :: eig(m)
    real(8) :: dtmp(n), etmp(n-1)
    real(8) :: w(n)          
    real(8) :: z(n,n)        
    real(8) :: work(5*n)     
    integer :: iwork(5*n)    
    integer :: ifail(n)
    integer :: info, m_out
    real(8), parameter :: vl = 0.0d0, vu = 0.0d0, abstol = -1.0d0
    dtmp(:) = alpha(:)
    etmp(:) = beta(:)
    call DSTEVX('V', 'I', n, dtmp, etmp, &
                vl, vu, 1, m, abstol, &
                m_out, eig, vec, n, &
                work, iwork, ifail, info)
    if ((info .ne. 0) .or. (m_out .ne. m)) then
        stop "Error: DSYEVX in calc_eigen_real_sym_tri failed!"
    end if
    return
end subroutine calc_eigen_real_sym_tri



subroutine calc_rho(nr, dr, nb, l_max, u_wf, occ, rho)
    implicit none
    real(8), intent(in) :: dr
    integer, intent(in) :: nr, nb, l_max
    real(8), intent(in) :: u_wf(1:nr, 1:nb, 0:l_max)
    real(8), intent(in) :: occ(1:nb, 0:l_max)
    real(8), intent(out) :: rho(1:nr)
    real(8) :: r
    integer :: l, ib, ir
    rho(:) = 0.0d0
    do l = 0, l_max
        do ib = 1, nb
            if (occ(ib, l) > 0) then
                do ir = 1, nr
                    r = ir * dr
                    rho(ir) = rho(ir) + occ(ib, l) * abs(u_wf(ir, ib, l) / r) ** 2
                end do
            end if
        end do
    end do
    return
end subroutine calc_rho


subroutine calc_v_h(nr, dr, rho, v_h)
    implicit none
    real(8), intent(in) :: dr
    integer, intent(in) :: nr
    real(8), intent(in) :: rho(1:nr)
    real(8), intent(out) :: v_h(1:nr)
    real(8) :: tmp, r
    integer :: i
    real(8), parameter :: pi = 4.0d0 * atan(1.0d0)
    v_h(:) = 0.0d0
    tmp = 0.0d0
    do i = 1, nr
        r = i * dr
        tmp = tmp + r ** 2 * rho(i) * dr
        v_h(i) = 4 * pi / r * tmp
    end do
    tmp = 0.0d0
    do i = nr, 1, -1
        r = i * dr
        tmp = tmp + r * rho(i) * dr
        v_h(i) = v_h(i) + 4 * pi * tmp
    end do
    return
end subroutine calc_v_h

subroutine calc_xc(x_func, c_func, nr, rho, exc, vxc)
    implicit none
    TYPE(xc_f03_func_t), intent(in) :: x_func
    TYPE(xc_f03_func_t), intent(in) :: c_func
    integer, intent(in) :: nr
    real(8), intent(in) :: rho(nr)
    real(8), intent(out) :: exc(nr), vxc(nr)
    real(8) :: ex(nr), vx(nr), ec(nr), vc(nr)
    TYPE(xc_f03_func_info_t) :: xc_info
    integer(8) :: nlen
    nlen = nr
    call xc_f03_lda_exc_vxc(x_func, nlen, rho(1), ex(1), vx(1))
    call xc_f03_lda_exc_vxc(c_func, nlen, rho(1), ec(1), vc(1))
    exc = ex + ec
    vxc = vx + vc
end subroutine calc_xc

subroutine solve_schrodinger(nr, dr, nb, l_max, v_eff, np, np_max, p_coeff, p_func, u_wf, eig)
    implicit none
    real(8), intent(in) :: dr
    integer, intent(in) :: nr, nb, l_max
    real(8), intent(in) :: v_eff(1:nr)
    integer, intent(in) :: np_max
    integer, intent(in) :: np(0:l_max)
    real(8), intent(in) :: p_coeff(1:np_max, 0:l_max)
    real(8), intent(in) :: p_func(1:nr, 1:np_max, 0:l_max)
    real(8), intent(inout) :: u_wf(1:nr, 1:nb, 0:l_max)
    real(8), intent(out) :: eig(1:nb, 0:l_max)
    integer :: l, ir, ib, ip, ir1, ir2
    real(8) :: r
    real(8), parameter :: pi = 4.0d0 * atan(1.0d0)
    
    real(8) :: h_mat(nr, nr)
    ! real(8) :: alpha(nr), beta(nr-1)
    
    ! beta(:) = -0.5d0 / dr**2
    do l = 0, l_max
        h_mat(:, :) = 0.0d0
        do ir = 1, nr
            r = ir * dr
            h_mat(ir, ir) = 1.0d0 / dr**2 + 0.5 * dble(l*(l+1)) / r**2 + v_eff(ir)
            ! alpha(ir) = 1.0d0 / dr**2 + 0.5 * dble(l*(l+1)) / r**2 + v_eff(ir)
        end do
        do ir = 1, nr-1
            h_mat(ir, ir+1) = -0.5d0 / dr**2
            h_mat(ir+1, ir) = -0.5d0 / dr**2
        end do
        do ip = 1, np(l)
            do ir1 = 1, nr
                do ir2 = 1, nr
                    h_mat(ir1, ir2) = h_mat(ir1, ir2) + p_coeff(ip, l) * p_func(ir1, ip, l) * p_func(ir2, ip, l) * dr
                end do
            end do
        end do
        call calc_eigen_real_sym(nr, nb, h_mat, u_wf(:, :, l), eig(:, l))
        ! call calc_eigen_real_sym_tri(nr, nb, alpha, beta, u_wf(:, :, l), eig(:, l))
        do ib = 1, nb
            call normalize(nr, dr, u_wf(:, ib, l))
        end do
    end do
end subroutine solve_schrodinger

subroutine calc_total_energy(nr, dr, rho, v_h, v_xc, e_xc, nb, l_max, eig, occ, e_tot)
    implicit none
    integer, intent(in) :: nr
    real(8), intent(in) :: dr
    real(8), intent(in) :: rho(nr), v_h(nr), v_xc(nr), e_xc(nr)
    integer, intent(in) :: nb, l_max
    real(8), intent(in) :: eig(nb, 0:l_max), occ(nb, 0:l_max)
    real(8), intent(out) :: e_tot
    real(8) :: r
    integer :: ir
    e_tot = sum(eig * occ)
    do ir = 1, nr
        r = ir * dr
        e_tot = e_tot + 4.0 * pi * r**2 * ( &
            & -0.5 * v_h(ir) - v_xc(ir) + e_xc(ir) &
            & ) * rho(ir) * dr
    end do
end subroutine calc_total_energy

subroutine calc_occ(n_elec, nb, l_max, eig, occ)
    implicit none
    integer, intent(in) :: n_elec
    integer, intent(in) :: nb, l_max
    real(8), intent(in) :: eig(nb, 0:l_max)
    real(8), intent(out) :: occ(nb, 0:l_max)
    real(8) :: e_tmp
    integer :: i, l, ib, l_tmp, ib_tmp
    occ(:, :) = 0.0
    do i = 1, n_elec
        e_tmp = maxval(eig)+1.0
        do l = 0, l_max
            do ib = 1, nb
                if (eig(ib, l) < e_tmp) then
                    if  (occ(ib, l) < 2*(2*l+1)) then
                        ib_tmp = ib
                        l_tmp = l
                        e_tmp = eig(ib, l)
                    end if
                end if
            end do
        end do
        occ(ib_tmp, l_tmp) = occ(ib_tmp, l_tmp) + 1
    end do
end subroutine calc_occ



subroutine calc_atom_wf(iz_num, n_elec, r_max, nr, nb, l_max, &
                         n_scf, a_mix, read_rho, update_occ, occ0, &
                         x_func_index, c_func_index, pp_file)
    implicit none
    integer :: iz_num 
    integer :: n_elec 
    real(8) :: r_max 
    integer :: nr 
    integer :: nb 
    integer :: l_max 
    integer :: n_scf 
    real(8) :: a_mix 
    logical :: read_rho
    logical :: update_occ
    real(8) :: occ0(1:9, 0:9)
    integer :: x_func_index, c_func_index
    real(8) :: dr
    real(8), allocatable :: u_wf(:, :, :), eig(:, :), v_eff(:), v_h(:)
    real(8), allocatable :: v_xc(:), e_xc(:)
    real(8), allocatable :: rho(:), rho_new(:),  occ(:, :)
    real(8) ::  r, e_tot, delta_rho
    integer :: ir, ib, l, i_scf, i_dummy, ip
    TYPE(xc_f03_func_t) :: x_func
    TYPE(xc_f03_func_t) :: c_func
    character(256) :: pp_file
    type(t_pp_data) :: pp_data
    real(8), allocatable :: p_func(:, :, :), v_loc(:)


    ! iz_num = -1
    n_elec = -1
    r_max = 10.0
    nr = 1000
    nb = 5
    l_max = -1
    n_scf = 50
    a_mix = 0.5
    read_rho = .false.
    update_occ = .true.
    occ0(:, :) = 0
    x_func_index = XC_LDA_X
    c_func_index = XC_LDA_C_PZ
    pp_file = ''

    ! open(unit=10, file="atom4_input.nml", status="old", action="read")
    ! close(10)


    allocate(u_wf(1:nr, 1:nb, 0:l_max), eig(1:nb, 0:l_max), v_eff(1:nr), v_h(1:nr))
    allocate(v_xc(1:nr), e_xc(1:nr))
    allocate(rho(1:nr),rho_new(1:nr),  occ(1:nb, 0:l_max))
    allocate(p_func(1:nr, 1:pp_data%np_max, 0:l_max), v_loc(1:nr))

    dr = r_max / nr
    occ(1:nb, 0:l_max) = occ0(1:nb, 0:l_max)
    if ((.not. update_occ) .and. (abs(sum(occ) - n_elec) > 1d-9)) &
        stop "Error: sum(occ) is not n_elec!"

    call xc_f03_func_init(x_func, x_func_index, XC_UNPOLARIZED)
    call xc_f03_func_init(c_func, c_func_index, XC_UNPOLARIZED)

    do l = 0, l_max
        do ip = 1, pp_data%np(l)
            do ir = 1, nr
                p_func(ir, ip, l) = interp_p_func(pp_data, ir*dr, ip, l)
            end do
        end do
    end do
    do ir = 1, nr
        v_loc(ir) = interp_v_loc(pp_data, ir*dr)
    end do



    if (read_rho) then
        open(100, file="rho.txt", action="read")
        do ir = 1, nr
            read(100, *) r, rho(ir)
        end do
        close(100)
    else
        rho(:) = dble(n_elec) / (4.0d0 / 3.0d0 * pi * r_max ** 3)
    end if

    do i_scf = 1, n_scf
        call calc_v_h(nr, dr, rho, v_h)
        call calc_xc(x_func, c_func, nr, rho, e_xc, v_xc)
        do ir = 1, nr
            r = dr * ir
            ! v_eff(ir) = -dble(iz_num) / r + v_h(ir) + v_xc(ir)
            v_eff(ir) = v_loc(ir) + v_h(ir) + v_xc(ir)
        end do
        call solve_schrodinger(nr, dr, nb, l_max, v_eff, &
            & pp_data%np, pp_data%np_max, pp_data%p_coeff, p_func, &
            & u_wf, eig)
        if (update_occ) call calc_occ(n_elec, nb, l_max, eig, occ)
        call calc_rho(nr, dr, nb, l_max, u_wf, occ, rho_new)
        call calc_total_energy(nr, dr, rho, v_h, v_xc, e_xc, nb, l_max, eig, occ, e_tot)
        delta_rho = maxval(abs(rho_new - rho))
        write(*, "(a,i3,a,es7.1e2,a,f14.6)") "# i_scf=", i_scf, ", delta_rho=", delta_rho , ", total_energy=  ", e_tot
        rho = a_mix * rho_new + (1.0d0 - a_mix) * rho
    end do

    ! write(*, "(a,a3,a3,2a14)") "#",  "ib", "l", "Eigenenergy", "Occupation"
    ! do l=0, l_max
    !     do ib = 1, nb
    !         write(*, "(a,i3,i3,2f14.6)") "#",  ib, l, eig(ib, l), occ(ib, l)
    !     end do
    ! end do

    ! open(100, file="rho.txt", action="write")
    ! do ir = 1, nr
    !     write(100, "(f12.6,es24.15e3)") ir*dr, rho(ir)
    ! end do
    ! close(100)
    ! open(100, file="v_h.txt", action="write")
    ! do ir = 1, nr
    !     write(100, "(f12.6,es24.15e3)") ir*dr, v_h(ir)
    ! end do
    ! close(100)
    ! open(100, file="v_xc.txt", action="write")
    ! do ir = 1, nr
    !     write(100, "(f12.6,es24.15e3)") ir*dr, v_xc(ir)
    ! end do
    ! open(100, file="v_loc.txt", action="write")
    ! do ir = 1, nr
    !     write(100, "(f12.6,es24.15e3)") ir*dr, v_loc(ir) 
    ! end do
    ! close(100)

end subroutine calc_atom_wf

end module

