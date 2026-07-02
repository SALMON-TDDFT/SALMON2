! Ground State Date Storage Module:

module gs_info_ssbe
    use math_constants, only: pi, zI
    implicit none

    type s_sbe_gs_info
        !Lattice information
        real(8) :: a_matrix(1:3, 1:3)
        real(8) :: b_matrix(1:3, 1:3)
        real(8) :: volume

        !Ground state (GS) electronic system information
        integer :: nk, nb, ne
        real(8), allocatable :: kpoint(:, :), kweight(:)
        real(8), allocatable :: eigen(:, :)
        real(8), allocatable :: occup(:, :)
        real(8), allocatable :: delta_omega(:, :, :)
        complex(8), allocatable :: p_mod_matrix(:, :, :, :)
        ! p_tm_matrix = <u|p|u>
        complex(8), allocatable :: p_tm_matrix(:, :, :, :)
        ! rvnl_tm_matrix = <u|-i[r, Vnl]|u>
        complex(8), allocatable :: rvnl_tm_matrix(:, :, :, :)
        complex(8), allocatable :: d_matrix(:, :, :, :)
        real(8), allocatable :: grad_k_eigen(:, :, :)

        ! GW collision-term inputs (Phase 2)
        real(8), allocatable :: gamma_gw(:, :)   ! (nb,nk) inelastic linewidth, a.u.
        real(8), allocatable :: f0_ref(:, :)     ! (nb,nk) cold occupation reference

        ! LG-SBE degeneracy inputs (prod_dk overlaps <u_{n,k}|u_{m,k+dk}>)
        integer :: nbvec                                ! number of dk-shift vectors = (2*ndk+1)**3
        integer, allocatable :: bvec(:, :)              ! (3,nbvec) integer dk shifts (jdk1,jdk2,jdk3)
        complex(8), allocatable :: prod_dk(:, :, :, :)  ! (nb,nb,nbvec,nk)

        ! LG-SBE Tier2 Pb3: non-Abelian Berry connection replacing the divergent
        ! interband dipole inside degenerate blocks (allocated only when 'gi').
        complex(8), allocatable :: xi(:, :, :, :)       ! (nb,nb,3,nk)  xi = s*i*logm(U)/dk
        logical,    allocatable :: xi_ok(:, :, :)       ! (nb,nb,nk)  .true. where xi is trustworthy

        !k-space grid and geometry information
        !NOTE: prepred for uniformally distributed k-grid....
        !integer :: num_kgrid(1:3)
    end type


contains


subroutine init_sbe_gs_info(gs, sysname, gs_directory, nk, nb, ne, a1, a2, a3, read_bin, icomm)
    use communication
    use filesystem, only: open_filehandle, get_filehandle
    use common_ssbe, only: grad_k_array_nb1d_double
    use salmon_global, only: gauge_sbe, file_sbe_prod_dk, sbe_lg_degen, num_kgrid, sbe_lg_degen_floor
    use degenerate_block_ssbe, only: build_xi, same_block, blend, theta_on, theta_off
    implicit none
    type(s_sbe_gs_info), intent(inout) :: gs
    character(*), intent(in) :: sysname
    character(*), intent(in) :: gs_directory
    integer, intent(in) :: nk
    integer, intent(in) :: nb
    integer, intent(in) :: ne
    real(8), intent(in) :: a1(1:3), a2(1:3), a3(1:3)
    logical, intent(in) :: read_bin
    integer, intent(in) :: icomm
    integer :: irank, nproc

    call comm_get_groupinfo(icomm, irank, nproc)

    gs%nk = nk
    gs%nb = nb
    gs%ne = ne
    !gs%num_kgrid(1:3) = num_kgrid(1:3)

    !Calculate b_matrix, volume_cell and volume_bz from a1..a3 vector.
    call calc_lattice_info()

    allocate(gs%kpoint(1:3, 1:nk))
    allocate(gs%kweight(1:nk))
    allocate(gs%eigen(1:nb, 1:nk))
    allocate(gs%occup(1:nb, 1:nk))
    allocate(gs%delta_omega(1:nb, 1:nb, 1:nk))
    allocate(gs%p_mod_matrix(1:nb, 1:nb, 1:3, 1:nk))
    allocate(gs%d_matrix(1:nb, 1:nb, 1:3, 1:nk))
    allocate(gs%p_tm_matrix(1:nb, 1:nb, 1:3, 1:nk))
    allocate(gs%rvnl_tm_matrix(1:nb, 1:nb, 1:3, 1:nk))
    allocate(gs%grad_k_eigen(1:nb, 1:3, 1:nk))

    if (irank == 0) then
        if (read_bin) then
            !Retrieve all data from binray
            write(*,*) "# read_sbe_gs_bin"
            call read_sbe_gs_bin()
        else
            !Retrieve eigenenergies from 'SYSNAME_eigen.data':
            write(*, '(a)') "# read_eigen_data"
            call read_eigen_data()
            !Retrieve k-points from 'SYSNAME_k.data':
            write(*, '(a)') "# read_k_data"
            call read_k_data()
            !Retrieve transition matrix from 'SYSNAME_tm.data':
            write(*, '(a)') "# read_tm_data"
            call read_tm_data()
            !Export all data from binray
            write(*, '(a)') "# save_sbe_gs_bin"
            call save_sbe_gs_bin()
        end if
    end if

    call comm_bcast(gs%kpoint, icomm, 0)
    call comm_bcast(gs%kweight, icomm, 0)
    call comm_bcast(gs%eigen, icomm, 0)
    call comm_bcast(gs%occup, icomm, 0)
    call comm_bcast(gs%p_tm_matrix, icomm, 0)
    call comm_bcast(gs%rvnl_tm_matrix, icomm, 0)

    !Retrieve k-space overlap products from 'file_sbe_prod_dk' (LG-SBE degeneracy):
    if (trim(sbe_lg_degen) == 'gi' .or. trim(sbe_lg_degen) == 'gifix') call read_prod_dk_data()

    !Calculate omega and d_matrix (neglecting diagonal part):
    if (irank == 0) write(*,"(a)") "# prepare_matrix"

    call prepare_matrix()
    call comm_bcast(gs%p_mod_matrix, icomm, 0)
    call comm_bcast(gs%delta_omega, icomm, 0)
    call comm_bcast(gs%d_matrix, icomm, 0) ! Experimental

    select case(trim(gauge_sbe))
    case ("length_gauge")
        call grad_k_array_nb1d_double(gs%nb, gs%nk, gs%b_matrix,  &
                                  &   gs%eigen, gs%grad_k_eigen)
    end select

    !Initial Occupation Number
    gs%occup(:,:) = 0d0 !!Experimental!!
    gs%occup(1:(ne/2),:) = 2d0 !!Experimental!!

contains

    ! Calculate lattice and reciprocal vectors
    subroutine calc_lattice_info()
        implicit none
        real(8) :: a12(1:3), a23(1:3), a31(1:3), volume
        real(8) :: b1(1:3), b2(1:3), b3(1:3)

        a12(1) = a1(2) * a2(3) - a1(3) * a2(2)
        a12(2) = a1(3) * a2(1) - a1(1) * a2(3)
        a12(3) = a1(1) * a2(2) - a1(2) * a2(1)
        a23(1) = a2(2) * a3(3) - a2(3) * a3(2)
        a23(2) = a2(3) * a3(1) - a2(1) * a3(3)
        a23(3) = a2(1) * a3(2) - a2(2) * a3(1)
        a31(1) = a3(2) * a1(3) - a3(3) * a1(2)
        a31(2) = a3(3) * a1(1) - a3(1) * a1(3)
        a31(3) = a3(1) * a1(2) - a3(2) * a1(1)
        volume = dot_product(a12, a3)
        b1(1:3) = (2d0 * pi / volume) * a23(1:3)
        b2(1:3) = (2d0 * pi / volume) * a31(1:3)
        b3(1:3) = (2d0 * pi / volume) * a12(1:3)

        gs%a_matrix(1:3, 1) = a1(1:3)
        gs%a_matrix(1:3, 2) = a2(1:3)
        gs%a_matrix(1:3, 3) = a3(1:3)
        gs%b_matrix(1, 1:3) = b1(1:3)
        gs%b_matrix(2, 1:3) = b2(1:3)
        gs%b_matrix(3, 1:3) = b3(1:3)
        gs%volume = volume
    end subroutine calc_lattice_info


    ! Read k-point coordinates from SALMON's output file
    subroutine read_k_data()
        implicit none
        character(256) :: dummy
        integer :: fh, ik, iik
        real(8) :: tmp(4)
        fh = open_filehandle(trim(gs_directory) // trim(sysname) // '_k.data', 'old')
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        do ik = 1, nk
            read(fh, *) iik, tmp(1:4)
            if (ik .ne. iik) stop "ik mismatch"
            gs%kpoint(1:3, ik) = tmp(1:3)
            gs%kweight(ik) = tmp(4)
        end do
        close(fh)
    end subroutine read_k_data


    ! Read eigenvalue data from SALMON's output file
    subroutine read_eigen_data()
        use inputoutput, only: au_energy_ev
        implicit none
        character(256) :: dummy
        integer :: fh, i, ik, iik, iib, ib
        real(8) :: tmp(2)
        real(8) :: efac

        fh = open_filehandle(trim(gs_directory) // trim(sysname) // '_eigen.data', 'old')
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        ! Third header line states the energy unit: "esp[a.u.]" or "esp[eV]".
        ! Convert eV files to a.u.; default to a.u. for tagless/unknown headers.
        if (index(dummy, '[eV]') > 0 .or. index(dummy, '[ev]') > 0) then
            efac = 1.0d0 / au_energy_ev
            write(*, "('#>',4x,a)") "read_eigen_data: detected eV units -> converting to a.u."
        else
            efac = 1.0d0
        end if
        do ik = 1, nk
            read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
            do ib = 1, nb
                read(fh, *) iib, tmp(1:2)
                if (ib .ne. iib) stop "ib mismatch"
                gs%eigen(ib, ik) = tmp(1) * efac
                ! gs%occup(ib, ik) = ctmp(2)
            end do
        end do
        close(fh)
    end subroutine read_eigen_data




    ! Read transition dipole moment from SALMON's output file
    subroutine read_tm_data()
        implicit none
        character(256) :: dummy
        integer :: fh, i, ik, ib, jb, iik, iib, jjb
        real(8) :: tmp(1:6)


        fh = open_filehandle(trim(gs_directory) // trim(sysname) // '_tm.data', 'old')
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        do ik = 1, nk
            do ib = 1, nb
                do jb = 1, nb
                    read(fh, *) iik, iib, jjb, tmp(1:6)
                    if (ik .ne. iik) stop "ik mismatch"
                    if (ib .ne. iib) stop "ib mismatch"
                    if (jb .ne. jjb) stop "jb mismatch"
                    gs%p_tm_matrix(ib, jb, 1, ik) = dcmplx(tmp(1), tmp(2))
                    gs%p_tm_matrix(ib, jb, 2, ik) = dcmplx(tmp(3), tmp(4))
                    gs%p_tm_matrix(ib, jb, 3, ik) = dcmplx(tmp(5), tmp(6))
                end do
            end do
        end do
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        do ik = 1, nk
            do ib = 1, nb
                do jb = 1, nb
                    read(fh, *) iik, iib, jjb, tmp(1:6)
                    if (ik .ne. iik) stop "ik mismatch"
                    if (ib .ne. iib) stop "ib mismatch"
                    if (jb .ne. jjb) stop "jb mismatch"
                    gs%rvnl_tm_matrix(ib, jb, 1, ik) = dcmplx(tmp(1), tmp(2))
                    gs%rvnl_tm_matrix(ib, jb, 2, ik) = dcmplx(tmp(3), tmp(4))
                    gs%rvnl_tm_matrix(ib, jb, 3, ik) = dcmplx(tmp(5), tmp(6))
                end do
            end do
        end do


        close(fh)
    end subroutine read_tm_data


    subroutine read_sbe_gs_bin()
        implicit none
        integer :: fh
        ! fh = get_filehandle()
        ! open(fh, file=trim(gs_directory) // trim(sysname) // '_sbe_gs.bin', form='unformatted', status='old')
        ! read(fh) gs%kpoint
        ! read(fh) gs%kweight
        ! read(fh) gs%eigen
        ! read(fh) gs%p_mod_matrix
        ! read(fh) gs%rvnl_tm_matrix
        ! ! read(fh) gs%prod_dk
        ! close(fh)
        ! return
    end subroutine read_sbe_gs_bin


    subroutine save_sbe_gs_bin()
        implicit none
        integer :: fh
        ! fh = get_filehandle()
        ! open(fh, file=trim(gs_directory) // trim(sysname) // '_sbe_gs.bin', form='unformatted', status='replace')
        ! write(fh) gs%kpoint
        ! write(fh) gs%kweight
        ! write(fh) gs%eigen
        ! write(fh) gs%p_mod_matrix
        ! write(fh) gs%rvnl_tm_matrix
        ! ! write(fh) gs%prod_dk
        ! close(fh)
        ! return
    end subroutine save_sbe_gs_bin


    ! Read the k-space overlap products <u_{io,ik}|u_{jo,ik+dk}> from the
    ! text file written by write_prod_dk_data (SALMON's '<sysname>_prod_dk.data').
    ! Root reads/validates, all ranks allocate, then the data are broadcast.
    subroutine read_prod_dk_data()
        implicit none
        integer :: fh, ios, ip, ierr
        logical :: file_exists
        character(256) :: cline
        integer :: file_no, file_nk, file_nk1, file_nk2, file_nk3, file_ndk
        integer :: ndk_loc, mdk, nrec, irec, ik_exp, ivec
        integer :: jdk1, jdk2, jdk3
        integer :: ik_r, ik1_r, ik2_r, ik3_r, jdk1_r, jdk2_r, jdk3_r, io_r, jo_r
        real(8) :: re_r, im_r

        ierr    = 0
        gs%nbvec = 0
        ndk_loc = 0
        file_no = 0

        ! --- root: open the file and parse the metadata header ---
        if (irank == 0) then
            write(*, '(a)') "# read_prod_dk_data"
            if (len_trim(file_sbe_prod_dk) == 0) then
                write(*, '(a)') "ERROR(read_prod_dk_data): 'file_sbe_prod_dk' is empty while 'sbe_lg_degen'='gi/gifix'."
                ierr = 1
            else
                inquire(file=trim(file_sbe_prod_dk), exist=file_exists)
                if (.not. file_exists) then
                    write(*, '(a)') "ERROR(read_prod_dk_data): file not found: " // trim(file_sbe_prod_dk)
                    ierr = 1
                end if
            end if

            if (ierr == 0) then
                fh = open_filehandle(trim(file_sbe_prod_dk), 'old')
                ! metadata line 1: "# no nk num_kgrid(1) num_kgrid(2) num_kgrid(3) ndk"
                read(fh, '(a)', iostat=ios) cline
                if (ios /= 0) then
                    write(*, '(a)') "ERROR(read_prod_dk_data): cannot read metadata header."
                    ierr = 1
                else
                    ip = index(cline, '#')
                    read(cline(ip+1:), *, iostat=ios) &
                        & file_no, file_nk, file_nk1, file_nk2, file_nk3, file_ndk
                    if (ios /= 0) then
                        write(*, '(a)') "ERROR(read_prod_dk_data): malformed metadata header."
                        ierr = 1
                    end if
                end if
            end if

            if (ierr == 0) then
                if (file_nk /= nk) then
                    write(*, '(a)') "ERROR(read_prod_dk_data): nk in file differs from SBE nk."
                    ierr = 1
                end if
                if (file_nk1 * file_nk2 * file_nk3 /= file_nk) then
                    write(*, '(a)') "ERROR(read_prod_dk_data): num_kgrid product differs from nk in file."
                    ierr = 1
                end if
                if (file_no < nb) then
                    write(*, '(a)') "ERROR(read_prod_dk_data): band window in file is smaller than SBE nb."
                    ierr = 1
                end if
                if (file_ndk < 0) then
                    write(*, '(a)') "ERROR(read_prod_dk_data): ndk in file is negative."
                    ierr = 1
                end if
            end if

            if (ierr == 0) then
                ndk_loc  = file_ndk
                gs%nbvec = (2 * ndk_loc + 1) ** 3
                ! metadata line 2: column legend (skip)
                read(fh, '(a)', iostat=ios) cline
            end if
        end if

        ! --- propagate metadata-stage error, then broadcast nbvec ---
        call comm_bcast(ierr, icomm, 0)
        if (ierr /= 0) stop 1
        call comm_bcast(gs%nbvec, icomm, 0)

        ! --- allocate on ALL ranks ---
        allocate(gs%bvec(1:3, 1:gs%nbvec))
        allocate(gs%prod_dk(1:nb, 1:nb, 1:gs%nbvec, 1:nk))
        gs%bvec(:, :)          = 0
        gs%prod_dk(:, :, :, :) = dcmplx(0d0, 0d0)

        ! --- root: build the dk-shift table and read all 11-column records ---
        if (irank == 0) then
            mdk  = 2 * ndk_loc + 1
            ivec = 0
            do jdk3 = -ndk_loc, ndk_loc
                do jdk2 = -ndk_loc, ndk_loc
                    do jdk1 = -ndk_loc, ndk_loc
                        ivec = ivec + 1
                        gs%bvec(1, ivec) = jdk1
                        gs%bvec(2, ivec) = jdk2
                        gs%bvec(3, ivec) = jdk3
                    end do
                end do
            end do

            ! writer emits nk * (2*ndk+1)**3 * no**2 records, ik outermost/slowest.
            nrec = nk * gs%nbvec * file_no * file_no
            do irec = 1, nrec
                read(fh, *, iostat=ios) &
                    & ik_r, ik1_r, ik2_r, ik3_r, jdk1_r, jdk2_r, jdk3_r, io_r, jo_r, re_r, im_r
                if (ios /= 0) then
                    write(*, '(a)') "ERROR(read_prod_dk_data): fewer records than expected."
                    ierr = 1
                    exit
                end if
                ! ik must run 1..nk in contiguous blocks (matches SBE k ordering)
                ik_exp = (irec - 1) / (gs%nbvec * file_no * file_no) + 1
                if (ik_r /= ik_exp) then
                    write(*, '(a)') "ERROR(read_prod_dk_data): ik ordering does not match SBE k-grid."
                    ierr = 1
                    exit
                end if
                if (abs(jdk1_r) > ndk_loc .or. abs(jdk2_r) > ndk_loc .or. abs(jdk3_r) > ndk_loc) then
                    write(*, '(a)') "ERROR(read_prod_dk_data): dk-shift index out of range."
                    ierr = 1
                    exit
                end if
                ! keep only the SBE band window (writer emits the full band window)
                if (io_r <= nb .and. jo_r <= nb) then
                    ivec = (jdk3_r + ndk_loc) * mdk * mdk &
                        & + (jdk2_r + ndk_loc) * mdk &
                        & + (jdk1_r + ndk_loc) + 1
                    gs%prod_dk(io_r, jo_r, ivec, ik_r) = dcmplx(re_r, im_r)
                end if
            end do

            ! record-count check: no surplus data lines beyond the expected block
            if (ierr == 0) then
                read(fh, '(a)', iostat=ios) cline
                if (ios == 0) then
                    if (len_trim(cline) > 0) then
                        write(*, '(a)') "ERROR(read_prod_dk_data): more records than expected."
                        ierr = 1
                    end if
                end if
            end if
            close(fh)
        end if

        ! --- propagate record-stage error, then broadcast the data ---
        call comm_bcast(ierr, icomm, 0)
        if (ierr /= 0) stop 1
        call comm_bcast(gs%bvec, icomm, 0)
        call comm_bcast(gs%prod_dk, icomm, 0)
    end subroutine read_prod_dk_data


    subroutine prepare_matrix()
        implicit none
        integer :: ik, ib, jb
        real(8) :: omega_eps
        real(8) :: x, w, resu, resp, resp_max
        integer :: nrej
        complex(8) :: dpdw(1:3)

        omega_eps = sbe_lg_degen_floor

        gs%p_mod_matrix = gs%p_tm_matrix + gs%rvnl_tm_matrix

        if (trim(sbe_lg_degen) == 'gi' .or. trim(sbe_lg_degen) == 'gifix') then
            ! ===== Pb3: non-Abelian xi inside degenerate blocks + smooth blend =====
            ! delta_omega first (needed by the blend), then build xi from prod_dk.
            do ik=1, nk
                do ib=1, nb
                    do jb=1, nb
                        gs%delta_omega(ib, jb, ik) = gs%eigen(ib, ik) - gs%eigen(jb, ik)
                    end do
                end do
            end do

            if (.not. allocated(gs%xi))    allocate(gs%xi(1:nb, 1:nb, 1:3, 1:nk))
            if (.not. allocated(gs%xi_ok)) allocate(gs%xi_ok(1:nb, 1:nb, 1:nk))
            call build_xi(nb, nk, gs%nbvec, gs%bvec, gs%prod_dk, gs%eigen, &
                        & gs%b_matrix, num_kgrid, gs%xi, gs%xi_ok, nrej, resu, &
                        & fixed_blocks=(trim(sbe_lg_degen) == 'gifix'))

            resp_max = 0d0
            do ik=1, nk
                do ib=1, nb
                    do jb=1, nb
                        x = abs(gs%delta_omega(ib, jb, ik))
                        if (same_block(ib, jb, ik) .and. ib /= jb .and. gs%xi_ok(ib, jb, ik)) then
                            ! degenerate-block pair: blend xi (x<=theta_on) into the
                            ! ordinary dipole i*p/delta_omega (x>=theta_off).
                            w = blend(x, theta_on, theta_off)
                            if (x > omega_eps) then
                                dpdw(1:3) = zi * gs%p_mod_matrix(ib, jb, 1:3, ik) &
                                          & / gs%delta_omega(ib, jb, ik)
                            else
                                dpdw(1:3) = (0d0, 0d0)   ! exact degeneracy: dipole undefined
                            end if
                            gs%d_matrix(ib, jb, 1:3, ik) = &
                                & w * gs%xi(ib, jb, 1:3, ik) + (1d0 - w) * dpdw(1:3)
                            ! diagnostic: xi vs i*p/dw where both are valid
                            if (x >= theta_on) then
                                resp = maxval(abs(gs%xi(ib, jb, 1:3, ik) - dpdw(1:3)))
                                if (resp > resp_max) resp_max = resp
                            end if
                        else if (omega_eps < x) then
                            gs%d_matrix(ib, jb, 1:3, ik) = &
                                & zi * (gs%p_mod_matrix(ib, jb, 1:3, ik)) &
                                & / gs%delta_omega(ib, jb, ik)
                        else
                            gs%d_matrix(ib, jb, 1:3, ik) = 0d0
                        end if
                    end do
                end do
            end do

            if (irank == 0) then
                write(*, '(a,i0)')     "# build_xi: rejected block-links          = ", nrej
                write(*, '(a,es12.4)') "# build_xi: max |U^H U - I| (polar health) = ", resu
                write(*, '(a,es12.4)') "# build_xi: max |xi - i p/dw| (both-valid) = ", resp_max
            end if
        else
            ! ===== default 'off': bit-identical to the pre-Pb3 dipole construction =====
            do ik=1, nk
                do ib=1, nb
                    do jb=1, nb
                        gs%delta_omega(ib, jb, ik) = gs%eigen(ib, ik) - gs%eigen(jb, ik)
                        if (omega_eps < abs(gs%delta_omega(ib, jb, ik))) then
                            ! gs%d_matrix(ib, jb, 1:3, ik) = &
                            !     & (zi * gs%p_mod_matrix(ib, jb, 1:3, ik) - gs%rvnl_tm_matrix(ib, jb, 1:3, ik)) &
                            !     & / gs%delta_omega(ib, jb, ik)
                            ! gs%d_matrix(ib, jb, 1:3, ik) = &
                            !     & zi * (gs%p_mod_matrix(ib, jb, 1:3, ik) +  gs%rvnl_tm_matrix(ib, jb, 1:3, ik)) &
                            !     & / gs%delta_omega(ib, jb, ik)
                            gs%d_matrix(ib, jb, 1:3, ik) = &
                                & zi * (gs%p_mod_matrix(ib, jb, 1:3, ik)) &
                                & / gs%delta_omega(ib, jb, ik)
                        else
                            gs%d_matrix(ib, jb, 1:3, ik) = 0d0
                        end if
                    end do
                end do
            end do
        end if
    end subroutine prepare_matrix


end subroutine init_sbe_gs_info

end module gs_info_ssbe

