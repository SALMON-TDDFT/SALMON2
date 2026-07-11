! Ground State Date Storage Module:

module gs_info_ssbe
    use math_constants, only: pi, zI
    use sbe_lg_mode_ssbe, only: uses_prod_dk, uses_xfull_links
    implicit none

    type s_sbe_gs_info
        !Lattice information
        real(8) :: a_matrix(1:3, 1:3)
        real(8) :: b_matrix(1:3, 1:3)
        real(8) :: volume

        !Ground state (GS) electronic system information
        integer :: nk, nb, ne
        ! SO-SBE spinor convention (spin='noncollinear', always with
        ! yn_spinorbit='y'): the export bands are SPINOR states -- 1 electron
        ! per occupied band instead of 2.  focc/nvb centralize the occupation
        ! convention so every consumer (occup init, trace bookkeeping, the
        ! norder_correction occupied-band loop) shares one source of truth.
        logical :: spinor = .false.  ! .true. = noncollinear/SOC spinor bands
        real(8) :: focc = 2d0        ! occupation per occupied band (2 / 1)
        integer :: nvb = 0           ! occupied bands inside the window (= ne/focc)
        ! metal detector (check_gicov_occupation, degenerate_block_ssbe): set
        ! once by init_sbe_gs_info alongside gs%occup, so every downstream
        ! consumer (occupation fill, fsum-rule deficiency tensor, the
        ! norder_correction occupied-band loop) can cheaply consult it without
        ! recomputing.  .false. (insulator/semiconductor) is the safe default.
        logical :: is_metal = .false.
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

        ! LG-SBE Phase 1 (gicov): block-diagonal Wilson-line transport + fixed partition
        complex(8), allocatable :: u_transport(:, :, :, :)  ! (nb,nb,3,nk) block-diagonal Wilson-line transport, gicov
        integer,    allocatable :: block_id(:, :)           ! (nb,nk) fixed-block partition, per-gs (NOT the module cache)

        ! VG completion (yn_sbe_vnl_exact='y'): all-order nonlocal V_nl(k+A) on a
        ! 1D kappa-stencil kappa_s = k + s*vnl_h*vnl_dir, read from
        ! file_sbe_vnl_kappa (written by write_sbe_vnl_kappa_data).  Stored as a
        ! per-rank k-SLICE [vnl_ik_min:vnl_ik_max] (split_range partition,
        ! asserted equal to the solver's slice at init_sbe_bloch_solver).
        ! On read, gs%rvnl_tm_matrix is OVERWRITTEN (full nk, all ranks) by the
        ! s=0 velocity matrices W_{i,0}, so fsum_D / p_mod / diagnostics and the
        ! exact readout share one binary-precision velocity convention.
        logical :: vnl_exact = .false.
        integer :: vnl_ns = 0
        real(8) :: vnl_h = 0d0
        real(8) :: vnl_dir(1:3) = 0d0
        integer :: vnl_ik_min = 0, vnl_ik_max = -1
        complex(8), allocatable :: vnl_V(:, :, :, :)    ! (nb,nb,-ns:ns, ik-slice)
        complex(8), allocatable :: vnl_W(:, :, :, :, :) ! (nb,nb,3,-ns:ns, ik-slice)

        !k-space grid and geometry information
        !NOTE: prepred for uniformally distributed k-grid....
        !integer :: num_kgrid(1:3)
    end type


contains


! Band-window selection (READER-side windowing).  Three distinct band counts:
!   nb     = the EXPORT band count (rows per k in SYSNAME_eigen/tm.data, as
!            written by the GS run) -- the read loops always consume ALL nb
!            records to keep the file alignment;
!   [nb_min, nb_hi] = the contiguous window actually stored into gs%*
!            (nb_hi = nstate_sbe upper cut, nb_min = nband_sbe_min lower cut);
!   nb_eff = nb_hi - nb_min + 1 = gs%nb (the propagated band count).
! Window index w <-> absolute band w + nb_min - 1.
! Bands 1..nb_min-1 are frozen as inert fully-occupied: they carry no dynamics
! and no current (symmetric filled bands), so they only enter the electron
! bookkeeping via gs%ne = ne - focc*(nb_min-1) (spinless: 2 electrons per
! band; spinor/noncollinear: 1 electron per band -- Kramers partners are
! separate bands).  nb_min = 1, nb_hi = nb reproduces the unwindowed behavior
! exactly.
subroutine init_sbe_gs_info(gs, sysname, gs_directory, nk, nb, nb_min, nb_hi, ne, a1, a2, a3, read_bin, icomm)
    use communication
    use filesystem, only: open_filehandle, get_filehandle
    use common_ssbe, only: grad_k_array_nb1d_double
    use salmon_global, only: gauge_sbe, file_sbe_prod_dk, sbe_lg_degen, num_kgrid, sbe_lg_degen_floor, spin, &
                           & yn_sbe_vnl_exact, file_sbe_vnl_kappa, nelec, natom
    use degenerate_block_ssbe, only: build_xi, same_block, blend, theta_on, theta_off, &
                                   & build_block_transport, &
                                   & prod_dk_axis_slot, gicov_prod_dk_nbvec, expected_prod_dk_nrec, &
                                   & check_gicov_occupation
    implicit none
    type(s_sbe_gs_info), intent(inout) :: gs
    character(*), intent(in) :: sysname
    character(*), intent(in) :: gs_directory
    integer, intent(in) :: nk
    integer, intent(in) :: nb
    integer, intent(in) :: nb_min
    integer, intent(in) :: nb_hi
    integer, intent(in) :: ne
    real(8), intent(in) :: a1(1:3), a2(1:3), a3(1:3)
    logical, intent(in) :: read_bin
    integer, intent(in) :: icomm
    integer :: irank, nproc
    integer :: nb_eff, ndeg
    logical :: is_metal

    call comm_get_groupinfo(icomm, irank, nproc)

    if (nb_min < 1 .or. nb_min > nb_hi .or. nb_hi > nb) then
        if (irank == 0) write(*, '(a,i0,a,i0,a,i0)') &
            & "ERROR(init_sbe_gs_info): band window out of range: [", nb_min, &
            & ", ", nb_hi, "] of export nb = ", nb
        stop 1
    end if
    nb_eff = nb_hi - (nb_min - 1)

    ! Occupation convention: spinless bands hold 2 electrons, spinor
    ! (noncollinear/SOC) bands hold 1 (the input checker enforces the
    ! supported SOC combinations before we get here).
    gs%spinor = (trim(spin) == 'noncollinear')
    ndeg = 2
    if (gs%spinor) ndeg = 1
    gs%focc = dble(ndeg)

    gs%nk = nk
    gs%nb = nb_eff
    gs%ne = ne - ndeg * (nb_min - 1)
    gs%nvb = gs%ne / ndeg
    !gs%num_kgrid(1:3) = num_kgrid(1:3)

    if (irank == 0 .and. (nb_min > 1 .or. nb_hi < nb)) then
        write(*, '(a,i0,a,i0,a,i0,a,i0,a,i0)') "# band window: bands ", &
            & nb_min, "..", nb_hi, " of the ", nb, "-band export -> nb_eff = ", &
            & nb_eff, ", nelec_eff = ", gs%ne
    end if

    !Calculate b_matrix, volume_cell and volume_bz from a1..a3 vector.
    call calc_lattice_info()

    allocate(gs%kpoint(1:3, 1:nk))
    allocate(gs%kweight(1:nk))
    allocate(gs%eigen(1:nb_eff, 1:nk))
    allocate(gs%occup(1:nb_eff, 1:nk))
    allocate(gs%delta_omega(1:nb_eff, 1:nb_eff, 1:nk))
    allocate(gs%p_mod_matrix(1:nb_eff, 1:nb_eff, 1:3, 1:nk))
    if (.not. uses_xfull_links(sbe_lg_degen)) then
        allocate(gs%d_matrix(1:nb_eff, 1:nb_eff, 1:3, 1:nk))
    end if
    allocate(gs%p_tm_matrix(1:nb_eff, 1:nb_eff, 1:3, 1:nk))
    allocate(gs%rvnl_tm_matrix(1:nb_eff, 1:nb_eff, 1:3, 1:nk))
    allocate(gs%grad_k_eigen(1:nb_eff, 1:3, 1:nk))

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

    ! Metal detector: gs%occup already holds the REAL per-(band,k) GS
    ! occupation just read from file and broadcast (every rank sees identical
    ! data). Run this BEFORE read_prod_dk_data/prepare_matrix/build_xi below
    ! so a metal misconfigured with sbe_lg_degen='gi'/'gifix' (architecturally
    ! insulator-only: both assume a gap-isolated manifold) is rejected up
    ! front rather than after paying for xi/block-transport construction on
    ! data those modes cannot correctly consume.
    call check_gicov_occupation(gs%nb, gs%nk, gs%occup, gs%focc, is_metal)
    gs%is_metal = is_metal
    if (is_metal .and. (trim(sbe_lg_degen) == 'gi' .or. trim(sbe_lg_degen) == 'gifix')) then
        if (irank == 0) write(*, '(a)') &
            & "ERROR(init_sbe_gs_info): metal-like fractional/k-varying occupation " // &
            & "detected; sbe_lg_degen='gi'/'gifix' assume a gap-isolated insulating " // &
            & "manifold -- use sbe_lg_degen='gicov' or gauge_sbe='velocity_gauge'."
        stop 1
    end if

    !Retrieve k-space overlap products from 'file_sbe_prod_dk' (LG-SBE degeneracy):
    if (uses_prod_dk(sbe_lg_degen)) call read_prod_dk_data()

    !VG completion: read the nonlocal kappa-stencil (V/W band matrices) and
    !overwrite rvnl_tm_matrix with the binary-precision W_{i,0}.  MUST run
    !before prepare_matrix() so p_mod_matrix inherits the overwrite.
    if (yn_sbe_vnl_exact == 'y') call read_vnl_kappa_data()

    !Calculate omega and d_matrix (neglecting diagonal part):
    if (irank == 0) write(*,"(a)") "# prepare_matrix"

    call prepare_matrix()
    call comm_bcast(gs%p_mod_matrix, icomm, 0)
    call comm_bcast(gs%delta_omega, icomm, 0)
    if (.not. uses_xfull_links(sbe_lg_degen)) then
        call comm_bcast(gs%d_matrix, icomm, 0) ! Experimental
    end if

    select case(trim(gauge_sbe))
    case ("length_gauge")
        call grad_k_array_nb1d_double(gs%nb, gs%nk, gs%b_matrix,  &
                                  &   gs%eigen, gs%grad_k_eigen)
    end select

    !Initial Occupation Number
    !Window bands 1..gs%nvb (= absolute nb_min..nb_min-1+nvb) are the occupied
    !bands inside the window (nvb = ne/2 spinless, = ne spinor); the frozen
    !bands 1..nb_min-1 are NOT stored (inert, fully occupied) and enter only
    !through gs%ne = ne - focc*(nb_min-1).
    !
    !gs%is_metal was already determined right after gs%occup was read+broadcast
    !(above, before read_prod_dk_data/prepare_matrix); reuse it here rather
    !than recomputing.
    if (.not. gs%is_metal) then
        ! insulator/semiconductor: keep the exact T=0 rigid fill -- bit-identical
        ! to today for Si/graphene (this is a REGRESSION GUARD, not a stylistic
        ! choice: it also shields against a not-fully-T=0 or slightly-
        ! unconverged insulator SCF perturbing gs%occup away from the exact
        ! 0/focc step the rest of the SBE machinery assumes for gap systems).
        gs%occup(:,:) = 0d0 !!Experimental!!
        if (gs%ne > 0) gs%occup(1:gs%nvb,:) = gs%focc !!Experimental!!
    else
        ! metal: keep the REAL per-(band,k) GS occupation just read from file
        ! (already in gs%occup from read_eigen_data) -- do NOT rigid-fill it.
        ! (sbe_lg_degen='gi'/'gifix' were already rejected above for a metal.)
        if (irank == 0) write(*, '(a)') &
            & "# metal-like occupation detected: using REAL per-(band,k) GS occupation " // &
            & "(Fermi-Dirac), NOT the insulator rigid fill."
    end if

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
        ! Noncollinear/SOC eigen exports (write_eigen) carry TWO spin blocks
        ! (is=1,2) that are DUPLICATE copies of the same spinor spectrum (see
        ! occupation.f90 mu2ne: jspin=2 duplicates jspin=1 in SO mode).  The
        ! loop below reads exactly the first nk k-blocks = the spin-1 block,
        ! then the file is closed -- the duplicate spin-2 block is skipped by
        ! construction, so the SAME reader serves spinless and spinor exports.
        do ik = 1, nk
            read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
            ! consume ALL nb export rows (keeps the per-k header alignment);
            ! store only the band window [nb_min : nb_hi] at window index ib-nb_min+1
            do ib = 1, nb
                read(fh, *) iib, tmp(1:2)
                if (ib .ne. iib) stop "ib mismatch"
                if (ib >= nb_min .and. ib <= nb_hi) gs%eigen(ib - nb_min + 1, ik) = tmp(1) * efac
                ! real per-(band,k) GS occupation (column 3, dimensionless --
                ! no efac unit conversion needed); same window index as eigen.
                if (ib >= nb_min .and. ib <= nb_hi) gs%occup(ib - nb_min + 1, ik) = tmp(2)
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
        ! consume ALL nb*nb export records per k (keeps the record alignment);
        ! store only pairs with both bands inside the window [nb_min : nb_hi]
        do ik = 1, nk
            do ib = 1, nb
                do jb = 1, nb
                    read(fh, *) iik, iib, jjb, tmp(1:6)
                    if (ik .ne. iik) stop "ik mismatch"
                    if (ib .ne. iib) stop "ib mismatch"
                    if (jb .ne. jjb) stop "jb mismatch"
                    if (ib >= nb_min .and. ib <= nb_hi .and. &
                        & jb >= nb_min .and. jb <= nb_hi) then
                        gs%p_tm_matrix(ib - nb_min + 1, jb - nb_min + 1, 1, ik) = dcmplx(tmp(1), tmp(2))
                        gs%p_tm_matrix(ib - nb_min + 1, jb - nb_min + 1, 2, ik) = dcmplx(tmp(3), tmp(4))
                        gs%p_tm_matrix(ib - nb_min + 1, jb - nb_min + 1, 3, ik) = dcmplx(tmp(5), tmp(6))
                    end if
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
                    if (ib >= nb_min .and. ib <= nb_hi .and. &
                        & jb >= nb_min .and. jb <= nb_hi) then
                        gs%rvnl_tm_matrix(ib - nb_min + 1, jb - nb_min + 1, 1, ik) = dcmplx(tmp(1), tmp(2))
                        gs%rvnl_tm_matrix(ib - nb_min + 1, jb - nb_min + 1, 2, ik) = dcmplx(tmp(3), tmp(4))
                        gs%rvnl_tm_matrix(ib - nb_min + 1, jb - nb_min + 1, 3, ik) = dcmplx(tmp(5), tmp(6))
                    end if
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


    ! Read the k-space overlap products <u_{io,ik}|u_{jo,ik+dk}> from
    ! 'file_sbe_prod_dk', auto-detecting the on-disk format:
    !   binary (magic 'SBEPRDK1', Task-5 write_prod_dk_data collective
    !     MPI-IO writer) -> collective MPI-IO reader below, mirroring
    !     read_vnl_kappa_data's structure (Task 4): collective open ->
    !     rank-0 header read + fail-closed validation + bcast ->
    !     split_range k-slice -> collective read_at_all -> finite +
    !     whole-file XOR-checksum validate (ALWAYS present in this header,
    !     unlike vnl_kappa's optional sidecar) -> slot-map into a local
    !     full-nk zero-fill buffer using the SAME reduced/gicov slot logic
    !     as the ASCII path -> chunked comm_summation reconstruct (disjoint
    !     per-rank k-slices tile 1..nk; N ranks wrote, any M ranks read the
    !     same canonical bytes = nproc-portable).
    !   ASCII (legacy list-directed, root-serial) -> the ORIGINAL parser,
    !     kept fully intact below as the rare pre-conversion safety net.
    ! Rank 0 decides the format from the first 8 bytes and broadcasts the
    ! verdict so every rank takes the identical branch.
    subroutine read_prod_dk_data()
        use util_ssbe, only: split_range
        use mpiio_export_ssbe, only: mpiio_open_read, mpiio_close, mpiio_read_at_all_z, &
                                   & mpiio_rd_c8, mpiio_rd_i, mpiio_rd_i8, mpiio_disp_k, &
                                   & mpiio_finite_z, mpiio_csum_local_z, mpiio_csum_reduce, &
                                   & mpiio_assert_dim, &
                                   & HDR_PRODK, MAGIC_PRODK, SBE_FMT_VERSION
        implicit none
        integer :: fh, ios, ip, ierr
        logical :: file_exists
        character(256) :: cline
        integer :: file_no, file_nk, file_nk1, file_nk2, file_nk3, file_ndk
        integer :: ndk_loc, mdk, ivec, nbvec_full, islot
        integer(8) :: nrec, irec, ik_exp, per_k_records
        logical :: reduced
        integer :: jdk1, jdk2, jdk3
        integer :: ik_r, ik1_r, ik2_r, ik3_r, jdk1_r, jdk2_r, jdk3_r, io_r, jo_r
        real(8) :: re_r, im_r

        ! ---- auto-detect + binary-path locals (Task 6) ----
        integer, parameter :: OFF = kind(HDR_PRODK)
        integer(OFF) :: disp
        integer :: fh_sniff, ierr_h
        logical :: is_bin, has_slice
        character(8) :: magic8
        integer :: f_ver, f_nk, f_no, f_ndk, f_nbvec_full, ihdr_ver(1), ihdr4(4)
        integer(8) :: fsize, expect, csum_expect, g0, nz, nz_k, csum_loc, csum_glob
        integer :: ifin_loc, ifin_any
        integer :: ik, jo, io, jd1, jd2, jd3, iv
        integer :: ik_lo, ik_hi, max_tile
        integer, allocatable :: itbl_min(:), itbl_max(:)
        complex(8), allocatable :: blk(:, :, :, :, :, :)
        complex(8), allocatable :: prod_dk_l(:, :, :, :)
        complex(8) :: dummy_z(1)
        ! Task-6 numeric round-trip gate: flag-gated (env SBE_PRODDK_DUMP=1) debug
        ! dump of gs%prod_dk + gs%bvec, shared by BOTH branches (placed after the
        ! is_bin block below) -- mirrors read_vnl_kappa_data's SBE_RVNL_DUMP gate.
        character(8) :: dbgval
        integer :: dbglen, dbgfh, npair, jk_d, ip_d
        integer :: idump_ik(3), ib_list(3), jb_list(3)

        ierr    = 0
        gs%nbvec = 0
        ndk_loc = 0
        file_no = 0
        nbvec_full = 0
        reduced = .false.
        is_bin  = .false.
        fsize   = 0_8

        ! --- rank-0: existence guard + 8-byte format sniff (plain sequential
        ! stream I/O -- a peek that only decides the branch, not a data
        ! read; a short/garbled peek falls through to the ASCII parser
        ! rather than hard-failing here) ---
        if (irank == 0) then
            if (len_trim(file_sbe_prod_dk) == 0) then
                write(*, '(a)') &
                    & "ERROR(read_prod_dk_data): 'file_sbe_prod_dk' is empty while 'sbe_lg_degen'='gi/gifix/gicov'."
                ierr = 1
            else
                inquire(file=trim(file_sbe_prod_dk), exist=file_exists, size=fsize)
                if (.not. file_exists) then
                    write(*, '(a)') "ERROR(read_prod_dk_data): file not found: " // trim(file_sbe_prod_dk)
                    ierr = 1
                end if
            end if
            if (ierr == 0) then
                fh_sniff = get_filehandle()
                open(fh_sniff, file=trim(file_sbe_prod_dk), form='unformatted', access='stream', &
                    & status='old', action='read', iostat=ios)
                if (ios /= 0) then
                    write(*, '(a)') "ERROR(read_prod_dk_data): cannot open file for format sniff."
                    ierr = 1
                else
                    magic8 = ' '
                    read(fh_sniff, iostat=ios) magic8
                    close(fh_sniff)
                    is_bin = (ios == 0 .and. magic8 == MAGIC_PRODK)
                end if
            end if
        end if
        call comm_bcast(ierr, icomm, 0)
        if (ierr /= 0) stop 1
        call comm_bcast(is_bin, icomm, 0)

        if (is_bin) then
            !=================================================================
            ! BINARY PATH: collective MPI-IO (SBEPRDK1 header, Task-5 writer)
            !=================================================================
            if (irank == 0) write(*, '(a)') "# read_prod_dk_data (binary, MPI-IO)"

            ! ---- collective open of the single shared file (all ranks) ----
            call mpiio_open_read(trim(file_sbe_prod_dk), icomm, fh, ierr)

            ! ---- rank-0 header read + fail-closed validation (fields that
            ! do not need host nk/nb_hi -- magic/version/ndk-sign/nbvec
            ! consistency/exact-size) ----
            if (irank == 0) then
                ierr_h = 0
                ! codex review (Task 6): default these to 0 BEFORE the reads so a
                ! partial-header read failure (ierr_h/=0, already fail-closed via
                ! the ierr==0 guards below) never unpacks undefined memory into
                ! f_ver/f_nk/f_no/f_ndk/f_nbvec_full/csum_expect, even transiently.
                magic8 = ' ';  ihdr_ver = 0;  ihdr4 = 0;  csum_expect = 0_8
                call mpiio_rd_c8(fh,  0_OFF, magic8, ierr);      ierr_h = ior(ierr_h, ierr)
                call mpiio_rd_i (fh,  8_OFF, ihdr_ver, 1, ierr); ierr_h = ior(ierr_h, ierr)
                call mpiio_rd_i (fh, 12_OFF, ihdr4, 4, ierr);    ierr_h = ior(ierr_h, ierr)
                call mpiio_rd_i8(fh, 28_OFF, csum_expect, ierr); ierr_h = ior(ierr_h, ierr)
                f_ver = ihdr_ver(1)
                f_nk = ihdr4(1);  f_no = ihdr4(2);  f_ndk = ihdr4(3);  f_nbvec_full = ihdr4(4)

                ierr = 0
                if (ierr_h /= 0) then
                    write(*, '(a)') "ERROR(read_prod_dk_data): header field read failed (MPI_File_read_at)."
                    ierr = 1
                end if
                if (ierr == 0 .and. magic8 /= MAGIC_PRODK) then
                    write(*, '(a)') "ERROR(read_prod_dk_data): bad magic (not a binary prod_dk export)."
                    ierr = 1
                end if
                if (ierr == 0 .and. f_ver /= SBE_FMT_VERSION) then
                    write(*, '(a,i0)') "ERROR(read_prod_dk_data): unsupported prod_dk format version = ", f_ver
                    ierr = 1
                end if
                if (ierr == 0 .and. f_ndk < 0) then
                    write(*, '(a)') "ERROR(read_prod_dk_data): ndk in file is negative."
                    ierr = 1
                end if
                if (ierr == 0 .and. f_ndk > 100) then
                    write(*, '(a)') "ERROR(read_prod_dk_data): implausible ndk in file header (corrupt?)."
                    ierr = 1
                end if
                if (ierr == 0 .and. f_nbvec_full /= (2*f_ndk+1)**3) then
                    write(*, '(a)') &
                        & "ERROR(read_prod_dk_data): nbvec_full in file is inconsistent with ndk (corrupt header)."
                    ierr = 1
                end if
                if (ierr == 0) then
                    ! exact size: header (36 B) + nk*nbvec_full*no^2 complex(8).
                    ! Catches truncation/surplus that a short-read check alone
                    ! would miss (mirrors read_vnl_kappa_data's expect gate).
                    expect = HDR_PRODK + int(f_nk,8) * int(f_nbvec_full,8) * int(f_no,8) * int(f_no,8) * 16_8
                    if (fsize /= expect) then
                        write(*, '(a,i0,a,i0,a)') "ERROR(read_prod_dk_data): file size ", fsize, &
                            & " differs from expected ", expect, " (truncated or surplus data)."
                        ierr = 1
                    end if
                end if
            end if
            call comm_bcast(ierr, icomm, 0)
            if (ierr /= 0) stop 1

            ! ---- broadcast header fields, then the two checks that need
            ! host nk/nb_hi, run identically (collective-safe) on every rank
            ! ----
            call comm_bcast(f_nk, icomm, 0)
            call comm_bcast(f_no, icomm, 0)
            call comm_bcast(f_ndk, icomm, 0)
            call comm_bcast(f_nbvec_full, icomm, 0)

            call mpiio_assert_dim(f_nk, nk, 'prod_dk nk', icomm)
            if (f_no < nb_hi) then
                if (irank == 0) write(*, '(a)') &
                    & "ERROR(read_prod_dk_data): band window in file is smaller than the SBE window top."
                stop 1
            end if

            ! ---- slot-map bookkeeping: SAME decision as the ASCII path
            ! (gicov_prod_dk_nbvec), now a pure function of already-
            ! broadcast header fields -> computed identically on every rank
            ! ----
            ndk_loc    = f_ndk
            nbvec_full = f_nbvec_full
            gs%nbvec   = gicov_prod_dk_nbvec(sbe_lg_degen, ndk_loc, nbvec_full)
            reduced    = (gs%nbvec /= nbvec_full)

            ! ---- allocate on ALL ranks (same shape as the ASCII path) ----
            allocate(gs%bvec(1:3, 1:gs%nbvec))
            allocate(gs%prod_dk(1:nb_eff, 1:nb_eff, 1:gs%nbvec, 1:nk))
            gs%bvec(:, :)          = 0
            gs%prod_dk(:, :, :, :) = dcmplx(0d0, 0d0)

            ! ---- dk-shift table: identical construction to the ASCII path
            ! ----
            mdk = 2 * ndk_loc + 1
            if (reduced) then
                gs%bvec(1:3, 1) = (/ 1, 0, 0 /)
                gs%bvec(1:3, 2) = (/ 0, 1, 0 /)
                gs%bvec(1:3, 3) = (/ 0, 0, 1 /)
            else
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
            end if

            ! ---- k-slice partition: identical split_range formula used by
            ! the vnl_kappa reader and the solver (nproc-portable: N ranks
            ! wrote, any M ranks read the same canonical bytes) ----
            allocate(itbl_min(0:nproc-1), itbl_max(0:nproc-1))
            call split_range(1, nk, nproc, itbl_min, itbl_max)
            ik_lo = itbl_min(irank);  ik_hi = itbl_max(irank)
            deallocate(itbl_min, itbl_max)
            has_slice = (ik_hi >= ik_lo)   ! .false. only when nproc > nk

            ! ---- collective read of this rank's contiguous k-slice (the
            ! FULL shell nbvec_full -- the on-disk layout does not shrink
            ! just because 'reduced' keeps fewer columns) ----
            nz_k = int(f_no,8) * int(f_no,8) * int(nbvec_full,8)
            max_tile = (nk + nproc - 1) / nproc
            if (nz_k * int(max_tile,8) > int(huge(0),8)) then
                if (irank == 0) write(*, '(a)') &
                    & "ERROR(read_prod_dk_data): per-slice element count exceeds MPI int count (chunked read needed)."
                stop 1
            end if

            if (has_slice) then
                allocate(blk(1:f_no, 1:f_no, -ndk_loc:ndk_loc, -ndk_loc:ndk_loc, -ndk_loc:ndk_loc, ik_lo:ik_hi))
                disp = mpiio_disp_k(HDR_PRODK, nz_k, ik_lo)
                nz   = nz_k * int(ik_hi - ik_lo + 1, 8)
                call mpiio_read_at_all_z(fh, disp, blk(1,1,-ndk_loc,-ndk_loc,-ndk_loc,ik_lo), nz, icomm, &
                    & 'prod_dk', ierr)
            else
                allocate(blk(1:f_no, 1:f_no, -ndk_loc:ndk_loc, -ndk_loc:ndk_loc, -ndk_loc:ndk_loc, 1:1))  ! placeholder
                nz = 0_8
                call mpiio_read_at_all_z(fh, HDR_PRODK, dummy_z, 0_8, icomm, 'prod_dk', ierr)  ! collective no-op
            end if

            ! ---- fail-closed validation: (3) finiteness, collective-safe
            ! (LAND via SUM of 0/1, matches read_vnl_kappa_data's idiom);
            ! (4) whole-file XOR checksum -- ALWAYS present in this header
            ! (unlike vnl_kappa's optional sidecar), so csum_expect stays
            ! rank-0-only (no int64 bcast needed: only rank 0 compares) ----
            ifin_loc = 0
            if (has_slice) then
                if (.not. mpiio_finite_z(blk(1,1,-ndk_loc,-ndk_loc,-ndk_loc,ik_lo), nz)) ifin_loc = 1
            end if
            call comm_summation(ifin_loc, ifin_any, icomm)
            if (ifin_any > 0) then
                if (irank == 0) write(*, '(a)') &
                    & "ERROR(read_prod_dk_data): non-finite element (NaN/Inf) in prod_dk data."
                stop 1
            end if

            csum_loc = 0_8
            if (has_slice) then
                g0 = nz_k * int(ik_lo - 1, 8)
                csum_loc = mpiio_csum_local_z(blk(1,1,-ndk_loc,-ndk_loc,-ndk_loc,ik_lo), nz, g0)
            end if
            csum_glob = mpiio_csum_reduce(csum_loc, icomm)
            if (irank == 0 .and. csum_glob /= csum_expect) then
                write(*, '(a)') "ERROR(read_prod_dk_data): checksum mismatch (data corruption / truncated write)."
                ierr = 1
            end if
            call comm_bcast(ierr, icomm, 0)
            if (ierr /= 0) stop 1

            ! ---- slot-map this rank's k-slice into a local full-nk
            ! zero-fill buffer, using the SAME reduced/gicov slot logic as
            ! the ASCII path (prod_dk_axis_slot / linear-pack), windowed to
            ! [nb_min:nb_hi] exactly like the ASCII path's io_r/jo_r window
            ! check ----
            allocate(prod_dk_l(1:nb_eff, 1:nb_eff, 1:gs%nbvec, 1:nk))
            prod_dk_l(:, :, :, :) = dcmplx(0d0, 0d0)
            do ik = ik_lo, ik_hi
                do jd3 = -ndk_loc, ndk_loc
                    do jd2 = -ndk_loc, ndk_loc
                        do jd1 = -ndk_loc, ndk_loc
                            if (reduced) then
                                iv = prod_dk_axis_slot(jd1, jd2, jd3)
                            else
                                iv = (jd3 + ndk_loc) * mdk * mdk + (jd2 + ndk_loc) * mdk + (jd1 + ndk_loc) + 1
                            end if
                            if (iv > 0) then
                                do jo = nb_min, nb_hi
                                    do io = nb_min, nb_hi
                                        prod_dk_l(io - nb_min + 1, jo - nb_min + 1, iv, ik) = &
                                            & blk(io, jo, jd1, jd2, jd3, ik)
                                    end do
                                end do
                            end if
                        end do
                    end do
                end do
            end do

            if (allocated(blk)) deallocate(blk)
            call mpiio_close(fh, ierr)   ! collective (all ranks)

            ! ---- chunked full-nk reconstruct: disjoint per-rank k-slices
            ! tile 1..nk, so comm_summation (allreduce SUM over the
            ! zero-fill buffers) IS the reconstruction -- reuses the ASCII
            ! path's max_elems_per_chunk=200000000_8 bound (comm_summation's
            ! count is default-integer, same 32-bit overflow risk as
            ! comm_bcast) ----
            block
                integer(8) :: elems_per_k, max_elems_per_chunk
                integer :: chunk_nk, ik0, ik1
                elems_per_k = int(nb_eff, 8) * int(nb_eff, 8) * int(gs%nbvec, 8)
                ! codex review (Task 6): chunk_nk=max(1,...) only bounds the chunk to
                ! ONE k-plane -- if a single k-plane's element count itself exceeded
                ! the default-integer count passed to comm_summation, the chunking
                ! loop could not shrink further and the cast below would wrap.  Not
                ! reachable at realistic nb_eff (needs nb_eff > ~46000 at nbvec=1, far
                ! beyond any SBE band window), but guard it explicitly and fail closed
                ! rather than silently wrap, matching this subroutine's own
                ! nz_k*max_tile guard above.
                if (elems_per_k > int(huge(0), 8)) then
                    if (irank == 0) write(*, '(a)') &
                        & "ERROR(read_prod_dk_data): per-k element count exceeds MPI int count (reconstruct needs finer chunking)."
                    stop 1
                end if
                max_elems_per_chunk = 200000000_8
                chunk_nk = max(1, int(max_elems_per_chunk / max(elems_per_k, 1_8)))
                ik0 = 1
                do while (ik0 <= nk)
                    ik1 = min(nk, ik0 + chunk_nk - 1)
                    call comm_summation(prod_dk_l(:,:,:,ik0:ik1), gs%prod_dk(:,:,:,ik0:ik1), &
                        & int(elems_per_k * int(ik1-ik0+1,8)), icomm)
                    ik0 = ik1 + 1
                end do
            end block
            deallocate(prod_dk_l)

            if (irank == 0) write(*, '(a)') "#   prod_dk: binary (MPI-IO) k-slice read, checksum verified."

        else
        !=====================================================================
        ! ASCII PATH (legacy, root-serial list-directed parser) -- kept
        ! UNCHANGED below as the rare pre-conversion safety net.
        !=====================================================================
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
                if (file_no < nb_hi) then
                    write(*, '(a)') "ERROR(read_prod_dk_data): band window in file is smaller than the SBE window top."
                    ierr = 1
                end if
                if (file_ndk < 0) then
                    write(*, '(a)') "ERROR(read_prod_dk_data): ndk in file is negative."
                    ierr = 1
                end if
            end if

            if (ierr == 0) then
                ndk_loc    = file_ndk
                nbvec_full = (2 * ndk_loc + 1) ** 3
                ! gicov memory-scale fix: keep only the +x/+y/+z unit-shift
                ! columns for 'gicov' (ndk>=1) -- build_block_transport,
                ! degenerate_block_ssbe.f90, is the only prod_dk consumer on
                ! this path and never reads anything else. gi/gifix/off (and
                ! gicov with ndk==0) are unchanged (see gicov_prod_dk_nbvec).
                gs%nbvec = gicov_prod_dk_nbvec(sbe_lg_degen, ndk_loc, nbvec_full)
                reduced  = (gs%nbvec /= nbvec_full)
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
        allocate(gs%prod_dk(1:nb_eff, 1:nb_eff, 1:gs%nbvec, 1:nk))
        gs%bvec(:, :)          = 0
        gs%prod_dk(:, :, :, :) = dcmplx(0d0, 0d0)

        ! --- root: build the dk-shift table and read all 11-column records ---
        if (irank == 0) then
            mdk  = 2 * ndk_loc + 1
            if (reduced) then
                ! gicov: bvec holds only the 3 kept +unit-shift directions, in
                ! the SAME fixed slot order prod_dk_axis_slot maps records to
                ! (so gs%bvec(:,islot) always names the column gs%prod_dk(:,
                ! :,islot,:) actually holds -- find_bvec's callers, e.g.
                ! build_block_transport, don't care about slot order, only
                ! about which shift a given column is).
                gs%bvec(1:3, 1) = (/ 1, 0, 0 /)
                gs%bvec(1:3, 2) = (/ 0, 1, 0 /)
                gs%bvec(1:3, 3) = (/ 0, 0, 1 /)
            else
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
            end if

            ! writer emits nk * (2*ndk+1)**3 * no**2 records, ik outermost/slowest
            ! -- ALWAYS the full shell (nbvec_full), regardless of what 'reduced'
            ! keeps: the file layout does not shrink just because the reader
            ! keeps fewer columns. int64 so it stays exact past 2**31 (hit
            ! already at a k20^3=8000 k-point mesh for a realistic file_no).
            nrec          = expected_prod_dk_nrec(nk, nbvec_full, file_no)
            per_k_records = nrec / max(int(nk, 8), 1_8)
            do irec = 1, nrec
                read(fh, *, iostat=ios) &
                    & ik_r, ik1_r, ik2_r, ik3_r, jdk1_r, jdk2_r, jdk3_r, io_r, jo_r, re_r, im_r
                if (ios /= 0) then
                    write(*, '(a)') "ERROR(read_prod_dk_data): fewer records than expected."
                    ierr = 1
                    exit
                end if
                ! ik must run 1..nk in contiguous blocks (matches SBE k ordering)
                ik_exp = (irec - 1_8) / per_k_records + 1_8
                if (int(ik_r, 8) /= ik_exp) then
                    write(*, '(a)') "ERROR(read_prod_dk_data): ik ordering does not match SBE k-grid."
                    ierr = 1
                    exit
                end if
                if (abs(jdk1_r) > ndk_loc .or. abs(jdk2_r) > ndk_loc .or. abs(jdk3_r) > ndk_loc) then
                    write(*, '(a)') "ERROR(read_prod_dk_data): dk-shift index out of range."
                    ierr = 1
                    exit
                end if
                ! keep only the SBE band window [nb_min : nb_hi] (writer emits
                ! the full band window); store at window indices io/jo - nb_min + 1
                if (io_r >= nb_min .and. io_r <= nb_hi .and. &
                    & jo_r >= nb_min .and. jo_r <= nb_hi) then
                    if (reduced) then
                        ! gicov: drop everything except the +x/+y/+z columns
                        ! (see the header comment above prod_dk_axis_slot,
                        ! degenerate_block_ssbe.f90, for why this is safe).
                        islot = prod_dk_axis_slot(jdk1_r, jdk2_r, jdk3_r)
                        if (islot > 0) then
                            gs%prod_dk(io_r - nb_min + 1, jo_r - nb_min + 1, islot, ik_r) = dcmplx(re_r, im_r)
                        end if
                    else
                        ivec = (jdk3_r + ndk_loc) * mdk * mdk &
                            & + (jdk2_r + ndk_loc) * mdk &
                            & + (jdk1_r + ndk_loc) + 1
                        gs%prod_dk(io_r - nb_min + 1, jo_r - nb_min + 1, ivec, ik_r) = dcmplx(re_r, im_r)
                    end if
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

        ! --- propagate record-stage error, then broadcast the data (prod_dk
        ! is bcast in k-CHUNKS so no single MPI_Bcast call's element COUNT
        ! risks the same 32-bit 'size(val)' class of overflow the nrec fix
        ! above targets -- comm_bcast_array4d_dcomplex, communication.f90,
        ! passes a default-integer size(val) to MPI_Bcast with no kind guard;
        ! reduced (gicov) storage already cuts this ~9x at ndk=1, this is
        ! defense-in-depth for the remaining gi/gifix/large-nb_eff cases) ---
        call comm_bcast(ierr, icomm, 0)
        if (ierr /= 0) stop 1
        call comm_bcast(gs%bvec, icomm, 0)
        block
            integer(8) :: elems_per_k, max_elems_per_chunk
            integer :: chunk_nk, ik0, ik1
            elems_per_k = int(nb_eff, 8) * int(nb_eff, 8) * int(gs%nbvec, 8)
            max_elems_per_chunk = 200000000_8  ! << 2**31-1, generous headroom
            chunk_nk = max(1, int(max_elems_per_chunk / max(elems_per_k, 1_8)))
            ik0 = 1
            do while (ik0 <= nk)
                ik1 = min(nk, ik0 + chunk_nk - 1)
                call comm_bcast(gs%prod_dk(:, :, :, ik0:ik1), icomm, 0)
                ik0 = ik1 + 1
            end do
        end block
        end if

        ! ---- Task-6 numeric round-trip gate: optional debug dump of the
        !      final gs%prod_dk + gs%bvec at a FIXED small element set (full
        !      precision es24.16).  OFF unless env SBE_PRODDK_DUMP=1
        !      (production-safe default); rank-0 only, deterministic
        !      indices.  Runs AFTER the is_bin branch above, so it exercises
        !      whichever path (binary MPI-IO or legacy ASCII) actually ran
        !      -- this is what the binary-vs-ASCII and M4-vs-M8 read-test
        !      dumps are diffed on. ----
        if (irank == 0) then
            call get_environment_variable('SBE_PRODDK_DUMP', dbgval, dbglen)
            if (dbglen > 0 .and. trim(dbgval) == '1') then
                idump_ik(1) = 1;  idump_ik(2) = (nk + 1) / 2;  idump_ik(3) = nk
                npair = 1;  ib_list(1) = 1;  jb_list(1) = 1
                if (nb_eff >= 2) then
                    ib_list(2) = 2;  jb_list(2) = 1
                    ib_list(3) = 1;  jb_list(3) = 2;  npair = 3
                end if
                dbgfh = get_filehandle()
                open(dbgfh, file='proddk_dump.txt', status='replace', action='write')
                write(dbgfh, '(a)') "# prod_dk dump: io jo iv ik  Re Im"
                write(dbgfh, '(a,i0,a,i0,a,i0)') "# nb_eff=", nb_eff, " nk=", nk, " nbvec=", gs%nbvec
                do iv = 1, gs%nbvec
                    do jk_d = 1, 3
                        if (idump_ik(jk_d) < 1 .or. idump_ik(jk_d) > nk) cycle
                        do ip_d = 1, npair
                            write(dbgfh, '(4(i6,1x),2es24.16)') ib_list(ip_d), jb_list(ip_d), iv, idump_ik(jk_d), &
                                & real (gs%prod_dk(ib_list(ip_d), jb_list(ip_d), iv, idump_ik(jk_d))), &
                                & aimag(gs%prod_dk(ib_list(ip_d), jb_list(ip_d), iv, idump_ik(jk_d)))
                        end do
                    end do
                end do
                write(dbgfh, '(a)') "# bvec: iv  jdk1 jdk2 jdk3"
                do iv = 1, gs%nbvec
                    write(dbgfh, '(i6,1x,3(i4,1x))') iv, gs%bvec(1,iv), gs%bvec(2,iv), gs%bvec(3,iv)
                end do
                close(dbgfh)
                write(*, '(a)') "#   prod_dk: SBE_PRODDK_DUMP=1 -> wrote proddk_dump.txt (Task-6 numeric gate)."
            end if
        end if
    end subroutine read_prod_dk_data


    ! VG completion reader (collective MPI-IO): nonlocal kappa-stencil V/W band
    ! matrices from the single shared file written by write_sbe_vnl_kappa_data
    ! (write.f90).  All ranks collectively MPI_File_open the file, rank 0 reads
    ! and fail-closed validates the 156-byte header (magic/ndim/nk/num_kgrid/nb/
    ! ns/h/dir/nelec/natom/lattice + exact file-size), then each rank
    ! collectively reads its own split_range k-slice (MPI_Get_count short-read is
    ! fail-closed).  Validation layers: finiteness always; the whole-file
    ! splitmix64 XOR checksum vs the '_vnl_kappa.chk' sidecar when present, else
    ! legacy count+finite only (byte-identical OLD root-serial files carry no
    ! sidecar, so the reader stays backward-compatible).  Every rank keeps (i)
    ! the full-nk s=0 velocity W_{i,0} by OVERWRITING gs%rvnl_tm_matrix
    ! (zero-fill + comm_summation reconstruct; binary precision; feeds
    ! p_mod/fsum_D) -- AFTER an all-reduced (LOR) staleness cross-check of
    ! W_{i,0} against the PRE-overwrite text-file rvnl_tm_matrix (catches a stale
    ! stencil file from another GS at 1e-3 rel-to-max; a wiring/staleness check,
    ! not a precision check) -- and (ii) the windowed V/W stencil for its own
    ! k-slice only.  nproc-portable: N ranks write, any M ranks read the same
    ! bytes (canonical global-index layout).
    subroutine read_vnl_kappa_data()
        use util_ssbe, only: split_range
        use mpiio_export_ssbe, only: mpiio_open_read, mpiio_close, mpiio_read_at_all_z, &
                                   & mpiio_rd_c8, mpiio_rd_i, mpiio_rd_d, mpiio_disp_k, &
                                   & mpiio_csum_local_z, mpiio_csum_reduce, mpiio_finite_z, &
                                   & HDR_VNL, MAGIC_VNL, MAGIC_CHK
        implicit none
        ! OFF == MPI_OFFSET_KIND, taken from the public HDR_VNL parameter so the
        ! header-offset literals get the right kind in BOTH the USE_MPI and the
        ! serial-fallback builds WITHOUT this module needing `use mpi`.
        integer, parameter :: OFF = kind(HDR_VNL)
        integer(OFF) :: disp
        integer :: fh, ios, ik, s, i, iw, jw, nlen
        integer :: f_ndim, f_nk, f_nb, f_ns, f_kgrid(3), f_nelec, f_natom
        integer :: ihdr4(4), ihdr2(2)
        real(8) :: dhdr5(5), avec9(9)
        real(8) :: f_h, f_dir(3), f_amax, f_avec(3,3)
        character(8) :: magic
        integer :: ierr, ierr_h, ik_lo, ik_hi, max_tile
        integer(8) :: fsize, expect
        integer(8) :: nz_k, nz, csum_loc, csum_glob, csum_expect
        real(8) :: wmax, dmax, denom, smax_loc
        complex(8) :: w0
        logical :: file_exists_vnl, has_sidecar, has_slice
        integer :: istale_loc, istale_any, ifin_loc, ifin_any, ik_worst
        integer :: chkfh, s_nk, s_nb, s_ns, s_ver
        character(256) :: file_vnl_chk
        character(8) :: cmagic
        ! Task-4 numeric round-trip gate: flag-gated (env SBE_RVNL_DUMP=1) debug dump
        character(8) :: dbgval
        integer :: dbglen, dbgfh, npair, idir_d, jk_d, ip_d
        integer :: idump_ik(3), ib_list(3), jb_list(3)
        integer, allocatable :: itbl_min(:), itbl_max(:)
        complex(8), allocatable :: tmp_slice(:, :, :, :, :)  ! (f_nb,f_nb,0:3,-ns:ns, ik_lo:ik_hi)
        complex(8), allocatable :: rvnl_full_l(:, :, :, :)   ! (nb_eff,nb_eff,3,nk) zero-fill reduce buffer
        complex(8) :: dummy_z(1)

        ierr = 0
        f_nb = 0;  f_ns = 0;  f_h = 0d0;  f_dir = 0d0
        has_sidecar = .false.;  csum_expect = 0_8

        ! ---- pre-open guards (rank 0): empty name / existence / size ----
        if (irank == 0) then
            write(*, '(a)') "# read_vnl_kappa_data"
            if (len_trim(file_sbe_vnl_kappa) == 0) then
                write(*, '(a)') "ERROR(read_vnl_kappa_data): 'file_sbe_vnl_kappa' is empty while 'yn_sbe_vnl_exact'='y'."
                ierr = 1
            end if
            if (ierr == 0) then
                inquire(file=trim(file_sbe_vnl_kappa), exist=file_exists_vnl, size=fsize)
                if (.not. file_exists_vnl) then
                    write(*, '(a)') "ERROR(read_vnl_kappa_data): file not found: " // trim(file_sbe_vnl_kappa)
                    ierr = 1
                end if
            end if
        end if
        ! collective consistency: all ranks decide together whether to open.
        call comm_bcast(ierr, icomm, 0)
        if (ierr /= 0) stop 1

        ! ---- collective open of the single shared file (all ranks of icomm) ----
        call mpiio_open_read(trim(file_sbe_vnl_kappa), icomm, fh, ierr)

        ! ---- rank-0 header read + fail-closed validation (156B, byte-identical
        !      to the legacy sequential-stream layout) ----
        if (irank == 0) then
            ierr_h = 0
            call mpiio_rd_c8(fh,  0_OFF, magic,   ierr);    ierr_h = ior(ierr_h, ierr)  ! 'SBEVNLK1'
            call mpiio_rd_i (fh,  8_OFF, ihdr4, 4, ierr);   ierr_h = ior(ierr_h, ierr)  ! ndim, nk, nb, ns
            call mpiio_rd_i (fh, 24_OFF, f_kgrid, 3, ierr); ierr_h = ior(ierr_h, ierr)  ! num_kgrid(1:3)
            call mpiio_rd_i (fh, 36_OFF, ihdr2, 2, ierr);   ierr_h = ior(ierr_h, ierr)  ! nelec, natom
            call mpiio_rd_d (fh, 44_OFF, dhdr5, 5, ierr);   ierr_h = ior(ierr_h, ierr)  ! h, edir(1:3), amax
            call mpiio_rd_d (fh, 84_OFF, avec9, 9, ierr);   ierr_h = ior(ierr_h, ierr)  ! primitive_a (col order)
            f_ndim  = ihdr4(1);  f_nk = ihdr4(2);  f_nb = ihdr4(3);  f_ns = ihdr4(4)
            f_nelec = ihdr2(1);  f_natom = ihdr2(2)
            f_h = dhdr5(1);  f_dir(1:3) = dhdr5(2:4);  f_amax = dhdr5(5)
            f_avec = reshape(avec9, [3, 3])

            ierr = 0
            ! fail-closed: an MPI_File_read_at error on any header field must not
            ! be silently discarded before the field-value validation below.
            if (ierr_h /= 0) then
                write(*, '(a)') "ERROR(read_vnl_kappa_data): header field read failed (MPI_File_read_at)."
                ierr = 1
            end if
            if (magic /= MAGIC_VNL) then
                write(*, '(a)') "ERROR(read_vnl_kappa_data): bad magic (not a vnl_kappa export)."
                ierr = 1
            end if
            if (f_ndim /= 1) then
                write(*, '(a,i0)') "ERROR(read_vnl_kappa_data): unsupported stencil ndim = ", f_ndim
                ierr = 1
            end if
            if (f_nk /= nk) then
                write(*, '(a)') "ERROR(read_vnl_kappa_data): nk in file differs from SBE nk."
                ierr = 1
            end if
            if (any(f_kgrid(1:3) /= num_kgrid(1:3))) then
                write(*, '(a)') "ERROR(read_vnl_kappa_data): num_kgrid in file differs from input."
                ierr = 1
            end if
            if (f_nb < nb_hi) then
                write(*, '(a)') "ERROR(read_vnl_kappa_data): band count in file is smaller than the SBE window top."
                ierr = 1
            end if
            if (f_ns < 4) then
                write(*, '(a)') "ERROR(read_vnl_kappa_data): ns in file must be >= 4."
                ierr = 1
            end if
            if (f_h <= 0d0) then
                write(*, '(a)') "ERROR(read_vnl_kappa_data): stencil spacing h must be > 0."
                ierr = 1
            end if
            if (abs(norm2(f_dir(1:3)) - 1d0) > 1d-8) then
                write(*, '(a)') "ERROR(read_vnl_kappa_data): stencil direction in file is not a unit vector."
                ierr = 1
            end if
            if (f_nelec /= ne) then
                ! ne here is the ABSOLUTE electron count passed by the caller
                ! (window reduction happens after this argument), so it must
                ! match the GS export fingerprint exactly.
                write(*, '(a,i0,a,i0)') "ERROR(read_vnl_kappa_data): nelec fingerprint mismatch: file ", &
                    & f_nelec, " vs input ", ne
                ierr = 1
            end if
            if (natom > 0 .and. f_natom /= natom) then
                write(*, '(a)') "ERROR(read_vnl_kappa_data): natom fingerprint mismatch."
                ierr = 1
            end if
            if (maxval(abs(f_avec(1:3,1) - gs%a_matrix(1:3,1))) > 1d-8 * maxval(abs(gs%a_matrix)) .or. &
              & maxval(abs(f_avec(1:3,2) - gs%a_matrix(1:3,2))) > 1d-8 * maxval(abs(gs%a_matrix)) .or. &
              & maxval(abs(f_avec(1:3,3) - gs%a_matrix(1:3,3))) > 1d-8 * maxval(abs(gs%a_matrix))) then
                write(*, '(a)') "ERROR(read_vnl_kappa_data): lattice-vector fingerprint mismatch."
                ierr = 1
            end if
            if (ierr == 0) then
                ! exact size: header (156 B) + nk*(2ns+1)*4*nb^2 complex(8).  The
                ! MPI_Get_count short-read guards truncation of the data region;
                ! this catches surplus data / gross corruption too.
                expect = 156_8 + int(f_nk,8) * int(2*f_ns+1,8) * 4_8 * int(f_nb,8)**2 * 16_8
                if (fsize /= expect) then
                    write(*, '(a,i0,a,i0,a)') "ERROR(read_vnl_kappa_data): file size ", fsize, &
                        & " differs from expected ", expect, " (truncated or surplus data)."
                    ierr = 1
                end if
            end if
        end if

        call comm_bcast(ierr, icomm, 0)
        if (ierr /= 0) stop 1
        call comm_bcast(f_nb, icomm, 0)
        call comm_bcast(f_ns, icomm, 0)
        call comm_bcast(f_h,  icomm, 0)
        call comm_bcast(f_dir, icomm, 0)

        ! ---- optional sidecar checksum (rank 0; legacy byte-identical files
        !      lack it -> fall back to count+finite legacy validation) ----
        if (irank == 0) then
            nlen = len_trim(file_sbe_vnl_kappa)
            if (nlen >= 4 .and. file_sbe_vnl_kappa(nlen-3:nlen) == '.bin') then
                file_vnl_chk = file_sbe_vnl_kappa(1:nlen-4) // '.chk'
            else
                file_vnl_chk = trim(file_sbe_vnl_kappa) // '.chk'
            end if
            inquire(file=trim(file_vnl_chk), exist=has_sidecar)
            if (has_sidecar) then
                chkfh = get_filehandle()
                open(chkfh, file=trim(file_vnl_chk), status='old', action='read', iostat=ios)
                if (ios /= 0) then
                    ! present-but-unreadable sidecar is suspicious -> fail closed
                    ! (only an ABSENT sidecar triggers the legacy fallback above).
                    write(*, '(a)') "ERROR(read_vnl_kappa_data): sidecar .chk exists but could not be opened."
                    ierr = 1
                else
                    read(chkfh, '(a)', iostat=ios) cmagic
                    if (ios == 0) read(chkfh, *, iostat=ios) s_nk, s_nb, s_ns, s_ver
                    ! csum is written as an integer(8) signed decimal by the writer
                    ! (write.f90: write(chkfh,'(i0)') csum) -> parse as integer(8).
                    if (ios == 0) read(chkfh, *, iostat=ios) csum_expect
                    close(chkfh)
                    if (ios /= 0 .or. cmagic /= MAGIC_CHK) then
                        write(*, '(a)') "ERROR(read_vnl_kappa_data): sidecar .chk is corrupt (magic/format)."
                        ierr = 1
                    else if (s_nk /= f_nk .or. s_nb /= f_nb .or. s_ns /= f_ns) then
                        write(*, '(a)') "ERROR(read_vnl_kappa_data): sidecar .chk dims disagree with the .bin header."
                        ierr = 1
                    end if
                end if
            end if
        end if
        call comm_bcast(ierr, icomm, 0)
        if (ierr /= 0) stop 1
        call comm_bcast(has_sidecar, icomm, 0)

        ! ---- k-slice partition (deterministic, identical to the solver's split;
        !      asserted against sbe%ik_min/ik_max in init_sbe_bloch_solver) ----
        allocate(itbl_min(0:nproc-1), itbl_max(0:nproc-1))
        call split_range(1, nk, nproc, itbl_min, itbl_max)
        ik_lo = itbl_min(irank);  ik_hi = itbl_max(irank)
        deallocate(itbl_min, itbl_max)
        has_slice = (ik_hi >= ik_lo)   ! .false. only when nproc > nk (empty rank)

        allocate(gs%vnl_V(1:nb_eff, 1:nb_eff, -f_ns:f_ns, ik_lo:ik_hi))
        allocate(gs%vnl_W(1:nb_eff, 1:nb_eff, 1:3, -f_ns:f_ns, ik_lo:ik_hi))

        ! ---- collective read of this rank's contiguous k-slice.  tmp_slice is
        !      allocated EXACTLY (f_nb,f_nb,0:3,-ns:ns, ik_lo:ik_hi): dim1=row
        !      (fastest) .. dim5=ik (slowest), so its Fortran array-element order
        !      IS the canonical on-disk byte order.  Empty ranks issue a count-0
        !      collective no-op so MPI_File_read_at_all stays well-formed. ----
        nz_k = int(f_nb,8) * int(f_nb,8) * 4_8 * (2_8*int(f_ns,8) + 1_8)
        ! fail-closed: mpiio_read_at_all_z passes a default-integer element count
        ! (int(nz)); guard the largest k-slice (= ceil(nk/nproc)) against 2^31-1
        ! wrap.  Deterministic on every rank -> collective-safe (all stop together,
        ! not one rank alone).  Realistic band/k counts never hit this; hard stop.
        max_tile = (nk + nproc - 1) / nproc
        if (nz_k * int(max_tile,8) > int(huge(0),8)) then
            if (irank == 0) write(*, '(a)') &
                & "ERROR(read_vnl_kappa_data): per-slice element count exceeds MPI int count (chunked read needed)."
            stop 1
        end if
        if (has_slice) then
            allocate(tmp_slice(1:f_nb, 1:f_nb, 0:3, -f_ns:f_ns, ik_lo:ik_hi))
            disp = mpiio_disp_k(HDR_VNL, nz_k, ik_lo)
            nz   = nz_k * int(ik_hi - ik_lo + 1, 8)
            call mpiio_read_at_all_z(fh, disp, tmp_slice(1,1,0,-f_ns,ik_lo), nz, icomm, 'vnl_kappa', ierr)
        else
            allocate(tmp_slice(1:f_nb, 1:f_nb, 0:3, -f_ns:f_ns, 1:1))  ! placeholder; never read as data
            nz = 0_8
            call mpiio_read_at_all_z(fh, HDR_VNL, dummy_z, 0_8, icomm, 'vnl_kappa', ierr)  ! collective no-op
        end if

        ! ---- validation layer (3): finiteness, always, collective (fail-closed) ----
        ifin_loc = 0
        if (has_slice) then
            if (.not. mpiio_finite_z(tmp_slice(1,1,0,-f_ns,ik_lo), nz)) ifin_loc = 1
        end if
        call comm_summation(ifin_loc, ifin_any, icomm)   ! LOR via SUM of 0/1
        if (ifin_any > 0) then
            if (irank == 0) write(*, '(a)') &
                & "ERROR(read_vnl_kappa_data): non-finite element (NaN/Inf) in vnl_kappa data."
            stop 1   ! all ranks reach here (ifin_any is allreduced) -> abort together
        end if

        ! ---- validation layer (4): whole-file XOR checksum, ONLY when the
        !      sidecar is present.  The writer checksummed the RAW on-disk complex
        !      data (buf_l, canonical order), so we must checksum the RAW read
        !      buffer tmp_slice with the same global-index offset (NOT the
        !      Hermitized W reconstruction).  XOR-fold => the M-rank split of the
        !      same nk gives the same global checksum as the N-rank writer. ----
        if (has_sidecar) then
            csum_loc = 0_8
            if (has_slice) csum_loc = mpiio_csum_local_z(tmp_slice(1,1,0,-f_ns,ik_lo), nz, nz_k*int(ik_lo-1,8))
            csum_glob = mpiio_csum_reduce(csum_loc, icomm)   ! collective (BXOR), value on all ranks
            if (irank == 0 .and. csum_glob /= csum_expect) then
                write(*, '(a)') "ERROR(read_vnl_kappa_data): sidecar checksum mismatch (data corruption / stale file)."
                ierr = 1
            end if
            call comm_bcast(ierr, icomm, 0)
            if (ierr /= 0) stop 1
        end if

        ! ---- (i) full-nk W_{i,0} reconstruct via zero-fill + comm_summation.
        !      Each rank fills ONLY its owned k-slice (disjoint slices tile
        !      1..nk); the all-reduce sums them into the full-nk matrix on every
        !      rank.  BEFORE the overwrite, cross-check this rank's slice against
        !      the pre-overwrite text-file rvnl_tm_matrix and LOR-reduce the fail
        !      flag so a mismatch on ANY rank aborts ALL. ----
        allocate(rvnl_full_l(1:nb_eff, 1:nb_eff, 1:3, 1:nk))
        rvnl_full_l = (0d0, 0d0)
        istale_loc = 0;  smax_loc = 0d0;  ik_worst = 0
        do ik = ik_lo, ik_hi
            dmax = 0d0
            do i = 1, 3
                do jw = nb_min, nb_hi
                    do iw = nb_min, nb_hi
                        w0 = 0.5d0 * (tmp_slice(iw, jw, i, 0, ik) + conjg(tmp_slice(jw, iw, i, 0, ik)))
                        rvnl_full_l(iw - nb_min + 1, jw - nb_min + 1, i, ik) = w0
                        dmax = max(dmax, abs(w0 - gs%rvnl_tm_matrix(iw - nb_min + 1, jw - nb_min + 1, i, ik)))
                    end do
                end do
            end do
            wmax  = maxval(abs(gs%rvnl_tm_matrix(:, :, :, ik)))
            denom = max(wmax, 1d-12)
            if (dmax / denom > 1d-3) then
                istale_loc = 1
                if (dmax / denom > smax_loc) then
                    smax_loc = dmax / denom;  ik_worst = ik
                end if
            end if
        end do
        call comm_summation(istale_loc, istale_any, icomm)   ! LOR via SUM of 0/1
        if (istale_any > 0) then
            if (istale_loc == 1) then
                write(*, '(a,i0,a,es12.4)') "ERROR(read_vnl_kappa_data): W(s=0) vs tm rvnl mismatch at ik=", &
                    & ik_worst, ", rel-to-max residual = ", smax_loc
                write(*, '(a)') "  (stale/mismatched 'file_sbe_vnl_kappa' for this GS export?)"
            end if
            stop 1   ! all ranks reach here (istale_any is allreduced) -> abort together
        end if
        call comm_summation(rvnl_full_l, gs%rvnl_tm_matrix, 3*nb_eff*nb_eff*nk, icomm)
        deallocate(rvnl_full_l)

        ! ---- Task-4 numeric round-trip gate: optional debug dump of the
        !      reconstructed full-nk rvnl_tm_matrix at a FIXED small element set
        !      (full precision es24.16).  OFF unless env SBE_RVNL_DUMP=1
        !      (production-safe default); rank-0 only, deterministic indices.  This
        !      is what the M4/M8/HEAD dumps are diffed on, so it exercises the
        !      offset + Hermitization + comm_summation reconstruct -- which the
        !      raw-tmp_slice sidecar checksum does NOT cover. ----
        if (irank == 0) then
            call get_environment_variable('SBE_RVNL_DUMP', dbgval, dbglen)
            if (dbglen > 0 .and. trim(dbgval) == '1') then
                idump_ik(1) = 1;  idump_ik(2) = (nk + 1) / 2;  idump_ik(3) = nk
                npair = 1;  ib_list(1) = 1;  jb_list(1) = 1
                if (nb_eff >= 2) then
                    ib_list(2) = 2;  jb_list(2) = 1
                    ib_list(3) = 1;  jb_list(3) = 2;  npair = 3
                end if
                dbgfh = get_filehandle()
                open(dbgfh, file='rvnl_dump.txt', status='replace', action='write')
                write(dbgfh, '(a)')         "# rvnl_tm_matrix reconstruct dump: ib jb idir ik  Re Im"
                write(dbgfh, '(a,i0,a,i0)') "# nb_eff=", nb_eff, " nk=", nk
                do idir_d = 1, 3
                    do jk_d = 1, 3
                        if (idump_ik(jk_d) < 1 .or. idump_ik(jk_d) > nk) cycle
                        do ip_d = 1, npair
                            write(dbgfh, '(4(i6,1x),2es24.16)') ib_list(ip_d), jb_list(ip_d), idir_d, idump_ik(jk_d), &
                                & real (gs%rvnl_tm_matrix(ib_list(ip_d), jb_list(ip_d), idir_d, idump_ik(jk_d))), &
                                & aimag(gs%rvnl_tm_matrix(ib_list(ip_d), jb_list(ip_d), idir_d, idump_ik(jk_d)))
                        end do
                    end do
                end do
                close(dbgfh)
                write(*, '(a)') "#   vnl_kappa: SBE_RVNL_DUMP=1 -> wrote rvnl_dump.txt (Task-4 numeric gate)."
            end if
        end if

        ! ---- (ii) windowed, Hermitized V/W stencil for this rank's k-slice ----
        do ik = ik_lo, ik_hi
            do s = -f_ns, f_ns
                do jw = nb_min, nb_hi
                    do iw = nb_min, nb_hi
                        gs%vnl_V(iw - nb_min + 1, jw - nb_min + 1, s, ik) = &
                            & 0.5d0 * (tmp_slice(iw, jw, 0, s, ik) + conjg(tmp_slice(jw, iw, 0, s, ik)))
                        do i = 1, 3
                            gs%vnl_W(iw - nb_min + 1, jw - nb_min + 1, i, s, ik) = &
                                & 0.5d0 * (tmp_slice(iw, jw, i, s, ik) + conjg(tmp_slice(jw, iw, i, s, ik)))
                        end do
                    end do
                end do
            end do
        end do

        if (allocated(tmp_slice)) deallocate(tmp_slice)
        call mpiio_close(fh, ierr)   ! collective MPI_File_close (all ranks)

        ! ---- header params -> gs (already bcast to all ranks above) ----
        gs%vnl_exact  = .true.
        gs%vnl_ns     = f_ns
        gs%vnl_h      = f_h
        gs%vnl_dir(1:3) = f_dir(1:3)
        gs%vnl_ik_min = ik_lo
        gs%vnl_ik_max = ik_hi

        if (irank == 0) then
            write(*, '(a,i0,a,es12.4,a,3f10.6)') "#   vnl_kappa: ns = ", f_ns, &
                & ", h = ", f_h, ", dir = ", f_dir(1:3)
            if (has_sidecar) then
                write(*, '(a)') "#   vnl_kappa: sidecar checksum verified."
            else
                write(*, '(a)') "#   vnl_kappa: no sidecar (.chk) -> legacy count+finite validation only."
            end if
        end if
    end subroutine read_vnl_kappa_data


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
                do ib=1, nb_eff
                    do jb=1, nb_eff
                        gs%delta_omega(ib, jb, ik) = gs%eigen(ib, ik) - gs%eigen(jb, ik)
                    end do
                end do
            end do

            if (.not. allocated(gs%xi))    allocate(gs%xi(1:nb_eff, 1:nb_eff, 1:3, 1:nk))
            if (.not. allocated(gs%xi_ok)) allocate(gs%xi_ok(1:nb_eff, 1:nb_eff, 1:nk))
            call build_xi(nb_eff, nk, gs%nbvec, gs%bvec, gs%prod_dk, gs%eigen, &
                        & gs%b_matrix, num_kgrid, gs%xi, gs%xi_ok, nrej, resu, &
                        & fixed_blocks=(trim(sbe_lg_degen) == 'gifix'))

            resp_max = 0d0
            do ik=1, nk
                do ib=1, nb_eff
                    do jb=1, nb_eff
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
        else if (uses_xfull_links(sbe_lg_degen)) then
            ! ===== X-full: single full-band block -> full-M Wilson transport =====
            ! block_id≡1 (one block spanning all nb bands) makes
            ! build_block_transport polar-factor the WHOLE nb×nb overlap M =
            ! the full-band Wilson-line transport. The full-band
            ! covariant_grad_block then supplies the WHOLE field term
            ! E·(∂_k ρ − i[ξ,ρ]) including the interband dipole, so d_matrix
            ! is unused by gicov_rhs (dropped there). No fixed blocks, no
            ! closure, no occupation guard.
            do ik=1, nk
                do ib=1, nb_eff
                    do jb=1, nb_eff
                        gs%delta_omega(ib, jb, ik) = gs%eigen(ib, ik) - gs%eigen(jb, ik)
                    end do
                end do
            end do

            ! X-full: ONE full-band block (block_id≡1) -> build_block_transport
            ! polar-factors the whole nb×nb overlap M = the full-band Wilson
            ! transport. No fixed blocks, no closure.
            if (.not. allocated(gs%block_id)) allocate(gs%block_id(1:nb_eff, 1:nk))
            gs%block_id(:, :) = 1                                   ! single full-band block
            if (.not. allocated(gs%u_transport)) allocate(gs%u_transport(1:nb_eff, 1:nb_eff, 1:3, 1:nk))
            call build_block_transport(nb_eff, nk, gs%nbvec, gs%bvec, gs%prod_dk, &
                                     & gs%block_id, num_kgrid, gs%u_transport, nrej)

            ! X-full: the full-band covariant transport supplies the WHOLE field
            ! term incl. the interband dipole (xi_inter = dipole), so d_matrix is
            ! unused by gicov_rhs and left unallocated.

            if (irank == 0) then
                write(*, '(a,i0)') "# build_block_transport: rejected blocks = ", nrej
            end if
        else
            ! ===== default 'off': bit-identical to the pre-Pb3 dipole construction =====
            do ik=1, nk
                do ib=1, nb_eff
                    do jb=1, nb_eff
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


! T2 Delta-omega dephasing gate weight (design: gw/gw_design/plans/
! 2026-07-08-t2-gate-shape.md).  MODULE-LEVEL (not nested in
! init_sbe_gs_info/prepare_matrix above, even though it lives next to the
! delta_omega construction those build) so bloch_solver_ssbe.f90's
! init_sbe_bloch_solver can precompute sbe%t2_gate_w from it, and so it is
! directly unit-testable.  Returns the real weight W in [0,1] multiplying the
! phenomenological -rho/t_2 dephasing term; the caller (gicov_rhs) applies
! 1/t_2 and the off-diagonal (n/=m) restriction itself.
!   W = 0                                     if |delta_omega| <= floor (exact-degeneracy clamp)
!   step : W = 1 if |delta_omega| >  theta  else 0     (strict >, matches the current hard gate)
!   gauss: W = 1 - exp( -(delta_omega/width)^2 )        (Thuemmler Eq.28; W(0)=0, W'(0)=0)
pure real(8) function t2_gate_weight(delta_omega, shape, theta, width, floor) result(w)
    implicit none
    real(8),      intent(in) :: delta_omega, theta, width, floor
    character(*), intent(in) :: shape
    real(8) :: adw
    adw = abs(delta_omega)
    if (adw <= floor) then
        w = 0d0
    else if (trim(shape) == 'gauss') then
        w = 1d0 - exp( -(delta_omega / width)**2 )
    else                                          ! 'step' (default)
        w = merge(1d0, 0d0, adw > theta)
    end if
end function t2_gate_weight

end module gs_info_ssbe

