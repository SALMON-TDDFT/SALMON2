!=======================================================================
  ! Calculate Hamiltonian matrix in fragment basis
  !=======================================================================
  logical function is_momentum_neighbor_axis(lg, s1, n1, s2, n2) result(ok)
    implicit none
    integer, intent(in) :: lg, s1, n1, s2, n2
    integer :: e1, e2, s1_next, s2_next

    e1 = s1 + n1 - 1
    e2 = s2 + n2 - 1
    s1_next = modulo(e1, lg) + 1
    s2_next = modulo(e2, lg) + 1
    ok = ((s1 == s2) .and. (n1 == n2)) .or. (s1 == s2_next) .or. (s2 == s1_next)
  end function is_momentum_neighbor_axis

  logical function is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col) result(is_pair)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag_row, ifrag_col
    integer :: axis
    logical :: axis_ok(3)

    is_pair = .false.
    if (ifrag_row == ifrag_col) then
      is_pair = .true.
      return
    end if

    do axis = 1, 3
      axis_ok(axis) = is_momentum_neighbor_axis(dg_frag%lgnum_total(axis), &
        dg_frag%ixyz_frag(axis, ifrag_row), dg_frag%nxyz_domain(axis, ifrag_row), &
        dg_frag%ixyz_frag(axis, ifrag_col), dg_frag%nxyz_domain(axis, ifrag_col))
    end do

    is_pair = all(axis_ok)
  end function is_momentum_neighbor_pair

  integer function find_momentum_block(dg_frag, ifrag_row, ifrag_col) result(iblk)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag_row, ifrag_col

    iblk = 0
    if (.not. allocated(dg_frag%momentum_block_map)) return
    if (ifrag_row < 1 .or. ifrag_row > size(dg_frag%momentum_block_map, 1)) return
    if (ifrag_col < 1 .or. ifrag_col > size(dg_frag%momentum_block_map, 2)) return
    iblk = dg_frag%momentum_block_map(ifrag_row, ifrag_col)
  end function find_momentum_block

  subroutine init_momentum_blocks(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer :: ifrag_row, ifrag_col, nblk, iblk
    integer :: nrow_max, ncol_max

    if (allocated(dg_frag%momentum_blocks)) then
      do iblk = 1, size(dg_frag%momentum_blocks)
        if (allocated(dg_frag%momentum_blocks(iblk)%val)) deallocate(dg_frag%momentum_blocks(iblk)%val)
      end do
      deallocate(dg_frag%momentum_blocks)
    end if
    if (allocated(dg_frag%momentum_block_map)) deallocate(dg_frag%momentum_block_map)

    nblk = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col)) nblk = nblk + 1
      end do
    end do

    dg_frag%n_momentum_blocks = nblk
    if (nblk <= 0) return
    allocate(dg_frag%momentum_blocks(nblk))
    allocate(dg_frag%momentum_block_map(dg_frag%n_frag, dg_frag%n_frag))
    dg_frag%momentum_block_map = 0

    iblk = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (.not. is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col)) cycle
        iblk = iblk + 1
        dg_frag%momentum_block_map(ifrag_row, ifrag_col) = iblk
        dg_frag%momentum_blocks(iblk)%ifrag_row = ifrag_row
        dg_frag%momentum_blocks(iblk)%ifrag_col = ifrag_col
        nrow_max = max(1, maxval(dg_frag%n_basis(ifrag_row, 1:dg_frag%nspin)))
        ncol_max = max(1, maxval(dg_frag%n_basis(ifrag_col, 1:dg_frag%nspin)))
        dg_frag%momentum_blocks(iblk)%nrow_max = nrow_max
        dg_frag%momentum_blocks(iblk)%ncol_max = ncol_max
        allocate(dg_frag%momentum_blocks(iblk)%val(3, nrow_max, ncol_max, dg_frag%nspin))
        dg_frag%momentum_blocks(iblk)%val = 0.0d0
      end do
    end do
  end subroutine init_momentum_blocks

  !=======================================================================
  ! Calculate initial Hamiltonian matrix from basis functions
  !
  ! Includes halo (ghost cell) exchange for accurate boundary treatment.
  ! System boundaries use PERIODIC boundary conditions (full system is periodic).
  ! Fragment boundaries are handled via MPI communication between neighboring fragments.
  ! The real-space fragment basis itself is shared across spin channels in the
  ! present non-SOI DG path; the nspin axis here labels spin-resolved projected
  ! matrices and basis indexing, not separate copies of phi_frag for each spin.
  !=======================================================================
  subroutine calculate_hamiltonian_matrix(dg_frag, system, lg, mg, stencil, &
                                         Vh, Vxc, Vpsl, pp, ppg)
    use structures
    use communication, only: comm_is_root, comm_summation
    use parallelization, only: nproc_size_global
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: lg, mg
    type(s_stencil),        intent(in)    :: stencil
    type(s_scalar),         intent(in)    :: Vh, Vxc(:), Vpsl
    type(s_pp_info),        intent(in)    :: pp
    type(s_pp_grid),        intent(in)    :: ppg
    
    integer :: ifrag, ispin, io, jo, i_local, nbf, ig_i, ig_j
    integer :: ndom(3)
    real(8) :: hvol
    real(8) :: max_p
    real(8) :: Ac_zero(3)
    integer :: is(3), ie(3)
    real(8), allocatable :: T_phi(:,:,:)  ! Kinetic energy operator applied to basis (fragment-local)
    real(8), allocatable :: H_phi(:,:,:)  ! Hamiltonian-applied field H|phi_j> = T|phi_j> + V|phi_j> (fragment-local)
    real(8), allocatable :: V_total(:,:,:)  ! Total potential V = Vpsl + Vh + Vxc
    real(8), allocatable :: partial_t(:), partial_h(:), reduced_t(:), reduced_h(:)
    
    if (.not. dg_frag%has_real_space_basis) then
      if (.not. allocated(dg_frag%H_mat)) then
        allocate(dg_frag%H_mat(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
      end if
      dg_frag%H_mat = 0.0d0
      return
    end if
    
    if (comm_is_root(dg_frag%id)) then
      write(*,*)
      write(*,*) "=== Preparing Hamiltonian Matrix ==="
    end if
    
    ! Step 1: Calculate momentum matrix elements (transition moments)
    ! Required for velocity gauge A·p coupling
    if (.not. allocated(dg_frag%momentum_blocks) .and. .not. allocated(dg_frag%momentum_mat)) then
      if (comm_is_root(dg_frag%id)) then
        write(*,*) "  [1/3] Calculating momentum matrix elements (p_ij)..."
        write(*,*) "        Using 4th-order finite difference stencil"
      end if
      call calculate_momentum_matrix(dg_frag, system, mg, stencil)
      call calculate_overlap_matrix(dg_frag, system, mg)
      write(*,'(1x,a,i0,a,i0,a,i0,a,a)') "        hamiltonian stage: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " stage=", "after-overlap-return"
      flush(6)
      if (comm_is_root(dg_frag%id)) then
        write(*,*) "        Momentum matrix calculated (for A·p coupling)"
        write(*,*) "        Overlap matrix S calculated (for generalized propagation)"
      end if
    else
      if (comm_is_root(dg_frag%id)) then
        write(*,*) "  [1/3] Momentum matrix already available"
      end if
      if (.not. allocated(dg_frag%S_mat)) then
        call calculate_overlap_matrix(dg_frag, system, mg)
        write(*,'(1x,a,i0,a,i0,a,i0,a,a)') "        hamiltonian stage: rank=", dg_frag%id, &
          " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " stage=", "after-overlap-return"
        flush(6)
      end if
    end if
    
    ! Step 2: Allocate Hamiltonian matrix
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "  [2/3] Constructing Hamiltonian matrix H = T + V..."
    end if
    
    ! Allocate only when needed to keep the momentum-matrix peak lower.
    if (.not. allocated(dg_frag%H_mat)) then
      allocate(dg_frag%H_mat(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
    end if
    if (.not. allocated(dg_frag%H_mat_kinetic)) then
      allocate(dg_frag%H_mat_kinetic(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
    end if
    dg_frag%H_mat = 0.0d0
    dg_frag%H_mat_kinetic = 0.0d0
    write(*,'(1x,a,i0,a,i0,a,i0,a,a)') "        hamiltonian stage: rank=", dg_frag%id, &
      " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " stage=", "after-hmat-alloc"
    flush(6)
    
    ! Exchange halo regions between fragments before stencil operations
    ! This ensures accurate Laplacian calculation at fragment boundaries
    call exchange_phi_frag_halo(dg_frag)
    write(*,'(1x,a,i0,a,i0,a,i0,a,a)') "        hamiltonian stage: rank=", dg_frag%id, &
      " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " stage=", "after-step2-halo"
    flush(6)
    
    hvol = system%hvol
    is = lg%is
    ie = lg%ie
    
    allocate(V_total(is(1):ie(1), is(2):ie(2), is(3):ie(3)))
    write(*,'(1x,a,i0,a,i0,a,i0,a,a)') "        hamiltonian stage: rank=", dg_frag%id, &
      " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " stage=", "after-vtotal-alloc"
    flush(6)
    
    ! Construct total potential: V = Vpsl + Vh + Vxc
    ! Note: This is used for initial H_mat calculation
    do ispin = 1, system%nspin
      call build_total_potential_grid(lg, Vh, Vxc(ispin), Vpsl, V_total)
      if (ispin == 1) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,a)') "        hamiltonian stage: rank=", dg_frag%id, &
          " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " stage=", "after-build-total-potential"
        flush(6)
      end if
      
      ! Loop over fragments assigned to this rank
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        ndom(:) = dg_frag%nxyz_domain(:, ifrag)
        allocate(T_phi(1:ndom(1), 1:ndom(2), 1:ndom(3)))
        allocate(H_phi(1:ndom(1), 1:ndom(2), 1:ndom(3)))
        
        ! Calculate Hamiltonian matrix elements for this fragment
        ! H_ij = <φ_i | T + V | φ_j> = T_ij + V_ij
        nbf = dg_frag%n_basis(ifrag, ispin)
        allocate(partial_t(nbf), partial_h(nbf), reduced_t(nbf), reduced_h(nbf))
        do jo = 1, nbf
          ig_j = dg_frag%index_basis(jo, ifrag, ispin)
          if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle

          call build_hpsi_for_basis(dg_frag, ifrag, i_local, jo, mg, stencil, V_total, T_phi, H_phi)

          ! Calculate matrix elements with all φ_i
          partial_t(:) = 0.0d0
          partial_h(:) = 0.0d0
          !$omp parallel do private(io, ig_i)
          do io = 1, nbf
            ig_i = dg_frag%index_basis(io, ifrag, ispin)
            if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle

            ! Kinetic energy matrix element: T_ij = ∫ φ_i (T|φ_j>) dr
            call integrate_basis_with_field(dg_frag, ifrag, i_local, io, mg, T_phi, hvol, partial_t(io))

            ! Store kinetic part
            call integrate_basis_with_field(dg_frag, ifrag, i_local, io, mg, H_phi, hvol, partial_h(io))

          end do
          !$omp end parallel do

          call comm_summation(partial_t, reduced_t, nbf, dg_frag%icomm_frag)
          call comm_summation(partial_h, reduced_h, nbf, dg_frag%icomm_frag)
          if (dg_frag%is_frag_root) then
            do io = 1, nbf
              ig_i = dg_frag%index_basis(io, ifrag, ispin)
              if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
              dg_frag%H_mat_kinetic(ig_i, ig_j, ispin) = reduced_t(io)
              dg_frag%H_mat(ig_i, ig_j, ispin) = reduced_h(io)
            end do
          end if

        end do  ! jo loop
        deallocate(partial_t, partial_h, reduced_t, reduced_h)
        deallocate(T_phi, H_phi)
          
        
      end do  ! ifrag loop
      
    end do  ! ispin loop
    
    ! CRITICAL: MPI aggregation of Hamiltonian matrix
    ! Each rank computed elements only for its assigned fragments
    ! Sum across all ranks to get complete global H_mat
    block
      real(8), allocatable :: H_mat_flat(:), H_mat_tmp_flat(:)
      integer :: mat_size
      mat_size = dg_frag%n_mat_max * dg_frag%n_mat_max * dg_frag%nspin
      allocate(H_mat_flat(mat_size), H_mat_tmp_flat(mat_size))
      
      H_mat_flat = reshape(dg_frag%H_mat, [mat_size])
      call comm_summation(H_mat_flat, H_mat_tmp_flat, mat_size, dg_frag%icomm)
      dg_frag%H_mat = reshape(H_mat_tmp_flat, [dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin])
      
      H_mat_flat = reshape(dg_frag%H_mat_kinetic, [mat_size])
      call comm_summation(H_mat_flat, H_mat_tmp_flat, mat_size, dg_frag%icomm)
      dg_frag%H_mat_kinetic = reshape(H_mat_tmp_flat, [dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin])
      
      deallocate(H_mat_flat, H_mat_tmp_flat)
    end block

    ! Enforce Hermiticity for the static Hamiltonian parts used in RT propagation.
    do ispin = 1, system%nspin
      do jo = 1, dg_frag%n_mat_max
        do io = jo + 1, dg_frag%n_mat_max
          dg_frag%H_mat_kinetic(io, jo, ispin) = 0.5d0 * (dg_frag%H_mat_kinetic(io, jo, ispin) + dg_frag%H_mat_kinetic(jo, io, ispin))
          dg_frag%H_mat_kinetic(jo, io, ispin) = dg_frag%H_mat_kinetic(io, jo, ispin)

          dg_frag%H_mat(io, jo, ispin) = 0.5d0 * (dg_frag%H_mat(io, jo, ispin) + dg_frag%H_mat(jo, io, ispin))
          dg_frag%H_mat(jo, io, ispin) = dg_frag%H_mat(io, jo, ispin)
        end do
      end do
    end do

    if (comm_is_root(dg_frag%id)) then
      write(*,*) "        Kinetic and potential terms computed"
    end if
    
    ! Step 3: Non-local pseudopotential is handled in time evolution
    ! with vector potential A(t), so it is not added to H_mat here.
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "  [3/3] Non-local PP handled in time evolution (A-dependent)"
    end if

    ! Build initial mixed basis once (fragment + orthogonalized PW) with A=0.
    if (dg_frag%use_plane_wave_basis .and. dg_frag%n_plane_waves > 0) then
      Ac_zero(:) = 0.0d0
      if (comm_is_root(dg_frag%id)) then
        write(*,*) "  [init] Building mixed basis at startup (A=0)"
      end if
      call diagonalize_mixed_basis_pw(dg_frag, system, Vh, Vxc, Vpsl, Ac_zero)
      dg_frag%coef_new(:, :, :) = dg_frag%coef(:, :, :)
    end if

    ! Initialize field-free reference Hamiltonian for adaptive-basis metric.
    if (allocated(dg_frag%H_mat_old)) then
      do ispin = 1, min(dg_frag%nspin, size(dg_frag%H_mat_old,3))
        do jo = 1, min(dg_frag%nstate_frag, size(dg_frag%H_mat_old,2))
          do io = 1, min(dg_frag%nstate_frag, size(dg_frag%H_mat_old,1))
            dg_frag%H_mat_old(io, jo, ispin) = cmplx(dg_frag%H_mat(io, jo, ispin), 0.0d0, kind=8)
          end do
        end do
      end do
    end if
    
    deallocate(V_total)
    
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "=== Hamiltonian Matrix Ready ==="
      write(*,*)
    end if
    
  end subroutine calculate_hamiltonian_matrix

  !=======================================================================
  ! Build total local potential on the given grid:
  !   V_total = Vpsl + Vh + Vxc(ispin)
  !=======================================================================
  subroutine build_total_potential_grid(grid, Vh, Vxc_spin, Vpsl, V_total)
    use structures
    implicit none
    type(s_rgrid), intent(in) :: grid
    type(s_scalar), intent(in) :: Vh, Vxc_spin, Vpsl
    real(8), intent(out) :: V_total(grid%is(1):grid%ie(1), grid%is(2):grid%ie(2), grid%is(3):grid%ie(3))
    integer :: ix, iy, iz

    do iz = grid%is(3), grid%ie(3)
      do iy = grid%is(2), grid%ie(2)
        do ix = grid%is(1), grid%ie(1)
          V_total(ix, iy, iz) = Vpsl%f(ix, iy, iz) + Vh%f(ix, iy, iz) + Vxc_spin%f(ix, iy, iz)
        end do
      end do
    end do
  end subroutine build_total_potential_grid

  !=======================================================================
  ! Build T|phi_j> and H|phi_j>=T|phi_j>+V|phi_j> for one fragment/basis state
  !=======================================================================
  subroutine build_hpsi_for_basis(dg_frag, ifrag, i_local, jo, mg, stencil, V_total, T_phi, H_phi)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: ifrag, i_local, jo
    type(s_rgrid), intent(in) :: mg
    type(s_stencil), intent(in) :: stencil
    real(8), intent(in) :: V_total(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    real(8), intent(out) :: T_phi(:,:,:)
    real(8), intent(out) :: H_phi(:,:,:)
    integer :: lx, ly, lz, gx, gy, gz
    integer :: iorg(3), ndom(3)
    integer :: loc_s(3), loc_e(3)

    call apply_kinetic_to_basis(dg_frag, i_local, jo, mg, stencil, T_phi)
    H_phi(:, :, :) = T_phi(:, :, :)

    iorg(:) = dg_frag%ixyz_frag(:, ifrag)
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    call get_fragment_local_range(dg_frag, ndom, loc_s, loc_e)
    do lz = loc_s(3), loc_e(3)
      gz = iorg(3) + lz - 1
      do ly = loc_s(2), loc_e(2)
        gy = iorg(2) + ly - 1
        do lx = loc_s(1), loc_e(1)
          gx = iorg(1) + lx - 1
          H_phi(lx, ly, lz) = H_phi(lx, ly, lz) + V_total(gx, gy, gz) * dg_frag%phi_frag(lx, ly, lz, jo, i_local)
        end do
      end do
    end do
  end subroutine build_hpsi_for_basis

  !=======================================================================
  ! Integrate one bra basis function against a real-space field
  !   integral = <phi_io | field>
  !=======================================================================
  subroutine integrate_basis_with_field(dg_frag, ifrag, i_local, io, mg, field, hvol, integral)
    use structures
    use communication, only: comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, i_local, io
    type(s_rgrid), intent(in) :: mg
    real(8), intent(in) :: field(:,:,:)
    real(8), intent(in) :: hvol
    real(8), intent(out) :: integral
    real(8) :: partial
    integer :: lx, ly, lz
    integer :: ndom(3), loc_s(3), loc_e(3)

    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    call get_fragment_local_range(dg_frag, ndom, loc_s, loc_e)
    partial = 0.0d0
    do lz = loc_s(3), loc_e(3)
      do ly = loc_s(2), loc_e(2)
        !$omp simd reduction(+:partial)
        do lx = loc_s(1), loc_e(1)
          partial = partial + dg_frag%phi_frag(lx, ly, lz, io, i_local) * field(lx, ly, lz) * hvol
        end do
      end do
    end do
    integral = partial
  end subroutine integrate_basis_with_field
  
  !=======================================================================
  ! Apply kinetic energy operator to a single basis function
  !
  ! T|φ> = -∇²/2 |φ> = -0.5 * Laplacian(φ)
  !
  ! Uses 4th-order finite difference stencil (requires ±4 grid points).
  ! With halo exchange, computation is valid over entire domain (1:nx, 1:ny, 1:nz).
  !
  ! System boundaries: PERIODIC (full system is periodic)
  ! Fragment boundaries: Halo exchange provides neighbor data via MPI
  !=======================================================================
  subroutine apply_kinetic_to_basis(dg_frag, i_local, jo, mg, stencil, T_phi)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer,                intent(in) :: i_local, jo
    type(s_rgrid),          intent(in) :: mg
    type(s_stencil),        intent(in) :: stencil
    real(8),                intent(out) :: T_phi(:,:,:)
    
    integer :: lx, ly, lz, ifrag
    real(8) :: v, lap0
    real(8) :: lapt(4,3)
    integer :: ndom(3), loc_s(3), loc_e(3)
    
    ! Extract stencil coefficients
    lap0 = stencil%coef_lap0
    lapt = stencil%coef_lap
    ifrag = dg_frag%ifrag_start + i_local - 1
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    call get_fragment_local_range(dg_frag, ndom, loc_s, loc_e)
    
    ! Note: phi_frag is allocated as (1-nb:nx+nb, 1-nb:ny+nb, 1-nb:nz+nb, ...)
    ! where nb = nxyz_buffer = 4 for 4th-order stencil
    ! The interior domain is (1:nx, 1:ny, 1:nz), and halo provides data for stencil
    ! operations near boundaries.
    !
    ! With halo exchange, stencil operations can access phi_frag(ix±4, iy±4, iz±4)
    ! for all interior points without boundary restrictions.
    
    ! Apply kinetic operator using finite difference stencil
    ! With halo regions available, we can compute over FULL interior domain
    !
    ! Note: exchange_phi_frag_halo() must be called before this routine
    
    T_phi = 0.0d0
    
    do lz = loc_s(3), loc_e(3)
      do ly = loc_s(2), loc_e(2)
        do lx = loc_s(1), loc_e(1)
          ! Compute Laplacian using 4th-order finite difference
          ! Stencil accesses phi_frag(ix±4, iy±4, iz±4) which now includes halo
          v = lapt(1,1) * (dg_frag%phi_frag(lx+1, ly, lz, jo, i_local) + &
                           dg_frag%phi_frag(lx-1, ly, lz, jo, i_local)) + &
              lapt(2,1) * (dg_frag%phi_frag(lx+2, ly, lz, jo, i_local) + &
                           dg_frag%phi_frag(lx-2, ly, lz, jo, i_local)) + &
              lapt(3,1) * (dg_frag%phi_frag(lx+3, ly, lz, jo, i_local) + &
                           dg_frag%phi_frag(lx-3, ly, lz, jo, i_local)) + &
              lapt(4,1) * (dg_frag%phi_frag(lx+4, ly, lz, jo, i_local) + &
                           dg_frag%phi_frag(lx-4, ly, lz, jo, i_local))
          
          v = v + &
              lapt(1,2) * (dg_frag%phi_frag(lx, ly+1, lz, jo, i_local) + &
                           dg_frag%phi_frag(lx, ly-1, lz, jo, i_local)) + &
              lapt(2,2) * (dg_frag%phi_frag(lx, ly+2, lz, jo, i_local) + &
                           dg_frag%phi_frag(lx, ly-2, lz, jo, i_local)) + &
              lapt(3,2) * (dg_frag%phi_frag(lx, ly+3, lz, jo, i_local) + &
                           dg_frag%phi_frag(lx, ly-3, lz, jo, i_local)) + &
              lapt(4,2) * (dg_frag%phi_frag(lx, ly+4, lz, jo, i_local) + &
                           dg_frag%phi_frag(lx, ly-4, lz, jo, i_local))
          
          v = v + &
              lapt(1,3) * (dg_frag%phi_frag(lx, ly, lz+1, jo, i_local) + &
                           dg_frag%phi_frag(lx, ly, lz-1, jo, i_local)) + &
              lapt(2,3) * (dg_frag%phi_frag(lx, ly, lz+2, jo, i_local) + &
                           dg_frag%phi_frag(lx, ly, lz-2, jo, i_local)) + &
              lapt(3,3) * (dg_frag%phi_frag(lx, ly, lz+3, jo, i_local) + &
                           dg_frag%phi_frag(lx, ly, lz-3, jo, i_local)) + &
              lapt(4,3) * (dg_frag%phi_frag(lx, ly, lz+4, jo, i_local) + &
                           dg_frag%phi_frag(lx, ly, lz-4, jo, i_local))
          
          ! T|φ> = (-∇²/2)|φ> = lap0*φ - 0.5 * (sum of neighbor terms)
          T_phi(lx, ly, lz) = lap0 * dg_frag%phi_frag(lx, ly, lz, jo, i_local) - 0.5d0 * v
          
        end do
      end do
    end do
    
  end subroutine apply_kinetic_to_basis

  subroutine get_fragment_local_range(dg_frag, ndom, loc_s, loc_e)
    use salmon_global, only: nproc_rgrid
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ndom(3)
    integer, intent(out) :: loc_s(3), loc_e(3)

    integer :: ipx, ipy, ipz, coords(3), nsize

    ipx = max(1, nproc_rgrid(1))
    ipy = max(1, nproc_rgrid(2))
    ipz = max(1, nproc_rgrid(3))

    if (dg_frag%id_frag < 0 .or. dg_frag%id_frag >= ipx * ipy * ipz) then
      stop "DG-Fragment RT: invalid fragment-local MPI rank in get_fragment_local_range"
    end if

    coords(1) = mod(dg_frag%id_frag, ipx)
    coords(2) = mod(dg_frag%id_frag / ipx, ipy)
    coords(3) = dg_frag%id_frag / max(1, ipx * ipy)

    nsize = (ndom(1) + ipx - 1) / ipx
    loc_s(1) = 1 + nsize * coords(1)
    loc_e(1) = min(ndom(1), loc_s(1) + nsize - 1)

    nsize = (ndom(2) + ipy - 1) / ipy
    loc_s(2) = 1 + nsize * coords(2)
    loc_e(2) = min(ndom(2), loc_s(2) + nsize - 1)

    nsize = (ndom(3) + ipz - 1) / ipz
    loc_s(3) = 1 + nsize * coords(3)
    loc_e(3) = min(ndom(3), loc_s(3) + nsize - 1)
  end subroutine get_fragment_local_range

  subroutine apply_gradient_to_basis_ops_local_2d(dg_frag, i_local, jo, mg, stencil, loc_s, loc_e, grad_phi, grad_local_2d)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer,                intent(in) :: i_local, jo
    type(s_rgrid),          intent(in) :: mg
    type(s_stencil),        intent(in) :: stencil
    integer,                intent(in) :: loc_s(3), loc_e(3)
    real(8),                intent(out) :: grad_phi(:,:,:,:)
    real(8),                intent(out) :: grad_local_2d(:,:)

    integer :: lx, ly, lz, ifrag, ndom(3), nloc1, nloc2, ipt
    real(8) :: nabt(4,3), gx, gy, gz

    nabt = stencil%coef_nab
    ifrag = dg_frag%ifrag_start + i_local - 1
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    nloc1 = loc_e(1) - loc_s(1) + 1
    nloc2 = loc_e(2) - loc_s(2) + 1

    !$omp parallel do collapse(3) private(lx, ly, lz, ipt, gx, gy, gz) schedule(static)
    do lz = 1, ndom(3)
      do ly = 1, ndom(2)
        do lx = 1, ndom(1)
          gx = &
              nabt(1,1) * (dg_frag%phi_frag(lx+1, ly, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx-1, ly, lz, jo, i_local)) + &
              nabt(2,1) * (dg_frag%phi_frag(lx+2, ly, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx-2, ly, lz, jo, i_local)) + &
              nabt(3,1) * (dg_frag%phi_frag(lx+3, ly, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx-3, ly, lz, jo, i_local)) + &
              nabt(4,1) * (dg_frag%phi_frag(lx+4, ly, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx-4, ly, lz, jo, i_local))

          gy = &
              nabt(1,2) * (dg_frag%phi_frag(lx, ly+1, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly-1, lz, jo, i_local)) + &
              nabt(2,2) * (dg_frag%phi_frag(lx, ly+2, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly-2, lz, jo, i_local)) + &
              nabt(3,2) * (dg_frag%phi_frag(lx, ly+3, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly-3, lz, jo, i_local)) + &
              nabt(4,2) * (dg_frag%phi_frag(lx, ly+4, lz, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly-4, lz, jo, i_local))

          gz = &
              nabt(1,3) * (dg_frag%phi_frag(lx, ly, lz+1, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly, lz-1, jo, i_local)) + &
              nabt(2,3) * (dg_frag%phi_frag(lx, ly, lz+2, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly, lz-2, jo, i_local)) + &
              nabt(3,3) * (dg_frag%phi_frag(lx, ly, lz+3, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly, lz-3, jo, i_local)) + &
              nabt(4,3) * (dg_frag%phi_frag(lx, ly, lz+4, jo, i_local) - &
                           dg_frag%phi_frag(lx, ly, lz-4, jo, i_local))

          grad_phi(lx, ly, lz, 1) = gx
          grad_phi(lx, ly, lz, 2) = gy
          grad_phi(lx, ly, lz, 3) = gz

          if (lx >= loc_s(1) .and. lx <= loc_e(1) .and. &
              ly >= loc_s(2) .and. ly <= loc_e(2) .and. &
              lz >= loc_s(3) .and. lz <= loc_e(3)) then
            ipt = ((lz - loc_s(3)) * nloc2 + (ly - loc_s(2))) * nloc1 + (lx - loc_s(1)) + 1
            grad_local_2d(ipt, 1) = gx
            grad_local_2d(ipt, 2) = gy
            grad_local_2d(ipt, 3) = gz
          end if
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine apply_gradient_to_basis_ops_local_2d

  !=======================================================================
  ! Add non-local pseudopotential contribution to Hamiltonian matrix
  !
  ! Calculates <φ_i|V_NL|φ_j> = Σ_ilma <φ_i|proj_ilma> V_ilma <proj_ilma|φ_j>
  ! where proj_ilma are the pseudopotential projector functions
  !
  ! NUMERICAL ACCURACY: Store unnormalized overlaps, apply rinv_uvu once
  ! This prevents rinv_uvu^2 error amplification and follows SALMON convention
  !=======================================================================
  subroutine add_nonlocal_pp_matrix(dg_frag, mg, ppg, nspin, hvol)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_rgrid),          intent(in) :: mg
    type(s_pp_grid),        intent(in) :: ppg
    integer,                intent(in)    :: nspin
    real(8),                intent(in)    :: hvol
    
    integer :: ifrag, ispin, io, jo, i_local, ilma, ia, j, ix, iy, iz, ig_i, ig_j, nbf
    integer :: is(3), ie(3), ifrag_count
    integer :: iorg(3), ndom(3), lx, ly, lz
    real(8), allocatable :: uVpsi(:,:,:)  ! Projector overlaps (unnormalized): (nstate_frag, Nlma, nspin)
    real(8) :: overlap_i, overlap_j, nlpp_contrib
    
    if (ppg%Nlma == 0) return  ! No non-local PP
    
    is = mg%is
    ie = mg%ie
    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    
    ! Allocate array for projector overlaps
    ! MEMORY OPTIMIZATION: Only store current fragment's data (removed ifrag_count dimension)
    allocate(uVpsi(dg_frag%nstate_frag, ppg%Nlma, nspin))
    uVpsi = 0.0d0
    
    ! Loop over fragments assigned to this rank
    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      iorg(:) = dg_frag%ixyz_frag(:, ifrag)
      ndom(:) = dg_frag%nxyz_domain(:, ifrag)
      
      ! Reset uVpsi for this fragment (memory reuse)
      uVpsi = 0.0d0
      
      do ispin = 1, nspin
        
        ! Calculate projector overlaps <φ_io|proj_ilma> for all basis and projectors
        ! OpenMP parallelization over basis functions
!$omp parallel do collapse(2) private(ilma, io, ia, j, ix, iy, iz, overlap_i)
        do ilma = 1, ppg%Nlma
          do io = 1, dg_frag%nstate_frag
            
            ia = ppg%ia_tbl(ilma)  ! Atom index for this projector
            
            ! Calculate <φ_io|proj_ilma> = Σ_j φ_io(r_j) * uV(r_j, ilma) * hvol
            ! NOTE: Store UNNORMALIZED overlap to avoid rinv_uvu^2 numerical error
            overlap_i = 0.0d0
            do j = 1, ppg%mps(ia)
              ix = ppg%jxyz(1, j, ia)
              iy = ppg%jxyz(2, j, ia)
              iz = ppg%jxyz(3, j, ia)
              
              ! Map global projector-grid index to fragment-local basis index.
              lx = ix - iorg(1) + 1
              ly = iy - iorg(2) + 1
              lz = iz - iorg(3) + 1
              if (lx >= 1 .and. lx <= ndom(1) .and. &
                  ly >= 1 .and. ly <= ndom(2) .and. &
                  lz >= 1 .and. lz <= ndom(3)) then
                overlap_i = overlap_i + &
                  dg_frag%phi_frag(lx, ly, lz, io, i_local) * ppg%uV(j, ilma) * hvol
              end if
            end do
            
            ! Store unnormalized overlap (normalization applied once in matrix calculation)
            uVpsi(io, ilma, ispin) = overlap_i
            
          end do  ! io loop
        end do  ! ilma loop
!$omp end parallel do
        
        ! Calculate matrix elements <φ_i|V_NL|φ_j> = Σ_ilma <φ_i|proj_ilma> V_ilma <proj_ilma|φ_j>
        ! where V_ilma is encoded in rinv_uvu (includes normalization and energy coefficient)
        ! OpenMP parallelization over matrix elements
        nbf = dg_frag%n_basis(ifrag, ispin)
!$omp parallel do collapse(2) private(jo, io, ilma, nlpp_contrib, overlap_i, overlap_j, ig_i, ig_j)
        do jo = 1, nbf
          do io = 1, nbf
            ig_i = dg_frag%index_basis(io, ifrag, ispin)
            ig_j = dg_frag%index_basis(jo, ifrag, ispin)
            if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
            if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle

! Sum over all projectors
            nlpp_contrib = 0.0d0
            do ilma = 1, ppg%Nlma
              
              ! Get unnormalized overlaps
              overlap_i = uVpsi(io, ilma, ispin)
              overlap_j = uVpsi(jo, ilma, ispin)
              
              ! V_NL matrix element contribution from this projector
              ! Physical formula: <i|V_NL|j> = Σ_ilma <i|proj> * V_coeff * <proj|j>
              ! where V_coeff = rinv_uvu contains normalization and energy
              ! NUMERICAL ACCURACY: Apply rinv_uvu ONCE to avoid error amplification
              nlpp_contrib = nlpp_contrib + overlap_i * overlap_j * ppg%rinv_uvu(ilma)
              
            end do  ! ilma loop
            
            ! Add non-local PP contribution to Hamiltonian matrix
!$omp atomic update
            dg_frag%H_mat(ig_i, ig_j, ispin) = dg_frag%H_mat(ig_i, ig_j, ispin) + nlpp_contrib
            
          end do  ! io loop
        end do  ! jo loop
!$omp end parallel do
        
      end do  ! ispin loop
      
    end do  ! ifrag loop
    
    deallocate(uVpsi)
    
  end subroutine add_nonlocal_pp_matrix

  !=======================================================================
  ! Calculate momentum matrix elements in fragment basis (velocity gauge)
  !=======================================================================
  subroutine calculate_momentum_matrix(dg_frag, system, mg, stencil)
    use structures
    use communication, only: comm_is_root, comm_summation, comm_get_groupinfo, comm_get_min, comm_get_max
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_stencil),        intent(in)    :: stencil
    
    integer :: ifrag, i_local, ispin, io, jo, idir
    integer :: ix, iy, iz, is(3), ie(3), i_halo, jfrag, n_basis_halo, ig_row, ig_col, ig_i, ig_j, l(3), d(3)
    integer :: lx, ly, lz, gx, gy, gz, iorg(3), ndom(3), loc_s(3), loc_e(3), halo_s(3), halo_e(3)
    integer :: phi_lb1, phi_ub1, phi_lb2, phi_ub2, phi_lb3, phi_ub3
    integer :: grad_lb1, grad_ub1, grad_lb2, grad_ub2, grad_lb3, grad_ub3
    integer :: iblk, iblk_rev, iblk_self, ii, jj, mat_size
    integer :: npts_local, npts_halo, ipt
    logical :: log_frag_progress
    real(8) :: hvol, integral
    real(8) :: momentum_gb
    real(8) :: max_p, pavg
    real(8) :: t0, t1
    real(8) :: time_halo_exchange, time_self_integral, time_halo_integral
    real(8) :: time_grad_total
    real(8) :: time_block_reduce, time_antisym
    real(8) :: time_reduce_pack, time_reduce_comm, time_reduce_unpack
    real(8) :: frag_grad_start, frag_self_start, frag_halo_start
    real(8), allocatable :: grad_phi(:,:,:,:)  ! gradient of basis function (x,y,z components, fragment-local)
    real(8), allocatable :: grad_local_2d(:,:), phi_local_2d(:,:), self_proj(:,:)
    real(8), allocatable :: grad_halo_2d(:,:), halo_buf_2d(:,:), halo_proj(:,:)
    
    if (.not. dg_frag%has_real_space_basis) return
    time_halo_exchange = 0.0d0
    time_self_integral = 0.0d0
    time_halo_integral = 0.0d0
    time_grad_total = 0.0d0
    time_block_reduce = 0.0d0
    time_antisym = 0.0d0
    time_reduce_pack = 0.0d0
    time_reduce_comm = 0.0d0
    time_reduce_unpack = 0.0d0
    
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "        Computing transition moments: <φ_i|∇|φ_j>"
      flush(6)
    end if
    momentum_gb = real(3_8 * int(dg_frag%n_mat_max, kind=8) * int(dg_frag%n_mat_max, kind=8) * &
      int(dg_frag%nspin, kind=8) * 8_8, 8) / 1.0d9
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,i0,a,f10.3,a)') "        n_mat_max=", dg_frag%n_mat_max, &
        " momentum_mat GB=", momentum_gb, " per rank"
      flush(6)
    end if
    
    ! Allocate momentum matrix: (3 directions, n_mat_max x n_mat_max, nspin)
    ! Momentum matrix elements for vector potential coupling: p_ij = <phi_i|p|phi_j>
    ! In velocity gauge: H(t) = H_0 - i*A(t)·∇ + A(t)^2/2
    ! The A·p term couples to momentum matrix elements
    ! The A^2/2 term is diagonal (diamagnetic contribution)
    if (allocated(dg_frag%momentum_mat)) deallocate(dg_frag%momentum_mat)
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "        momentum alloc begin"
      flush(6)
    end if
    call init_momentum_blocks(dg_frag)
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "        momentum real alloc done"
      flush(6)
    end if
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "        momentum zero real done"
      flush(6)
    end if
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "        momentum alloc done"
      flush(6)
    end if
    is = mg%is
    ie = mg%ie
    hvol = system%hvol
    
    ! Exchange halo regions before stencil operations
    call cpu_time(t0)
    call exchange_phi_frag_halo(dg_frag)
    call cpu_time(t1)
    time_halo_exchange = time_halo_exchange + (t1 - t0)
    write(*,'(1x,a,i0,a,i0,a,i0,a,a)') "        momentum stage: rank=", dg_frag%id, &
      " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " stage=", "after-halo-exchange"
    flush(6)
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "        momentum halo exchange done"
      flush(6)
    end if
    
    ! Loop over spin
    do ispin = 1, system%nspin
      ! Loop over local fragments
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        log_frag_progress = dg_frag%is_frag_root .and. ifrag == dg_frag%ifrag_start
        frag_grad_start = time_grad_total
        frag_self_start = time_self_integral
        frag_halo_start = time_halo_integral
        if (ifrag == dg_frag%ifrag_start) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        momentum stage: rank=", dg_frag%id, &
            " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, &
            " ifrag=", ifrag, " stage=", "fragment-begin"
          flush(6)
        end if
        if (log_frag_progress) then
          write(*,'(1x,a,i0,a,i0)') "        momentum fragment begin: ifrag=", ifrag, " i_local=", i_local
          flush(6)
        end if
        iorg(:) = dg_frag%ixyz_frag(:, ifrag)
        ndom(:) = dg_frag%nxyz_domain(:, ifrag)
        call get_fragment_local_range(dg_frag, ndom, loc_s, loc_e)
        
        ! Cache the local basis matrix once per fragment; it does not depend on jo.
        if (i_local < 1 .or. i_local > size(dg_frag%phi_frag, 5)) then
          write(*,'(1x,a,i0,a,i0)') "DG-Fragment RT invalid i_local=", i_local, " phi_frag dim5=", size(dg_frag%phi_frag, 5)
          stop "DG-Fragment RT: invalid fragment-local basis index"
        end if
        phi_lb1 = lbound(dg_frag%phi_frag, 1)
        phi_ub1 = ubound(dg_frag%phi_frag, 1)
        phi_lb2 = lbound(dg_frag%phi_frag, 2)
        phi_ub2 = ubound(dg_frag%phi_frag, 2)
        phi_lb3 = lbound(dg_frag%phi_frag, 3)
        phi_ub3 = ubound(dg_frag%phi_frag, 3)
        if (loc_s(1) < phi_lb1 .or. loc_e(1) > phi_ub1 .or. &
            loc_s(2) < phi_lb2 .or. loc_e(2) > phi_ub2 .or. &
            loc_s(3) < phi_lb3 .or. loc_e(3) > phi_ub3) then
          write(*,'(1x,a,1x,3(i0,1x),a,1x,3(i0,1x),a,1x,3(i0,1x),a,1x,3(i0,1x))') &
            "DG-Fragment RT momentum phi_frag local range out of bounds: loc_s=", &
            loc_s(1), loc_s(2), loc_s(3), "loc_e=", loc_e(1), loc_e(2), loc_e(3), &
            "lb=", phi_lb1, phi_lb2, phi_lb3, "ub=", phi_ub1, phi_ub2, phi_ub3
          stop "DG-Fragment RT: momentum phi_frag local range out of bounds"
        end if
        npts_local = (loc_e(1) - loc_s(1) + 1) * (loc_e(2) - loc_s(2) + 1) * (loc_e(3) - loc_s(3) + 1)
        allocate(phi_local_2d(npts_local, dg_frag%n_basis(ifrag, ispin)), &
                 grad_local_2d(npts_local, 3), self_proj(dg_frag%n_basis(ifrag, ispin), 3))

        do io = 1, dg_frag%n_basis(ifrag, ispin)
          if (log_frag_progress .and. io == 1) then
            write(*,'(1x,a,i0,a,i0,a,i0)') "        momentum first io: io=", io, " n_basis=", dg_frag%n_basis(ifrag, ispin), &
              " phi_dim4=", size(dg_frag%phi_frag, 4)
            write(*,'(1x,a,3(i0,1x),a,3(i0,1x))') &
              "        momentum first local range: loc_s=", loc_s(1), loc_s(2), loc_s(3), &
              " loc_e=", loc_e(1), loc_e(2), loc_e(3)
            flush(6)
          end if
          if (io < 1 .or. io > size(dg_frag%phi_frag, 4)) then
            write(*,'(1x,a,i0,a,i0)') "DG-Fragment RT invalid io=", io, " phi_frag dim4=", size(dg_frag%phi_frag, 4)
            stop "DG-Fragment RT: invalid basis-function index"
          end if
          ipt = 0
          do lz = loc_s(3), loc_e(3)
            do ly = loc_s(2), loc_e(2)
              do lx = loc_s(1), loc_e(1)
                ipt = ipt + 1
                phi_local_2d(ipt, io) = dg_frag%phi_frag(lx, ly, lz, io, i_local)
              end do
            end do
          end do
        end do
        if (ifrag == dg_frag%ifrag_start) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        momentum stage: rank=", dg_frag%id, &
            " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, &
            " ifrag=", ifrag, " stage=", "after-local-pack"
          flush(6)
        end if
        iblk_self = find_momentum_block(dg_frag, ifrag, ifrag)

        ! Loop over basis functions in fragment j (ket side)
        ! Keep this loop serial to avoid per-thread duplication of large grad_phi buffers.
        ! Parallelism is still provided inside apply_gradient_to_basis and SIMD in accumulations.
        do jo = 1, dg_frag%n_basis(ifrag, ispin)
          if (log_frag_progress .and. jo == 1) then
            write(*,*) "        momentum first jo begin"
            flush(6)
          end if

          allocate(grad_phi(1:ndom(1), 1:ndom(2), 1:ndom(3), 3))
          call cpu_time(t0)
          call apply_gradient_to_basis_ops_local_2d(dg_frag, i_local, jo, mg, stencil, loc_s, loc_e, grad_phi, grad_local_2d)
          call cpu_time(t1)
          time_grad_total = time_grad_total + (t1 - t0)
          if (log_frag_progress .and. jo == 1) then
            write(*,*) "        momentum first gradient done"
            flush(6)
          end if
          if (log_frag_progress) then
            if (jo == 1 .or. jo == dg_frag%n_basis(ifrag, ispin) .or. &
                mod(jo, max(1, dg_frag%n_basis(ifrag, ispin) / 4)) == 0) then
              write(*,'(1x,a,i0,a,i0,a,1pe12.4)') "        momentum basis progress: jo=", jo, "/", &
                dg_frag%n_basis(ifrag, ispin), " grad=", time_grad_total
              flush(6)
            end if
          end if

          grad_lb1 = lbound(grad_phi, 1)
          grad_ub1 = ubound(grad_phi, 1)
          grad_lb2 = lbound(grad_phi, 2)
          grad_ub2 = ubound(grad_phi, 2)
          grad_lb3 = lbound(grad_phi, 3)
          grad_ub3 = ubound(grad_phi, 3)
          if (loc_s(1) < grad_lb1 .or. loc_e(1) > grad_ub1 .or. &
              loc_s(2) < grad_lb2 .or. loc_e(2) > grad_ub2 .or. &
              loc_s(3) < grad_lb3 .or. loc_e(3) > grad_ub3) then
            write(*,'(1x,a,1x,3(i0,1x),a,1x,3(i0,1x),a,1x,3(i0,1x),a,1x,3(i0,1x))') &
              "DG-Fragment RT momentum grad_phi local range out of bounds: loc_s=", &
              loc_s(1), loc_s(2), loc_s(3), "loc_e=", loc_e(1), loc_e(2), loc_e(3), &
              "lb=", grad_lb1, grad_lb2, grad_lb3, "ub=", grad_ub1, grad_ub2, grad_ub3
            stop "DG-Fragment RT: momentum grad_phi local range out of bounds"
          end if

          ig_j = dg_frag%index_basis(jo, ifrag, ispin)
          call cpu_time(t0)
          call dgemm('T', 'N', dg_frag%n_basis(ifrag, ispin), 3, npts_local, hvol, phi_local_2d, npts_local, &
            grad_local_2d, npts_local, 0.0d0, self_proj, dg_frag%n_basis(ifrag, ispin))
          call cpu_time(t1)
          time_self_integral = time_self_integral + (t1 - t0)

          do io = 1, dg_frag%n_basis(ifrag, ispin)
            ig_i = dg_frag%index_basis(io, ifrag, ispin)
            if (log_frag_progress .and. jo == 1 .and. io == 1) then
              write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "        momentum first index: ig_i=", ig_i, " ig_j=", ig_j, &
                " n_mat_max=", dg_frag%n_mat_max, " ispin=", ispin
              write(*,*) "        momentum first idir begin"
              flush(6)
            end if
            if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
            if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle
            do idir = 1, 3
              integral = self_proj(io, idir)
              if (iblk_self > 0) dg_frag%momentum_blocks(iblk_self)%val(idir, io, jo, ispin) = integral
              if (log_frag_progress .and. jo == 1 .and. io == 1 .and. idir == 1) then
                write(*,'(1x,a,1x,es12.4)') "        momentum first idir done integral=", integral
                flush(6)
              end if
            end do
          end do

          ig_col = ig_j
          if (ig_col >= 1 .and. ig_col <= dg_frag%n_mat_max) then
            do i_halo = 1, dg_frag%n_halo
              if (dg_frag%halo(i_halo)%ifrag_dst /= ifrag) cycle
              jfrag = dg_frag%halo(i_halo)%ifrag_src
              n_basis_halo = dg_frag%n_basis(jfrag, ispin)
              l = dg_frag%halo(i_halo)%length
              d = dg_frag%halo(i_halo)%dsp_send
              halo_s(:) = max(loc_s(:), d(:) + 1)
              halo_e(:) = min(loc_e(:), d(:) + l(:))
              if (any(halo_s(:) > halo_e(:))) cycle
              if (log_frag_progress) then
                if (jo == 1 .or. jo == dg_frag%n_basis(ifrag, ispin) .or. &
                    mod(jo, max(1, dg_frag%n_basis(ifrag, ispin) / 4)) == 0) then
                  write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,3(i0,1x),a,1pe12.4)') &
                    "        momentum halo detail: jo=", jo, " i_halo=", i_halo, &
                    " jfrag=", jfrag, " n_basis_halo=", n_basis_halo, " halo_len=", &
                    l(1), l(2), l(3), " time_halo=", time_halo_integral
                  flush(6)
                end if
              end if

              npts_halo = (halo_e(1) - halo_s(1) + 1) * (halo_e(2) - halo_s(2) + 1) * (halo_e(3) - halo_s(3) + 1)
              if (npts_halo <= 0 .or. n_basis_halo <= 0) cycle

              allocate(grad_halo_2d(npts_halo, 3), halo_buf_2d(npts_halo, n_basis_halo), halo_proj(n_basis_halo, 3))

              ipt = 0
              do lz = halo_s(3), halo_e(3)
                do ly = halo_s(2), halo_e(2)
                  do lx = halo_s(1), halo_e(1)
                    ipt = ipt + 1
                    grad_halo_2d(ipt, 1:3) = grad_phi(lx, ly, lz, 1:3)
                  end do
                end do
              end do

              do io = 1, n_basis_halo
                ipt = 0
                do lz = halo_s(3), halo_e(3)
                  iz = lz - d(3)
                  do ly = halo_s(2), halo_e(2)
                    iy = ly - d(2)
                    do lx = halo_s(1), halo_e(1)
                      ix = lx - d(1)
                      ipt = ipt + 1
                      halo_buf_2d(ipt, io) = dg_frag%halo(i_halo)%buf_recv(ix, iy, iz, io, 1)
                    end do
                  end do
                end do
              end do

              call cpu_time(t0)
              call dgemm('T', 'N', n_basis_halo, 3, npts_halo, hvol, halo_buf_2d, npts_halo, &
                grad_halo_2d, npts_halo, 0.0d0, halo_proj, n_basis_halo)
              call cpu_time(t1)
              time_halo_integral = time_halo_integral + (t1 - t0)

              iblk = find_momentum_block(dg_frag, jfrag, ifrag)
              iblk_rev = find_momentum_block(dg_frag, ifrag, jfrag)
              do io = 1, n_basis_halo
                ig_row = dg_frag%index_basis(io, jfrag, ispin)
                if (ig_row < 1 .or. ig_row > dg_frag%n_mat_max) cycle
                do idir = 1, 3
                  integral = halo_proj(io, idir)
                  if (iblk > 0) then
                    dg_frag%momentum_blocks(iblk)%val(idir, io, jo, ispin) = &
                      dg_frag%momentum_blocks(iblk)%val(idir, io, jo, ispin) + 0.5d0 * integral
                  end if
                  if (iblk_rev > 0) then
                    dg_frag%momentum_blocks(iblk_rev)%val(idir, jo, io, ispin) = &
                      dg_frag%momentum_blocks(iblk_rev)%val(idir, jo, io, ispin) - 0.5d0 * integral
                  end if
                end do
              end do

              deallocate(grad_halo_2d, halo_buf_2d, halo_proj)
            end do
          end if

          deallocate(grad_phi)
        end do  ! jo
        if (ifrag == dg_frag%ifrag_start) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        momentum stage: rank=", dg_frag%id, &
            " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, &
            " ifrag=", ifrag, " stage=", "after-jo-loop"
          flush(6)
        end if

        deallocate(phi_local_2d, grad_local_2d, self_proj)

        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,l1,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
            "        momentum fragment done: rank=", dg_frag%id, " id_frag=", dg_frag%id_frag, &
            " ifrag=", ifrag, " ispin=", ispin, " root=", dg_frag%is_frag_root, &
            " grad=", time_grad_total - frag_grad_start, &
            " self=", time_self_integral - frag_self_start, &
            " halo=", time_halo_integral - frag_halo_start
        flush(6)
      end do  ! ifrag
    end do  ! ispin
    
    ! MPI aggregation of fragment-neighbor momentum blocks.
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "        momentum reduce begin"
      flush(6)
    end if
    block
      real(8), allocatable :: send_flat(:), recv_flat(:)
      real(8) :: t_reduce_start, t_reduce_end
      integer :: nrow, ncol
      integer :: total_size, offset_flat, chunk_size, chunk_begin, chunk_count
      integer :: total_size_min, total_size_max, nblk_min, nblk_max, ifrag_chk
      real(8) :: meta_sig_blocks, meta_sig_basis
      real(8) :: meta_sig_blocks_min(1), meta_sig_blocks_max(1)
      real(8) :: meta_sig_basis_min(1), meta_sig_basis_max(1)
      call cpu_time(t_reduce_start)
      total_size = 0
      do iblk = 1, dg_frag%n_momentum_blocks
        nrow = dg_frag%momentum_blocks(iblk)%nrow_max
        ncol = dg_frag%momentum_blocks(iblk)%ncol_max
        total_size = total_size + 3 * nrow * ncol * dg_frag%nspin
      end do
      if (comm_is_root(dg_frag%id)) then
        write(*,*) "        momentum reduce metadata begin"
        flush(6)
      end if
      write(*,'(1x,a,i0,a,i0,a,i0,a,l1,a,i0)') "        momentum reduce reach: rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, &
        " root=", dg_frag%is_frag_root, " nblk=", dg_frag%n_momentum_blocks
      flush(6)

      nblk_max = dg_frag%n_momentum_blocks
      call comm_get_max(nblk_max, dg_frag%icomm)
      nblk_min = -dg_frag%n_momentum_blocks
      call comm_get_max(nblk_min, dg_frag%icomm)
      nblk_min = -nblk_min
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,i0,a,i0)') "        momentum reduce metadata nblk done: min=", nblk_min, " max=", nblk_max
        flush(6)
      end if
      if (nblk_min /= nblk_max) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "        [FATAL] momentum reduce metadata mismatch: rank=", &
          dg_frag%id, " nblk=", dg_frag%n_momentum_blocks, " min=", nblk_min, " max=", nblk_max
        flush(6)
        stop "DG-Fragment RT: momentum reduce n_momentum_blocks mismatch across MPI ranks"
      end if

      meta_sig_blocks = 0.0d0
      do iblk = 1, dg_frag%n_momentum_blocks
        meta_sig_blocks = meta_sig_blocks + &
          dble(iblk) * 1.0d0 + &
          dble(dg_frag%momentum_blocks(iblk)%ifrag_row) * 1.0d3 + &
          dble(dg_frag%momentum_blocks(iblk)%ifrag_col) * 1.0d6 + &
          dble(dg_frag%momentum_blocks(iblk)%nrow_max) * 1.0d9 + &
          dble(dg_frag%momentum_blocks(iblk)%ncol_max) * 1.0d12
      end do
      meta_sig_blocks_min(1) = meta_sig_blocks
      meta_sig_blocks_max(1) = meta_sig_blocks
      call comm_get_min(meta_sig_blocks_min, meta_sig_blocks_min, 1, dg_frag%icomm)
      call comm_get_max(meta_sig_blocks_max, meta_sig_blocks_max, 1, dg_frag%icomm)
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,1pe14.6,a,1pe14.6)') "        momentum reduce metadata block-sig done: min=", &
          meta_sig_blocks_min(1), " max=", meta_sig_blocks_max(1)
        flush(6)
      end if

      meta_sig_basis = 0.0d0
      do ispin = 1, dg_frag%nspin
        do ifrag_chk = 1, dg_frag%n_frag
          meta_sig_basis = meta_sig_basis + &
            dble(ifrag_chk) * 1.0d0 + &
            dble(ispin) * 1.0d3 + &
            dble(dg_frag%n_basis(ifrag_chk, ispin)) * 1.0d6
        end do
      end do
      meta_sig_basis_min(1) = meta_sig_basis
      meta_sig_basis_max(1) = meta_sig_basis
      call comm_get_min(meta_sig_basis_min, meta_sig_basis_min, 1, dg_frag%icomm)
      call comm_get_max(meta_sig_basis_max, meta_sig_basis_max, 1, dg_frag%icomm)
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,1pe14.6,a,1pe14.6)') "        momentum reduce metadata basis-sig done: min=", &
          meta_sig_basis_min(1), " max=", meta_sig_basis_max(1)
        flush(6)
      end if

      if (abs(meta_sig_blocks_max(1) - meta_sig_blocks_min(1)) > 1.0d-9 .or. &
          abs(meta_sig_basis_max(1) - meta_sig_basis_min(1)) > 1.0d-9) then
        write(*,'(1x,a,i0,a,1pe14.6,a,1pe14.6,a,1pe14.6,a,1pe14.6)') &
          "        [FATAL] momentum reduce metadata mismatch: rank=", dg_frag%id, &
          " block_min=", meta_sig_blocks_min(1), " block_max=", meta_sig_blocks_max(1), &
          " basis_min=", meta_sig_basis_min(1), " basis_max=", meta_sig_basis_max(1)
        flush(6)
        stop "DG-Fragment RT: momentum reduce metadata signature mismatch across MPI ranks"
      end if

      ! All ranks must call the same number of collective operations with the
      ! same count. If total_size diverges, chunked allreduce can deadlock.
      total_size_max = total_size
      call comm_get_max(total_size_max, dg_frag%icomm)
      total_size_min = -total_size
      call comm_get_max(total_size_min, dg_frag%icomm)
      total_size_min = -total_size_min
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,i0,a,i0)') "        momentum reduce metadata total-size done: min=", total_size_min, &
          " max=", total_size_max
        flush(6)
      end if
      if (total_size_min /= total_size_max) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "        [FATAL] momentum reduce size mismatch: rank=", &
          dg_frag%id, " local=", total_size, " min=", total_size_min, " max=", total_size_max
        flush(6)
        stop "DG-Fragment RT: momentum reduce total_size mismatch across MPI ranks"
      end if

      if (total_size > 0) then
        allocate(send_flat(total_size), recv_flat(total_size))
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,i0)') "        momentum reduce pack begin total_size=", total_size
          flush(6)
        end if
        call cpu_time(t0)
        offset_flat = 1
        do iblk = 1, dg_frag%n_momentum_blocks
          nrow = dg_frag%momentum_blocks(iblk)%nrow_max
          ncol = dg_frag%momentum_blocks(iblk)%ncol_max
          do ispin = 1, dg_frag%nspin
            do jj = 1, ncol
              do ii = 1, nrow
                do idir = 1, 3
                  send_flat(offset_flat) = dg_frag%momentum_blocks(iblk)%val(idir, ii, jj, ispin)
                  offset_flat = offset_flat + 1
                end do
              end do
            end do
          end do
        end do
        call cpu_time(t1)
        time_reduce_pack = time_reduce_pack + (t1 - t0)
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,1pe12.4)') "        momentum reduce pack done time=", time_reduce_pack
          flush(6)
        end if
        chunk_size = 262144
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,i0)') "        momentum reduce comm begin chunk_size=", chunk_size
          flush(6)
        end if
        call cpu_time(t0)
        chunk_begin = 1
        do while (chunk_begin <= total_size)
          chunk_count = min(chunk_size, total_size - chunk_begin + 1)
          call comm_summation(send_flat(chunk_begin:chunk_begin + chunk_count - 1), &
                              recv_flat(chunk_begin:chunk_begin + chunk_count - 1), chunk_count, dg_frag%icomm)
          if (comm_is_root(dg_frag%id)) then
            write(*,'(1x,a,i0,a,i0)') "        momentum reduce comm chunk done: begin=", chunk_begin, &
              " count=", chunk_count
            flush(6)
          end if
          chunk_begin = chunk_begin + chunk_count
        end do
        call cpu_time(t1)
        time_reduce_comm = time_reduce_comm + (t1 - t0)
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,1pe12.4)') "        momentum reduce comm done time=", time_reduce_comm
          flush(6)
        end if
        call cpu_time(t0)
        offset_flat = 1
        do iblk = 1, dg_frag%n_momentum_blocks
          nrow = dg_frag%momentum_blocks(iblk)%nrow_max
          ncol = dg_frag%momentum_blocks(iblk)%ncol_max
          do ispin = 1, dg_frag%nspin
            do jj = 1, ncol
              do ii = 1, nrow
                do idir = 1, 3
                  dg_frag%momentum_blocks(iblk)%val(idir, ii, jj, ispin) = recv_flat(offset_flat)
                  offset_flat = offset_flat + 1
                end do
              end do
            end do
          end do
        end do
        call cpu_time(t1)
        time_reduce_unpack = time_reduce_unpack + (t1 - t0)
        if (comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,1pe12.4)') "        momentum reduce unpack done time=", time_reduce_unpack
          flush(6)
        end if
        deallocate(send_flat, recv_flat)
      end if
      call cpu_time(t_reduce_end)
      time_block_reduce = time_block_reduce + (t_reduce_end - t_reduce_start)
    end block
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') "        momentum reduce done time=", time_block_reduce, &
        " pack=", time_reduce_pack, " comm=", time_reduce_comm, " unpack=", time_reduce_unpack
      write(*,*) "        momentum antisym begin"
      flush(6)
    end if

    ! Enforce anti-symmetry blockwise: self blocks against themselves, off-diagonal
    ! blocks against the reverse ordered fragment pair.
    call cpu_time(t0)
    do ispin = 1, system%nspin
      do iblk = 1, dg_frag%n_momentum_blocks
        ifrag = dg_frag%momentum_blocks(iblk)%ifrag_row
        jfrag = dg_frag%momentum_blocks(iblk)%ifrag_col
        if (ifrag == jfrag) then
          do idir = 1, 3
            do ii = 1, dg_frag%n_basis(ifrag, ispin)
              dg_frag%momentum_blocks(iblk)%val(idir, ii, ii, ispin) = 0.0d0
              do jj = ii + 1, dg_frag%n_basis(ifrag, ispin)
                pavg = 0.5d0 * (dg_frag%momentum_blocks(iblk)%val(idir, ii, jj, ispin) - &
                                dg_frag%momentum_blocks(iblk)%val(idir, jj, ii, ispin))
                dg_frag%momentum_blocks(iblk)%val(idir, ii, jj, ispin) = pavg
                dg_frag%momentum_blocks(iblk)%val(idir, jj, ii, ispin) = -pavg
              end do
            end do
          end do
        else
          iblk_rev = find_momentum_block(dg_frag, jfrag, ifrag)
          if (iblk_rev <= 0 .or. iblk >= iblk_rev) cycle
          do idir = 1, 3
            do ii = 1, dg_frag%n_basis(ifrag, ispin)
              do jj = 1, dg_frag%n_basis(jfrag, ispin)
                pavg = 0.5d0 * (dg_frag%momentum_blocks(iblk)%val(idir, ii, jj, ispin) - &
                                dg_frag%momentum_blocks(iblk_rev)%val(idir, jj, ii, ispin))
                dg_frag%momentum_blocks(iblk)%val(idir, ii, jj, ispin) = pavg
                dg_frag%momentum_blocks(iblk_rev)%val(idir, jj, ii, ispin) = -pavg
              end do
            end do
          end do
        end if
      end do
    end do
    call cpu_time(t1)
    time_antisym = time_antisym + (t1 - t0)
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,1pe12.4)') "        momentum antisym done time=", time_antisym
      flush(6)
    end if

    max_p = 0.0d0
    do iblk = 1, dg_frag%n_momentum_blocks
      max_p = max(max_p, maxval(abs(dg_frag%momentum_blocks(iblk)%val)))
    end do
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
        "        momentum timing: halo_exchange=", time_halo_exchange, &
        " grad=", time_grad_total, " self=", time_self_integral, &
        " halo=", time_halo_integral, &
        " reduce=", time_block_reduce, " antisym=", time_antisym
      write(*,'(a,1pe12.4)') "        Max |momentum_mat|: ", max_p
      write(*,'(a,i0,a,i0,a)') "        Total matrix elements: ", &
                               3 * dg_frag%n_mat_max * dg_frag%n_mat_max * system%nspin, &
                               " (", 3, " directions × basis states × spin)"
    end if
    
  end subroutine calculate_momentum_matrix

  !=======================================================================
  ! Calculate overlap matrix in DG basis (S_ij = <phi_i|phi_j>)
  !=======================================================================
  subroutine ensure_overlap_prop_available(dg_frag, n_use)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: n_use

    if (n_use <= 0) return
    if (.not. dg_frag%overlap_prop_root_authoritative) return
    if (dg_frag%is_frag_root) return

    stop "ensure_overlap_prop_available not implemented yet"
  end subroutine ensure_overlap_prop_available

  subroutine calculate_overlap_matrix(dg_frag, system, mg)
    use structures
    use communication, only: comm_summation, comm_is_root
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg

    integer :: ifrag, i_local, ispin, io, jo
    integer :: ix, iy, iz, is(3), ie(3), i_halo, jfrag, n_basis_halo
    integer :: ig_row, ig_col, l(3), d(3), ii, jj
    integer :: lx, ly, lz, iorg(3), ndom(3), loc_s(3), loc_e(3), halo_s(3), halo_e(3)
    integer :: phi_lb1, phi_lb2, phi_lb3, phi_ub1, phi_ub2, phi_ub3
    integer :: buf_lb1, buf_lb2, buf_lb3, buf_ub1, buf_ub2, buf_ub3
    integer :: offset_flat, sample_n, nan_count
    integer :: mat_size, n_eval, lwork, info_eig
    logical :: log_frag_progress, log_reduce_diag
    real(8) :: hvol, integral, savg, s_min, s_max, cond_est
    real(8) :: t0, t1, time_self_integral, time_halo_integral, time_reduce_total
    real(8) :: frag_self_start, frag_halo_start
    real(8) :: work_query(1)
    real(8) :: sum_abs_flat, max_abs_flat
    real(8), allocatable :: S_flat(:), S_tmp_flat(:)
    real(8), allocatable :: S_eval(:,:), eigvals(:), eig_work(:)

    if (.not. dg_frag%has_real_space_basis) return
    if (.not. allocated(dg_frag%index_basis) .or. .not. allocated(dg_frag%n_mat)) return

    if (allocated(dg_frag%S_mat)) deallocate(dg_frag%S_mat)
    if (allocated(dg_frag%S_mat_prop)) deallocate(dg_frag%S_mat_prop)
    if (.not. allocated(dg_frag%S_mat)) then
      allocate(dg_frag%S_mat(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
    end if
    if (.not. allocated(dg_frag%S_mat_prop)) then
      allocate(dg_frag%S_mat_prop(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
    end if
    dg_frag%S_mat = 0.0d0
    dg_frag%S_mat_prop = 0.0d0
    dg_frag%has_global_overlap_copy = .true.
    dg_frag%overlap_prop_root_authoritative = .false.

    is = mg%is
    ie = mg%ie
    hvol = system%hvol
    time_self_integral = 0.0d0
    time_halo_integral = 0.0d0
    time_reduce_total = 0.0d0
    phi_lb1 = lbound(dg_frag%phi_frag, 1)
    phi_lb2 = lbound(dg_frag%phi_frag, 2)
    phi_lb3 = lbound(dg_frag%phi_frag, 3)
    phi_ub1 = ubound(dg_frag%phi_frag, 1)
    phi_ub2 = ubound(dg_frag%phi_frag, 2)
    phi_ub3 = ubound(dg_frag%phi_frag, 3)

    call exchange_phi_frag_halo(dg_frag)
    write(*,'(1x,a,i0,a,i0,a,i0,a,a)') "        overlap stage: rank=", dg_frag%id, &
      " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " stage=", "after-halo-exchange"
    flush(6)

    do ispin = 1, system%nspin
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        log_frag_progress = .true.
        frag_self_start = time_self_integral
        frag_halo_start = time_halo_integral
        if (i_local < 1 .or. i_local > size(dg_frag%phi_frag, 5)) then
          write(*,*) "[FATAL] overlap invalid i_local: rank=", dg_frag%id, " id_frag=", dg_frag%id_frag, &
            " ifrag=", ifrag, " i_local=", i_local, " phi_dim5=", size(dg_frag%phi_frag, 5)
          stop 1
        end if
        iorg(:) = dg_frag%ixyz_frag(:, ifrag)
        ndom(:) = dg_frag%nxyz_domain(:, ifrag)
        call get_fragment_local_range(dg_frag, ndom, loc_s, loc_e)
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        overlap stage: rank=", dg_frag%id, &
          " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, &
          " ifrag=", ifrag, " stage=", "fragment-begin"
        flush(6)
        if (log_frag_progress) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,3(i0,1x),a,3(i0,1x),a,i0)') &
            "        overlap fragment begin: ifrag=", ifrag, " i_local=", i_local, " ispin=", ispin, &
            " loc_s=", loc_s(1), loc_s(2), loc_s(3), " loc_e=", loc_e(1), loc_e(2), loc_e(3), &
            " n_basis=", dg_frag%n_basis(ifrag, ispin)
          flush(6)
        end if
        if (loc_s(1) < phi_lb1 .or. loc_e(1) > phi_ub1 .or. &
            loc_s(2) < phi_lb2 .or. loc_e(2) > phi_ub2 .or. &
            loc_s(3) < phi_lb3 .or. loc_e(3) > phi_ub3) then
          write(*,*) "[FATAL] overlap local range out of bounds: ifrag=", ifrag, "ispin=", ispin, &
            "loc_s=", loc_s(1), loc_s(2), loc_s(3), "loc_e=", loc_e(1), loc_e(2), loc_e(3), &
            "phi_lb=", phi_lb1, phi_lb2, phi_lb3, "phi_ub=", phi_ub1, phi_ub2, phi_ub3
          stop 1
        end if

        do jo = 1, dg_frag%n_basis(ifrag, ispin)
          if (jo > size(dg_frag%phi_frag, 4)) then
            write(*,*) "[FATAL] overlap jo exceeds phi dim4: rank=", dg_frag%id, " ifrag=", ifrag, &
              " jo=", jo, " phi_dim4=", size(dg_frag%phi_frag, 4)
            stop 1
          end if
          ig_col = dg_frag%index_basis(jo, ifrag, ispin)
          if (ig_col < 1 .or. ig_col > dg_frag%n_mat_max) cycle

          do io = 1, dg_frag%n_basis(ifrag, ispin)
            if (io > size(dg_frag%phi_frag, 4)) then
              write(*,*) "[FATAL] overlap io exceeds phi dim4: rank=", dg_frag%id, " ifrag=", ifrag, &
                " io=", io, " phi_dim4=", size(dg_frag%phi_frag, 4)
              stop 1
            end if
            ig_row = dg_frag%index_basis(io, ifrag, ispin)
            if (ig_row < 1 .or. ig_row > dg_frag%n_mat_max) cycle
            integral = 0.0d0
            call cpu_time(t0)
            do lz = loc_s(3), loc_e(3)
              do ly = loc_s(2), loc_e(2)
                do lx = loc_s(1), loc_e(1)
                  integral = integral + dg_frag%phi_frag(lx, ly, lz, io, i_local) * &
                             dg_frag%phi_frag(lx, ly, lz, jo, i_local) * hvol
                end do
              end do
            end do
            call cpu_time(t1)
            time_self_integral = time_self_integral + (t1 - t0)
            dg_frag%S_mat(ig_row, ig_col, ispin) = integral
          end do
          if (log_frag_progress) then
            if (jo == 1 .or. jo == dg_frag%n_basis(ifrag, ispin) .or. &
                mod(jo, max(1, dg_frag%n_basis(ifrag, ispin) / 4)) == 0) then
              write(*,'(1x,a,i0,a,i0,a,1pe12.4)') "        overlap self progress: jo=", jo, "/", &
                dg_frag%n_basis(ifrag, ispin), " self=", time_self_integral - frag_self_start
              flush(6)
            end if
          end if

          do i_halo = 1, dg_frag%n_halo
            if (dg_frag%halo(i_halo)%ifrag_dst /= ifrag) cycle
            jfrag = dg_frag%halo(i_halo)%ifrag_src
            if (jfrag < 1) cycle
            n_basis_halo = dg_frag%n_basis(jfrag, ispin)
            l = dg_frag%halo(i_halo)%length
            d = dg_frag%halo(i_halo)%dsp_send
            buf_lb1 = lbound(dg_frag%halo(i_halo)%buf_recv, 1)
            buf_lb2 = lbound(dg_frag%halo(i_halo)%buf_recv, 2)
            buf_lb3 = lbound(dg_frag%halo(i_halo)%buf_recv, 3)
            buf_ub1 = ubound(dg_frag%halo(i_halo)%buf_recv, 1)
            buf_ub2 = ubound(dg_frag%halo(i_halo)%buf_recv, 2)
            buf_ub3 = ubound(dg_frag%halo(i_halo)%buf_recv, 3)
            if (size(dg_frag%halo(i_halo)%buf_recv, 5) < 1) then
              write(*,*) "[FATAL] overlap halo buf dim5 invalid: rank=", dg_frag%id, " i_halo=", i_halo
              stop 1
            end if
            if (n_basis_halo > size(dg_frag%halo(i_halo)%buf_recv, 4)) then
              write(*,*) "[FATAL] overlap halo basis exceeds buf dim4: rank=", dg_frag%id, &
                " i_halo=", i_halo, " n_basis_halo=", n_basis_halo, " buf_dim4=", size(dg_frag%halo(i_halo)%buf_recv, 4)
              stop 1
            end if
            halo_s(:) = max(loc_s(:), d(:) + 1)
            halo_e(:) = min(loc_e(:), d(:) + l(:))
            if (any(halo_s(:) > halo_e(:))) cycle

            do io = 1, n_basis_halo
              ig_row = dg_frag%index_basis(io, jfrag, ispin)
              if (ig_row < 1 .or. ig_row > dg_frag%n_mat_max) cycle
              integral = 0.0d0
              call cpu_time(t0)
              do lz = halo_s(3), halo_e(3)
                iz = lz - d(3)
                do ly = halo_s(2), halo_e(2)
                  iy = ly - d(2)
                  do lx = halo_s(1), halo_e(1)
                    ix = lx - d(1)
                    if (ix < buf_lb1 .or. ix > buf_ub1 .or. &
                        iy < buf_lb2 .or. iy > buf_ub2 .or. &
                        iz < buf_lb3 .or. iz > buf_ub3) then
                      write(*,*) "[FATAL] overlap halo index out of bounds: ifrag=", ifrag, &
                        "jfrag=", jfrag, "i_halo=", i_halo, "idx=", ix, iy, iz, &
                        "buf_lb=", buf_lb1, buf_lb2, buf_lb3, "buf_ub=", buf_ub1, buf_ub2, buf_ub3
                      stop 1
                    end if
                    integral = integral + dg_frag%halo(i_halo)%buf_recv(ix, iy, iz, io, 1) * &
                               dg_frag%phi_frag(lx, ly, lz, jo, i_local) * hvol
                  end do
                end do
              end do
              call cpu_time(t1)
              time_halo_integral = time_halo_integral + (t1 - t0)
              dg_frag%S_mat(ig_row, ig_col, ispin) = dg_frag%S_mat(ig_row, ig_col, ispin) + 0.5d0 * integral
              dg_frag%S_mat(ig_col, ig_row, ispin) = dg_frag%S_mat(ig_col, ig_row, ispin) + 0.5d0 * integral
            end do
            if (log_frag_progress) then
              if (jo == 1 .or. jo == dg_frag%n_basis(ifrag, ispin) .or. &
                  mod(jo, max(1, dg_frag%n_basis(ifrag, ispin) / 4)) == 0) then
                write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,3(i0,1x),a,1pe12.4)') &
                  "        overlap halo detail: jo=", jo, " i_halo=", i_halo, " jfrag=", jfrag, &
                  " n_basis_halo=", n_basis_halo, " halo_len=", l(1), l(2), l(3), &
                  " halo=", time_halo_integral - frag_halo_start
                flush(6)
              end if
            end if
          end do

        end do
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') "        overlap stage: rank=", dg_frag%id, &
          " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, &
          " ifrag=", ifrag, " stage=", "after-jo-loop"
        flush(6)
        write(*,'(1x,a,i0,a,i0,a,i0,a,l1,a,1pe12.4,a,1pe12.4)') &
          "        overlap fragment done: rank=", dg_frag%id, " id_frag=", dg_frag%id_frag, &
          " ifrag=", ifrag, " root=", dg_frag%is_frag_root, " self=", time_self_integral - frag_self_start, &
          " halo=", time_halo_integral - frag_halo_start
        flush(6)
      end do
    end do

    write(*,*) "        overlap reduce begin"
    flush(6)
    call cpu_time(t0)
    mat_size = dg_frag%n_mat_max * dg_frag%n_mat_max * dg_frag%nspin
    allocate(S_flat(mat_size), S_tmp_flat(mat_size))
    log_reduce_diag = (dg_frag%ifrag_group == 42 .or. dg_frag%ifrag_group == 43)
    sample_n = min(8, mat_size)
    offset_flat = 1
    do ispin = 1, dg_frag%nspin
      do jj = 1, dg_frag%n_mat_max
        do ii = 1, dg_frag%n_mat_max
          S_flat(offset_flat) = dg_frag%S_mat(ii, jj, ispin)
          offset_flat = offset_flat + 1
        end do
      end do
    end do
    if (offset_flat /= mat_size + 1) then
      write(*,*) "[FATAL] overlap reduce pack mismatch: rank=", dg_frag%id, " id_frag=", dg_frag%id_frag, &
        " ifrag_group=", dg_frag%ifrag_group, " offset=", offset_flat, " mat_size=", mat_size
      stop 1
    end if
    if (log_reduce_diag) then
      nan_count = count(S_flat /= S_flat)
      sum_abs_flat = sum(abs(S_flat))
      max_abs_flat = maxval(abs(S_flat))
      write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,l1,a,i0,a,i0,a,i0,a,i0,a,1pe12.4,a,1pe12.4,a,i0)') &
        "        overlap reduce stats: stage=", "before-frag-comm", " rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, " root=", dg_frag%is_frag_root, &
        " isize_frag=", dg_frag%isize_frag, " isize=", dg_frag%isize, " mat_size=", mat_size, &
        " sample_n=", sample_n, " sumabs=", sum_abs_flat, " maxabs=", max_abs_flat, " nan=", nan_count
      if (sample_n > 0) then
        write(*,'(1x,a,8(1pe12.4,1x))') "        overlap reduce sample: ", S_flat(1:sample_n)
      end if
      flush(6)
    end if
    call comm_summation(S_flat, S_tmp_flat, mat_size, dg_frag%icomm_frag)
    if (log_reduce_diag) then
      nan_count = count(S_tmp_flat /= S_tmp_flat)
      sum_abs_flat = sum(abs(S_tmp_flat))
      max_abs_flat = maxval(abs(S_tmp_flat))
      write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,1pe12.4,a,1pe12.4,a,i0)') &
        "        overlap reduce stats: stage=", "after-frag-comm", " rank=", dg_frag%id, &
        " id_frag=", dg_frag%id_frag, " ifrag_group=", dg_frag%ifrag_group, &
        " sumabs=", sum_abs_flat, " maxabs=", max_abs_flat, " nan=", nan_count
      if (sample_n > 0) then
        write(*,'(1x,a,8(1pe12.4,1x))') "        overlap reduce sample: ", S_tmp_flat(1:sample_n)
      end if
      flush(6)
    end if
    offset_flat = 1
    do ispin = 1, dg_frag%nspin
      do jj = 1, dg_frag%n_mat_max
        do ii = 1, dg_frag%n_mat_max
          dg_frag%S_mat(ii, jj, ispin) = S_tmp_flat(offset_flat)
          offset_flat = offset_flat + 1
        end do
      end do
    end do
    if (offset_flat /= mat_size + 1) then
      write(*,*) "[FATAL] overlap reduce unpack mismatch (frag): rank=", dg_frag%id, " id_frag=", dg_frag%id_frag, &
        " ifrag_group=", dg_frag%ifrag_group, " offset=", offset_flat, " mat_size=", mat_size
      stop 1
    end if
    if (.not. dg_frag%is_frag_root) dg_frag%S_mat = 0.0d0
    deallocate(S_flat, S_tmp_flat)
    call cpu_time(t1)
    time_reduce_total = time_reduce_total + (t1 - t0)
    write(*,'(1x,a,1pe12.4)') "        overlap reduce done time=", time_reduce_total
    flush(6)

    do ispin = 1, dg_frag%nspin
      do ii = 1, dg_frag%n_mat_max
        if (dg_frag%S_mat(ii, ii, ispin) < 1.0d-12) dg_frag%S_mat(ii, ii, ispin) = 1.0d0
        do jj = ii + 1, dg_frag%n_mat_max
          savg = 0.5d0 * (dg_frag%S_mat(ii, jj, ispin) + dg_frag%S_mat(jj, ii, ispin))
          dg_frag%S_mat(ii, jj, ispin) = savg
          dg_frag%S_mat(jj, ii, ispin) = savg
        end do
      end do
    end do

    do ispin = 1, dg_frag%nspin
      n_eval = dg_frag%n_mat(ispin)
      if (n_eval <= 1) cycle
      if (comm_is_root(dg_frag%id)) then
        allocate(S_eval(n_eval, n_eval), eigvals(n_eval))
        S_eval(:, :) = dg_frag%S_mat(1:n_eval, 1:n_eval, ispin)
        lwork = -1
        call dsyev('N', 'U', n_eval, S_eval, n_eval, eigvals, work_query, lwork, info_eig)
        if (info_eig == 0) then
          lwork = max(1, int(work_query(1)))
          allocate(eig_work(lwork))
          call dsyev('N', 'U', n_eval, S_eval, n_eval, eigvals, eig_work, lwork, info_eig)
          deallocate(eig_work)
        end if
        if (info_eig == 0) then
          s_min = eigvals(1)
          s_max = eigvals(n_eval)
          if (abs(s_min) > 1.0d-14) then
            cond_est = abs(s_max / s_min)
          else
            cond_est = huge(1.0d0)
          end if
          write(*,'(a,i0,a,1pe12.4,a,1pe12.4,a,1pe12.4)') &
            '        [S-eig] spin=', ispin, ' min=', s_min, ' max=', s_max, ' cond~', cond_est
        else
          write(*,'(a,i0,a,i0)') '        [S-eig] spin=', ispin, ' dsyev failed, info=', info_eig
        end if
        deallocate(S_eval, eigvals)
      end if
    end do

    ! Default propagation overlap equals the raw fragment overlap.
    dg_frag%S_mat_prop(:, :, :) = dg_frag%S_mat(:, :, :)
    dg_frag%has_global_overlap_copy = .false.
    dg_frag%overlap_prop_root_authoritative = .true.

  end subroutine calculate_overlap_matrix
