!=======================================================================
  ! Calculate Hamiltonian matrix in fragment basis
  !=======================================================================
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
    real(8) :: hvol
    real(8) :: max_p
    real(8) :: Ac_zero(3)
    integer :: is(3), ie(3)
    real(8), allocatable :: T_phi(:,:,:)  ! Kinetic energy operator applied to basis
    real(8), allocatable :: H_phi(:,:,:)  ! Hamiltonian-applied field H|phi_j> = T|phi_j> + V|phi_j>
    real(8), allocatable :: V_total(:,:,:)  ! Total potential V = Vpsl + Vh + Vxc
    real(8), allocatable :: partial_t(:), partial_h(:), reduced_t(:), reduced_h(:)
    
    if (.not. dg_frag%has_real_space_basis) then
      if (.not. allocated(dg_frag%H_mat)) then
        allocate(dg_frag%H_mat(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
      end if
      if (.not. allocated(dg_frag%H_mat_c)) then
        allocate(dg_frag%H_mat_c(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
      end if
      dg_frag%H_mat = 0.0d0
      dg_frag%H_mat_c = (0.0d0, 0.0d0)
      return
    end if
    
    if (comm_is_root(dg_frag%id)) then
      write(*,*)
      write(*,*) "=== Preparing Hamiltonian Matrix ==="
    end if
    
    ! Step 1: Calculate momentum matrix elements (transition moments)
    ! Required for velocity gauge A·p coupling
    if (.not. allocated(dg_frag%momentum_mat) .or. .not. allocated(dg_frag%momentum_mat_c)) then
      if (comm_is_root(dg_frag%id)) then
        write(*,*) "  [1/3] Calculating momentum matrix elements (p_ij)..."
        write(*,*) "        Using 4th-order finite difference stencil"
      end if
      call calculate_momentum_matrix(dg_frag, system, mg, stencil)
      call calculate_overlap_matrix(dg_frag, system, mg)
      if (comm_is_root(dg_frag%id)) then
        write(*,*) "        Momentum matrix calculated (for A·p coupling)"
        write(*,*) "        Overlap matrix S calculated (for generalized propagation)"
      end if
    else
      if (comm_is_root(dg_frag%id)) then
        write(*,*) "  [1/3] Momentum matrix already available"
      end if
      if (.not. allocated(dg_frag%S_mat) .or. .not. allocated(dg_frag%S_mat_c)) then
        call calculate_overlap_matrix(dg_frag, system, mg)
      end if
    end if
    
    ! Step 2: Allocate Hamiltonian matrix
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "  [2/3] Constructing Hamiltonian matrix H = T + V..."
    end if
    
    ! Only allocate if not already allocated (may be pre-allocated in initialization)
    if (.not. allocated(dg_frag%H_mat)) then
      allocate(dg_frag%H_mat(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
    end if
    if (.not. allocated(dg_frag%H_mat_c)) then
      allocate(dg_frag%H_mat_c(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
    end if
    dg_frag%H_mat = 0.0d0
    dg_frag%H_mat_c = (0.0d0, 0.0d0)
    
    ! Exchange halo regions between fragments before stencil operations
    ! This ensures accurate Laplacian calculation at fragment boundaries
    call exchange_phi_frag_halo(dg_frag)
    
    hvol = system%hvol
    is = lg%is
    ie = lg%ie
    
    ! Allocate work arrays
    allocate(T_phi(is(1):ie(1), is(2):ie(2), is(3):ie(3)))
    allocate(H_phi(is(1):ie(1), is(2):ie(2), is(3):ie(3)))
    allocate(V_total(is(1):ie(1), is(2):ie(2), is(3):ie(3)))
    
    ! Construct total potential: V = Vpsl + Vh + Vxc
    ! Note: This is used for initial H_mat calculation
    do ispin = 1, system%nspin
      call build_total_potential_grid(lg, Vh, Vxc(ispin), Vpsl, V_total)
      
      ! Loop over fragments assigned to this rank
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        
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

    dg_frag%H_mat_c(:, :, :) = cmplx(dg_frag%H_mat(:, :, :), 0.0d0, kind=8)
    
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
    
    deallocate(T_phi, H_phi, V_total)
    
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
    real(8), intent(out) :: T_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    real(8), intent(out) :: H_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
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
          H_phi(gx, gy, gz) = H_phi(gx, gy, gz) + V_total(gx, gy, gz) * dg_frag%phi_frag(lx, ly, lz, jo, i_local)
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
    real(8), intent(in) :: field(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    real(8), intent(in) :: hvol
    real(8), intent(out) :: integral
    real(8) :: partial
    integer :: lx, ly, lz, gx, gy, gz
    integer :: iorg(3), ndom(3), loc_s(3), loc_e(3)

    iorg(:) = dg_frag%ixyz_frag(:, ifrag)
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    call get_fragment_local_range(dg_frag, ndom, loc_s, loc_e)
    partial = 0.0d0
    do lz = loc_s(3), loc_e(3)
      gz = iorg(3) + lz - 1
      do ly = loc_s(2), loc_e(2)
        gy = iorg(2) + ly - 1
        !$omp simd reduction(+:partial)
        do lx = loc_s(1), loc_e(1)
          gx = iorg(1) + lx - 1
          partial = partial + dg_frag%phi_frag(lx, ly, lz, io, i_local) * field(gx, gy, gz) * hvol
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
    real(8),                intent(out) :: T_phi(mg%is(1):mg%ie(1), &
                                                 mg%is(2):mg%ie(2), &
                                                 mg%is(3):mg%ie(3))
    
    integer :: ix, iy, iz, lx, ly, lz, gx, gy, gz, ifrag
    real(8) :: v, lap0
    real(8) :: lapt(4,3)
    integer :: is(3), ie(3), iorg(3), ndom(3), loc_s(3), loc_e(3)
    
    ! Extract stencil coefficients
    lap0 = stencil%coef_lap0
    lapt = stencil%coef_lap
    is = mg%is
    ie = mg%ie
    ifrag = dg_frag%ifrag_start + i_local - 1
    iorg(:) = dg_frag%ixyz_frag(:, ifrag)
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
      gz = iorg(3) + lz - 1
      do ly = loc_s(2), loc_e(2)
        gy = iorg(2) + ly - 1
        do lx = loc_s(1), loc_e(1)
          gx = iorg(1) + lx - 1
          
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
          T_phi(gx, gy, gz) = lap0 * dg_frag%phi_frag(lx, ly, lz, jo, i_local) - 0.5d0 * v
          
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
    use communication, only: comm_is_root, comm_summation, comm_get_groupinfo
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_stencil),        intent(in)    :: stencil
    
    integer :: ifrag, i_local, ispin, io, jo, idir
    integer :: ix, iy, iz, is(3), ie(3), i_halo, jfrag, n_basis_halo, ig_row, ig_col, ig_i, ig_j, l(3), d(3)
    integer :: lx, ly, lz, gx, gy, gz, iorg(3), ndom(3), loc_s(3), loc_e(3), halo_s(3), halo_e(3)
    real(8) :: hvol, integral
    real(8) :: max_p
    real(8), allocatable :: grad_phi(:,:,:,:)  ! gradient of basis function (x,y,z components)
    
    if (.not. dg_frag%has_real_space_basis) return
    
    if (comm_is_root(dg_frag%id)) then
      write(*,*) "        Computing transition moments: <φ_i|∇|φ_j>"
    end if
    
    ! Allocate momentum matrix: (3 directions, n_mat_max x n_mat_max, nspin)
    ! Momentum matrix elements for vector potential coupling: p_ij = <phi_i|p|phi_j>
    ! In velocity gauge: H(t) = H_0 - i*A(t)·∇ + A(t)^2/2
    ! The A·p term couples to momentum matrix elements
    ! The A^2/2 term is diagonal (diamagnetic contribution)
    if (allocated(dg_frag%momentum_mat)) deallocate(dg_frag%momentum_mat)
    if (allocated(dg_frag%momentum_mat_c)) deallocate(dg_frag%momentum_mat_c)
    allocate(dg_frag%momentum_mat(3, dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
    allocate(dg_frag%momentum_mat_c(3, dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
    dg_frag%momentum_mat = 0.0d0    
    dg_frag%momentum_mat_c = (0.0d0, 0.0d0)
    is = mg%is
    ie = mg%ie
    hvol = system%hvol
    
    ! Exchange halo regions before stencil operations
    call exchange_phi_frag_halo(dg_frag)
    
    ! Loop over spin
    do ispin = 1, system%nspin
      ! Loop over local fragments
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        iorg(:) = dg_frag%ixyz_frag(:, ifrag)
        ndom(:) = dg_frag%nxyz_domain(:, ifrag)
        call get_fragment_local_range(dg_frag, ndom, loc_s, loc_e)
        
        ! Loop over basis functions in fragment j (ket side)
        ! Note: Each thread allocates its own grad_phi to avoid race conditions
        !$omp parallel do private(jo,io,idir,ix,iy,iz,integral,grad_phi) collapse(1)
        do jo = 1, dg_frag%n_basis(ifrag, ispin)
          
          ! Allocate thread-local workspace for gradient
          allocate(grad_phi(is(1):ie(1), is(2):ie(2), is(3):ie(3), 3))
          
          ! Calculate gradient of phi_j using stencil
          call apply_gradient_to_basis_ops(dg_frag, i_local, jo, mg, stencil, grad_phi)
          
          ! Loop over basis functions in fragment i (bra side)
          do io = 1, dg_frag%n_basis(ifrag, ispin)
            ig_i = dg_frag%index_basis(io, ifrag, ispin)
            ig_j = dg_frag%index_basis(jo, ifrag, ispin)
            if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
            if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle
            
            ! Loop over x, y, z directions
            do idir = 1, 3
              
              ! Compute matrix element: p_ij = ∫ φ_i(r) * (∂φ_j/∂dir) dr
              ! Note: momentum operator uses p = -i∇; the -i is applied in time evolution
              integral = 0.0d0
              do lz = loc_s(3), loc_e(3)
                gz = iorg(3) + lz - 1
                do ly = loc_s(2), loc_e(2)
                  gy = iorg(2) + ly - 1
                  !$omp simd reduction(+:integral)
                  do lx = loc_s(1), loc_e(1)
                    gx = iorg(1) + lx - 1
                    integral = integral + &
                      dg_frag%phi_frag(lx, ly, lz, io, i_local) * &
                      grad_phi(gx, gy, gz, idir) * hvol
                  end do
                end do
              end do
              
              ! Store in global momentum matrix
              dg_frag%momentum_mat(idir, ig_i, ig_j, ispin) = integral
              
            end do  ! idir
          end do  ! io
          
          ! Deallocate thread-local workspace
          deallocate(grad_phi)
          
        end do  ! jo
        !$omp end parallel do
        
        ! Inter-fragment momentum blocks from halo regions.
        ! Use half-weight updates to avoid double counting since each neighbor pair
        ! appears from both destination fragments across the full communicator.
        !$omp parallel do private(jo,i_halo,jfrag,n_basis_halo,l,d,ig_col,io,ig_row,idir,ix,iy,iz,integral,grad_phi) collapse(1)
        do jo = 1, dg_frag%n_basis(ifrag, ispin)
          allocate(grad_phi(is(1):ie(1), is(2):ie(2), is(3):ie(3), 3))
          call apply_gradient_to_basis_ops(dg_frag, i_local, jo, mg, stencil, grad_phi)
          ig_col = dg_frag%index_basis(jo, ifrag, ispin)
          if (ig_col < 1 .or. ig_col > dg_frag%n_mat_max) then
            deallocate(grad_phi)
            cycle
          end if

          do i_halo = 1, dg_frag%n_halo
            if (dg_frag%halo(i_halo)%ifrag_dst /= ifrag) cycle
            jfrag = dg_frag%halo(i_halo)%ifrag_src
            n_basis_halo = dg_frag%n_basis(jfrag, ispin)
            l = dg_frag%halo(i_halo)%length
            d = dg_frag%halo(i_halo)%dsp_send
            halo_s(:) = max(loc_s(:), d(:) + 1)
            halo_e(:) = min(loc_e(:), d(:) + l(:))
            if (any(halo_s(:) > halo_e(:))) cycle

            do io = 1, n_basis_halo
              ig_row = dg_frag%index_basis(io, jfrag, ispin)
              if (ig_row < 1 .or. ig_row > dg_frag%n_mat_max) cycle
              do idir = 1, 3
                integral = 0.0d0
                do lz = halo_s(3), halo_e(3)
                  gz = iorg(3) + lz - 1
                  iz = lz - d(3)
                  do ly = halo_s(2), halo_e(2)
                    gy = iorg(2) + ly - 1
                    iy = ly - d(2)
                    !$omp simd reduction(+:integral)
                    do lx = halo_s(1), halo_e(1)
                      gx = iorg(1) + lx - 1
                      ix = lx - d(1)
                      integral = integral + &
                        dg_frag%halo(i_halo)%buf_recv(ix, iy, iz, io, 1) * &
                        grad_phi(gx, gy, gz, idir) * hvol
                    end do
                  end do
                end do

                dg_frag%momentum_mat(idir, ig_row, ig_col, ispin) = &
                  dg_frag%momentum_mat(idir, ig_row, ig_col, ispin) + 0.5d0 * integral
                dg_frag%momentum_mat(idir, ig_col, ig_row, ispin) = &
                  dg_frag%momentum_mat(idir, ig_col, ig_row, ispin) - 0.5d0 * integral
              end do
            end do
          end do

          deallocate(grad_phi)
        end do
        !$omp end parallel do

      end do  ! ifrag
    end do  ! ispin
    
    ! CRITICAL: MPI aggregation of momentum matrix
    ! Each rank computed elements only for its assigned fragments
    ! Sum across all ranks to get complete global momentum matrix
    block
      real(8), allocatable :: momentum_mat_flat(:), momentum_mat_tmp_flat(:)
      integer :: mat_size
      mat_size = 3 * dg_frag%n_mat_max * dg_frag%n_mat_max * dg_frag%nspin
      allocate(momentum_mat_flat(mat_size), momentum_mat_tmp_flat(mat_size))

      momentum_mat_flat = reshape(dg_frag%momentum_mat, [mat_size])
      call comm_summation(momentum_mat_flat, momentum_mat_tmp_flat, mat_size, dg_frag%icomm_frag)
      dg_frag%momentum_mat = reshape(momentum_mat_tmp_flat, [3, dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin])
      if (.not. dg_frag%is_frag_root) dg_frag%momentum_mat = 0.0d0
      
      momentum_mat_flat = reshape(dg_frag%momentum_mat, [mat_size])
      call comm_summation(momentum_mat_flat, momentum_mat_tmp_flat, mat_size, dg_frag%icomm)
      dg_frag%momentum_mat = reshape(momentum_mat_tmp_flat, [3, dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin])
      
      deallocate(momentum_mat_flat, momentum_mat_tmp_flat)
    end block

    ! Enforce anti-symmetry for real-basis gradient operator: P = -P^T.
    ! This suppresses accumulation of discretization asymmetry at fragment boundaries.
    block
      integer :: ii, jj
      real(8) :: pavg
      do ispin = 1, system%nspin
        do idir = 1, 3
          do ii = 1, dg_frag%n_mat_max
            dg_frag%momentum_mat(idir, ii, ii, ispin) = 0.0d0
            do jj = ii + 1, dg_frag%n_mat_max
              pavg = 0.5d0 * (dg_frag%momentum_mat(idir, ii, jj, ispin) - dg_frag%momentum_mat(idir, jj, ii, ispin))
              dg_frag%momentum_mat(idir, ii, jj, ispin) = pavg
              dg_frag%momentum_mat(idir, jj, ii, ispin) = -pavg
            end do
          end do
        end do
      end do
    end block

    max_p = maxval(abs(dg_frag%momentum_mat))
    dg_frag%momentum_mat_c = cmplx(dg_frag%momentum_mat, 0.0d0, kind=8)
    if (comm_is_root(dg_frag%id)) then
      write(*,'(a,1pe12.4)') "        Max |momentum_mat|: ", max_p
      write(*,'(a,i0,a,i0,a)') "        Total matrix elements: ", &
                               3 * dg_frag%n_mat_max * dg_frag%n_mat_max * system%nspin, &
                               " (", 3, " directions × basis states × spin)"
    end if
    
  end subroutine calculate_momentum_matrix

  !=======================================================================
  ! Calculate overlap matrix in DG basis (S_ij = <phi_i|phi_j>)
  !=======================================================================
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
    integer :: mat_size, n_eval, lwork, info_eig
    real(8) :: hvol, integral, savg, s_min, s_max, cond_est
    real(8) :: work_query(1)
    real(8), allocatable :: S_flat(:), S_tmp_flat(:)
    real(8), allocatable :: S_eval(:,:), eigvals(:), eig_work(:)

    if (.not. dg_frag%has_real_space_basis) return
    if (.not. allocated(dg_frag%index_basis) .or. .not. allocated(dg_frag%n_mat)) return

    if (allocated(dg_frag%S_mat)) deallocate(dg_frag%S_mat)
    if (allocated(dg_frag%S_mat_prop)) deallocate(dg_frag%S_mat_prop)
    if (allocated(dg_frag%S_mat_c)) deallocate(dg_frag%S_mat_c)
    if (allocated(dg_frag%S_mat_prop_c)) deallocate(dg_frag%S_mat_prop_c)
    if (.not. allocated(dg_frag%S_mat)) then
      allocate(dg_frag%S_mat(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
    end if
    if (.not. allocated(dg_frag%S_mat_prop)) then
      allocate(dg_frag%S_mat_prop(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
    end if
    if (.not. allocated(dg_frag%S_mat_c)) then
      allocate(dg_frag%S_mat_c(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
    end if
    if (.not. allocated(dg_frag%S_mat_prop_c)) then
      allocate(dg_frag%S_mat_prop_c(dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin))
    end if
    dg_frag%S_mat = 0.0d0
    dg_frag%S_mat_prop = 0.0d0
    dg_frag%S_mat_c = (0.0d0, 0.0d0)
    dg_frag%S_mat_prop_c = (0.0d0, 0.0d0)

    is = mg%is
    ie = mg%ie
    hvol = system%hvol

    call exchange_phi_frag_halo(dg_frag)

    do ispin = 1, system%nspin
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        iorg(:) = dg_frag%ixyz_frag(:, ifrag)
        ndom(:) = dg_frag%nxyz_domain(:, ifrag)
        call get_fragment_local_range(dg_frag, ndom, loc_s, loc_e)

        do jo = 1, dg_frag%n_basis(ifrag, ispin)
          ig_col = dg_frag%index_basis(jo, ifrag, ispin)
          if (ig_col < 1 .or. ig_col > dg_frag%n_mat_max) cycle

          do io = 1, dg_frag%n_basis(ifrag, ispin)
            ig_row = dg_frag%index_basis(io, ifrag, ispin)
            if (ig_row < 1 .or. ig_row > dg_frag%n_mat_max) cycle
            integral = 0.0d0
            do lz = loc_s(3), loc_e(3)
              do ly = loc_s(2), loc_e(2)
                do lx = loc_s(1), loc_e(1)
                  integral = integral + dg_frag%phi_frag(lx, ly, lz, io, i_local) * &
                             dg_frag%phi_frag(lx, ly, lz, jo, i_local) * hvol
                end do
              end do
            end do
            dg_frag%S_mat(ig_row, ig_col, ispin) = integral
          end do

          do i_halo = 1, dg_frag%n_halo
            if (dg_frag%halo(i_halo)%ifrag_dst /= ifrag) cycle
            jfrag = dg_frag%halo(i_halo)%ifrag_src
            if (jfrag < 1) cycle
            n_basis_halo = dg_frag%n_basis(jfrag, ispin)
            l = dg_frag%halo(i_halo)%length
            d = dg_frag%halo(i_halo)%dsp_send
            halo_s(:) = max(loc_s(:), d(:) + 1)
            halo_e(:) = min(loc_e(:), d(:) + l(:))
            if (any(halo_s(:) > halo_e(:))) cycle

            do io = 1, n_basis_halo
              ig_row = dg_frag%index_basis(io, jfrag, ispin)
              if (ig_row < 1 .or. ig_row > dg_frag%n_mat_max) cycle
              integral = 0.0d0
              do lz = halo_s(3), halo_e(3)
                iz = lz - d(3)
                do ly = halo_s(2), halo_e(2)
                  iy = ly - d(2)
                  do lx = halo_s(1), halo_e(1)
                    ix = lx - d(1)
                    integral = integral + dg_frag%halo(i_halo)%buf_recv(ix, iy, iz, io, 1) * &
                               dg_frag%phi_frag(lx, ly, lz, jo, i_local) * hvol
                  end do
                end do
              end do
              dg_frag%S_mat(ig_row, ig_col, ispin) = dg_frag%S_mat(ig_row, ig_col, ispin) + 0.5d0 * integral
              dg_frag%S_mat(ig_col, ig_row, ispin) = dg_frag%S_mat(ig_col, ig_row, ispin) + 0.5d0 * integral
            end do
          end do

        end do
      end do
    end do

    mat_size = dg_frag%n_mat_max * dg_frag%n_mat_max * dg_frag%nspin
    allocate(S_flat(mat_size), S_tmp_flat(mat_size))
    S_flat = reshape(dg_frag%S_mat, [mat_size])
    call comm_summation(S_flat, S_tmp_flat, mat_size, dg_frag%icomm_frag)
    dg_frag%S_mat = reshape(S_tmp_flat, [dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin])
    if (.not. dg_frag%is_frag_root) dg_frag%S_mat = 0.0d0

    S_flat = reshape(dg_frag%S_mat, [mat_size])
    call comm_summation(S_flat, S_tmp_flat, mat_size, dg_frag%icomm)
    dg_frag%S_mat = reshape(S_tmp_flat, [dg_frag%n_mat_max, dg_frag%n_mat_max, dg_frag%nspin])
    deallocate(S_flat, S_tmp_flat)

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
      if (comm_is_root(dg_frag%id)) then
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
      end if
      deallocate(S_eval, eigvals)
    end do

    ! Default propagation overlap equals the raw fragment overlap.
    dg_frag%S_mat_prop(:, :, :) = dg_frag%S_mat(:, :, :)
    dg_frag%S_mat_c(:, :, :) = cmplx(dg_frag%S_mat(:, :, :), 0.0d0, kind=8)
    dg_frag%S_mat_prop_c(:, :, :) = dg_frag%S_mat_c(:, :, :)

  end subroutine calculate_overlap_matrix
