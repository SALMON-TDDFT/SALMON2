!=======================================================================
  ! Calculate Hamiltonian matrix in fragment basis
  !=======================================================================
  !=======================================================================
  ! Calculate initial Hamiltonian matrix from basis functions
  !
  ! Includes halo (ghost cell) exchange for accurate boundary treatment.
  ! System boundaries use PERIODIC boundary conditions (full system is periodic).
  ! Fragment boundaries are handled via MPI communication between neighboring fragments.
  !=======================================================================
  subroutine calculate_hamiltonian_matrix(dg_frag, system, lg, mg, stencil, &
                                         Vh, Vxc, Vpsl, pp, ppg)
    use structures
    use communication, only: comm_is_root, comm_summation
    use parallelization, only: nproc_size_global
    use rt_dg_plane_wave, only: prepare_mixed_basis_startup
    use rt_dg_fragment_ops, only: copy_matrix_blocks_to_complex_dense, copy_matrix_blocks_metric_to_complex_dense, &
      symmetrize_real_matrix_blocks
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
    complex(8) :: integral_t, integral_h
    real(8) :: max_p
    real(8) :: Ac_zero(3)
    integer :: n_metric
    integer :: is(3), ie(3)
    complex(8), allocatable :: T_phi(:,:,:)   ! Kinetic energy operator applied to basis
    complex(8), allocatable :: H_phi(:,:,:)   ! Hamiltonian-applied field H|phi_j> = T|phi_j> + V|phi_j>
    real(8), allocatable :: V_total(:,:,:)  ! Total potential V = Vpsl + Vh + Vxc
    complex(8), allocatable :: H_metric_ref(:,:)
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
      if (.not. allocated(dg_frag%S_mat)) call calculate_overlap_matrix(dg_frag, system, mg)
    end if
    
    ! Step 2: Allocate Hamiltonian matrix
    write(*,*) "  [2/3] Constructing Hamiltonian matrix H = T + V..."
    
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
    is = mg%is
    ie = mg%ie
    
    ! Allocate work arrays
    allocate(T_phi(is(1):ie(1), is(2):ie(2), is(3):ie(3)))
    allocate(H_phi(is(1):ie(1), is(2):ie(2), is(3):ie(3)))
    allocate(V_total(is(1):ie(1), is(2):ie(2), is(3):ie(3)))
    
    ! Construct total potential: V = Vpsl + Vh + Vxc
    ! Note: This is used for initial H_mat calculation
    do ispin = 1, system%nspin
      call build_total_potential_grid(mg, Vh, Vxc(ispin), Vpsl, V_total)
      
      ! Loop over fragments assigned to this rank
      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1
        
        ! Calculate Hamiltonian matrix elements for this fragment
        ! H_ij = <φ_i | T + V | φ_j> = T_ij + V_ij
        nbf = dg_frag%n_basis(ifrag, ispin)
        do jo = 1, nbf
          ig_j = dg_frag%index_basis(jo, ifrag, ispin)
          if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle

          call build_hpsi_for_basis(dg_frag, ifrag, i_local, jo, mg, stencil, V_total, T_phi, H_phi)

          ! Calculate matrix elements with all φ_i
          !$omp parallel do private(io, integral_t, integral_h, ig_i)
          do io = 1, nbf
            ig_i = dg_frag%index_basis(io, ifrag, ispin)
            if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle

            ! Kinetic energy matrix element: T_ij = ∫ φ_i (T|φ_j>) dr
            call integrate_basis_with_field(dg_frag, ifrag, i_local, io, mg, T_phi, hvol, integral_t)

            ! Store kinetic part
            dg_frag%H_mat_kinetic(ig_i, ig_j, ispin) = real(integral_t, kind=8)
            call integrate_basis_with_field(dg_frag, ifrag, i_local, io, mg, H_phi, hvol, integral_h)
            dg_frag%H_mat(ig_i, ig_j, ispin) = real(integral_h, kind=8)

          end do
          !$omp end parallel do

        end do  ! jo loop
          
        
      end do  ! ifrag loop
      
    end do  ! ispin loop
    
    call init_matrix_blocks(dg_frag, dg_frag%H_mat_blocks, dg_frag%H_block_map, dg_frag%n_H_blocks)
    call sync_dense_matrix_to_blocks(dg_frag, dg_frag%H_mat, dg_frag%H_mat_blocks, dg_frag%H_block_map)
    call init_matrix_blocks(dg_frag, dg_frag%H_mat_kinetic_blocks, dg_frag%H_block_map, dg_frag%n_H_blocks)
    call sync_dense_matrix_to_blocks(dg_frag, dg_frag%H_mat_kinetic, dg_frag%H_mat_kinetic_blocks, dg_frag%H_block_map)

    ! Global Hamiltonian aggregation via fragment blocks to avoid a monolithic dense allreduce.
    call reduce_matrix_blocks(dg_frag, dg_frag%H_mat_blocks, "hmat-soi", dg_frag%icomm)
    call reduce_matrix_blocks(dg_frag, dg_frag%H_mat_kinetic_blocks, "hmat-kinetic-soi", dg_frag%icomm)
    call symmetrize_real_matrix_blocks(dg_frag, dg_frag%H_mat_blocks)
    call symmetrize_real_matrix_blocks(dg_frag, dg_frag%H_mat_kinetic_blocks)
    dg_frag%H_mat_c(:, :, :) = (0.0d0, 0.0d0)
    do ispin = 1, dg_frag%nspin
      call copy_matrix_blocks_to_complex_dense(dg_frag, dg_frag%H_mat_blocks, ispin, dg_frag%H_mat_c(1:dg_frag%n_mat_max, 1:dg_frag%n_mat_max, ispin))
    end do
    if (allocated(dg_frag%H_mat)) deallocate(dg_frag%H_mat)
    if (allocated(dg_frag%H_mat_kinetic)) deallocate(dg_frag%H_mat_kinetic)
    
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
      call prepare_mixed_basis_startup(dg_frag, system, Vh, Vxc, Vpsl, Ac_zero)
      dg_frag%coef_new(:, :, :) = dg_frag%coef(:, :, :)
    end if

    ! Initialize field-free reference Hamiltonian for adaptive-basis metric.
    if (allocated(dg_frag%H_mat_old)) then
      dg_frag%H_mat_old(:, :, :) = (0.0d0, 0.0d0)
      n_metric = min(dg_frag%nstate_frag, size(dg_frag%H_mat_old, 1), size(dg_frag%H_mat_old, 2))
      if (n_metric > 0) then
        allocate(H_metric_ref(n_metric, n_metric))
      end if
      do ispin = 1, min(dg_frag%nspin, size(dg_frag%H_mat_old,3))
        if (n_metric <= 0) cycle
        H_metric_ref(:, :) = (0.0d0, 0.0d0)
        call copy_matrix_blocks_metric_to_complex_dense(dg_frag, dg_frag%H_mat_blocks, ispin, n_metric, H_metric_ref)
        dg_frag%H_mat_old(1:n_metric, 1:n_metric, ispin) = H_metric_ref(:, :)
      end do
      if (allocated(H_metric_ref)) deallocate(H_metric_ref)
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
    complex(8), intent(out) :: T_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    complex(8), intent(out) :: H_phi(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    integer :: lx, ly, lz, gx, gy, gz
    integer :: iorg(3), ndom(3)

    call apply_kinetic_to_basis(dg_frag, i_local, jo, mg, stencil, T_phi)
    H_phi(:, :, :) = T_phi(:, :, :)

    iorg(:) = dg_frag%ixyz_frag(:, ifrag)
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
!$omp parallel do private(lz, ly, lx, gx, gy, gz) schedule(static)
    do lz = 1, ndom(3)
      gz = iorg(3) + lz - 1
      do ly = 1, ndom(2)
        gy = iorg(2) + ly - 1
!$omp simd private(gx)
        do lx = 1, ndom(1)
          gx = iorg(1) + lx - 1
          if (allocated(dg_frag%phi_frag_c)) then
            H_phi(gx, gy, gz) = H_phi(gx, gy, gz) + V_total(gx, gy, gz) * dg_frag%phi_frag_c(lx, ly, lz, jo, i_local)
          else
            H_phi(gx, gy, gz) = H_phi(gx, gy, gz) + V_total(gx, gy, gz) * cmplx(dg_frag%phi_frag(lx, ly, lz, jo, i_local), 0.0d0, kind=8)
          end if
        end do
      end do
    end do
!$omp end parallel do
  end subroutine build_hpsi_for_basis

  !=======================================================================
  ! Integrate one bra basis function against a real-space field
  !   integral = <phi_io | field>
  !=======================================================================
  subroutine integrate_basis_with_field(dg_frag, ifrag, i_local, io, mg, field, hvol, integral)
    use structures
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: ifrag, i_local, io
    type(s_rgrid), intent(in) :: mg
    complex(8), intent(in) :: field(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3))
    real(8), intent(in) :: hvol
    complex(8), intent(out) :: integral
    integer :: lx, ly, lz, gx, gy, gz
    integer :: iorg(3), ndom(3)

    iorg(:) = dg_frag%ixyz_frag(:, ifrag)
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    integral = (0.0d0, 0.0d0)
    do lz = 1, ndom(3)
      gz = iorg(3) + lz - 1
      do ly = 1, ndom(2)
        gy = iorg(2) + ly - 1
        do lx = 1, ndom(1)
          gx = iorg(1) + lx - 1
          if (allocated(dg_frag%phi_frag_c)) then
            integral = integral + conjg(dg_frag%phi_frag_c(lx, ly, lz, io, i_local)) * field(gx, gy, gz) * hvol
          else
            integral = integral + cmplx(dg_frag%phi_frag(lx, ly, lz, io, i_local), 0.0d0, kind=8) * field(gx, gy, gz) * hvol
          end if
        end do
      end do
    end do
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
    complex(8),             intent(out) :: T_phi(mg%is(1):mg%ie(1), &
                                                 mg%is(2):mg%ie(2), &
                                                 mg%is(3):mg%ie(3))
    
    integer :: lx, ly, lz, gx, gy, gz, ifrag
    complex(8) :: v
    real(8) :: lap0
    real(8) :: lapt(4,3)
    integer :: is(3), ie(3), iorg(3), ndom(3)
    
    ! Extract stencil coefficients
    lap0 = stencil%coef_lap0
    lapt = stencil%coef_lap
    is = mg%is
    ie = mg%ie
    ifrag = dg_frag%ifrag_start + i_local - 1
    iorg(:) = dg_frag%ixyz_frag(:, ifrag)
    ndom(:) = dg_frag%nxyz_domain(:, ifrag)
    
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
    
    T_phi = (0.0d0, 0.0d0)
    
!$omp parallel do private(lz, ly, lx, gx, gy, gz, v) schedule(static)
    do lz = 1, ndom(3)
      gz = iorg(3) + lz - 1
      do ly = 1, ndom(2)
        gy = iorg(2) + ly - 1
!$omp simd private(gx, v)
        do lx = 1, ndom(1)
          gx = iorg(1) + lx - 1
          
          ! Compute Laplacian using 4th-order finite difference
          ! Stencil accesses phi_frag(ix±4, iy±4, iz±4) which now includes halo
          if (allocated(dg_frag%phi_frag_c)) then
            v = lapt(1,1) * (dg_frag%phi_frag_c(lx+1, ly, lz, jo, i_local) + &
                             dg_frag%phi_frag_c(lx-1, ly, lz, jo, i_local)) + &
                lapt(2,1) * (dg_frag%phi_frag_c(lx+2, ly, lz, jo, i_local) + &
                             dg_frag%phi_frag_c(lx-2, ly, lz, jo, i_local)) + &
                lapt(3,1) * (dg_frag%phi_frag_c(lx+3, ly, lz, jo, i_local) + &
                             dg_frag%phi_frag_c(lx-3, ly, lz, jo, i_local)) + &
                lapt(4,1) * (dg_frag%phi_frag_c(lx+4, ly, lz, jo, i_local) + &
                             dg_frag%phi_frag_c(lx-4, ly, lz, jo, i_local))
            v = v + &
                lapt(1,2) * (dg_frag%phi_frag_c(lx, ly+1, lz, jo, i_local) + &
                             dg_frag%phi_frag_c(lx, ly-1, lz, jo, i_local)) + &
                lapt(2,2) * (dg_frag%phi_frag_c(lx, ly+2, lz, jo, i_local) + &
                             dg_frag%phi_frag_c(lx, ly-2, lz, jo, i_local)) + &
                lapt(3,2) * (dg_frag%phi_frag_c(lx, ly+3, lz, jo, i_local) + &
                             dg_frag%phi_frag_c(lx, ly-3, lz, jo, i_local)) + &
                lapt(4,2) * (dg_frag%phi_frag_c(lx, ly+4, lz, jo, i_local) + &
                             dg_frag%phi_frag_c(lx, ly-4, lz, jo, i_local))
            v = v + &
                lapt(1,3) * (dg_frag%phi_frag_c(lx, ly, lz+1, jo, i_local) + &
                             dg_frag%phi_frag_c(lx, ly, lz-1, jo, i_local)) + &
                lapt(2,3) * (dg_frag%phi_frag_c(lx, ly, lz+2, jo, i_local) + &
                             dg_frag%phi_frag_c(lx, ly, lz-2, jo, i_local)) + &
                lapt(3,3) * (dg_frag%phi_frag_c(lx, ly, lz+3, jo, i_local) + &
                             dg_frag%phi_frag_c(lx, ly, lz-3, jo, i_local)) + &
                lapt(4,3) * (dg_frag%phi_frag_c(lx, ly, lz+4, jo, i_local) + &
                             dg_frag%phi_frag_c(lx, ly, lz-4, jo, i_local))
            T_phi(gx, gy, gz) = lap0 * dg_frag%phi_frag_c(lx, ly, lz, jo, i_local) - 0.5d0 * v
          else
            v = cmplx(lapt(1,1) * (dg_frag%phi_frag(lx+1, ly, lz, jo, i_local) + &
                                   dg_frag%phi_frag(lx-1, ly, lz, jo, i_local)) + &
                      lapt(2,1) * (dg_frag%phi_frag(lx+2, ly, lz, jo, i_local) + &
                                   dg_frag%phi_frag(lx-2, ly, lz, jo, i_local)) + &
                      lapt(3,1) * (dg_frag%phi_frag(lx+3, ly, lz, jo, i_local) + &
                                   dg_frag%phi_frag(lx-3, ly, lz, jo, i_local)) + &
                      lapt(4,1) * (dg_frag%phi_frag(lx+4, ly, lz, jo, i_local) + &
                                   dg_frag%phi_frag(lx-4, ly, lz, jo, i_local)), 0.0d0, kind=8)
            v = v + cmplx(lapt(1,2) * (dg_frag%phi_frag(lx, ly+1, lz, jo, i_local) + &
                                       dg_frag%phi_frag(lx, ly-1, lz, jo, i_local)) + &
                          lapt(2,2) * (dg_frag%phi_frag(lx, ly+2, lz, jo, i_local) + &
                                       dg_frag%phi_frag(lx, ly-2, lz, jo, i_local)) + &
                          lapt(3,2) * (dg_frag%phi_frag(lx, ly+3, lz, jo, i_local) + &
                                       dg_frag%phi_frag(lx, ly-3, lz, jo, i_local)) + &
                          lapt(4,2) * (dg_frag%phi_frag(lx, ly+4, lz, jo, i_local) + &
                                       dg_frag%phi_frag(lx, ly-4, lz, jo, i_local)), 0.0d0, kind=8)
            v = v + cmplx(lapt(1,3) * (dg_frag%phi_frag(lx, ly, lz+1, jo, i_local) + &
                                       dg_frag%phi_frag(lx, ly, lz-1, jo, i_local)) + &
                          lapt(2,3) * (dg_frag%phi_frag(lx, ly, lz+2, jo, i_local) + &
                                       dg_frag%phi_frag(lx, ly, lz-2, jo, i_local)) + &
                          lapt(3,3) * (dg_frag%phi_frag(lx, ly, lz+3, jo, i_local) + &
                                       dg_frag%phi_frag(lx, ly, lz-3, jo, i_local)) + &
                          lapt(4,3) * (dg_frag%phi_frag(lx, ly, lz+4, jo, i_local) + &
                                       dg_frag%phi_frag(lx, ly, lz-4, jo, i_local)), 0.0d0, kind=8)
            T_phi(gx, gy, gz) = cmplx(lap0 * dg_frag%phi_frag(lx, ly, lz, jo, i_local), 0.0d0, kind=8) - 0.5d0 * v
          end if
          
        end do
      end do
    end do
!$omp end parallel do
    
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
    complex(8), allocatable :: uVpsi(:,:,:)  ! Projector overlaps (unnormalized): (nstate_frag, Nlma, nspin)
    complex(8) :: overlap_i, overlap_j, nlpp_contrib
    
    if (ppg%Nlma == 0) return  ! No non-local PP
    
    is = mg%is
    ie = mg%ie
    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    
    ! Allocate array for projector overlaps
    ! MEMORY OPTIMIZATION: Only store current fragment's data (removed ifrag_count dimension)
    allocate(uVpsi(dg_frag%nstate_frag, ppg%Nlma, nspin))
    uVpsi = (0.0d0, 0.0d0)
    
    ! Loop over fragments assigned to this rank
    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1
      iorg(:) = dg_frag%ixyz_frag(:, ifrag)
      ndom(:) = dg_frag%nxyz_domain(:, ifrag)
      
      ! Reset uVpsi for this fragment (memory reuse)
      uVpsi = (0.0d0, 0.0d0)
      
      do ispin = 1, nspin
        
        ! Calculate projector overlaps <φ_io|proj_ilma> for all basis and projectors
        ! OpenMP parallelization over basis functions
!$omp parallel do collapse(2) private(ilma, io, ia, j, ix, iy, iz, overlap_i)
        do ilma = 1, ppg%Nlma
          do io = 1, dg_frag%nstate_frag
            
            ia = ppg%ia_tbl(ilma)  ! Atom index for this projector
            
            ! Calculate <φ_io|proj_ilma> = Σ_j φ_io(r_j) * uV(r_j, ilma) * hvol
            ! NOTE: Store UNNORMALIZED overlap to avoid rinv_uvu^2 numerical error
            overlap_i = (0.0d0, 0.0d0)
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
                if (allocated(dg_frag%phi_frag_c)) then
                  overlap_i = overlap_i + &
                    conjg(dg_frag%phi_frag_c(lx, ly, lz, io, i_local)) * cmplx(ppg%uV(j, ilma), 0.0d0, kind=8) * hvol
                else
                  overlap_i = overlap_i + &
                    cmplx(dg_frag%phi_frag(lx, ly, lz, io, i_local), 0.0d0, kind=8) * cmplx(ppg%uV(j, ilma), 0.0d0, kind=8) * hvol
                end if
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
            nlpp_contrib = (0.0d0, 0.0d0)
            do ilma = 1, ppg%Nlma
              
              ! Get unnormalized overlaps
              overlap_i = uVpsi(io, ilma, ispin)
              overlap_j = uVpsi(jo, ilma, ispin)
              
              ! V_NL matrix element contribution from this projector
              ! Physical formula: <i|V_NL|j> = Σ_ilma <i|proj> * V_coeff * <proj|j>
              ! where V_coeff = rinv_uvu contains normalization and energy
              ! NUMERICAL ACCURACY: Apply rinv_uvu ONCE to avoid error amplification
              nlpp_contrib = nlpp_contrib + conjg(overlap_i) * overlap_j * cmplx(ppg%rinv_uvu(ilma), 0.0d0, kind=8)
              
            end do  ! ilma loop
            
            ! Add non-local PP contribution to Hamiltonian matrix
!$omp atomic update
            dg_frag%H_mat(ig_i, ig_j, ispin) = dg_frag%H_mat(ig_i, ig_j, ispin) + real(nlpp_contrib, kind=8)
            
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
    use rt_dg_fragment_types, only: momentum_block_info
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_stencil),        intent(in)    :: stencil
    
    integer :: ifrag, i_local, ispin, io, jo, idir, iblk
    integer :: ix, iy, iz, is(3), ie(3), i_halo, jfrag, n_basis_halo, ig_row, ig_col, ig_i, ig_j, l(3), d(3)
    integer :: lx, ly, lz, gx, gy, gz, iorg(3), ndom(3), loc_s(3), loc_e(3), halo_s(3), halo_e(3)
    integer :: nrow, ncol, ii, jj
    real(8) :: hvol
    complex(8) :: integral
    real(8) :: max_p
    real(8), allocatable :: grad_phi(:,:,:,:)  ! gradient of basis function (x,y,z components)
    type(momentum_block_info), allocatable :: momentum_blocks_re(:), momentum_blocks_im(:)
    integer, allocatable :: momentum_block_map_local(:,:)
    integer :: n_momentum_blocks_local
    
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
              integral = (0.0d0, 0.0d0)
              do lz = loc_s(3), loc_e(3)
                gz = iorg(3) + lz - 1
                do ly = loc_s(2), loc_e(2)
                  gy = iorg(2) + ly - 1
                  do lx = loc_s(1), loc_e(1)
                    gx = iorg(1) + lx - 1
                    if (allocated(dg_frag%phi_frag_c)) then
                      integral = integral + &
                        conjg(dg_frag%phi_frag_c(lx, ly, lz, io, i_local)) * &
                        cmplx(grad_phi(gx, gy, gz, idir), 0.0d0, kind=8) * hvol
                    else
                      integral = integral + &
                        cmplx(dg_frag%phi_frag(lx, ly, lz, io, i_local), 0.0d0, kind=8) * &
                        cmplx(grad_phi(gx, gy, gz, idir), 0.0d0, kind=8) * hvol
                    end if
                  end do
                end do
              end do
              
              ! Store in global momentum matrix
              dg_frag%momentum_mat_c(idir, ig_i, ig_j, ispin) = integral
              
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
                integral = (0.0d0, 0.0d0)
                do lz = halo_s(3), halo_e(3)
                  gz = iorg(3) + lz - 1
                  iz = lz - d(3)
                  do ly = halo_s(2), halo_e(2)
                    gy = iorg(2) + ly - 1
                    iy = ly - d(2)
                    do lx = halo_s(1), halo_e(1)
                      gx = iorg(1) + lx - 1
                      ix = lx - d(1)
                      if (allocated(dg_frag%halo(i_halo)%buf_recv_c)) then
                        integral = integral + &
                          conjg(dg_frag%halo(i_halo)%buf_recv_c(ix, iy, iz, io, 1)) * &
                          cmplx(grad_phi(gx, gy, gz, idir), 0.0d0, kind=8) * hvol
                      else
                        integral = integral + &
                          cmplx(dg_frag%halo(i_halo)%buf_recv(ix, iy, iz, io, 1), 0.0d0, kind=8) * &
                          cmplx(grad_phi(gx, gy, gz, idir), 0.0d0, kind=8) * hvol
                      end if
                    end do
                  end do
                end do

                dg_frag%momentum_mat_c(idir, ig_row, ig_col, ispin) = &
                  dg_frag%momentum_mat_c(idir, ig_row, ig_col, ispin) + 0.5d0 * integral
                dg_frag%momentum_mat_c(idir, ig_col, ig_row, ispin) = &
                  dg_frag%momentum_mat_c(idir, ig_col, ig_row, ispin) - 0.5d0 * conjg(integral)
              end do
            end do
          end do

          deallocate(grad_phi)
        end do
        !$omp end parallel do

      end do  ! ifrag
    end do  ! ispin
    
    call init_complex_momentum_blocks(dg_frag, momentum_blocks_re, momentum_blocks_im, &
                                      momentum_block_map_local, n_momentum_blocks_local)
    call sync_complex_dense_momentum_to_blocks(dg_frag, dg_frag%momentum_mat_c, momentum_blocks_re, momentum_blocks_im, &
                                               momentum_block_map_local)
    call reduce_complex_momentum_blocks(dg_frag, momentum_blocks_re, momentum_blocks_im, "momentum-soi", dg_frag%icomm)
    call sync_blocks_to_complex_dense_momentum(dg_frag, momentum_blocks_re, momentum_blocks_im, momentum_block_map_local, &
                                               dg_frag%momentum_mat_c)

    if (allocated(momentum_blocks_re)) then
      do iblk = 1, size(momentum_blocks_re)
        if (allocated(momentum_blocks_re(iblk)%val)) deallocate(momentum_blocks_re(iblk)%val)
      end do
      deallocate(momentum_blocks_re)
    end if
    if (allocated(momentum_blocks_im)) then
      do iblk = 1, size(momentum_blocks_im)
        if (allocated(momentum_blocks_im(iblk)%val)) deallocate(momentum_blocks_im(iblk)%val)
      end do
      deallocate(momentum_blocks_im)
    end if
    if (allocated(momentum_block_map_local)) deallocate(momentum_block_map_local)

    ! Enforce skew-Hermitian structure for gradient operator in complex basis:
    !   P = -P^\dagger.
    block
      integer :: ii, jj
      complex(8) :: pavg
      do ispin = 1, system%nspin
        do idir = 1, 3
          do ii = 1, dg_frag%n_mat_max
            dg_frag%momentum_mat_c(idir, ii, ii, ispin) = (0.0d0, 0.0d0)
            do jj = ii + 1, dg_frag%n_mat_max
              pavg = 0.5d0 * (dg_frag%momentum_mat_c(idir, ii, jj, ispin) - conjg(dg_frag%momentum_mat_c(idir, jj, ii, ispin)))
              dg_frag%momentum_mat_c(idir, ii, jj, ispin) = pavg
              dg_frag%momentum_mat_c(idir, jj, ii, ispin) = -conjg(pavg)
            end do
          end do
        end do
      end do
    end block

    dg_frag%momentum_mat = real(dg_frag%momentum_mat_c, kind=8)
    call init_momentum_blocks(dg_frag)
    do iblk = 1, dg_frag%n_momentum_blocks
      ifrag = dg_frag%momentum_blocks(iblk)%ifrag_row
      jfrag = dg_frag%momentum_blocks(iblk)%ifrag_col
      do ispin = 1, dg_frag%nspin
        nrow = dg_frag%n_basis(ifrag, ispin)
        ncol = dg_frag%n_basis(jfrag, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle
        do idir = 1, 3
          do jj = 1, ncol
            ig_j = dg_frag%index_basis(jj, jfrag, ispin)
            if (ig_j < 1 .or. ig_j > dg_frag%n_mat_max) cycle
            do ii = 1, nrow
              ig_i = dg_frag%index_basis(ii, ifrag, ispin)
              if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
              dg_frag%momentum_blocks(iblk)%val(idir, ii, jj, ispin) = dg_frag%momentum_mat(idir, ig_i, ig_j, ispin)
            end do
          end do
        end do
      end do
    end do

    max_p = maxval(abs(dg_frag%momentum_mat_c))
    if (comm_is_root(dg_frag%id)) then
      write(*,'(a,1pe12.4)') "        Max |momentum_mat|: ", max_p
      write(*,'(a,i0,a,i0,a)') "        Total matrix elements: ", &
                               3 * dg_frag%n_mat_max * dg_frag%n_mat_max * system%nspin, &
                               " (", 3, " directions × basis states × spin)"
    end if
    
  end subroutine calculate_momentum_matrix

  integer function find_matrix_block(block_map, ifrag_row, ifrag_col) result(iblk)
    implicit none
    integer, intent(in) :: block_map(:, :)
    integer, intent(in) :: ifrag_row, ifrag_col

    iblk = 0
    if (ifrag_row < 1 .or. ifrag_row > size(block_map, 1)) return
    if (ifrag_col < 1 .or. ifrag_col > size(block_map, 2)) return
    iblk = block_map(ifrag_row, ifrag_col)
  end function find_matrix_block

  logical function matrix_row_is_locally_owned(dg_frag, row_idx, ispin) result(is_local)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    integer, intent(in) :: row_idx, ispin

    is_local = .false.
    if (.not. allocated(dg_frag%coef_owner)) return
    if (ispin < 1 .or. ispin > size(dg_frag%coef_owner, 2)) return
    if (row_idx < 1 .or. row_idx > size(dg_frag%coef_owner, 1)) return
    is_local = (dg_frag%coef_owner(row_idx, ispin) == dg_frag%id)
  end function matrix_row_is_locally_owned

  subroutine init_matrix_blocks(dg_frag, blocks, block_map, n_blocks)
    use rt_dg_fragment_types, only: matrix_block_info
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(matrix_block_info), allocatable, intent(inout) :: blocks(:)
    integer, allocatable, intent(inout) :: block_map(:, :)
    integer, intent(out) :: n_blocks
    integer :: ifrag_row, ifrag_col, iblk, nrow_max, ncol_max

    if (allocated(blocks)) then
      do iblk = 1, size(blocks)
        if (allocated(blocks(iblk)%val)) deallocate(blocks(iblk)%val)
      end do
      deallocate(blocks)
    end if
    if (allocated(block_map)) deallocate(block_map)
    call ensure_momentum_neighbor_pair_cache(dg_frag)

    n_blocks = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col)) n_blocks = n_blocks + 1
      end do
    end do
    if (n_blocks <= 0) return

    allocate(blocks(n_blocks))
    allocate(block_map(dg_frag%n_frag, dg_frag%n_frag))
    block_map = 0

    iblk = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (.not. is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col)) cycle
        iblk = iblk + 1
        nrow_max = max(1, maxval(dg_frag%n_basis(ifrag_row, 1:dg_frag%nspin)))
        ncol_max = max(1, maxval(dg_frag%n_basis(ifrag_col, 1:dg_frag%nspin)))
        block_map(ifrag_row, ifrag_col) = iblk
        blocks(iblk)%ifrag_row = ifrag_row
        blocks(iblk)%ifrag_col = ifrag_col
        blocks(iblk)%nrow_max = nrow_max
        blocks(iblk)%ncol_max = ncol_max
        allocate(blocks(iblk)%val(nrow_max, ncol_max, dg_frag%nspin))
        blocks(iblk)%val = 0.0d0
      end do
    end do
  end subroutine init_matrix_blocks

  subroutine sync_dense_matrix_to_blocks(dg_frag, mat, blocks, block_map)
    use rt_dg_fragment_types, only: matrix_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    real(8), intent(in) :: mat(:, :, :)
    type(matrix_block_info), intent(inout) :: blocks(:)
    integer, intent(in) :: block_map(:, :)
    integer :: ifrag_row, ifrag_col, iblk, ispin, ii, jj, ig_i, ig_j
    integer :: nrow, ncol

    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        iblk = find_matrix_block(block_map, ifrag_row, ifrag_col)
        if (iblk <= 0) cycle
        blocks(iblk)%val(:, :, :) = 0.0d0
        do ispin = 1, dg_frag%nspin
          nrow = dg_frag%n_basis(ifrag_row, ispin)
          ncol = dg_frag%n_basis(ifrag_col, ispin)
          if (nrow <= 0 .or. ncol <= 0) cycle
          do jj = 1, ncol
            ig_j = dg_frag%index_basis(jj, ifrag_col, ispin)
            if (ig_j < 1 .or. ig_j > size(mat, 2)) cycle
            do ii = 1, nrow
              ig_i = dg_frag%index_basis(ii, ifrag_row, ispin)
              if (ig_i < 1 .or. ig_i > size(mat, 1)) cycle
              blocks(iblk)%val(ii, jj, ispin) = mat(ig_i, ig_j, ispin)
            end do
          end do
        end do
      end do
    end do
  end subroutine sync_dense_matrix_to_blocks

  subroutine sync_blocks_to_dense_matrix(dg_frag, blocks, block_map, mat)
    use rt_dg_fragment_types, only: matrix_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), intent(in) :: blocks(:)
    integer, intent(in) :: block_map(:, :)
    real(8), intent(inout) :: mat(:, :, :)
    integer :: ifrag_row, ifrag_col, iblk, ispin, ii, jj, ig_i, ig_j
    integer :: nrow, ncol, idx_ii, idx_jj, valid_row_count, valid_col_count
    integer, allocatable :: row_gid(:), col_gid(:), valid_row_ids(:), valid_col_ids(:)

    mat(:, :, :) = 0.0d0
    allocate(row_gid(size(dg_frag%index_basis, 1)), col_gid(size(dg_frag%index_basis, 1)))
    allocate(valid_row_ids(size(dg_frag%index_basis, 1)), valid_col_ids(size(dg_frag%index_basis, 1)))
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        iblk = find_matrix_block(block_map, ifrag_row, ifrag_col)
        if (iblk <= 0) cycle
        do ispin = 1, dg_frag%nspin
          nrow = dg_frag%n_basis(ifrag_row, ispin)
          ncol = dg_frag%n_basis(ifrag_col, ispin)
          if (nrow <= 0 .or. ncol <= 0) cycle
          valid_row_count = 0
          do ii = 1, nrow
            row_gid(ii) = dg_frag%index_basis(ii, ifrag_row, ispin)
            if (row_gid(ii) < 1 .or. row_gid(ii) > size(mat, 1)) cycle
            valid_row_count = valid_row_count + 1
            valid_row_ids(valid_row_count) = ii
          end do
          valid_col_count = 0
          do jj = 1, ncol
            col_gid(jj) = dg_frag%index_basis(jj, ifrag_col, ispin)
            if (col_gid(jj) < 1 .or. col_gid(jj) > size(mat, 2)) cycle
            valid_col_count = valid_col_count + 1
            valid_col_ids(valid_col_count) = jj
          end do
          do idx_jj = 1, valid_col_count
            jj = valid_col_ids(idx_jj)
            ig_j = col_gid(jj)
            do idx_ii = 1, valid_row_count
              ii = valid_row_ids(idx_ii)
              ig_i = row_gid(ii)
              mat(ig_i, ig_j, ispin) = blocks(iblk)%val(ii, jj, ispin)
            end do
          end do
        end do
      end do
    end do
    deallocate(row_gid, col_gid, valid_row_ids, valid_col_ids)
  end subroutine sync_blocks_to_dense_matrix

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
    if (allocated(dg_frag%momentum_neighbor_pair_cache)) then
      if (ifrag_row >= 1 .and. ifrag_row <= size(dg_frag%momentum_neighbor_pair_cache, 1) .and. &
          ifrag_col >= 1 .and. ifrag_col <= size(dg_frag%momentum_neighbor_pair_cache, 2)) then
        is_pair = dg_frag%momentum_neighbor_pair_cache(ifrag_row, ifrag_col)
        return
      end if
    end if

    do axis = 1, 3
      axis_ok(axis) = is_momentum_neighbor_axis(dg_frag%lgnum_total(axis), &
        dg_frag%ixyz_frag(axis, ifrag_row), dg_frag%nxyz_domain(axis, ifrag_row), &
        dg_frag%ixyz_frag(axis, ifrag_col), dg_frag%nxyz_domain(axis, ifrag_col))
    end do

    is_pair = all(axis_ok)
  end function is_momentum_neighbor_pair

  subroutine ensure_momentum_neighbor_pair_cache(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer :: ifrag_row, ifrag_col, axis
    logical :: axis_ok(3)

    if (allocated(dg_frag%momentum_neighbor_pair_cache)) return
    allocate(dg_frag%momentum_neighbor_pair_cache(dg_frag%n_frag, dg_frag%n_frag))
    dg_frag%momentum_neighbor_pair_cache(:, :) = .false.
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (ifrag_row == ifrag_col) then
          dg_frag%momentum_neighbor_pair_cache(ifrag_row, ifrag_col) = .true.
        else
          do axis = 1, 3
            axis_ok(axis) = is_momentum_neighbor_axis(dg_frag%lgnum_total(axis), &
              dg_frag%ixyz_frag(axis, ifrag_row), dg_frag%nxyz_domain(axis, ifrag_row), &
              dg_frag%ixyz_frag(axis, ifrag_col), dg_frag%nxyz_domain(axis, ifrag_col))
          end do
          dg_frag%momentum_neighbor_pair_cache(ifrag_row, ifrag_col) = all(axis_ok)
        end if
      end do
    end do
  end subroutine ensure_momentum_neighbor_pair_cache

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
    call ensure_momentum_neighbor_pair_cache(dg_frag)

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

  subroutine init_complex_momentum_blocks(dg_frag, blocks_re, blocks_im, block_map, n_blocks)
    use rt_dg_fragment_types, only: momentum_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(momentum_block_info), allocatable, intent(inout) :: blocks_re(:), blocks_im(:)
    integer, allocatable, intent(inout) :: block_map(:, :)
    integer, intent(out) :: n_blocks
    integer :: ifrag_row, ifrag_col, nblk, iblk
    integer :: nrow_max, ncol_max

    if (allocated(blocks_re)) then
      do iblk = 1, size(blocks_re)
        if (allocated(blocks_re(iblk)%val)) deallocate(blocks_re(iblk)%val)
      end do
      deallocate(blocks_re)
    end if
    if (allocated(blocks_im)) then
      do iblk = 1, size(blocks_im)
        if (allocated(blocks_im(iblk)%val)) deallocate(blocks_im(iblk)%val)
      end do
      deallocate(blocks_im)
    end if
    if (allocated(block_map)) deallocate(block_map)

    nblk = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col)) nblk = nblk + 1
      end do
    end do

    n_blocks = nblk
    if (nblk <= 0) return
    allocate(blocks_re(nblk), blocks_im(nblk))
    allocate(block_map(dg_frag%n_frag, dg_frag%n_frag))
    block_map = 0

    iblk = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        if (.not. is_momentum_neighbor_pair(dg_frag, ifrag_row, ifrag_col)) cycle
        iblk = iblk + 1
        block_map(ifrag_row, ifrag_col) = iblk
        nrow_max = max(1, maxval(dg_frag%n_basis(ifrag_row, 1:dg_frag%nspin)))
        ncol_max = max(1, maxval(dg_frag%n_basis(ifrag_col, 1:dg_frag%nspin)))

        blocks_re(iblk)%ifrag_row = ifrag_row
        blocks_re(iblk)%ifrag_col = ifrag_col
        blocks_re(iblk)%nrow_max = nrow_max
        blocks_re(iblk)%ncol_max = ncol_max
        allocate(blocks_re(iblk)%val(3, nrow_max, ncol_max, dg_frag%nspin))
        blocks_re(iblk)%val = 0.0d0

        blocks_im(iblk)%ifrag_row = ifrag_row
        blocks_im(iblk)%ifrag_col = ifrag_col
        blocks_im(iblk)%nrow_max = nrow_max
        blocks_im(iblk)%ncol_max = ncol_max
        allocate(blocks_im(iblk)%val(3, nrow_max, ncol_max, dg_frag%nspin))
        blocks_im(iblk)%val = 0.0d0
      end do
    end do
  end subroutine init_complex_momentum_blocks

  integer function find_momentum_block_runtime(block_map, ifrag_row, ifrag_col) result(iblk)
    implicit none
    integer, intent(in) :: block_map(:, :)
    integer, intent(in) :: ifrag_row, ifrag_col

    iblk = 0
    if (ifrag_row < 1 .or. ifrag_row > size(block_map, 1)) return
    if (ifrag_col < 1 .or. ifrag_col > size(block_map, 2)) return
    iblk = block_map(ifrag_row, ifrag_col)
  end function find_momentum_block_runtime

  subroutine sync_complex_dense_momentum_to_blocks(dg_frag, mat, blocks_re, blocks_im, block_map)
    use rt_dg_fragment_types, only: momentum_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    complex(8), intent(in) :: mat(:, :, :, :)
    type(momentum_block_info), intent(inout) :: blocks_re(:), blocks_im(:)
    integer, intent(in) :: block_map(:, :)
    integer :: ifrag_row, ifrag_col, iblk, idir, ispin, ii, jj, ig_i, ig_j
    integer :: nrow, ncol

    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        iblk = find_momentum_block_runtime(block_map, ifrag_row, ifrag_col)
        if (iblk <= 0) cycle
        blocks_re(iblk)%val(:, :, :, :) = 0.0d0
        blocks_im(iblk)%val(:, :, :, :) = 0.0d0
        do ispin = 1, dg_frag%nspin
          nrow = dg_frag%n_basis(ifrag_row, ispin)
          ncol = dg_frag%n_basis(ifrag_col, ispin)
          if (nrow <= 0 .or. ncol <= 0) cycle
          do jj = 1, ncol
            ig_j = dg_frag%index_basis(jj, ifrag_col, ispin)
            if (ig_j < 1 .or. ig_j > size(mat, 3)) cycle
            do ii = 1, nrow
              ig_i = dg_frag%index_basis(ii, ifrag_row, ispin)
              if (ig_i < 1 .or. ig_i > size(mat, 2)) cycle
              do idir = 1, 3
                blocks_re(iblk)%val(idir, ii, jj, ispin) = real(mat(idir, ig_i, ig_j, ispin), kind=8)
                blocks_im(iblk)%val(idir, ii, jj, ispin) = aimag(mat(idir, ig_i, ig_j, ispin))
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine sync_complex_dense_momentum_to_blocks

  subroutine sync_blocks_to_complex_dense_momentum(dg_frag, blocks_re, blocks_im, block_map, mat)
    use rt_dg_fragment_types, only: momentum_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(momentum_block_info), intent(in) :: blocks_re(:), blocks_im(:)
    integer, intent(in) :: block_map(:, :)
    complex(8), intent(inout) :: mat(:, :, :, :)
    integer :: ifrag_row, ifrag_col, iblk, idir, ispin, ii, jj, ig_i, ig_j
    integer :: nrow, ncol

    mat(:, :, :, :) = (0.0d0, 0.0d0)
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        iblk = find_momentum_block_runtime(block_map, ifrag_row, ifrag_col)
        if (iblk <= 0) cycle
        do ispin = 1, dg_frag%nspin
          nrow = dg_frag%n_basis(ifrag_row, ispin)
          ncol = dg_frag%n_basis(ifrag_col, ispin)
          if (nrow <= 0 .or. ncol <= 0) cycle
          do jj = 1, ncol
            ig_j = dg_frag%index_basis(jj, ifrag_col, ispin)
            if (ig_j < 1 .or. ig_j > size(mat, 3)) cycle
            do ii = 1, nrow
              ig_i = dg_frag%index_basis(ii, ifrag_row, ispin)
              if (ig_i < 1 .or. ig_i > size(mat, 2)) cycle
              do idir = 1, 3
                mat(idir, ig_i, ig_j, ispin) = cmplx(blocks_re(iblk)%val(idir, ii, jj, ispin), &
                                                     blocks_im(iblk)%val(idir, ii, jj, ispin), kind=8)
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine sync_blocks_to_complex_dense_momentum

  subroutine reduce_complex_momentum_blocks(dg_frag, blocks_re, blocks_im, label, icomm_reduce)
    use communication, only: comm_get_max, comm_is_root, comm_summation
    use rt_dg_fragment_types, only: momentum_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(momentum_block_info), intent(inout) :: blocks_re(:), blocks_im(:)
    character(*), intent(in) :: label
    integer, intent(in) :: icomm_reduce
    integer, parameter :: reduce_chunk_size = 262144
    real(8), allocatable :: send_block(:), recv_block(:)
    integer :: iblk, idir, ispin, ii, jj
    integer :: nrow, ncol, block_size, max_block_size, total_active_size
    integer :: total_active_min, total_active_max, max_block_size_global
    integer :: chunk_begin, chunk_count, offset_flat

    max_block_size = 0
    total_active_size = 0
    do iblk = 1, size(blocks_re)
      do ispin = 1, dg_frag%nspin
        nrow = dg_frag%n_basis(blocks_re(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(blocks_re(iblk)%ifrag_col, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle
        block_size = 3 * nrow * ncol
        max_block_size = max(max_block_size, block_size)
        total_active_size = total_active_size + block_size
      end do
    end do

    max_block_size_global = max_block_size
    call comm_get_max(max_block_size_global, icomm_reduce)
    total_active_max = total_active_size
    call comm_get_max(total_active_max, icomm_reduce)
    total_active_min = -total_active_size
    call comm_get_max(total_active_min, icomm_reduce)
    total_active_min = -total_active_min

    if (total_active_min /= total_active_max) then
      write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,i0)') "        [FATAL] Hamiltonian block size mismatch: label=", &
        trim(label), " rank=", dg_frag%id, " local=", total_active_size, &
        " min=", total_active_min, " max=", total_active_max
      flush(6)
      stop 1
    end if

    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,a,a,i0,a,i0,a,i0)') "        hamiltonian block reduce begin: label=", trim(label), &
        " total_active=", total_active_size, " max_block=", max_block_size_global, &
        " chunk_size=", reduce_chunk_size
      flush(6)
    end if

    if (max_block_size_global <= 0) return
    allocate(send_block(max_block_size_global), recv_block(max_block_size_global))

    do iblk = 1, size(blocks_re)
      do ispin = 1, dg_frag%nspin
        nrow = dg_frag%n_basis(blocks_re(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(blocks_re(iblk)%ifrag_col, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle
        block_size = 3 * nrow * ncol

        offset_flat = 1
        do idir = 1, 3
          do jj = 1, ncol
            do ii = 1, nrow
              send_block(offset_flat) = blocks_re(iblk)%val(idir, ii, jj, ispin)
              offset_flat = offset_flat + 1
            end do
          end do
        end do
        chunk_begin = 1
        do while (chunk_begin <= block_size)
          chunk_count = min(reduce_chunk_size, block_size - chunk_begin + 1)
          call comm_summation(send_block(chunk_begin:chunk_begin + chunk_count - 1), &
                              recv_block(chunk_begin:chunk_begin + chunk_count - 1), chunk_count, icomm_reduce)
          chunk_begin = chunk_begin + chunk_count
        end do
        offset_flat = 1
        do idir = 1, 3
          do jj = 1, ncol
            do ii = 1, nrow
              blocks_re(iblk)%val(idir, ii, jj, ispin) = recv_block(offset_flat)
              offset_flat = offset_flat + 1
            end do
          end do
        end do

        offset_flat = 1
        do idir = 1, 3
          do jj = 1, ncol
            do ii = 1, nrow
              send_block(offset_flat) = blocks_im(iblk)%val(idir, ii, jj, ispin)
              offset_flat = offset_flat + 1
            end do
          end do
        end do
        chunk_begin = 1
        do while (chunk_begin <= block_size)
          chunk_count = min(reduce_chunk_size, block_size - chunk_begin + 1)
          call comm_summation(send_block(chunk_begin:chunk_begin + chunk_count - 1), &
                              recv_block(chunk_begin:chunk_begin + chunk_count - 1), chunk_count, icomm_reduce)
          chunk_begin = chunk_begin + chunk_count
        end do
        offset_flat = 1
        do idir = 1, 3
          do jj = 1, ncol
            do ii = 1, nrow
              blocks_im(iblk)%val(idir, ii, jj, ispin) = recv_block(offset_flat)
              offset_flat = offset_flat + 1
            end do
          end do
        end do
      end do
    end do

    deallocate(send_block, recv_block)
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,a,a,i0)') "        hamiltonian block reduce done: label=", trim(label), &
        " total_active=", total_active_size
      flush(6)
    end if
  end subroutine reduce_complex_momentum_blocks

  subroutine reduce_matrix_blocks(dg_frag, blocks, label, icomm_reduce)
    use communication, only: comm_get_max, comm_is_root, comm_summation
    use rt_dg_fragment_types, only: matrix_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), intent(inout) :: blocks(:)
    character(*), intent(in) :: label
    integer, intent(in) :: icomm_reduce
    integer, parameter :: reduce_chunk_size = 262144
    real(8), allocatable :: send_block(:), recv_block(:)
    integer :: iblk, ispin, ii, jj
    integer :: nrow, ncol, block_size, max_block_size, total_active_size
    integer :: total_active_min, total_active_max, max_block_size_global
    integer :: chunk_begin, chunk_count, offset_flat

    max_block_size = 0
    total_active_size = 0
    do iblk = 1, size(blocks)
      do ispin = 1, dg_frag%nspin
        nrow = dg_frag%n_basis(blocks(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(blocks(iblk)%ifrag_col, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle
        block_size = nrow * ncol
        max_block_size = max(max_block_size, block_size)
        total_active_size = total_active_size + block_size
      end do
    end do

    max_block_size_global = max_block_size
    call comm_get_max(max_block_size_global, icomm_reduce)
    total_active_max = total_active_size
    call comm_get_max(total_active_max, icomm_reduce)
    total_active_min = -total_active_size
    call comm_get_max(total_active_min, icomm_reduce)
    total_active_min = -total_active_min

    if (total_active_min /= total_active_max) then
      write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,i0)') "        [FATAL] Hamiltonian block size mismatch: label=", &
        trim(label), " rank=", dg_frag%id, " local=", total_active_size, &
        " min=", total_active_min, " max=", total_active_max
      flush(6)
      stop 1
    end if

    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,a,a,i0,a,i0,a,i0)') "        hamiltonian block reduce begin: label=", trim(label), &
        " total_active=", total_active_size, " max_block=", max_block_size_global, &
        " chunk_size=", reduce_chunk_size
      flush(6)
    end if

    if (max_block_size_global <= 0) return
    allocate(send_block(max_block_size_global), recv_block(max_block_size_global))

    do iblk = 1, size(blocks)
      do ispin = 1, dg_frag%nspin
        nrow = dg_frag%n_basis(blocks(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(blocks(iblk)%ifrag_col, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle
        block_size = nrow * ncol
        offset_flat = 1
        do jj = 1, ncol
          do ii = 1, nrow
            send_block(offset_flat) = blocks(iblk)%val(ii, jj, ispin)
            offset_flat = offset_flat + 1
          end do
        end do

        chunk_begin = 1
        do while (chunk_begin <= block_size)
          chunk_count = min(reduce_chunk_size, block_size - chunk_begin + 1)
          call comm_summation(send_block(chunk_begin:chunk_begin + chunk_count - 1), &
                              recv_block(chunk_begin:chunk_begin + chunk_count - 1), chunk_count, icomm_reduce)
          chunk_begin = chunk_begin + chunk_count
        end do

        offset_flat = 1
        do jj = 1, ncol
          do ii = 1, nrow
            blocks(iblk)%val(ii, jj, ispin) = recv_block(offset_flat)
            offset_flat = offset_flat + 1
          end do
        end do
      end do
    end do

    deallocate(send_block, recv_block)
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,a,a,i0)') "        hamiltonian block reduce done: label=", trim(label), &
        " total_active=", total_active_size
      flush(6)
    end if
  end subroutine reduce_matrix_blocks

  subroutine init_complex_matrix_blocks(dg_frag, blocks_re, blocks_im, block_map, n_blocks)
    use rt_dg_fragment_types, only: matrix_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), allocatable, intent(inout) :: blocks_re(:), blocks_im(:)
    integer, allocatable, intent(inout) :: block_map(:, :)
    integer, intent(out) :: n_blocks
    integer :: ifrag_row, ifrag_col, iblk, nrow_max, ncol_max

    if (allocated(blocks_re)) then
      do iblk = 1, size(blocks_re)
        if (allocated(blocks_re(iblk)%val)) deallocate(blocks_re(iblk)%val)
      end do
      deallocate(blocks_re)
    end if
    if (allocated(blocks_im)) then
      do iblk = 1, size(blocks_im)
        if (allocated(blocks_im(iblk)%val)) deallocate(blocks_im(iblk)%val)
      end do
      deallocate(blocks_im)
    end if
    if (allocated(block_map)) deallocate(block_map)

    n_blocks = dg_frag%n_frag * dg_frag%n_frag
    if (n_blocks <= 0) return

    allocate(blocks_re(n_blocks), blocks_im(n_blocks))
    allocate(block_map(dg_frag%n_frag, dg_frag%n_frag))
    block_map = 0

    iblk = 0
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        iblk = iblk + 1
        nrow_max = max(1, maxval(dg_frag%n_basis(ifrag_row, 1:dg_frag%nspin)))
        ncol_max = max(1, maxval(dg_frag%n_basis(ifrag_col, 1:dg_frag%nspin)))
        block_map(ifrag_row, ifrag_col) = iblk

        blocks_re(iblk)%ifrag_row = ifrag_row
        blocks_re(iblk)%ifrag_col = ifrag_col
        blocks_re(iblk)%nrow_max = nrow_max
        blocks_re(iblk)%ncol_max = ncol_max
        allocate(blocks_re(iblk)%val(nrow_max, ncol_max, dg_frag%nspin))
        blocks_re(iblk)%val = 0.0d0

        blocks_im(iblk)%ifrag_row = ifrag_row
        blocks_im(iblk)%ifrag_col = ifrag_col
        blocks_im(iblk)%nrow_max = nrow_max
        blocks_im(iblk)%ncol_max = ncol_max
        allocate(blocks_im(iblk)%val(nrow_max, ncol_max, dg_frag%nspin))
        blocks_im(iblk)%val = 0.0d0
      end do
    end do
  end subroutine init_complex_matrix_blocks

  subroutine sync_complex_dense_matrix_to_blocks(dg_frag, mat, blocks_re, blocks_im, block_map)
    use rt_dg_fragment_types, only: matrix_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    complex(8), intent(in) :: mat(:, :, :)
    type(matrix_block_info), intent(inout) :: blocks_re(:), blocks_im(:)
    integer, intent(in) :: block_map(:, :)
    integer :: ifrag_row, ifrag_col, iblk, ispin, ii, jj, ig_i, ig_j
    integer :: nrow, ncol

    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        iblk = find_complex_matrix_block(block_map, ifrag_row, ifrag_col)
        if (iblk <= 0) cycle
        blocks_re(iblk)%val(:, :, :) = 0.0d0
        blocks_im(iblk)%val(:, :, :) = 0.0d0
        do ispin = 1, dg_frag%nspin
          nrow = dg_frag%n_basis(ifrag_row, ispin)
          ncol = dg_frag%n_basis(ifrag_col, ispin)
          if (nrow <= 0 .or. ncol <= 0) cycle
          do jj = 1, ncol
            ig_j = dg_frag%index_basis(jj, ifrag_col, ispin)
            if (ig_j < 1 .or. ig_j > size(mat, 2)) cycle
            do ii = 1, nrow
              ig_i = dg_frag%index_basis(ii, ifrag_row, ispin)
              if (ig_i < 1 .or. ig_i > size(mat, 1)) cycle
              blocks_re(iblk)%val(ii, jj, ispin) = real(mat(ig_i, ig_j, ispin), kind=8)
              blocks_im(iblk)%val(ii, jj, ispin) = aimag(mat(ig_i, ig_j, ispin))
            end do
          end do
        end do
      end do
    end do
  end subroutine sync_complex_dense_matrix_to_blocks

  subroutine sync_blocks_to_complex_dense_matrix(dg_frag, blocks_re, blocks_im, block_map, mat)
    use rt_dg_fragment_types, only: matrix_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), intent(in) :: blocks_re(:), blocks_im(:)
    integer, intent(in) :: block_map(:, :)
    complex(8), intent(inout) :: mat(:, :, :)
    integer :: ifrag_row, ifrag_col, iblk, ispin, ii, jj, ig_i, ig_j
    integer :: nrow, ncol

    mat(:, :, :) = (0.0d0, 0.0d0)
    do ifrag_col = 1, dg_frag%n_frag
      do ifrag_row = 1, dg_frag%n_frag
        iblk = find_complex_matrix_block(block_map, ifrag_row, ifrag_col)
        if (iblk <= 0) cycle
        do ispin = 1, dg_frag%nspin
          nrow = dg_frag%n_basis(ifrag_row, ispin)
          ncol = dg_frag%n_basis(ifrag_col, ispin)
          if (nrow <= 0 .or. ncol <= 0) cycle
          do jj = 1, ncol
            ig_j = dg_frag%index_basis(jj, ifrag_col, ispin)
            if (ig_j < 1 .or. ig_j > size(mat, 2)) cycle
            do ii = 1, nrow
              ig_i = dg_frag%index_basis(ii, ifrag_row, ispin)
              if (ig_i < 1 .or. ig_i > size(mat, 1)) cycle
              mat(ig_i, ig_j, ispin) = cmplx(blocks_re(iblk)%val(ii, jj, ispin), &
                                             blocks_im(iblk)%val(ii, jj, ispin), kind=8)
            end do
          end do
        end do
      end do
    end do
  end subroutine sync_blocks_to_complex_dense_matrix

  subroutine reduce_complex_matrix_blocks(dg_frag, blocks_re, blocks_im, label, icomm_reduce)
    use communication, only: comm_get_max, comm_is_root, comm_summation
    use rt_dg_fragment_types, only: matrix_block_info
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(matrix_block_info), intent(inout) :: blocks_re(:), blocks_im(:)
    character(*), intent(in) :: label
    integer, intent(in) :: icomm_reduce
    integer, parameter :: reduce_chunk_size = 262144
    real(8), allocatable :: send_block(:), recv_block(:)
    integer :: iblk, ispin, ii, jj
    integer :: nrow, ncol, block_size, max_block_size, total_active_size
    integer :: total_active_min, total_active_max, max_block_size_global
    integer :: chunk_begin, chunk_count, offset_flat

    max_block_size = 0
    total_active_size = 0
    do iblk = 1, size(blocks_re)
      do ispin = 1, dg_frag%nspin
        nrow = dg_frag%n_basis(blocks_re(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(blocks_re(iblk)%ifrag_col, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle
        block_size = nrow * ncol
        max_block_size = max(max_block_size, block_size)
        total_active_size = total_active_size + block_size
      end do
    end do

    max_block_size_global = max_block_size
    call comm_get_max(max_block_size_global, icomm_reduce)
    total_active_max = total_active_size
    call comm_get_max(total_active_max, icomm_reduce)
    total_active_min = -total_active_size
    call comm_get_max(total_active_min, icomm_reduce)
    total_active_min = -total_active_min

    if (total_active_min /= total_active_max) then
      write(*,'(1x,a,a,a,i0,a,i0,a,i0,a,i0)') "        [FATAL] Hamiltonian block size mismatch: label=", &
        trim(label), " rank=", dg_frag%id, " local=", total_active_size, &
        " min=", total_active_min, " max=", total_active_max
      flush(6)
      stop 1
    end if

    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,a,a,i0,a,i0,a,i0)') "        hamiltonian block reduce begin: label=", trim(label), &
        " total_active=", total_active_size, " max_block=", max_block_size_global, &
        " chunk_size=", reduce_chunk_size
      flush(6)
    end if

    if (max_block_size_global <= 0) return
    allocate(send_block(max_block_size_global), recv_block(max_block_size_global))

    do iblk = 1, size(blocks_re)
      do ispin = 1, dg_frag%nspin
        nrow = dg_frag%n_basis(blocks_re(iblk)%ifrag_row, ispin)
        ncol = dg_frag%n_basis(blocks_re(iblk)%ifrag_col, ispin)
        if (nrow <= 0 .or. ncol <= 0) cycle
        block_size = nrow * ncol

        offset_flat = 1
        do jj = 1, ncol
          do ii = 1, nrow
            send_block(offset_flat) = blocks_re(iblk)%val(ii, jj, ispin)
            offset_flat = offset_flat + 1
          end do
        end do
        chunk_begin = 1
        do while (chunk_begin <= block_size)
          chunk_count = min(reduce_chunk_size, block_size - chunk_begin + 1)
          call comm_summation(send_block(chunk_begin:chunk_begin + chunk_count - 1), &
                              recv_block(chunk_begin:chunk_begin + chunk_count - 1), chunk_count, icomm_reduce)
          chunk_begin = chunk_begin + chunk_count
        end do
        offset_flat = 1
        do jj = 1, ncol
          do ii = 1, nrow
            blocks_re(iblk)%val(ii, jj, ispin) = recv_block(offset_flat)
            offset_flat = offset_flat + 1
          end do
        end do

        offset_flat = 1
        do jj = 1, ncol
          do ii = 1, nrow
            send_block(offset_flat) = blocks_im(iblk)%val(ii, jj, ispin)
            offset_flat = offset_flat + 1
          end do
        end do
        chunk_begin = 1
        do while (chunk_begin <= block_size)
          chunk_count = min(reduce_chunk_size, block_size - chunk_begin + 1)
          call comm_summation(send_block(chunk_begin:chunk_begin + chunk_count - 1), &
                              recv_block(chunk_begin:chunk_begin + chunk_count - 1), chunk_count, icomm_reduce)
          chunk_begin = chunk_begin + chunk_count
        end do
        offset_flat = 1
        do jj = 1, ncol
          do ii = 1, nrow
            blocks_im(iblk)%val(ii, jj, ispin) = recv_block(offset_flat)
            offset_flat = offset_flat + 1
          end do
        end do
      end do
    end do

    deallocate(send_block, recv_block)
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,a,a,i0)') "        hamiltonian block reduce done: label=", trim(label), &
        " total_active=", total_active_size
      flush(6)
    end if
  end subroutine reduce_complex_matrix_blocks

  integer function find_complex_matrix_block(block_map, ifrag_row, ifrag_col) result(iblk)
    implicit none
    integer, intent(in) :: block_map(:, :)
    integer, intent(in) :: ifrag_row, ifrag_col

    iblk = 0
    if (ifrag_row < 1 .or. ifrag_row > size(block_map, 1)) return
    if (ifrag_col < 1 .or. ifrag_col > size(block_map, 2)) return
    iblk = block_map(ifrag_row, ifrag_col)
  end function find_complex_matrix_block

  !=======================================================================
  ! Calculate overlap matrix in DG basis (S_ij = <phi_i|phi_j>)
  !=======================================================================
  subroutine calculate_overlap_matrix(dg_frag, system, mg)
    use structures
    use rt_dg_fragment_types, only: matrix_block_info
    use communication, only: comm_is_root, COMM_GROUP_NULL
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg

    integer :: ifrag, i_local, ispin, io, jo
    integer :: ix, iy, iz, is(3), ie(3), i_halo, jfrag, n_basis_halo
    integer :: ig_row, ig_col, l(3), d(3), ii, jj
    integer :: lx, ly, lz, iorg(3), ndom(3)
    integer :: n_eval, lwork, info_eig, n_blocks, icomm_reduce
    real(8) :: hvol, s_min, s_max, cond_est
    complex(8) :: integral
    complex(8) :: cwork_query(1)
    complex(8), allocatable :: S_eval(:,:), eig_work(:)
    real(8), allocatable :: eigvals(:), rwork(:)
    type(matrix_block_info), allocatable :: S_blocks_re(:), S_blocks_im(:)
    integer, allocatable :: S_block_map_local(:,:)

    if (.not. dg_frag%has_real_space_basis) return
    if (.not. allocated(dg_frag%index_basis) .or. .not. allocated(dg_frag%n_mat)) return

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

        do jo = 1, dg_frag%n_basis(ifrag, ispin)
          ig_col = dg_frag%index_basis(jo, ifrag, ispin)
          if (ig_col < 1 .or. ig_col > dg_frag%n_mat_max) cycle

          do io = 1, dg_frag%n_basis(ifrag, ispin)
            ig_row = dg_frag%index_basis(io, ifrag, ispin)
            if (ig_row < 1 .or. ig_row > dg_frag%n_mat_max) cycle
            integral = (0.0d0, 0.0d0)
            do lz = 1, ndom(3)
              do ly = 1, ndom(2)
                do lx = 1, ndom(1)
                  if (allocated(dg_frag%phi_frag_c)) then
                    integral = integral + conjg(dg_frag%phi_frag_c(lx, ly, lz, io, i_local)) * &
                               dg_frag%phi_frag_c(lx, ly, lz, jo, i_local) * hvol
                  else
                    integral = integral + cmplx(dg_frag%phi_frag(lx, ly, lz, io, i_local), 0.0d0, kind=8) * &
                               cmplx(dg_frag%phi_frag(lx, ly, lz, jo, i_local), 0.0d0, kind=8) * hvol
                  end if
                end do
              end do
            end do
            dg_frag%S_mat_c(ig_row, ig_col, ispin) = integral
          end do

          do i_halo = 1, dg_frag%n_halo
            if (dg_frag%halo(i_halo)%ifrag_dst /= ifrag) cycle
            jfrag = dg_frag%halo(i_halo)%ifrag_src
            if (jfrag < 1) cycle
            n_basis_halo = dg_frag%n_basis(jfrag, ispin)
            l = dg_frag%halo(i_halo)%length
            d = dg_frag%halo(i_halo)%dsp_send

            do io = 1, n_basis_halo
              ig_row = dg_frag%index_basis(io, jfrag, ispin)
              if (ig_row < 1 .or. ig_row > dg_frag%n_mat_max) cycle
              integral = (0.0d0, 0.0d0)
              do iz = 1, l(3)
                do iy = 1, l(2)
                  do ix = 1, l(1)
                    if (allocated(dg_frag%halo(i_halo)%buf_recv_c) .and. allocated(dg_frag%phi_frag_c)) then
                      integral = integral + conjg(dg_frag%halo(i_halo)%buf_recv_c(ix, iy, iz, io, 1)) * &
                                 dg_frag%phi_frag_c(d(1) + ix, d(2) + iy, d(3) + iz, jo, i_local) * hvol
                    else
                      integral = integral + cmplx(dg_frag%halo(i_halo)%buf_recv(ix, iy, iz, io, 1), 0.0d0, kind=8) * &
                                 cmplx(dg_frag%phi_frag(d(1) + ix, d(2) + iy, d(3) + iz, jo, i_local), 0.0d0, kind=8) * hvol
                    end if
                  end do
                end do
              end do
              dg_frag%S_mat_c(ig_row, ig_col, ispin) = dg_frag%S_mat_c(ig_row, ig_col, ispin) + 0.5d0 * integral
              dg_frag%S_mat_c(ig_col, ig_row, ispin) = dg_frag%S_mat_c(ig_col, ig_row, ispin) + 0.5d0 * conjg(integral)
            end do
          end do

        end do
      end do
    end do

    icomm_reduce = dg_frag%icomm
    if (dg_frag%icomm_frag /= COMM_GROUP_NULL) icomm_reduce = dg_frag%icomm_frag

    call init_complex_matrix_blocks(dg_frag, S_blocks_re, S_blocks_im, S_block_map_local, n_blocks)
    call sync_complex_dense_matrix_to_blocks(dg_frag, dg_frag%S_mat_c, S_blocks_re, S_blocks_im, S_block_map_local)
    call reduce_complex_matrix_blocks(dg_frag, S_blocks_re, S_blocks_im, "smat-soi", icomm_reduce)
    call sync_blocks_to_complex_dense_matrix(dg_frag, S_blocks_re, S_blocks_im, S_block_map_local, dg_frag%S_mat_c)
    if (icomm_reduce == dg_frag%icomm_frag .and. .not. dg_frag%is_frag_root) then
      do ispin = 1, dg_frag%nspin
        do ii = 1, dg_frag%n_mat_max
          if (matrix_row_is_locally_owned(dg_frag, ii, ispin)) cycle
          dg_frag%S_mat_c(ii, :, ispin) = (0.0d0, 0.0d0)
        end do
      end do
    end if

    do ispin = 1, dg_frag%nspin
      do ii = 1, dg_frag%n_mat_max
        if (real(dg_frag%S_mat_c(ii, ii, ispin), kind=8) < 1.0d-12) dg_frag%S_mat_c(ii, ii, ispin) = (1.0d0, 0.0d0)
        do jj = ii + 1, dg_frag%n_mat_max
          dg_frag%S_mat_c(ii, jj, ispin) = 0.5d0 * (dg_frag%S_mat_c(ii, jj, ispin) + conjg(dg_frag%S_mat_c(jj, ii, ispin)))
          dg_frag%S_mat_c(jj, ii, ispin) = conjg(dg_frag%S_mat_c(ii, jj, ispin))
        end do
      end do
    end do

    dg_frag%S_mat = real(dg_frag%S_mat_c, kind=8)

    ! Default propagation overlap equals the raw fragment overlap.
    dg_frag%S_mat_prop(:, :, :) = dg_frag%S_mat(:, :, :)
    dg_frag%S_mat_prop_c(:, :, :) = dg_frag%S_mat_c(:, :, :)
    call init_matrix_blocks(dg_frag, dg_frag%S_mat_blocks, dg_frag%S_block_map, dg_frag%n_S_blocks)
    call sync_dense_matrix_to_blocks(dg_frag, dg_frag%S_mat, dg_frag%S_mat_blocks, dg_frag%S_block_map)
    call init_matrix_blocks(dg_frag, dg_frag%S_mat_prop_blocks, dg_frag%S_block_map, dg_frag%n_S_blocks)
    call sync_dense_matrix_to_blocks(dg_frag, dg_frag%S_mat_prop, dg_frag%S_mat_prop_blocks, dg_frag%S_block_map)
    dg_frag%has_global_overlap_copy = .false.
    dg_frag%overlap_prop_root_authoritative = .false.

    if (allocated(S_blocks_re)) then
      do ii = 1, size(S_blocks_re)
        if (allocated(S_blocks_re(ii)%val)) deallocate(S_blocks_re(ii)%val)
      end do
      deallocate(S_blocks_re)
    end if
    if (allocated(S_blocks_im)) then
      do ii = 1, size(S_blocks_im)
        if (allocated(S_blocks_im(ii)%val)) deallocate(S_blocks_im(ii)%val)
      end do
      deallocate(S_blocks_im)
    end if
    if (allocated(S_block_map_local)) deallocate(S_block_map_local)

  end subroutine calculate_overlap_matrix
