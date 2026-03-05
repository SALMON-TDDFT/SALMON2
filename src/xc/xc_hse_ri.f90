!========================================================================================!
! Module: xc_hse_ri
! 
! Purpose: RI/DF approximation for HSE exchange in DG-Fragment RT-TDDFT
!          Plan C implementation with 3-index integrals
!
! Author: Plan C Implementation Team
! Date: 2026-02-23
!========================================================================================!

module xc_hse_ri
  use auxiliary_basis
  use structures, only: s_rgrid
  implicit none
  
  private
  public :: init_hse_ri_fragment
  public :: calc_exact_exchange_hse_ri
  public :: deallocate_hse_ri_fragment
  
  ! RI data structure for HSE exchange
  type, public :: hse_ri_data_t
    logical :: initialized = .false.
    integer :: n_basis                          ! Number of basis functions
    integer :: n_aux                            ! Number of auxiliary functions
    integer :: n_chol                           ! Number of Cholesky vectors (CD-RI)
    logical :: use_cd_ri = .false.              ! Use Cholesky decomposition RI
    real(8), allocatable :: B_ijP(:,:,:)       ! 3-index integrals (n_basis, n_basis, n_aux)
    real(8), allocatable :: V_inv_PQ(:,:)      ! Inverse of Coulomb matrix (n_aux, n_aux) - standard RI
    real(8), allocatable :: L_PK(:,:)          ! Cholesky vectors (n_aux, n_chol) - CD-RI
    type(auxiliary_basis_t) :: aux_basis       ! Auxiliary basis set
  end type hse_ri_data_t
  
contains

  !--------------------------------------------------------------------------------------!
  ! Initialize RI data for a fragment (called once at the beginning)
  ! This is the expensive O(N²_basis × N_aux × L⁶) step
  !--------------------------------------------------------------------------------------!
  subroutine init_hse_ri_fragment(ri_data, phi_frag, lg, mg, ng, hvol, &
                                  natom, atom_coords, atom_types, hse_omega, &
                                  use_cd_ri, cd_ri_threshold)
    implicit none
    type(hse_ri_data_t), intent(inout) :: ri_data
    type(s_rgrid), intent(in) :: lg                         ! Local grid info
    type(s_rgrid), intent(in) :: mg                         ! Global grid info
    integer, intent(in) :: ng                               ! Number of basis functions
    real(8), intent(in) :: hvol                             ! Grid volume element
    integer, intent(in) :: natom                            ! Number of atoms in fragment
    real(8), intent(in) :: atom_coords(3, natom)          ! Atomic coordinates
    integer, intent(in) :: atom_types(natom)              ! Atomic types
    real(8), intent(in) :: hse_omega                        ! HSE screening parameter
    logical, intent(in), optional :: use_cd_ri              ! Use CD-RI (default: .false.)
    real(8), intent(in), optional :: cd_ri_threshold        ! CD-RI threshold (default: 1.0d-8)
    real(8), intent(in) :: phi_frag(lg%is(1):lg%ie(1), &
                                    lg%is(2):lg%ie(2), &
                                    lg%is(3):lg%ie(3), ng)  ! Basis functions on grid
    
    integer :: i, j, P
    real(8) :: threshold
    
    ! Store dimensions
    ri_data%n_basis = ng
    
    ! Set CD-RI options
    if (present(use_cd_ri)) then
      ri_data%use_cd_ri = use_cd_ri
    else
      ri_data%use_cd_ri = .false.
    end if
    
    if (present(cd_ri_threshold)) then
      threshold = cd_ri_threshold
    else
      threshold = 1.0d-8  ! Default threshold
    end if
    
    ! Initialize auxiliary basis
    call init_auxiliary_basis(ri_data%aux_basis, natom, atom_coords, atom_types)
    ri_data%n_aux = get_n_auxiliary(ri_data%aux_basis)
    
    write(*, '(A, I0, A, I0, A, F0.2)') &
      'Initializing RI-HSE: N_basis=', ri_data%n_basis, &
      ', N_aux=', ri_data%n_aux, &
      ', Ratio=', real(ri_data%n_aux) / real(ri_data%n_basis)
    
    ! Allocate 3-index integrals
    allocate(ri_data%B_ijP(ri_data%n_basis, ri_data%n_basis, ri_data%n_aux))
    ri_data%B_ijP = 0.0d0
    
    ! Calculate 3-index integrals: B_ijP = (ij|P)
    write(*, '(A)') 'Computing 3-index integrals B_ijP...'
    call compute_3index_integrals(ri_data%B_ijP, phi_frag, lg, mg, &
                                  ri_data%n_basis, ri_data%aux_basis, &
                                  ri_data%n_aux, hvol, hse_omega)
    
    ! Calculate and invert Coulomb matrix: V_PQ = (P|Q)
    if (ri_data%use_cd_ri) then
      ! CD-RI: Cholesky decomposition with threshold
      write(*, '(A)') 'Computing Coulomb matrix with CD-RI...'
      call compute_coulomb_cholesky(ri_data%L_PK, ri_data%n_chol, &
                                    ri_data%aux_basis, lg, mg, hvol, &
                                    hse_omega, threshold)
      write(*, '(A, I0, A, I0, A, F6.2, A)') &
        '  CD-RI: N_chol=', ri_data%n_chol, ' / N_aux=', ri_data%n_aux, &
        ' (', 100.0d0 * (1.0d0 - real(ri_data%n_chol)/real(ri_data%n_aux)), '% reduction)'
    else
      ! Standard RI: Full matrix inversion
      allocate(ri_data%V_inv_PQ(ri_data%n_aux, ri_data%n_aux))
      ri_data%V_inv_PQ = 0.0d0
      write(*, '(A)') 'Computing Coulomb matrix V_PQ and its inverse...'
      call compute_coulomb_matrix_inverse(ri_data%V_inv_PQ, ri_data%aux_basis, &
                                          lg, mg, hvol, hse_omega)
    end if
    
    ri_data%initialized = .true.
    write(*, '(A)') 'RI-HSE initialization completed!'
    
  end subroutine init_hse_ri_fragment

  !--------------------------------------------------------------------------------------!
  ! Compute 3-index integrals: B_ijP = (ij|P)
  ! = ∫∫ φ_i(r1) φ_j(r1) × [1/|r1-r2|] × χ_P(r2) dr1 dr2
  !--------------------------------------------------------------------------------------!
  subroutine compute_3index_integrals(B_ijP, phi_frag, lg, mg, n_basis, &
                                      aux_basis, n_aux, hvol, hse_omega)
    implicit none
    real(8), intent(out) :: B_ijP(:,:,:)           ! (n_basis, n_basis, n_aux)
    real(8), intent(in) :: phi_frag(:,:,:,:)       ! Basis functions on grid
    type(s_rgrid), intent(in) :: lg, mg
    integer, intent(in) :: n_basis, n_aux
    type(auxiliary_basis_t), intent(in) :: aux_basis
    real(8), intent(in) :: hvol
    real(8), intent(in) :: hse_omega
    
    integer :: i, j, P
    integer :: ix1, iy1, iz1, ix2, iy2, iz2
    real(8) :: r1(3), r2(3), distance, phi_i_r1, phi_j_r1, chi_P_r2
    real(8) :: coulomb_kernel
    real(8) :: integral_val
    integer :: progress_count, total_count
    real(8) :: aux_center(3), r1_to_aux_center
    real(8), parameter :: cutoff_distance = 15.0d0  ! Distance cutoff in bohr (HSE06 screening)
    integer :: n_skipped
    
    total_count = n_basis * n_basis * n_aux
    progress_count = 0
    n_skipped = 0
    
    !$omp parallel do collapse(3) private(i,j,P,ix1,iy1,iz1,ix2,iy2,iz2,r1,r2, &
    !$omp& distance,phi_i_r1,phi_j_r1,chi_P_r2,coulomb_kernel,integral_val, &
    !$omp& aux_center,r1_to_aux_center) reduction(+:n_skipped) &
    !$omp& schedule(dynamic)
    do P = 1, n_aux
      do j = 1, n_basis
        do i = 1, n_basis
          
          ! Get auxiliary basis center for distance-based screening
          aux_center(1) = aux_basis%center(1, P)
          aux_center(2) = aux_basis%center(2, P)
          aux_center(3) = aux_basis%center(3, P)
          
          integral_val = 0.0d0
          
          ! Double loop over grid points
          do iz1 = lg%is(3), lg%ie(3)
            do iy1 = lg%is(2), lg%ie(2)
              do ix1 = lg%is(1), lg%ie(1)
                
                ! Get r1 and basis function values
                r1(1) = mg%coordinate(ix1, 1)
                r1(2) = mg%coordinate(iy1, 2)
                r1(3) = mg%coordinate(iz1, 3)
                phi_i_r1 = phi_frag(ix1, iy1, iz1, i)
                phi_j_r1 = phi_frag(ix1, iy1, iz1, j)
                
                ! Skip if basis product is negligible
                if (abs(phi_i_r1 * phi_j_r1) < 1.0d-12) cycle
                
                ! Distance-based screening: check distance from r1 to auxiliary center
                ! If r1 is far from auxiliary center, skip the r2 integration
                r1_to_aux_center = sqrt((r1(1)-aux_center(1))**2 + &
                                       (r1(2)-aux_center(2))**2 + &
                                       (r1(3)-aux_center(3))**2)
                if (r1_to_aux_center > cutoff_distance) then
                  n_skipped = n_skipped + 1
                  cycle
                end if
                
                do iz2 = lg%is(3), lg%ie(3)
                  do iy2 = lg%is(2), lg%ie(2)
                    do ix2 = lg%is(1), lg%ie(1)
                      
                      ! Get r2 and auxiliary function value
                      r2(1) = mg%coordinate(ix2, 1)
                      r2(2) = mg%coordinate(iy2, 2)
                      r2(3) = mg%coordinate(iz2, 3)
                      
                      distance = sqrt((r1(1)-r2(1))**2 + (r1(2)-r2(2))**2 + (r1(3)-r2(3))**2)
                      
                      ! HSE short-range screened Coulomb kernel
                      if (distance < 1.0d-10) then
                        cycle
                      else
                        coulomb_kernel = erfc(hse_omega * distance) / distance
                      end if
                      
                      chi_P_r2 = calc_auxiliary_function(aux_basis, P, r2)
                      
                      ! Accumulate integral
                      integral_val = integral_val + phi_i_r1 * phi_j_r1 * coulomb_kernel * chi_P_r2
                      
                    end do
                  end do
                end do
                
              end do
            end do
          end do
          
          ! Multiply by volume elements
          B_ijP(i, j, P) = integral_val * hvol * hvol
          
          !$omp atomic
          progress_count = progress_count + 1
          
          ! Progress report (every 1%)
          if (mod(progress_count, total_count/100 + 1) == 0) then
            !$omp critical
            write(*, '(A, F6.2, A)', advance='no') char(13), 100.0d0 * progress_count / total_count, '% '
            !$omp end critical
          end if
          
        end do
      end do
    end do
    !$omp end parallel do
    
    write(*, '(A)') char(13) // '100.00% completed!'
    
    ! Report screening statistics
    if (n_skipped > 0) then
      write(*,'(A,I0,A,I0,A,F6.2,A)') '  Distance-based screening: skipped ', n_skipped, &
                ' / ', total_count, ' grid points (', &
                (100.0d0 * n_skipped) / max(total_count, 1), '% reduction)'
    end if
    
  end subroutine compute_3index_integrals

  !--------------------------------------------------------------------------------------!
  ! Compute Coulomb matrix and its inverse: V_PQ = (P|Q)
  ! Uses Cholesky decomposition for numerical stability
  !--------------------------------------------------------------------------------------!
  subroutine compute_coulomb_matrix_inverse(V_inv_PQ, aux_basis, lg, mg, hvol, hse_omega)
    implicit none
    real(8), intent(out) :: V_inv_PQ(:,:)          ! (n_aux, n_aux)
    type(auxiliary_basis_t), intent(in) :: aux_basis
    type(s_rgrid), intent(in) :: lg, mg
    real(8), intent(in) :: hvol
    real(8), intent(in) :: hse_omega
    
    integer :: n_aux, P, Q
    real(8) :: center_distance
    real(8), allocatable :: V_PQ(:,:)
    integer :: info
    
    n_aux = aux_basis%n_aux
    allocate(V_PQ(n_aux, n_aux))
    V_PQ = 0.0d0
    
    ! Compute diagonal elements (self-interaction, analytical)
    do P = 1, n_aux
      ! Analytical formula for Gaussian self-interaction
      ! (P|P) = (2α/π)^(3/2) × integral of exp(-2α r²) / r
      ! For HSE with screening ω, this is more complex, so we use numerical integration
      V_PQ(P, P) = compute_aux_self_interaction(aux_basis, P, lg, mg, hvol, hse_omega)
    end do
    
    ! Compute off-diagonal elements with center-based SR kernel approximation
    !$omp parallel do collapse(2) private(P,Q,center_distance)
    do Q = 1, n_aux
      do P = 1, Q-1
        center_distance = sqrt((aux_basis%center(1,P) - aux_basis%center(1,Q))**2 + &
                               (aux_basis%center(2,P) - aux_basis%center(2,Q))**2 + &
                               (aux_basis%center(3,P) - aux_basis%center(3,Q))**2)
        if (center_distance < 1.0d-12) then
          V_PQ(P, Q) = 0.5d0 * (V_PQ(P, P) + V_PQ(Q, Q))
        else
          V_PQ(P, Q) = erfc(hse_omega * center_distance) / center_distance
        end if
        
        V_PQ(Q, P) = V_PQ(P, Q)  ! Symmetry
        
      end do
    end do
    !$omp end parallel do
    
    ! Cholesky decomposition: V = L × L^T
    call dpotrf('L', n_aux, V_PQ, n_aux, info)
    if (info /= 0) then
      write(*, '(A, I0)') 'ERROR: Cholesky decomposition failed, info=', info
      stop
    end if
    
    ! Invert using Cholesky factor: V^(-1) = (L^T)^(-1) × L^(-1)
    call dpotri('L', n_aux, V_PQ, n_aux, info)
    if (info /= 0) then
      write(*, '(A, I0)') 'ERROR: Matrix inversion failed, info=', info
      stop
    end if
    
    ! Copy lower triangle to upper (dpotri only fills lower)
    do Q = 1, n_aux
      do P = 1, Q-1
        V_PQ(P, Q) = V_PQ(Q, P)
      end do
    end do
    
    V_inv_PQ = V_PQ
    
    deallocate(V_PQ)
    
  end subroutine compute_coulomb_matrix_inverse

  !--------------------------------------------------------------------------------------!
  ! Compute Coulomb matrix with Cholesky Decomposition (CD-RI)
  ! V_PQ ≈ Σ_K L_PK × L_QK with threshold-based truncation
  !--------------------------------------------------------------------------------------!
  subroutine compute_coulomb_cholesky(L_PK, n_chol, aux_basis, lg, mg, hvol, &
                                      hse_omega, threshold)
    implicit none
    real(8), allocatable, intent(out) :: L_PK(:,:)     ! Cholesky vectors (n_aux, n_chol)
    integer, intent(out) :: n_chol                      ! Number of Cholesky vectors retained
    type(auxiliary_basis_t), intent(in) :: aux_basis
    type(s_rgrid), intent(in) :: lg, mg
    real(8), intent(in) :: hvol, hse_omega, threshold
    
    integer :: n_aux, P, Q, K
    real(8) :: center_distance
    real(8), allocatable :: V_PQ(:,:), L_full(:,:)
    integer :: info
    real(8) :: max_diag, current_diag
    
    n_aux = aux_basis%n_aux
    allocate(V_PQ(n_aux, n_aux))
    V_PQ = 0.0d0
    
    ! Step 1: Compute Coulomb matrix V_PQ (same as standard RI)
    write(*, '(A)') '  Computing Coulomb matrix V_PQ...'
    
    ! Diagonal elements
    do P = 1, n_aux
      V_PQ(P, P) = compute_aux_self_interaction(aux_basis, P, lg, mg, hvol, hse_omega)
    end do
    
    ! Off-diagonal elements with center-based SR kernel approximation
    !$omp parallel do collapse(2) private(P,Q,center_distance)
    do Q = 1, n_aux
      do P = 1, Q-1
        center_distance = sqrt((aux_basis%center(1,P) - aux_basis%center(1,Q))**2 + &
                               (aux_basis%center(2,P) - aux_basis%center(2,Q))**2 + &
                               (aux_basis%center(3,P) - aux_basis%center(3,Q))**2)
        if (center_distance < 1.0d-12) then
          V_PQ(P, Q) = 0.5d0 * (V_PQ(P, P) + V_PQ(Q, Q))
        else
          V_PQ(P, Q) = erfc(hse_omega * center_distance) / center_distance
        end if
        V_PQ(Q, P) = V_PQ(P, Q)
      end do
    end do
    !$omp end parallel do
    
    ! Step 2: Perform Cholesky decomposition
    write(*, '(A)') '  Performing Cholesky decomposition...'
    allocate(L_full(n_aux, n_aux))
    L_full = V_PQ
    
    call dpotrf('L', n_aux, L_full, n_aux, info)
    if (info /= 0) then
      write(*, '(A, I0)') 'ERROR: Cholesky decomposition failed with info=', info
      stop
    end if
    
    ! Step 3: Determine number of significant Cholesky vectors
    max_diag = 0.0d0
    do P = 1, n_aux
      max_diag = max(max_diag, abs(L_full(P, P)))
    end do
    
    n_chol = 0
    do K = 1, n_aux
      current_diag = abs(L_full(K, K))
      if (current_diag > threshold * max_diag) then
        n_chol = K
      else
        exit  ! Remaining vectors are below threshold
      end if
    end do
    
    ! Step 4: Extract significant Cholesky vectors
    if (n_chol == 0) then
      write(*, '(A)') 'WARNING: No Cholesky vectors above threshold, using at least 1'
      n_chol = 1
    end if
    
    allocate(L_PK(n_aux, n_chol))
    do K = 1, n_chol
      do P = K, n_aux  ! Lower triangular
        L_PK(P, K) = L_full(P, K)
      end do
      do P = 1, K-1  ! Upper part is zero
        L_PK(P, K) = 0.0d0
      end do
    end do
    
    ! Report memory savings
    write(*, '(A, F8.3, A, F8.3, A, F6.2, A)') &
      '  Memory: Full=', (n_aux*n_aux*8.0d0)/(1024.0d0**2), ' MB, CD-RI=', &
      (n_aux*n_chol*8.0d0)/(1024.0d0**2), ' MB (', &
      100.0d0 * (1.0d0 - real(n_chol)/real(n_aux)), '% reduction)'
    
    deallocate(V_PQ, L_full)
    
  end subroutine compute_coulomb_cholesky

  !--------------------------------------------------------------------------------------!
  ! Compute self-interaction for auxiliary function
  !--------------------------------------------------------------------------------------!
  function compute_aux_self_interaction(aux_basis, P, lg, mg, hvol, hse_omega) result(value)
    implicit none
    type(auxiliary_basis_t), intent(in) :: aux_basis
    integer, intent(in) :: P
    type(s_rgrid), intent(in) :: lg, mg
    real(8), intent(in) :: hvol, hse_omega
    real(8) :: value
    
    real(8) :: alpha
    real(8), parameter :: pi = 3.14159265358979323846d0
    
    ! Analytical formula for Gaussian self-interaction with HSE screening
    alpha = aux_basis%alpha(P)
    
    ! Simplified: (P|P) ≈ (2α/π)^(3/2) × sqrt(π) / sqrt(2α + ω²)
    value = ((2.0d0 * alpha / pi)**1.5d0) * sqrt(pi) / sqrt(2.0d0 * alpha + hse_omega**2)
    
  end function compute_aux_self_interaction

  !--------------------------------------------------------------------------------------!
  ! Calculate HSE exchange matrix using RI approximation
  ! This is the fast O(N²_basis × N_aux × N_occ) step at each timestep
  !--------------------------------------------------------------------------------------!
  subroutine calc_exact_exchange_hse_ri(v_x, ri_data, density_matrix, hse_alpha, n_occ)
    implicit none
    real(8), intent(inout) :: v_x(:,:)                  ! Exchange matrix (n_basis, n_basis)
    type(hse_ri_data_t), intent(in) :: ri_data
    real(8), intent(in) :: density_matrix(:,:)          ! Density matrix (n_basis, n_basis)
    real(8), intent(in) :: hse_alpha
    integer, intent(in) :: n_occ
    
    integer :: n_basis, n_aux
    integer :: i, j, k, l, P, Q
    real(8), allocatable:: C_klQ(:,:,:)                ! Intermediate array
    real(8) :: v_x_ij
    
    if (.not. ri_data%initialized) then
      write(*, '(A)') 'ERROR: RI data not initialized!'
      stop
    end if
    
    n_basis = ri_data%n_basis
    n_aux = ri_data%n_aux
    
    ! Branch: CD-RI or standard RI
    if (ri_data%use_cd_ri) then
      ! CD-RI: Use Cholesky vectors L_PK
      ! C_klK = Σ_P B_klP × L_PK
      allocate(C_klQ(n_basis, n_basis, ri_data%n_chol))
      C_klQ = 0.0d0
      
      call dgemm('N', 'N', n_basis*n_basis, ri_data%n_chol, n_aux, &
                 1.0d0, ri_data%B_ijP, n_basis*n_basis, &
                 ri_data%L_PK, n_aux, &
                 0.0d0, C_klQ, n_basis*n_basis)
      
      ! Exchange matrix: V_x(i,j) = -α Σ_kl Σ_K B_ijK × C_klK × D_kl
      ! where B_ijK = Σ_P B_ijP × L_PK (computed on-the-fly)
      !$omp parallel do collapse(2) private(i,j,k,l,Q,v_x_ij)
      do j = 1, n_basis
        do i = 1, n_basis
          v_x_ij = 0.0d0
          do Q = 1, ri_data%n_chol  ! Loop over Cholesky vectors
            do l = 1, n_basis
              do k = 1, n_basis
                v_x_ij = v_x_ij + C_klQ(i, j, Q) * C_klQ(k, l, Q) * density_matrix(k, l)
              end do
            end do
          end do
          v_x(i, j) = v_x(i, j) - hse_alpha * v_x_ij
        end do
      end do
      !$omp end parallel do
      
    else
      ! Standard RI: Use full V^(-1)_PQ
      ! C_klQ = Σ_P B_klP × V^(-1)_PQ
      allocate(C_klQ(n_basis, n_basis, n_aux))
      C_klQ = 0.0d0
      
      call dgemm('N', 'N', n_basis*n_basis, n_aux, n_aux, &
                 1.0d0, ri_data%B_ijP, n_basis*n_basis, &
                 ri_data%V_inv_PQ, n_aux, &
                 0.0d0, C_klQ, n_basis*n_basis)
      
      ! Exchange matrix: V_x(i,j) = -α Σ_kl Σ_Q B_ijQ × C_klQ × D_kl
      !$omp parallel do collapse(2) private(i,j,k,l,Q,v_x_ij)
      do j = 1, n_basis
        do i = 1, n_basis
          v_x_ij = 0.0d0
          do Q = 1, n_aux
            do l = 1, n_basis
              do k = 1, n_basis
                v_x_ij = v_x_ij + ri_data%B_ijP(i, j, Q) * C_klQ(k, l, Q) * density_matrix(k, l)
              end do
            end do
          end do
          v_x(i, j) = v_x(i, j) - hse_alpha * v_x_ij
        end do
      end do
      !$omp end parallel do
    end if
    
    deallocate(C_klQ)
    
  end subroutine calc_exact_exchange_hse_ri

  !--------------------------------------------------------------------------------------!
  ! Deallocate RI data
  !--------------------------------------------------------------------------------------!
  subroutine deallocate_hse_ri_fragment(ri_data)
    implicit none
    type(hse_ri_data_t), intent(inout) :: ri_data
    
    if (allocated(ri_data%B_ijP)) deallocate(ri_data%B_ijP)
    if (allocated(ri_data%V_inv_PQ)) deallocate(ri_data%V_inv_PQ)
    if (allocated(ri_data%L_PK)) deallocate(ri_data%L_PK)
    call deallocate_auxiliary_basis(ri_data%aux_basis)
    
    ri_data%initialized = .false.
    
  end subroutine deallocate_hse_ri_fragment

end module xc_hse_ri
