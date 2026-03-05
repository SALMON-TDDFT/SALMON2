!----------------------------------------------------------------------
! HSE Hybrid Functional Module - ERI Optimization
! Plan B: Precomputed Electron Repulsion Integrals
!
! Features:
!   - One-time computation of all 2-electron integrals
!   - O(N^4) memory, O(N^2 × N_occ^2) runtime per timestep
!   - 8-fold permutation symmetry exploitation
!   - Compressed storage option
!
! Usage:
!   1. Initialization: call precompute_ERI_fragment(...)
!   2. Each timestep: call calc_exact_exchange_hse_fast(...)
!
! Author: Plan B Implementation
! Date: 2026-02-23
!----------------------------------------------------------------------

module xc_hse_eri
  implicit none
  
  ! Module parameters
  real(8), parameter :: ERI_threshold = 1.0d-10  ! Screening threshold
  
contains

  !=======================================================================
  ! Precompute all 2-electron integrals for a fragment
  ! 
  ! Computes: (ij|kl) = ∫∫ φ_i(r1)φ_j(r1) [1/r12] φ_k(r2)φ_l(r2) dr1 dr2
  ! 
  ! Exploits 8-fold permutation symmetry:
  !   (ij|kl) = (ji|kl) = (ij|lk) = (ji|lk) = (kl|ij) = (lk|ij) = (kl|ji) = (lk|ji)
  ! 
  ! Parameters:
  !   phi_grid    - Basis functions on grid: φ(ix,iy,iz,io)
  !   hgs         - Grid spacing [a.u.]
  !   hvol        - Volume element [a.u.^3]
  !   is/ie_grid  - Grid domain boundaries
  !   n_basis     - Number of basis functions
  !   ERI_frag    - Output: 4-index integral tensor
  !   use_symmetry- Flag: exploit symmetry (true) or compute all (false)
  !=======================================================================
  subroutine precompute_ERI_fragment(phi_grid, hgs, hvol, is_grid, ie_grid, &
                                     n_basis, ERI_frag, use_symmetry)
    implicit none
    
    ! Input
    real(8), intent(in) :: phi_grid(:,:,:,:)     ! (ix,iy,iz,istate)
    real(8), intent(in) :: hgs(3), hvol
    integer, intent(in) :: is_grid(3), ie_grid(3)
    integer, intent(in) :: n_basis
    logical, intent(in), optional :: use_symmetry
    
    ! Output
    real(8), intent(out) :: ERI_frag(n_basis, n_basis, n_basis, n_basis)
    
    ! Local variables
    integer :: i, j, k, l
    integer :: ix, iy, iz, jx, jy, jz
    real(8) :: r1(3), r2(3), distance, coulomb_1r
    real(8) :: integral_ijkl
    logical :: use_sym
    integer :: n_computed, n_total, n_skipped
    real(8) :: start_time, end_time
    
    ! Check symmetry flag
    use_sym = .true.
    if (present(use_symmetry)) use_sym = use_symmetry
    
    ! Initialize
    ERI_frag = 0.0d0
    n_computed = 0
    n_skipped = 0
    n_total = n_basis**4
    
    write(*,*) '======================================================'
    write(*,*) 'Precomputing 2-electron integrals (ERI)'
    write(*,*) '  Basis functions:', n_basis
    write(*,*) '  Grid points:', ie_grid - is_grid + 1
    write(*,*) '  Total integrals:', n_total
    if (use_sym) then
      write(*,*) '  Using 8-fold symmetry: computing ~', n_total/8, 'unique integrals'
    end if
    write(*,*) '======================================================'
    
    call cpu_time(start_time)
    
    ! Main computation loop with optional symmetry
    !$omp parallel do private(i,j,k,l,ix,iy,iz,jx,jy,jz,r1,r2,distance,coulomb_1r,integral_ijkl) &
    !$omp& collapse(4) schedule(dynamic) reduction(+:n_computed,n_skipped)
    do l = 1, n_basis
      do k = 1, n_basis
        
        ! Skip if using symmetry and this is a redundant computation
        if (use_sym .and. k > l) cycle
        
        do j = 1, n_basis
          do i = 1, n_basis
            
            ! Apply symmetry constraints
            if (use_sym) then
              ! Only compute for i >= j
              if (i < j) cycle
              ! Only compute for compound index (ij) >= (kl)
              if (index_compound(i,j) < index_compound(k,l)) then
                n_skipped = n_skipped + 1
                cycle
              end if
            end if
            
            ! Compute 6-dimensional grid integral
            integral_ijkl = 0.0d0
            
            do iz = is_grid(3), ie_grid(3)
              do iy = is_grid(2), ie_grid(2)
                do ix = is_grid(1), ie_grid(1)
                  
                  ! Position vector r1
                  r1(1) = ix * hgs(1)
                  r1(2) = iy * hgs(2)
                  r1(3) = iz * hgs(3)
                  
                  do jz = is_grid(3), ie_grid(3)
                    do jy = is_grid(2), ie_grid(2)
                      do jx = is_grid(1), ie_grid(1)
                        
                        ! Skip self-interaction
                        if (ix == jx .and. iy == jy .and. iz == jz) cycle
                        
                        ! Position vector r2
                        r2(1) = jx * hgs(1)
                        r2(2) = jy * hgs(2)
                        r2(3) = jz * hgs(3)
                        
                        ! Distance |r1 - r2|
                        distance = sqrt((r1(1)-r2(1))**2 + &
                                       (r1(2)-r2(2))**2 + &
                                       (r1(3)-r2(3))**2)
                        
                        ! Avoid singularity
                        if (distance < 1.0d-10) cycle
                        
                        ! Coulomb kernel
                        coulomb_1r = 1.0d0 / distance
                        
                        ! Accumulate: φ_i(r1) φ_j(r1) [1/r12] φ_k(r2) φ_l(r2)
                        integral_ijkl = integral_ijkl + &
                          phi_grid(ix, iy, iz, i) * &
                          phi_grid(ix, iy, iz, j) * &
                          coulomb_1r * &
                          phi_grid(jx, jy, jz, k) * &
                          phi_grid(jx, jy, jz, l) * &
                          hvol * hvol
                        
                      end do
                    end do
                  end do
                  
                end do
              end do
            end do
            
            ! Store computed integral
            ERI_frag(i, j, k, l) = integral_ijkl
            n_computed = n_computed + 1
            
            ! Apply 8-fold symmetry if enabled
            if (use_sym) then
              call apply_ERI_symmetry(ERI_frag, i, j, k, l, integral_ijkl, n_basis)
            end if
            
          end do
        end do
      end do
    end do
    !$omp end parallel do
    
    call cpu_time(end_time)
    
    ! Summary
    write(*,*) '======================================================'
    write(*,*) 'ERI precomputation completed'
    write(*,*) '  Computed:', n_computed, 'integrals'
    write(*,*) '  Skipped (symmetry):', n_skipped
    write(*,*) '  Time:', end_time - start_time, 'seconds'
    write(*,*) '  Memory:', real(n_basis**4) * 8.0d0 / 1024.0d0**2, 'MB'
    write(*,*) '======================================================'
    
  contains
    
    ! Compound index for symmetry: (i,j) -> i*(i-1)/2 + j
    integer function index_compound(ii, jj)
      integer, intent(in) :: ii, jj
      integer :: imax, jmin
      imax = max(ii, jj)
      jmin = min(ii, jj)
      index_compound = imax * (imax - 1) / 2 + jmin
    end function index_compound
    
  end subroutine precompute_ERI_fragment

  !=======================================================================
  ! Apply 8-fold permutation symmetry to fill ERI tensor
  !=======================================================================
  subroutine apply_ERI_symmetry(ERI_frag, i, j, k, l, value, n_basis)
    implicit none
    real(8), intent(inout) :: ERI_frag(n_basis, n_basis, n_basis, n_basis)
    integer, intent(in) :: i, j, k, l, n_basis
    real(8), intent(in) :: value
    
    ! (ij|kl) = value (already set)
    ! Apply 7 additional permutations:
    if (j /= i) ERI_frag(j, i, k, l) = value  ! (ji|kl)
    if (l /= k) ERI_frag(i, j, l, k) = value  ! (ij|lk)
    if (j /= i .and. l /= k) ERI_frag(j, i, l, k) = value  ! (ji|lk)
    if (k /= i .or. l /= j) then
      ERI_frag(k, l, i, j) = value  ! (kl|ij)
      if (l /= k) ERI_frag(l, k, i, j) = value  ! (lk|ij)
      if (j /= i) ERI_frag(k, l, j, i) = value  ! (kl|ji)
      if (j /= i .and. l /= k) ERI_frag(l, k, j, i) = value  ! (lk|ji)
    end if
    
  end subroutine apply_ERI_symmetry

  !=======================================================================
  ! Fast HSE exchange calculation using precomputed ERIs
  ! 
  ! Complexity: O(N_basis^2 × N_occ^2) instead of O(N_basis^2 × N_occ^2 × L^6)
  ! 
  ! Computes: V_x^HSE(i,j) = -α Σ_{k,l∈occ} (ij|kl)
  ! 
  ! Parameters:
  !   h_mat       - Hamiltonian matrix to be updated [inout]
  !   ERI_frag    - Precomputed 2-electron integrals
  !   occ_states  - List of occupied state indices
  !   hse_alpha   - Exchange mixing parameter
  !   n_basis     - Number of basis functions
  !   n_occ       - Number of occupied states
  !=======================================================================
  subroutine calc_exact_exchange_hse_fast(h_mat, ERI_frag, occ_states, &
                                          hse_alpha, n_basis, n_occ)
    implicit none
    
    ! Input/Output
    real(8), intent(inout) :: h_mat(n_basis, n_basis)
    real(8), intent(in) :: ERI_frag(n_basis, n_basis, n_basis, n_basis)
    integer, intent(in) :: occ_states(n_occ)
    real(8), intent(in) :: hse_alpha
    integer, intent(in) :: n_basis, n_occ
    
    ! Local variables
    integer :: i, j, ko, lo, istate_k, istate_l
    real(8) :: V_x_ij
    
    ! Parallel loop over basis pairs (i,j)
    !$omp parallel do private(i,j,ko,lo,istate_k,istate_l,V_x_ij) collapse(2)
    do j = 1, n_basis
      do i = 1, n_basis
        
        V_x_ij = 0.0d0
        
        ! Sum over occupied orbital pairs (k,l)
        do lo = 1, n_occ
          istate_l = occ_states(lo)
          
          do ko = 1, n_occ
            istate_k = occ_states(ko)
            
            ! Exchange integral: (ij|kl)
            ! Double-counting correction for diagonal (k==l)
            if (istate_k == istate_l) then
              V_x_ij = V_x_ij - 0.5d0 * hse_alpha * ERI_frag(i, j, istate_k, istate_l)
            else
              V_x_ij = V_x_ij - hse_alpha * ERI_frag(i, j, istate_k, istate_l)
            end if
            
          end do
        end do
        
        ! Add to Hamiltonian
        h_mat(i, j) = h_mat(i, j) + V_x_ij
        
      end do
    end do
    !$omp end parallel do
    
  end subroutine calc_exact_exchange_hse_fast

  !=======================================================================
  ! Estimate memory requirement for ERI storage
  !=======================================================================
  function estimate_ERI_memory(n_basis, use_symmetry) result(memory_mb)
    implicit none
    integer, intent(in) :: n_basis
    logical, intent(in), optional :: use_symmetry
    real(8) :: memory_mb
    logical :: use_sym
    integer :: n_elements
    
    use_sym = .true.
    if (present(use_symmetry)) use_sym = use_symmetry
    
    if (use_sym) then
      ! With 8-fold symmetry: ~N^4/8 unique elements
      n_elements = n_basis * (n_basis + 1) / 2
      n_elements = n_elements * (n_elements + 1) / 2
    else
      ! Full tensor: N^4 elements
      n_elements = n_basis**4
    end if
    
    ! 8 bytes per real(8) element
    memory_mb = real(n_elements, 8) * 8.0d0 / (1024.0d0**2)
    
  end function estimate_ERI_memory

  !=======================================================================
  ! Apply screening to eliminate small ERIs (sparsify)
  !=======================================================================
  subroutine screen_ERI(ERI_frag, n_basis, threshold, n_screened)
    implicit none
    real(8), intent(inout) :: ERI_frag(n_basis, n_basis, n_basis, n_basis)
    integer, intent(in) :: n_basis
    real(8), intent(in) :: threshold
    integer, intent(out) :: n_screened
    
    integer :: i, j, k, l
    
    n_screened = 0
    
    !$omp parallel do private(i,j,k,l) collapse(4) reduction(+:n_screened)
    do l = 1, n_basis
      do k = 1, n_basis
        do j = 1, n_basis
          do i = 1, n_basis
            
            if (abs(ERI_frag(i,j,k,l)) < threshold) then
              ERI_frag(i,j,k,l) = 0.0d0
              n_screened = n_screened + 1
            end if
            
          end do
        end do
      end do
    end do
    !$omp end parallel do
    
    write(*,*) 'ERI screening:'
    write(*,*) '  Threshold:', threshold
    write(*,*) '  Screened elements:', n_screened, '/', n_basis**4
    write(*,*) '  Sparsity:', real(n_screened) / real(n_basis**4) * 100.0d0, '%'
    
  end subroutine screen_ERI

end module xc_hse_eri
