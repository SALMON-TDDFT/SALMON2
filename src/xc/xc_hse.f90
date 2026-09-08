!
!  Copyright 2024-2025 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
!=======================================================================
! HSE (Heyd-Scuseria-Ernzerhof) Hybrid Functional Module
! Exact exchange computation in real-space basis
!=======================================================================

module xc_hse

  implicit none

  character(256), parameter :: hse_module_version = "1.0.0"

contains

  !=======================================================================
  ! Compute HSE exact exchange contribution to Hamiltonian matrix
  ! Plan A: Density matrix method (direct grid integration)
  ! 
  ! Parameters:
  !   h_mat        - Hamiltonian matrix to be updated [inout]
  !   phi_grid     - Real-space basis functions on grid: phi(ix,iy,iz,io)
  !   occ_states   - Occupied state indices (1:n_occ)
  !   hgs          - Grid spacing vector (3)
  !   hvol         - Grid volume element
  !   hse_alpha    - Exact exchange mixing parameter
  !   is_ij, ie_ij - Grid domain indices for basis evaluation
  !   n_basis      - Number of basis states
  !   n_occ        - Number of occupied states
  !=======================================================================
  subroutine calc_exact_exchange_hse(h_mat, phi_grid, occ_states, hgs, hvol, &
                                      hse_alpha, is_grid, ie_grid, &
                                      n_basis, n_occ)
    implicit none
    
    ! Input/Output
    real(8), intent(inout) :: h_mat(:, :)        ! H_mat(io, jo): matrix to be updated
    real(8), intent(in) :: phi_grid(:, :, :, :)  ! phi_grid(ix,iy,iz,io): basis functions on grid
    integer, intent(in) :: occ_states(:)         ! occ_states(n_occ): occupied state indices
    real(8), intent(in) :: hgs(3)                ! Grid spacing
    real(8), intent(in) :: hvol                  ! Grid volume element
    real(8), intent(in) :: hse_alpha             ! Exchange mixing parameter
    integer, intent(in) :: is_grid(3), ie_grid(3) ! Grid domain
    integer, intent(in) :: n_basis, n_occ
    
    ! Local variables
    integer :: io, jo, ko, lo
    integer :: ix, iy, iz, jx, jy, jz
    real(8) :: V_x_ij, coulomb_integral_ijkl
    real(8) :: distance_ij, coulomb_1r
    real(8) :: r1_rel(3), r2_rel(3)
    integer :: is(3), ie(3)
    integer :: istate_k, istate_l
    
    is = is_grid
    ie = ie_grid
    
    ! ===================================================================
    ! Loop over basis-state pairs (i,j) and occupied pairs (k,l)
    ! Compute exact exchange contribution to Hamiltonian
    ! ===================================================================
    do jo = 1, n_basis
      do io = 1, n_basis
        
        V_x_ij = 0.0d0
        
        ! Loop over occupied orbital pairs (k,l)
        do lo = 1, n_occ
          istate_l = occ_states(lo)
          
          do ko = 1, n_occ
            istate_k = occ_states(ko)
            
            ! Compute two-electron integral <io,jo | 1/r | istate_k,istate_l>
            coulomb_integral_ijkl = 0.0d0
            
            ! Grid integration: 6-fold loop over r1=(ix,iy,iz) and r2=(jx,jy,jz)
            do iz = is(3), ie(3)
              do iy = is(2), ie(2)
                do ix = is(1), ie(1)
                  
                  do jz = is(3), ie(3)
                    do jy = is(2), ie(2)
                      do jx = is(1), ie(1)
                        
                        ! Skip self-interaction
                        if (jx == ix .and. jy == iy .and. jz == iz) cycle
                        
                        ! Compute distance between grid points r1 and r2
                        r1_rel(1) = (ix - 1) * hgs(1)
                        r1_rel(2) = (iy - 1) * hgs(2)
                        r1_rel(3) = (iz - 1) * hgs(3)
                        
                        r2_rel(1) = (jx - 1) * hgs(1)
                        r2_rel(2) = (jy - 1) * hgs(2)
                        r2_rel(3) = (jz - 1) * hgs(3)
                        
                        distance_ij = sqrt((r1_rel(1)-r2_rel(1))**2 + &
                                          (r1_rel(2)-r2_rel(2))**2 + &
                                          (r1_rel(3)-r2_rel(3))**2)
                        
                        ! Avoid division by zero
                        if (distance_ij < 1.0d-10) cycle
                        
                        ! Coulomb kernel: 1/|r1-r2|
                        coulomb_1r = 1.0d0 / distance_ij
                        
                        ! Add contribution: φ_io(r1) φ_istate_k(r1) [1/r12] φ_istate_l(r2) φ_jo(r2)
                        coulomb_integral_ijkl = coulomb_integral_ijkl + &
                          phi_grid(ix, iy, iz, io) * &
                          phi_grid(ix, iy, iz, istate_k) * &
                          coulomb_1r * &
                          phi_grid(jx, jy, jz, istate_l) * &
                          phi_grid(jx, jy, jz, jo) * &
                          hvol * hvol
                        
                      end do
                    end do
                  end do
                  
                end do
              end do
            end do
            
            ! Accumulate exchange contribution to V_x_ij
            ! Double-counting correction for diagonal k==l terms
            if (istate_k == istate_l) then
              V_x_ij = V_x_ij - 0.5d0 * hse_alpha * coulomb_integral_ijkl
            else
              V_x_ij = V_x_ij - hse_alpha * coulomb_integral_ijkl
            end if
            
          end do
        end do
        
        ! Add exchange energy contribution to Hamiltonian matrix
        h_mat(io, jo) = h_mat(io, jo) + V_x_ij
        
      end do
    end do

  end subroutine calc_exact_exchange_hse

  !=======================================================================
  ! Version of HSE exchange for DG-Fragment real-time calculation
  ! Handles fragment-local occupied state indices
  !=======================================================================
  subroutine calc_exact_exchange_hse_fragment(h_mat, phi_frag, ifrag_local, &
                                              hgs, hvol, hse_alpha, &
                                              is_grid, ie_grid, &
                                              n_basis, n_occ, nelec)
    implicit none
    
    real(8), intent(inout) :: h_mat(:, :)
    real(8), intent(in) :: phi_frag(:, :, :, :, :)  ! phi_frag(ix,iy,iz,io,ifrag)
    integer, intent(in) :: ifrag_local              ! Local fragment index
    real(8), intent(in) :: hgs(3)
    real(8), intent(in) :: hvol
    real(8), intent(in) :: hse_alpha
    integer, intent(in) :: is_grid(3), ie_grid(3)
    integer, intent(in) :: n_basis, n_occ, nelec
    
    real(8), allocatable :: phi_grid(:, :, :, :)
    integer, allocatable :: occ_states(:)
    integer :: io, is(3), ie(3)
    
    is = is_grid
    ie = ie_grid
    
    ! Create local phi_grid array from phi_frag interior (1:nxyz_domain).
    ! phi_frag includes halo regions (indices 1-nb:nxyz+nb); use explicit
    ! is:ie indices to extract only the interior without shape mismatch.
    allocate(phi_grid(is(1):ie(1), is(2):ie(2), is(3):ie(3), n_basis))

    do io = 1, n_basis
      phi_grid(is(1):ie(1), is(2):ie(2), is(3):ie(3), io) = &
        phi_frag(is(1):ie(1), is(2):ie(2), is(3):ie(3), io, ifrag_local)
    end do
    
    ! Assume first n_occ states are occupied
    allocate(occ_states(n_occ))
    do io = 1, n_occ
      occ_states(io) = io
    end do
    
    ! Call generic routine
    call calc_exact_exchange_hse(h_mat, phi_grid, occ_states, hgs, hvol, &
                                  hse_alpha, is_grid, ie_grid, &
                                  n_basis, n_occ)
    
    deallocate(phi_grid, occ_states)
    
  end subroutine calc_exact_exchange_hse_fragment

  !=======================================================================
  ! Compute HSE exact exchange energy and potential (Vxc contribution)
  ! Density matrix approach: suitable for ground-state (GS) calculations
  ! 
  ! This version takes density matrix elements instead of orbital wavefunctions
  ! and computes both the energy and the potential contribution
  !=======================================================================
  subroutine calc_xc_hse_from_rho(rho_grid, vxc_grid, e_xc_hse, &
                                   hgs, hvol, hse_alpha, &
                                   is_grid, ie_grid, &
                                   mg)
    use structures, only: s_rgrid
    implicit none
    
    ! Charge density on grid (input)
    real(8), intent(in)  :: rho_grid(:, :, :)
    ! XC potential to be updated (output)
    real(8), intent(inout) :: vxc_grid(:, :, :)
    ! XC energy (output)
    real(8), intent(out) :: e_xc_hse
    ! Grid parameters
    real(8), intent(in) :: hgs(3)
    real(8), intent(in) :: hvol
    real(8), intent(in) :: hse_alpha
    integer, intent(in) :: is_grid(3), ie_grid(3)
    type(s_rgrid), intent(in) :: mg
    
    ! Local variables
    integer :: ix, iy, iz, jx, jy, jz
    real(8) :: distance, coulomb_1r, r1_rel(3), r2_rel(3)
    real(8) :: rho_ij, rho_ji, coulomb_integral
    integer :: is(3), ie(3)
    
    is = is_grid
    ie = ie_grid
    e_xc_hse = 0.0d0
    
    ! ===================================================================
    ! Compute HSE exchange energy and potential from charge density
    ! V_x^HSE(r) = - hse_alpha * ∫ ρ(r') / |r - r'| d³r'
    ! E_x^HSE = - hse_alpha/2 * ∫∫ ρ(r) ρ(r') / |r - r'| dr dr'
    ! ===================================================================
    
    ! Initialize Vxc contribution
    vxc_grid(:, :, :) = 0.0d0
    
    ! Loop over grid points r
    do iz = is(3), ie(3)
      do iy = is(2), ie(2)
        do ix = is(1), ie(1)
          
          rho_ij = rho_grid(ix, iy, iz)
          
          ! Loop over grid points r'
          do jz = is(3), ie(3)
            do jy = is(2), ie(2)
              do jx = is(1), ie(1)
                
                ! Skip self-interaction (or Small distance)
                if (jx == ix .and. jy == iy .and. jz == iz) cycle
                
                ! Compute distance |r - r'|
                r1_rel(1) = (ix - 1) * hgs(1)
                r1_rel(2) = (iy - 1) * hgs(2)
                r1_rel(3) = (iz - 1) * hgs(3)
                
                r2_rel(1) = (jx - 1) * hgs(1)
                r2_rel(2) = (jy - 1) * hgs(2)
                r2_rel(3) = (jz - 1) * hgs(3)
                
                distance = sqrt((r1_rel(1)-r2_rel(1))**2 + &
                               (r1_rel(2)-r2_rel(2))**2 + &
                               (r1_rel(3)-r2_rel(3))**2)
                
                if (distance < 1.0d-10) cycle
                
                ! Coulomb kernel
                coulomb_1r = 1.0d0 / distance
                
                ! Accumulate potential at point r due to charge at r'
                rho_ji = rho_grid(jx, jy, jz)
                vxc_grid(ix, iy, iz) = vxc_grid(ix, iy, iz) - &
                                       hse_alpha * rho_ji * coulomb_1r * hvol
                
                ! Accumulate energy: E_x ∝ ρ(r) * [∫ ρ(r')/|r-r'| dr']
                e_xc_hse = e_xc_hse - 0.5d0 * hse_alpha * &
                          rho_ij * rho_ji * coulomb_1r * hvol * hvol
                
              end do
            end do
          end do
          
        end do
      end do
    end do
    
  end subroutine calc_xc_hse_from_rho

  !=======================================================================
  ! FFT-based HSE exact exchange (using existing Hartree solver)
  ! For ground state calculations 
  ! 
  ! Reuses Hartree/Poisson infrastructure for Coulomb convolution:
  ! V_x^HSE(r) = -hse_alpha * hartree(-ρ_exchange, r)
  ! 
  ! Parameters:
  !   rho_grid   - Input charge density (for exchange)
  !   vxc_grid   - Output exchange potential (to be updated)
  !   e_xc_hse   - Output exchange energy
  !   system     - DFT system structure (contains grid info)
  !   lg, mg     - Local and global grid structures
  !   info       - Parallel info
  !   fg         - Reciprocal grid
  !   poisson    - Poisson solver structure
  !   srg_scalar - Send/receive grid
  !   stencil    - Stencil for finite difference
  !   hse_alpha  - Exchange mixing parameter
  !=======================================================================
  subroutine calc_xc_hse_fft(rho_grid, vxc_grid, e_xc_hse, &
                             system, lg, mg, info, fg, poisson, &
                             srg_scalar, stencil, hse_alpha)
    use structures, only: s_dft_system, s_rgrid, s_parallel_info, &
                         s_reciprocal_grid, s_poisson, s_sendrecv_grid, &
                         s_stencil, s_scalar, allocate_scalar, deallocate_scalar
    use hartree_sub, only: hartree
    use communication, only: comm_summation
    implicit none
    
    real(8), intent(in) :: rho_grid(:, :, :)
    real(8), intent(inout) :: vxc_grid(:, :, :)
    real(8), intent(out) :: e_xc_hse
    type(s_dft_system), intent(in) :: system
    type(s_rgrid), intent(in) :: lg, mg
    type(s_parallel_info), intent(in) :: info
    type(s_reciprocal_grid), intent(in) :: fg
    type(s_poisson), intent(inout) :: poisson
    type(s_sendrecv_grid), intent(inout) :: srg_scalar
    type(s_stencil), intent(in) :: stencil
    real(8), intent(in) :: hse_alpha
    
    ! Local variables
    type(s_scalar) :: rho_hse, vh_hse
    integer :: ix, iy, iz
    real(8) :: hvol, e_xc_hse_local
    
    ! ===================================================================
    ! Implementation:
    ! 1. Create temporary density field for Hartree solver
    ! 2. Call existing Hartree/Poisson solver (FFT-based)
    ! 3. Extract result and add to Vxc with -hse_alpha scaling
    ! ===================================================================
    
    hvol = system%hvol
    
    ! Allocate temporary scalar fields
    call allocate_scalar(mg, rho_hse)
    call allocate_scalar(mg, vh_hse)
    
    ! Copy input density to Hartree field
    rho_hse%f(:, :, :) = rho_grid(:, :, :)
    
    ! Compute Coulomb potential via existing Hartree solver (FFT)
    ! This solves: ∇²Vh = 4π*rho  (in atomic units)
    call hartree(lg, mg, info, system, fg, poisson, srg_scalar, &
                 stencil, rho_hse, vh_hse)
    
    ! Add to XC potential with exchange mixing factor
    ! V_x^HSE += -hse_alpha * Vh_coulomb
    vxc_grid(:, :, :) = vxc_grid(:, :, :) - hse_alpha * vh_hse%f(:, :, :)
    
    ! Compute exchange energy: E_x = -hse_alpha/2 * ∫ ρ * Vh_coulomb d³r
    e_xc_hse_local = 0.0d0
    do iz = mg%is(3), mg%ie(3)
      do iy = mg%is(2), mg%ie(2)
        do ix = mg%is(1), mg%ie(1)
          e_xc_hse_local = e_xc_hse_local - 0.5d0 * hse_alpha * &
                     rho_grid(ix, iy, iz) * vh_hse%f(ix, iy, iz) * hvol
        end do
      end do
    end do
    
    ! Sum over all MPI processes
    call comm_summation(e_xc_hse_local, e_xc_hse, info%icomm_r)
    
    ! Clean up
    call deallocate_scalar(rho_hse)
    call deallocate_scalar(vh_hse)
    
  end subroutine calc_xc_hse_fft

end module xc_hse
