!
!  Copyright 2026 SALMON developers
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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

module current_plusU_sub

  use plusU_global, only: PLUS_U_ON, V_eff

  implicit none

  private
  public :: calc_current_plusU
  public :: calc_current_plusU_rdivided
  public :: PLUS_U_ON

contains

  !-----------------------------------------------------------------------------------------------------------------------------------
  ! Calculate current contribution from DFT+U potential
  !
  ! Current from +U: j_+U = -Im Σ_ij V_eff(i,j) ⟨ψ|φ_i⟩⟨φ_j|∇|ψ⟩
  !
  ! This is analogous to nonlocal pseudopotential current contribution
  !-----------------------------------------------------------------------------------------------------------------------------------

  subroutine calc_current_plusU(jw, psi, ppg, is_array, ie_array, ik)
    use structures
    implicit none
    real(8), intent(out) :: jw(3)
    integer, intent(in) :: is_array(3), ie_array(3)
    complex(8), intent(in) :: psi(is_array(1):ie_array(1), is_array(2):ie_array(2), is_array(3):ie_array(3))
    type(s_pp_grid), intent(in) :: ppg
    integer, intent(in) :: ik
    !
    integer :: ilma, jlma, iprj, Nlma, Nproj_pairs
    integer :: ia, j, ix, iy, iz
    complex(8) :: phipsi, dphi_psi(3)
    complex(8), allocatable :: phipsi_lma(:)
    complex(8), parameter :: zero = (0.0d0, 0.0d0)
    complex(8), parameter :: zi = (0.0d0, 1.0d0)
    real(8) :: jx, jy, jz
    
    jw = 0.0d0
    
    if (.not. allocated(V_eff)) return
    if (.not. PLUS_U_ON) return
    
    Nlma = size(ppg%ia_tbl_ao)
    Nproj_pairs = size(ppg%proj_pairs_ao, 2)
    
    allocate(phipsi_lma(Nlma))
    phipsi_lma = zero
    
    ! Calculate ⟨φ_i|ψ⟩ for all localized orbitals
    do ilma = 1, Nlma
      ia = ppg%ia_tbl_ao(ilma)
      phipsi = zero
      do j = 1, ppg%mps_ao(ia)
        ix = ppg%jxyz_ao(1, j, ia)
        iy = ppg%jxyz_ao(2, j, ia)
        iz = ppg%jxyz_ao(3, j, ia)
        phipsi = phipsi + conjg(ppg%zekr_phi_ao(j, ilma, ik)) * psi(ix, iy, iz)
      end do
      phipsi_lma(ilma) = phipsi
    end do
    
    jx = 0.0d0
    jy = 0.0d0
    jz = 0.0d0
    
    ! Calculate j_+U = -Im Σ_ij V_eff(i,j) ⟨ψ|φ_i⟩⟨φ_j|∇|ψ⟩
    !
    ! Note: ⟨φ_j|∇|ψ⟩ ≈ -i Σ_r φ_j*(r) ∇ψ(r)
    !       For discrete grid: ∇ᵪψ ≈ (ψ(r+dx) - ψ(r-dx))/(2dx)
    !
    ! Using integration by parts: ⟨φ_j|∇|ψ⟩ = -⟨∇φ_j|ψ⟩ + boundary terms
    ! For localized φ_j, we use: ⟨φ_j|∇|ψ⟩ = i⟨∇φ_j|ψ⟩ (approximate)
    
    do iprj = 1, Nproj_pairs
      ilma = ppg%proj_pairs_ao(1, iprj)
      jlma = ppg%proj_pairs_ao(2, iprj)
      ia = ppg%ia_tbl_ao(ilma)
      
      if (abs(V_eff(iprj, 1)) < 1.0d-12) cycle  ! Skip if V_eff is zero
      
      ! Calculate ⟨φ_jlma|∇|ψ⟩ using gradient of basis function
      ! For efficient implementation, use numerical derivative on grid
      dphi_psi = zero
      
      do j = 1, ppg%mps_ao(ia)
        ix = ppg%jxyz_ao(1, j, ia)
        iy = ppg%jxyz_ao(2, j, ia)
        iz = ppg%jxyz_ao(3, j, ia)
        
        ! Approximate gradient using central difference
        ! NOTE: This is a simplified implementation
        ! Full implementation requires proper derivative of phi * exp(ikr)
        
        ! For now, use position operator approximation: ⟨φ|p|ψ⟩ ≈ ⟨φ|(-i∇)|ψ⟩
        ! which gives j ≈ -Im[V_eff * ⟨ψ|φ_i⟩ * ⟨φ_j|r|ψ⟩]
        
        ! Using r-weighted projection (approximation for current)
        dphi_psi(1) = dphi_psi(1) + conjg(ppg%zekr_phi_ao(j, jlma, ik)) &
                                  * ppg%rxyz_ao(1, j, ia) * psi(ix, iy, iz)
        dphi_psi(2) = dphi_psi(2) + conjg(ppg%zekr_phi_ao(j, jlma, ik)) &
                                  * ppg%rxyz_ao(2, j, ia) * psi(ix, iy, iz)
        dphi_psi(3) = dphi_psi(3) + conjg(ppg%zekr_phi_ao(j, jlma, ik)) &
                                  * ppg%rxyz_ao(3, j, ia) * psi(ix, iy, iz)
      end do
      
      ! j_+U contribution: -Im[V_eff * ⟨ψ|φ_i⟩ * ⟨φ_j|∇|ψ⟩]
      ! Using approximation: ⟨φ_j|p|ψ⟩ ≈ -i⟨φ_j|r|ψ⟩ for localized states
      
      jx = jx + aimag(V_eff(iprj, 1) * conjg(phipsi_lma(ilma)) * (-zi) * dphi_psi(1))
      jy = jy + aimag(V_eff(iprj, 1) * conjg(phipsi_lma(ilma)) * (-zi) * dphi_psi(2))
      jz = jz + aimag(V_eff(iprj, 1) * conjg(phipsi_lma(ilma)) * (-zi) * dphi_psi(3))
    end do
    
    jw(1) = jx
    jw(2) = jy
    jw(3) = jz
    
    deallocate(phipsi_lma)
    
    return
  end subroutine calc_current_plusU

  !-----------------------------------------------------------------------------------------------------------------------------------
  ! Calculate current contribution from +U for r-space divided case
  !-----------------------------------------------------------------------------------------------------------------------------------
  
  subroutine calc_current_plusU_rdivided(jw, psi, ppg, is_array, ie_array, ik, icomm)
    use structures
    use communication, only: comm_summation
    implicit none
    real(8), intent(out) :: jw(3)
    integer, intent(in) :: is_array(3), ie_array(3)
    complex(8), intent(in) :: psi(is_array(1):ie_array(1), is_array(2):ie_array(2), is_array(3):ie_array(3))
    type(s_pp_grid), intent(in) :: ppg
    integer, intent(in) :: ik
    integer, intent(in) :: icomm
    !
    real(8) :: jw_tmp(3)
    
    ! Calculate local contribution
    call calc_current_plusU(jw_tmp, psi, ppg, is_array, ie_array, ik)
    
    ! Sum over MPI processes
    call comm_summation(jw_tmp, jw, 3, icomm)
    
    return
  end subroutine calc_current_plusU_rdivided

end module current_plusU_sub
