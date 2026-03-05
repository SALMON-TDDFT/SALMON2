!=======================================================================
! DFT+U support routines for DG-Fragment RT-TDDFT
!=======================================================================
! This module provides DFT+U correction for fragment basis time evolution
! to prevent spurious metallic response in Mott insulators.
!
! Key features:
! - +U density matrix calculation from fragment coefficients
! - Hamiltonian correction with +U term
! - Self-consistent +U update during time evolution
!=======================================================================

module plusu_fragment_support
  use structures
  use plusU_global, only: PLUS_U_ON, U_eff, V_eff, dm_mms_nla
  implicit none
  
  private
  public :: init_plusu_fragment_framework
  public :: update_plusu_hamiltonian_fragment
  public :: calculate_plusu_density_matrix_fragment
  
contains

  !=======================================================================
  ! Initialize DFT+U framework for fragment basis
  !=======================================================================
  subroutine init_plusu_fragment_framework(n_mat_max, nspin, U_mat)
    implicit none
    integer, intent(in) :: n_mat_max, nspin
    real(8), allocatable, intent(out) :: U_mat(:,:,:)
    
    if (.not. PLUS_U_ON) return
    
    ! Allocate +U matrix for Hamiltonian correction
    allocate(U_mat(n_mat_max, n_mat_max, nspin))
    U_mat = 0.0d0
    
    write(*,*)
    write(*,*) "  *** DFT+U initialized for fragment basis ***"
    write(*,*) "  Note: Full +U implementation requires:"
    write(*,*) "    1. Orbital character analysis from DC-LCFO"
    write(*,*) "    2. Density matrix projection onto localized orbitals"
    write(*,*) "    3. Self-consistent +U update during time evolution"
    write(*,*)
    write(*,*) "  Current status: Framework ready"
    write(*,*) "  Action required: Implement orbital projection in DC-LCFO output"
    write(*,*)
    
  end subroutine init_plusu_fragment_framework

  !=======================================================================
  ! Calculate +U density matrix from fragment coefficients
  !=======================================================================
  subroutine calculate_plusu_density_matrix_fragment(coef, rocc, nspin, nstate_tot, nstate_frag, dm)
    implicit none
    complex(8), intent(in) :: coef(:,:,:)  ! (nstate_frag, nstate_tot, nspin)
    real(8), intent(in) :: rocc(:,:)       ! occupation numbers
    integer, intent(in) :: nspin, nstate_tot, nstate_frag
    complex(8), allocatable, intent(out) :: dm(:,:,:)  ! density matrix
    
    integer :: io, jo, ko, ispin
    real(8) :: occ_factor
    
    ! Allocate density matrix
    if (.not. allocated(dm)) allocate(dm(nstate_frag, nstate_frag, nspin))
    dm = (0.0d0, 0.0d0)
    
    ! Calculate density matrix: dm(j,k) = sum_i occ(i) * c_j^*(i) * c_k(i)
    ! This gives occupation of fragment basis states
    !$omp parallel do collapse(2) private(io,jo,ko,occ_factor)
    do ispin = 1, nspin
      do io = 1, nstate_tot
        if (rocc(io, ispin) < 1.0d-10) cycle
        occ_factor = rocc(io, ispin)
        
        do jo = 1, nstate_frag
          do ko = 1, nstate_frag
            dm(jo, ko, ispin) = dm(jo, ko, ispin) + &
                                occ_factor * conjg(coef(jo, io, ispin)) * coef(ko, io, ispin)
          end do
        end do
      end do
    end do
    !$omp end parallel do
    
  end subroutine calculate_plusu_density_matrix_fragment

  !=======================================================================
  ! Update Hamiltonian with +U correction
  !=======================================================================
  subroutine update_plusu_hamiltonian_fragment(H_mat, U_mat, dm, U_eff_param, n_mat_max, nspin)
    implicit none
    real(8), intent(inout) :: H_mat(:,:,:)
    real(8), intent(inout) :: U_mat(:,:,:)
    complex(8), intent(in) :: dm(:,:,:)
    real(8), intent(in) :: U_eff_param  ! Effective U parameter
    integer, intent(in) :: n_mat_max, nspin
    
    integer :: io, jo, ispin
    real(8) :: U_correction
    
    if (.not. PLUS_U_ON) return
    
    ! DFT+U correction: V_U = U * (1/2 - n)
    ! where n is the density matrix
    ! Full expression: H_ij += U * (delta_ij * (1/2 - sum_k n_kk) - n_ij)
    !
    ! Simplified for fragment basis:
    ! +U opens a gap by splitting occupied and unoccupied states
    
    do ispin = 1, nspin
      ! +U potential matrix element
      ! This is a placeholder - full implementation requires:
      ! 1. Projection of dm onto localized d/f orbitals
      ! 2. Rotationally invariant +U formulation
      ! 3. Self-consistent iteration

      ! Diagonal: shift energy based on occupation
      ! More occupied => higher energy (Coulomb repulsion)
      do io = 1, n_mat_max
        U_correction = U_eff_param * (0.5d0 - real(dm(io, io, ispin)))
        U_mat(io, io, ispin) = U_correction
        H_mat(io, io, ispin) = H_mat(io, io, ispin) + U_correction
      end do

      ! Off-diagonal: suppress hybridization
      do jo = 1, n_mat_max
        do io = 1, jo - 1
          U_correction = -U_eff_param * real(dm(io, jo, ispin))
          U_mat(io, jo, ispin) = U_correction
          H_mat(io, jo, ispin) = H_mat(io, jo, ispin) + U_correction
        end do
        do io = jo + 1, n_mat_max
          U_correction = -U_eff_param * real(dm(io, jo, ispin))
          U_mat(io, jo, ispin) = U_correction
          H_mat(io, jo, ispin) = H_mat(io, jo, ispin) + U_correction
        end do
      end do
    end do
    
  end subroutine update_plusu_hamiltonian_fragment

end module plusu_fragment_support
