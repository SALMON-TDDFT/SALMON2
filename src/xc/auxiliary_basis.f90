!========================================================================================!
! Module: auxiliary_basis
! 
! Purpose: Generate and manage auxiliary basis functions for RI/DF approximation
!          Used in Plan C (RI/DF) implementation of HSE exchange
!
! Author: Plan C Implementation Team
! Date: 2026-02-23
!========================================================================================!

module auxiliary_basis
  implicit none
  
  private
  public :: init_auxiliary_basis
  public :: calc_auxiliary_function
  public :: get_n_auxiliary
  public :: deallocate_auxiliary_basis
  
  ! Auxiliary basis data structure
  type, public :: auxiliary_basis_t
    integer :: n_aux                    ! Total number of auxiliary functions
    integer :: n_atom                   ! Number of atoms in fragment
    integer, allocatable :: atom_index(:)    ! Atom index for each aux function
    integer, allocatable :: l_quantum(:)     ! Angular momentum l for each aux function
    integer, allocatable :: m_quantum(:)     ! Magnetic quantum number m
    real(8), allocatable :: alpha(:)         ! Exponent for each Gaussian
    real(8), allocatable :: center(:,:)      ! Center position (3, n_aux)
  end type auxiliary_basis_t
  
contains

  !--------------------------------------------------------------------------------------!
  ! Initialize auxiliary basis set for a fragment
  ! Strategy: even-tempered Gaussian basis with s, p, d functions
  !--------------------------------------------------------------------------------------!
  subroutine init_auxiliary_basis(aux_basis, natom, atom_coords, atom_types, ratio)
    implicit none
    type(auxiliary_basis_t), intent(inout) :: aux_basis
    integer, intent(in) :: natom
    real(8), intent(in) :: atom_coords(3, natom)
    integer, intent(in) :: atom_types(natom)
    real(8), intent(in), optional :: ratio  ! Expansion ratio (default: 2.5)
    
    integer :: iatom, l, m, i_exp, idx
    integer :: n_exp_s, n_exp_p, n_exp_d
    real(8) :: alpha_start_s, alpha_start_p, alpha_start_d
    real(8) :: alpha_ratio
    real(8) :: alpha_val
    
    ! Parameters for even-tempered basis
    if (present(ratio)) then
      alpha_ratio = ratio
    else
      alpha_ratio = 2.5d0  ! Default expansion ratio
    end if
    
    ! Number of exponents per angular momentum
    ! For N_aux ~ 3 × N_basis, use: s(5), p(4), d(3) per atom
    n_exp_s = 5
    n_exp_p = 4
    n_exp_d = 3
    
    ! Starting exponents (optimized for typical molecular systems)
    alpha_start_s = 0.05d0
    alpha_start_p = 0.10d0
    alpha_start_d = 0.20d0
    
    ! Calculate total number of auxiliary functions
    ! s: 1 function, p: 3 functions, d: 5 functions per exponent
    aux_basis%n_atom = natom
    if (natom <= 0) then
      aux_basis%n_aux = 0
      allocate(aux_basis%atom_index(0))
      allocate(aux_basis%l_quantum(0))
      allocate(aux_basis%m_quantum(0))
      allocate(aux_basis%alpha(0))
      allocate(aux_basis%center(3, 0))
      return
    end if
    aux_basis%n_aux = natom * (n_exp_s * 1 + n_exp_p * 3 + n_exp_d * 5)
    
    ! Allocate arrays
    allocate(aux_basis%atom_index(aux_basis%n_aux))
    allocate(aux_basis%l_quantum(aux_basis%n_aux))
    allocate(aux_basis%m_quantum(aux_basis%n_aux))
    allocate(aux_basis%alpha(aux_basis%n_aux))
    allocate(aux_basis%center(3, aux_basis%n_aux))
    
    ! Generate auxiliary functions
    idx = 0
    
    do iatom = 1, natom
      ! s functions (l=0, m=0)
      do i_exp = 1, n_exp_s
        alpha_val = alpha_start_s * (alpha_ratio ** (i_exp - 1))
        idx = idx + 1
        aux_basis%atom_index(idx) = iatom
        aux_basis%l_quantum(idx) = 0
        aux_basis%m_quantum(idx) = 0
        aux_basis%alpha(idx) = alpha_val
        aux_basis%center(:, idx) = atom_coords(:, iatom)
      end do
      
      ! p functions (l=1, m=-1,0,+1)
      do i_exp = 1, n_exp_p
        alpha_val = alpha_start_p * (alpha_ratio ** (i_exp - 1))
        do m = -1, 1
          idx = idx + 1
          aux_basis%atom_index(idx) = iatom
          aux_basis%l_quantum(idx) = 1
          aux_basis%m_quantum(idx) = m
          aux_basis%alpha(idx) = alpha_val
          aux_basis%center(:, idx) = atom_coords(:, iatom)
        end do
      end do
      
      ! d functions (l=2, m=-2,-1,0,+1,+2)
      do i_exp = 1, n_exp_d
        alpha_val = alpha_start_d * (alpha_ratio ** (i_exp - 1))
        do m = -2, 2
          idx = idx + 1
          aux_basis%atom_index(idx) = iatom
          aux_basis%l_quantum(idx) = 2
          aux_basis%m_quantum(idx) = m
          aux_basis%alpha(idx) = alpha_val
          aux_basis%center(:, idx) = atom_coords(:, iatom)
        end do
      end do
    end do
    
    if (idx /= aux_basis%n_aux) then
      write(*, '(A, I0, A, I0)') 'ERROR: aux basis count mismatch: ', idx, ' vs ', aux_basis%n_aux
      stop
    end if
    
  end subroutine init_auxiliary_basis

  !--------------------------------------------------------------------------------------!
  ! Calculate auxiliary basis function value at a point
  ! χ_P(r) = N × r^l × exp(-α |r-R|²) × Y_lm(θ,φ)
  !--------------------------------------------------------------------------------------!
  function calc_auxiliary_function(aux_basis, idx, position) result(value)
    implicit none
    type(auxiliary_basis_t), intent(in) :: aux_basis
    integer, intent(in) :: idx
    real(8), intent(in) :: position(3)
    real(8) :: value
    
    real(8) :: r_vec(3), r2, r, theta, phi
    real(8) :: norm_factor, radial_part, angular_part
    integer :: l, m
    real(8) :: alpha
    real(8), parameter :: pi = 3.14159265358979323846d0
    
    ! Relative position
    r_vec = position - aux_basis%center(:, idx)
    r2 = sum(r_vec**2)
    r = sqrt(r2)
    
    ! Get quantum numbers
    l = aux_basis%l_quantum(idx)
    m = aux_basis%m_quantum(idx)
    alpha = aux_basis%alpha(idx)
    
    ! Normalization factor for Gaussian
    ! N = (2α/π)^(3/4) × (4α)^(l/2)
    norm_factor = (2.0d0 * alpha / pi)**(0.75d0) * (4.0d0 * alpha)**(0.5d0 * l)
    
    ! Radial part: r^l × exp(-α r²)
    if (r < 1.0d-10) then
      if (l == 0) then
        radial_part = norm_factor * exp(-alpha * r2)
      else
        radial_part = 0.0d0
      end if
    else
      radial_part = norm_factor * (r**l) * exp(-alpha * r2)
    end if
    
    ! Angular part: Y_lm(θ, φ) (simplified for s, p, d)
    angular_part = calc_spherical_harmonic(l, m, r_vec, r)
    
    value = radial_part * angular_part
    
  end function calc_auxiliary_function

  !--------------------------------------------------------------------------------------!
  ! Calculate spherical harmonic Y_lm (simplified for l=0,1,2)
  !--------------------------------------------------------------------------------------!
  function calc_spherical_harmonic(l, m, r_vec, r) result(ylm)
    implicit none
    integer, intent(in) :: l, m
    real(8), intent(in) :: r_vec(3)
    real(8), intent(in) :: r
    real(8) :: ylm
    
    real(8) :: x, y, z
    real(8), parameter :: pi = 3.14159265358979323846d0
    real(8) :: norm
    
    x = r_vec(1)
    y = r_vec(2)
    z = r_vec(3)
    
    select case(l)
    case(0)  ! s function
      ylm = 1.0d0 / sqrt(4.0d0 * pi)
      
    case(1)  ! p functions
      norm = sqrt(3.0d0 / (4.0d0 * pi))
      if (r < 1.0d-10) then
        ylm = 0.0d0
      else
        select case(m)
        case(-1)
          ylm = norm * y / r  ! p_y
        case(0)
          ylm = norm * z / r  ! p_z
        case(1)
          ylm = norm * x / r  ! p_x
        end select
      end if
      
    case(2)  ! d functions
      norm = sqrt(15.0d0 / (4.0d0 * pi))
      if (r < 1.0d-10) then
        ylm = 0.0d0
      else
        select case(m)
        case(-2)
          ylm = 0.5d0 * norm * x * y / r**2  ! d_xy
        case(-1)
          ylm = norm * y * z / r**2  ! d_yz
        case(0)
          ylm = sqrt(5.0d0/(16.0d0*pi)) * (3.0d0*z**2 - r**2) / r**2  ! d_3z²-r²
        case(1)
          ylm = norm * x * z / r**2  ! d_xz
        case(2)
          ylm = 0.5d0 * norm * (x**2 - y**2) / r**2  ! d_x²-y²
        end select
      end if
      
    case default
      ylm = 0.0d0
    end select
    
  end function calc_spherical_harmonic

  !--------------------------------------------------------------------------------------!
  ! Get number of auxiliary functions
  !--------------------------------------------------------------------------------------!
  function get_n_auxiliary(aux_basis) result(n_aux)
    implicit none
    type(auxiliary_basis_t), intent(in) :: aux_basis
    integer :: n_aux
    
    n_aux = aux_basis%n_aux
    
  end function get_n_auxiliary

  !--------------------------------------------------------------------------------------!
  ! Deallocate auxiliary basis
  !--------------------------------------------------------------------------------------!
  subroutine deallocate_auxiliary_basis(aux_basis)
    implicit none
    type(auxiliary_basis_t), intent(inout) :: aux_basis
    
    if (allocated(aux_basis%atom_index)) deallocate(aux_basis%atom_index)
    if (allocated(aux_basis%l_quantum)) deallocate(aux_basis%l_quantum)
    if (allocated(aux_basis%m_quantum)) deallocate(aux_basis%m_quantum)
    if (allocated(aux_basis%alpha)) deallocate(aux_basis%alpha)
    if (allocated(aux_basis%center)) deallocate(aux_basis%center)
    
  end subroutine deallocate_auxiliary_basis

end module auxiliary_basis
