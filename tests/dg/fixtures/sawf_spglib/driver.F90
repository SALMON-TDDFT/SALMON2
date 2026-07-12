#include "config.h"
program check_sawf_spglib_backend_driver
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use lcfo_wannier_sawf, only: t_sawf_symop, load_sawf_symmetry_auto
  implicit none

#ifndef HAVE_SPGLIB
  call check_fallback
#else
  call check_diamond
  call check_translated_order
  call check_noncentrosymmetric
  call check_failure(4, 0.04d0, 'operation count')
  call check_failure(5, 0.05d0, 'C API')
  call check_failure(6, 0.00d0, 'inverse')
  call check_failure(7, 0.07d0, 'operation limit')
  call check_invalid_inputs
#endif
  write(*,'(a)') 'PASS focused SAWF Spglib backend'

contains

  subroutine check_fallback
    real(8) :: lattice(3,3), frac_pos(3,1)
    integer :: species(1)
    type(t_sawf_symop), allocatable :: operations(:)
    logical :: ok
    character(512) :: message

    call identity_lattice(lattice)
    frac_pos = 0.0d0
    species = 14
    call load_sawf_symmetry_auto(lattice, frac_pos, species, 1.0d-8, &
      operations, ok, message)
    call require(.not. ok .and. index(message, 'USE_SPGLIB') > 0, &
      'OFF fallback must request USE_SPGLIB')
  end subroutine check_fallback


  subroutine check_diamond
    real(8) :: lattice(3,3), frac_pos(3,2)
    integer :: species(2)
    type(t_sawf_symop), allocatable :: operations(:)
    logical :: ok
    character(512) :: message

    lattice(:,1) = [0.0d0, 0.5d0, 0.5d0]
    lattice(:,2) = [0.5d0, 0.0d0, 0.5d0]
    lattice(:,3) = [0.5d0, 0.5d0, 0.0d0]
    frac_pos(:,1) = 0.0d0
    frac_pos(:,2) = 0.25d0
    species = 6
    call load_sawf_symmetry_auto(lattice, frac_pos, species, 1.0d-8, &
      operations, ok, message)
    call require(ok, 'diamond mock failed: '//trim(message))
    call require(size(operations) == 2, 'diamond operation count')
    call require(all(operations(2)%W == reshape( &
      [-1,0,0,0,-1,0,0,0,-1],[3,3])), 'diamond inversion')
    call require(all(operations(2)%atom_map == [2,1]), 'diamond atom map')
  end subroutine check_diamond


  subroutine check_translated_order
    real(8) :: lattice(3,3), frac_pos(3,2)
    integer :: species(2)
    type(t_sawf_symop), allocatable :: operations(:)
    logical :: ok
    character(512) :: message

    lattice(:,1) = [2.0d0, 0.0d0, 0.0d0]
    lattice(:,2) = [0.4d0, 1.5d0, 0.0d0]
    lattice(:,3) = [0.2d0, 0.3d0, 1.2d0]
    frac_pos = 0.0d0
    frac_pos(1,:) = [0.5d0, 0.0d0]
    species = 14
    call load_sawf_symmetry_auto(lattice, frac_pos, species, 1.0d-8, &
      operations, ok, message)
    call require(ok, 'translated mock failed: '//trim(message))
    call require(size(operations) == 2, 'translated operation count')
    call require(all(operations(2)%atom_map == [2,1]), &
      'translated atom ordering')
  end subroutine check_translated_order


  subroutine check_noncentrosymmetric
    real(8) :: lattice(3,3), frac_pos(3,3)
    integer :: species(3)
    type(t_sawf_symop), allocatable :: operations(:)
    logical :: ok
    character(512) :: message

    lattice(:,1) = [1.7d0, 0.0d0, 0.0d0]
    lattice(:,2) = [0.2d0, 1.4d0, 0.0d0]
    lattice(:,3) = [0.1d0, 0.3d0, 1.1d0]
    frac_pos(:,1) = [0.11d0, 0.23d0, 0.37d0]
    frac_pos(:,2) = [0.41d0, 0.07d0, 0.62d0]
    frac_pos(:,3) = [0.79d0, 0.53d0, 0.19d0]
    species = [14, 8, 1]
    call load_sawf_symmetry_auto(lattice, frac_pos, species, 1.0d-8, &
      operations, ok, message)
    call require(ok, 'noncentrosymmetric mock failed: '//trim(message))
    call require(size(operations) == 1, &
      'noncentrosymmetric identity-only operation count')
  end subroutine check_noncentrosymmetric


  subroutine check_failure(species_id, marker, expected)
    integer, intent(in) :: species_id
    real(8), intent(in) :: marker
    character(*), intent(in) :: expected
    real(8) :: lattice(3,3), frac_pos(3,1)
    integer :: species(1)
    type(t_sawf_symop), allocatable :: operations(:)
    logical :: ok
    character(512) :: message

    call identity_lattice(lattice)
    frac_pos = 0.0d0
    frac_pos(1,1) = marker
    species = species_id
    call load_sawf_symmetry_auto(lattice, frac_pos, species, 1.0d-8, &
      operations, ok, message)
    call require(.not. ok .and. index(message, expected) > 0, &
      'missing failure diagnostic: '//trim(expected)//' actual='//trim(message))
  end subroutine check_failure


  subroutine check_invalid_inputs
    real(8), allocatable :: frac_pos(:,:)
    real(8) :: lattice(3,3), tolerance
    integer, allocatable :: species(:)
    type(t_sawf_symop), allocatable :: operations(:)
    logical :: ok
    character(512) :: message

    call identity_lattice(lattice)
    tolerance = 1.0d-8
    allocate(frac_pos(2,1), species(1))
    frac_pos = 0.0d0
    species = 6
    call load_sawf_symmetry_auto(lattice, frac_pos, species, tolerance, &
      operations, ok, message)
    call require(.not. ok .and. index(message, 'dimensions') > 0, &
      'dimension validation')
    deallocate(frac_pos, species)

    allocate(frac_pos(3,1), species(1))
    frac_pos = 0.0d0
    species = 0
    call load_sawf_symmetry_auto(lattice, frac_pos, species, tolerance, &
      operations, ok, message)
    call require(.not. ok .and. index(message, 'species') > 0, &
      'species validation')

    species = 6
    frac_pos(1,1) = ieee_value(0.0d0, ieee_quiet_nan)
    call load_sawf_symmetry_auto(lattice, frac_pos, species, tolerance, &
      operations, ok, message)
    call require(.not. ok .and. index(message, 'finite') > 0, &
      'finite position validation')
  end subroutine check_invalid_inputs


  subroutine identity_lattice(lattice)
    real(8), intent(out) :: lattice(3,3)
    lattice = 0.0d0
    lattice(1,1) = 1.0d0
    lattice(2,2) = 1.0d0
    lattice(3,3) = 1.0d0
  end subroutine identity_lattice


  subroutine require(condition, description)
    logical, intent(in) :: condition
    character(*), intent(in) :: description
    if (.not. condition) error stop description
  end subroutine require
end program check_sawf_spglib_backend_driver
