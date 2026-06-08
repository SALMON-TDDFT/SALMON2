!
!  Empirical Pseudopotential Method (EPM): Cohen-Bergstresser local form factors
!
!  Reference: M. L. Cohen and T. K. Bergstresser, Phys. Rev. 141, 789 (1966).
!  Local pseudopotential form factors for GaAs (zincblende), tabulated at the
!  squared reciprocal-lattice-vector magnitudes G^2 = 3, 4, 8, 11 in units of
!  (2*pi/a)^2. Values are quoted by the original authors in Rydberg; this
!  module returns them already converted to Hartree atomic units (the
!  convention used throughout the rest of SALMON), i.e. divided by 2.
!
module epm_cohen_bergstresser
    implicit none

    private
    public :: cb_get_form_factors, cb_tau_zincblende, cb_lattice_vectors_fcc

contains

    ! Symmetric/antisymmetric local pseudopotential form factors V^S(G^2), V^A(G^2)
    ! in Hartree atomic units. G2 is the squared length of the reciprocal lattice
    ! vector in units of (2*pi/a)^2 (an integer for the fcc/zincblende lattice).
    subroutine cb_get_form_factors(material, G2, VS_ha, VA_ha)
        implicit none
        character(*), intent(in)  :: material
        integer,      intent(in)  :: G2
        real(8),      intent(out) :: VS_ha, VA_ha
        real(8), parameter :: ry_to_ha = 0.5d0
        real(8) :: VS_ry, VA_ry

        VS_ry = 0d0
        VA_ry = 0d0

        select case (trim(material))
        case ('GaAs')
            ! Cohen-Bergstresser (1966), Table 2:
            !   V^S(3)  = -0.23 Ry,  V^A(3)  = +0.07 Ry
            !   V^S(4)  =  0.00 Ry,  V^A(4)  = +0.05 Ry   (V^S(4)=0: structure factor vanishes)
            !   V^S(8)  = +0.01 Ry,  V^A(8)  =  0.00 Ry
            !   V^S(11) = +0.06 Ry,  V^A(11) = +0.01 Ry
            select case (G2)
            case (3)
                VS_ry = -0.23d0;  VA_ry =  0.07d0
            case (4)
                VS_ry =  0.00d0;  VA_ry =  0.05d0
            case (8)
                VS_ry =  0.01d0;  VA_ry =  0.00d0
            case (11)
                VS_ry =  0.06d0;  VA_ry =  0.01d0
            end select
        case default
            stop 'epm_cohen_bergstresser: unsupported epm_material (only "GaAs" form factors are tabulated)'
        end select

        VS_ha = VS_ry * ry_to_ha
        VA_ha = VA_ry * ry_to_ha
    end subroutine cb_get_form_factors


    ! Zincblende internal displacement between the two-atom basis, with the
    ! origin placed midway between the cation and the anion (Cohen-Bergstresser
    ! convention): tau = (a/8)*(1,1,1). This fixes the sign of V^A consistently
    ! with the tabulated form factors above.
    function cb_tau_zincblende(a_lattice) result(tau)
        implicit none
        real(8), intent(in) :: a_lattice
        real(8) :: tau(3)
        tau(1:3) = a_lattice / 8.0d0
    end function cb_tau_zincblende


    ! Conventional fcc primitive lattice vectors for the zincblende structure
    ! (Cohen-Bergstresser convention):
    !   a1 = (a/2)(0,1,1),  a2 = (a/2)(1,0,1),  a3 = (a/2)(1,1,0)
    subroutine cb_lattice_vectors_fcc(a_lattice, a1, a2, a3)
        implicit none
        real(8), intent(in)  :: a_lattice
        real(8), intent(out) :: a1(3), a2(3), a3(3)
        real(8) :: h
        h = 0.5d0 * a_lattice
        a1 = (/ 0d0, h,   h   /)
        a2 = (/ h,   0d0, h   /)
        a3 = (/ h,   h,   0d0 /)
    end subroutine cb_lattice_vectors_fcc

end module epm_cohen_bergstresser
