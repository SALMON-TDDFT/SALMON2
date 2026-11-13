! Shared state layout; no global calculation state.
module sbe_types
 implicit none
    type s_sbe_gs_info
        !Lattice information
        real(8) :: a_matrix(1:3, 1:3)
        real(8) :: b_matrix(1:3, 1:3)
        real(8) :: volume

        !Ground state (GS) electronic system information
        integer :: nk, nb, ne
        real(8), allocatable :: kpoint(:, :), kweight(:)
        real(8), allocatable :: eigen(:, :)
        real(8), allocatable :: occup(:, :)
        real(8), allocatable :: delta_omega(:, :, :)
        complex(8), allocatable :: p_mod_matrix(:, :, :, :)
        ! p_tm_matrix = <u|p|u>
        complex(8), allocatable :: p_tm_matrix(:, :, :, :)
        ! rvnl_tm_matrix = <u|-i[r, Vnl]|u>
        complex(8), allocatable :: rvnl_tm_matrix(:, :, :, :)
        complex(8), allocatable :: d_matrix(:, :, :, :)
        real(8), allocatable :: grad_k_eigen(:, :, :)

        !k-space grid and geometry information
        !NOTE: prepred for uniformally distributed k-grid....
        !integer :: num_kgrid(1:3)
    end type
    type s_sbe_bloch_solver
        !k-points for real-time SBE calculation
        integer :: nk, nb
        integer :: ik_max, ik_min
        complex(8), allocatable :: rho(:, :, :)
        logical :: flag_vnl_correction
        complex(8), allocatable :: qnm(:, :, :)
        complex(8), allocatable :: grad_qnm(:, :, :, :)
        complex(8), allocatable :: qnm_new(:, :, :)
        complex(8), allocatable :: dqnm_stock(:, :, :, :)
        complex(8), allocatable :: dnm_i(:, :, :, :)
        real(8), allocatable    :: abs_dnm(:, :, :)
        complex(8), allocatable :: exp_iphi(:, :, :)
    end type
end module
