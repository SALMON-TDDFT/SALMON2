! src/ssbe/sbe_collision_gw_core.f90
! GW-derived collision-term arithmetic for the semiconductor Bloch solver.
! Pure (no module dependencies) so it is unit-testable on its own.
! Linearized-RTA quasiparticle-lifetime sink: a diagonal occupation loss toward
! the cold reference f0, and an off-diagonal non-Hermitian linewidth dephasing.
module sbe_collision_gw_core
  implicit none
  private
  public :: add_collision_diag, add_collision_offdiag, add_collision_vg

contains

  ! Diagonal occupation loss into a right-hand-side accumulator (length-gauge).
  subroutine add_collision_diag(dqnm, qnm, gamma, f0ref, nb, nk, ik_min, ik_max)
    implicit none
    integer, intent(in)       :: nb, nk, ik_min, ik_max
    complex(8), intent(inout) :: dqnm(nb, nb, ik_min:ik_max)
    complex(8), intent(in)    :: qnm (nb, nb, ik_min:ik_max)
    real(8), intent(in)       :: gamma(nb, nk)
    real(8), intent(in)       :: f0ref(nb, nk)
    integer :: ik, ib
    !$omp parallel do default(shared) private(ik, ib) collapse(2)
    do ik = ik_min, ik_max
      do ib = 1, nb
        dqnm(ib, ib, ik) = dqnm(ib, ib, ik) &
          & - gamma(ib, ik) * ( qnm(ib, ib, ik) - f0ref(ib, ik) )
      end do
    end do
  end subroutine add_collision_diag

  ! Off-diagonal dephasing into the RHS accumulator (length-gauge).
  ! deph_mode 'gw'/'both' add the GW linewidth; 't2' adds nothing.
  subroutine add_collision_offdiag(dqnm, qnm, gamma, nb, nk, ik_min, ik_max, deph_mode)
    implicit none
    integer, intent(in)       :: nb, nk, ik_min, ik_max
    complex(8), intent(inout) :: dqnm(nb, nb, ik_min:ik_max)
    complex(8), intent(in)    :: qnm (nb, nb, ik_min:ik_max)
    real(8), intent(in)       :: gamma(nb, nk)
    character(*), intent(in)  :: deph_mode
    integer :: ik, ib, jb
    real(8) :: g
    select case (trim(deph_mode))
    case ('gw', 'both')
      !$omp parallel do default(shared) private(ik, ib, jb, g) collapse(2)
      do ik = ik_min, ik_max
        do ib = 1, nb
          do jb = 1, nb
            if (ib == jb) cycle
            g = 0.5d0 * ( gamma(ib, ik) + gamma(jb, ik) )
            dqnm(ib, jb, ik) = dqnm(ib, jb, ik) - g * qnm(ib, jb, ik)
          end do
        end do
      end do
    case default
      ! 't2': keep only the solver's existing t_2 dephasing.
    end select
  end subroutine add_collision_offdiag

  ! Velocity-gauge forward-Euler collision step applied directly to rho.
  ! VG carries no t_2 baseline, so 'gw' and 'both' behave identically here.
  subroutine add_collision_vg(rho, gamma, f0ref, dt, nb, nk, ik_min, ik_max, deph_mode)
    implicit none
    integer, intent(in)       :: nb, nk, ik_min, ik_max
    complex(8), intent(inout) :: rho(nb, nb, ik_min:ik_max)
    real(8), intent(in)       :: gamma(nb, nk)
    real(8), intent(in)       :: f0ref(nb, nk)
    real(8), intent(in)       :: dt
    character(*), intent(in)  :: deph_mode
    integer :: ik, ib, jb
    real(8) :: g
    logical :: do_off
    do_off = (trim(deph_mode) == 'gw' .or. trim(deph_mode) == 'both')
    !$omp parallel do default(shared) private(ik, ib, jb, g) collapse(1)
    do ik = ik_min, ik_max
      do ib = 1, nb
        rho(ib, ib, ik) = rho(ib, ib, ik) &
          & - gamma(ib, ik) * ( rho(ib, ib, ik) - f0ref(ib, ik) ) * dt
      end do
      if (do_off) then
        do ib = 1, nb
          do jb = 1, nb
            if (ib == jb) cycle
            g = 0.5d0 * ( gamma(ib, ik) + gamma(jb, ik) )
            rho(ib, jb, ik) = rho(ib, jb, ik) - g * rho(ib, jb, ik) * dt
          end do
        end do
      end if
    end do
  end subroutine add_collision_vg

end module sbe_collision_gw_core
