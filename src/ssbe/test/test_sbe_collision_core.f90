! src/ssbe/test/test_sbe_collision_core.f90
program test_sbe_collision_core
  use sbe_collision_gw_core
  implicit none
  integer, parameter :: nb = 2, nk = 1
  complex(8) :: dqnm(nb,nb,1:1), qnm(nb,nb,1:1), rho(nb,nb,1:1)
  real(8)    :: gamma(nb,nk), f0ref(nb,nk)
  real(8)    :: dt, tol
  integer    :: nfail
  complex(8) :: expect

  nfail = 0
  tol   = 1.0d-13

  ! --- diagonal loss: dqnm(1,1) += -gamma*(rho-f0) ---
  gamma  = 0.0d0; gamma(1,1) = 0.5d0
  f0ref  = 0.0d0
  qnm    = 0.0d0; qnm(1,1,1) = dcmplx(0.3d0, 0.0d0)
  dqnm   = dcmplx(7.0d0, 0.0d0)            ! pre-existing RHS must be preserved
  call add_collision_diag(dqnm, qnm, gamma, f0ref, nb, nk, 1, 1)
  expect = dcmplx(7.0d0, 0.0d0) - 0.5d0*(dcmplx(0.3d0,0.0d0) - 0.0d0)
  call check(dqnm(1,1,1), expect, tol, "diag", nfail)

  ! --- off-diagonal dephasing, mode 'gw': dqnm(1,2) += -(g1+g2)/2 * rho(1,2) ---
  gamma(1,1) = 0.5d0; gamma(2,1) = 1.5d0
  qnm  = 0.0d0; qnm(1,2,1) = dcmplx(0.0d0, 0.4d0)
  dqnm = 0.0d0
  call add_collision_offdiag(dqnm, qnm, gamma, nb, nk, 1, 1, "gw")
  expect = -0.5d0*(0.5d0+1.5d0)*dcmplx(0.0d0,0.4d0)
  call check(dqnm(1,2,1), expect, tol, "offdiag-gw", nfail)

  ! --- off-diagonal mode 't2' is a no-op ---
  dqnm = dcmplx(2.0d0, 0.0d0)
  call add_collision_offdiag(dqnm, qnm, gamma, nb, nk, 1, 1, "t2")
  call check(dqnm(1,2,1), dcmplx(2.0d0,0.0d0), tol, "offdiag-t2-noop", nfail)

  ! --- VG forward-Euler: rho(1,1) += -gamma*(rho-f0)*dt ---
  dt = 0.1d0
  gamma = 0.0d0; gamma(1,1) = 0.5d0
  f0ref = 0.0d0
  rho   = 0.0d0; rho(1,1,1) = dcmplx(0.3d0, 0.0d0)
  call add_collision_vg(rho, gamma, f0ref, dt, nb, nk, 1, 1, "gw")
  expect = dcmplx(0.3d0,0.0d0) - 0.5d0*(dcmplx(0.3d0,0.0d0)-0.0d0)*dt
  call check(rho(1,1,1), expect, tol, "vg-diag", nfail)

  if (nfail == 0) then
    write(*,*) "ALL CORE TESTS PASSED"
  else
    write(*,*) "FAILED:", nfail
    stop 1
  end if
contains
  subroutine check(got, want, tol, label, nfail)
    complex(8), intent(in) :: got, want
    real(8), intent(in)    :: tol
    character(*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (abs(got - want) > tol) then
      write(*,'(a,a,2(2es14.5))') "MISMATCH ", label, got, want
      nfail = nfail + 1
    else
      write(*,'(a,a)') "ok ", label
    end if
  end subroutine check
end program test_sbe_collision_core
