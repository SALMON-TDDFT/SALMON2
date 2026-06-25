program test_load_gw_rate
  use sbe_collision_gw, only: load_gw_rate
  implicit none
  integer, parameter :: nb = 2, nk = 1
  real(8) :: gamma_nk(nb, nk)
  real(8) :: au_time_fs, tol, g_au_expect
  integer :: nfail
  au_time_fs = 0.02418884326505d0
  tol   = 1.0d-12
  nfail = 0

  ! fixture A: unit_system = au, Gamma already in 1/au_time
  open(31, file='/tmp/rate_au.data', status='replace')
  write(31,'(a)') '# on-shell Im Sigma_c per (n,k)'
  write(31,'(a)') '# unit_system = au'
  write(31,'(a)') '# 1:n 2:k 3:Eks 4:Eqp 5:ReSigc 6:ImSigc 7:Gamma 8:occ'
  write(31,'(2i5,4es15.6,es13.4,f6.1)') 1,1, 0d0,0d0,0d0,0d0, 0.2000d0, 2.0d0
  write(31,'(2i5,4es15.6,es13.4,f6.1)') 2,1, 0d0,0d0,0d0,0d0, 0.5000d0, 0.0d0
  close(31)
  gamma_nk = -1d0
  call load_gw_rate('/tmp/rate_au.data', nb, nk, gamma_nk, -1)
  call check(gamma_nk(1,1), 0.2d0, tol, "au-band1", nfail)
  call check(gamma_nk(2,1), 0.5d0, tol, "au-band2", nfail)

  ! fixture B: unit_system = A_eV_fs, Gamma in 1/fs -> *au_time_fs to a.u.
  open(32, file='/tmp/rate_ev.data', status='replace')
  write(32,'(a)') '# on-shell Im Sigma_c per (n,k)'
  write(32,'(a)') '# unit_system = A_eV_fs'
  write(32,'(a)') '# 1:n 2:k 3:Eks 4:Eqp 5:ReSigc 6:ImSigc 7:Gamma 8:occ'
  write(32,'(2i5,4es15.6,es13.4,f6.1)') 1,1, 0d0,0d0,0d0,0d0, 10.0000d0, 2.0d0
  close(32)
  gamma_nk = -1d0
  call load_gw_rate('/tmp/rate_ev.data', nb, nk, gamma_nk, -1)
  g_au_expect = 10.0d0 * au_time_fs
  call check(gamma_nk(1,1), g_au_expect, tol, "evfs-band1", nfail)
  call check(gamma_nk(2,1), 0.0d0,       tol, "evfs-missing-band2-zero", nfail)

  ! fixture C: unit_system = a.u. (alternate canonical tag)
  open(33, file='/tmp/rate_au_dot.data', status='replace')
  write(33,'(a)') '# on-shell Im Sigma_c per (n,k)'
  write(33,'(a)') '# unit_system = a.u.'
  write(33,'(a)') '# 1:n 2:k 3:Eks 4:Eqp 5:ReSigc 6:ImSigc 7:Gamma 8:occ'
  write(33,'(2i5,4es15.6,es13.4,f6.1)') 1,1, 0d0,0d0,0d0,0d0, 0.3000d0, 2.0d0
  write(33,'(2i5,4es15.6,es13.4,f6.1)') 2,1, 0d0,0d0,0d0,0d0, 0.7000d0, 0.0d0
  close(33)
  gamma_nk = -1d0
  call load_gw_rate('/tmp/rate_au_dot.data', nb, nk, gamma_nk, -1)
  call check(gamma_nk(1,1), 0.3d0, tol, "a.u.-band1", nfail)
  call check(gamma_nk(2,1), 0.7d0, tol, "a.u.-band2", nfail)

  ! fixture D: no unit_system line (tagless) -> a.u. fallback, no warning
  open(34, file='/tmp/rate_tagless.data', status='replace')
  write(34,'(a)') '# on-shell Im Sigma_c per (n,k)'
  write(34,'(a)') '# 1:n 2:k 3:Eks 4:Eqp 5:ReSigc 6:ImSigc 7:Gamma 8:occ'
  write(34,'(2i5,4es15.6,es13.4,f6.1)') 1,1, 0d0,0d0,0d0,0d0, 0.4000d0, 2.0d0
  write(34,'(2i5,4es15.6,es13.4,f6.1)') 2,1, 0d0,0d0,0d0,0d0, 0.6000d0, 0.0d0
  close(34)
  gamma_nk = -1d0
  call load_gw_rate('/tmp/rate_tagless.data', nb, nk, gamma_nk, -1)
  call check(gamma_nk(1,1), 0.4d0, tol, "tagless-band1", nfail)
  call check(gamma_nk(2,1), 0.6d0, tol, "tagless-band2", nfail)

  ! fixture E: unknown tag -> a.u. fallback + warning printed (not checked as error)
  open(35, file='/tmp/rate_unknown.data', status='replace')
  write(35,'(a)') '# on-shell Im Sigma_c per (n,k)'
  write(35,'(a)') '# unit_system = bohr_Ry'
  write(35,'(a)') '# 1:n 2:k 3:Eks 4:Eqp 5:ReSigc 6:ImSigc 7:Gamma 8:occ'
  write(35,'(2i5,4es15.6,es13.4,f6.1)') 1,1, 0d0,0d0,0d0,0d0, 0.1000d0, 2.0d0
  write(35,'(2i5,4es15.6,es13.4,f6.1)') 2,1, 0d0,0d0,0d0,0d0, 0.9000d0, 0.0d0
  close(35)
  gamma_nk = -1d0
  call load_gw_rate('/tmp/rate_unknown.data', nb, nk, gamma_nk, -1)
  call check(gamma_nk(1,1), 0.1d0, tol, "unknown-tag-fallback-band1", nfail)
  call check(gamma_nk(2,1), 0.9d0, tol, "unknown-tag-fallback-band2", nfail)

  if (nfail == 0) then
    write(*,*) "ALL LOADER TESTS PASSED"
  else
    write(*,*) "FAILED:", nfail; stop 1
  end if
contains
  subroutine check(got, want, tol, label, nfail)
    real(8), intent(in) :: got, want, tol
    character(*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (abs(got - want) > tol) then
      write(*,'(a,a,2es16.7)') "MISMATCH ", label, got, want; nfail = nfail + 1
    else
      write(*,'(a,a)') "ok ", label
    end if
  end subroutine check
end program test_load_gw_rate
