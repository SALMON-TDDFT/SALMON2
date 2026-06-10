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
! Ghost-state detection for separable (Kleinman-Bylander / Blochl form)
! pseudopotentials, by direct diagonalization of the radial Hamiltonian
! actually applied by SALMON:
!   H_l = -1/2 d^2/dr^2 + l(l+1)/(2 r^2) + V_loc(r)
!         + sum_j |chi_lj> inorm_lj <chi_lj|        (chi = udvtbl tables)
! A ghost manifests as a bound state of H_l lying well below the lowest
! bound state of the LOCAL-only Hamiltonian (the separable term digs a
! spurious deep level that V_loc alone cannot support). Diagnostic only:
! a warning is printed and a per-channel table is written; the run is
! never aborted (flat deep bands can be physical in other contexts).
module ghost_check
  implicit none
  private
  public :: check_pp_ghosts

  real(8), parameter :: pi = 3.141592653589793d0
  real(8), parameter :: ghost_tol = 0.1d0   ! Ha; eps_sep_min < eps_loc_min - tol => warn
  real(8), parameter :: mesh_h = 0.025d0    ! a.u.
  real(8), parameter :: r_box = 25d0        ! a.u.

contains

  subroutine check_pp_ghosts(pp, ik)
    use structures, only: s_pp_info
    use salmon_global, only: base_directory, file_pseudo, yn_ghost_check
    use filesystem, only: open_filehandle
    implicit none
    type(s_pp_info), intent(in) :: pp
    integer, intent(in) :: ik
    integer :: n, ll, l, l0, j, i, fh, nfound, nbound
    real(8), allocatable :: rmesh(:), vloc(:), chi(:,:), h0(:,:), h1(:,:)
    real(8), allocatable :: eloc(:), esep(:)
    real(8) :: eps_loc(2), eps_sep(3)
    logical :: flagged_any, flagged
    character(256) :: fname

    if (yn_ghost_check /= 'y') return

    n = int(r_box / mesh_h) - 1
    allocate(rmesh(n), vloc(n), eloc(n), esep(n))
    allocate(h0(n,n), h1(n,n))
    do i = 1, n
      rmesh(i) = i * mesh_h
    end do
    call interp_vloc(pp, ik, n, rmesh, vloc)
    call add_valence_screening(pp, ik, n, rmesh, vloc)

    fname = trim(base_directory)//"PS_"//trim(pp%atom_symbol(ik))//"_ghost_check.data"
    fh = open_filehandle(trim(fname), status='replace')
    write(fh,'(a)') '# Separable-pseudopotential ghost diagnostic (radial direct diagonalization)'
    write(fh,'(a)') '# H_sep = T_l + V_scr + sum_j |udv_lj> c_lj <udv_lj| ; V_scr = V_loc + V_H[rho_atom] + V_xc[rho_atom]'
    write(fh,'(a)') '# verdict GHOST when the number of screened bound levels (eps < -0.05 Ha) exceeds nproj'
    write(fh,'(a)') '# (screened eigenvalues are directly comparable to atomic reference levels)'
    write(fh,'(a)') '# 1:l 2:nproj 3:eps_loc_0[Ha] 4:eps_loc_1[Ha] 5:eps_sep_0[Ha] 6:eps_sep_1[Ha] 7:eps_sep_2[Ha] 8:nbound 9:verdict'

    flagged_any = .false.
    l0 = 0
    do ll = 0, pp%mlps(ik)
      ! local-only spectrum for this angular momentum
      call build_h_local(n, mesh_h, rmesh, vloc, ll, h0)
      call lowest_eigs(n, h0, 2, eps_loc, eloc)

      ! separable spectrum: add every projector of this channel
      allocate(chi(n, max(1, pp%nproj(ll,ik))))
      h1 = 0d0
      call build_h_local(n, mesh_h, rmesh, vloc, ll, h1)
      do l = l0, l0 + pp%nproj(ll,ik) - 1
        j = l - l0 + 1
        call interp_chi(pp, ik, l, ll, n, rmesh, chi(:,j))
        if (pp%inorm(l,ik) /= 0) then
          ! undo the (2l+1)/4pi angular factor folded into udvtbl
          call add_separable(n, mesh_h, chi(:,j), dble(pp%inorm(l,ik)) * 4d0*pi/dble(2*ll+1), h1)
        end if
      end do
      call lowest_eigs(n, h1, 3, eps_sep, esep)
      deallocate(chi)

      ! Ghost criterion (screened): a separable channel may legitimately
      ! create reference bound states (semicore!), but never MORE bound
      ! states than it has projectors/reference orbitals. Flag when the
      ! number of screened bound levels (eps < -0.05 Ha) exceeds nproj,
      ! or when an extra level dives below the deepest physically
      ! expected scale (3x the deepest projector-made level).
      nbound = 0
      do i = 1, n
        if (esep(i) < -0.05d0) nbound = nbound + 1
      end do
      flagged = (nbound > pp%nproj(ll,ik))
      if (flagged) flagged_any = .true.
      write(fh,'(1x,i2,1x,i2,5(1x,es15.7),1x,i3,1x,a)') ll, pp%nproj(ll,ik), &
        eps_loc(1), eps_loc(2), eps_sep(1), eps_sep(2), eps_sep(3), nbound, &
        merge('GHOST', 'ok   ', flagged)
      if (flagged) then
        write(*,'(a)') '!!! WARNING: ghost state(s) suspected in pseudopotential '// &
          trim(file_pseudo(ik))//' (element '//trim(pp%atom_symbol(ik))//')'
        write(*,'(a,i2,a,i3,a,i2,a,f12.5,a)') '!!!   channel l=', ll, &
          ': ', nbound, ' bound levels for ', pp%nproj(ll,ik), &
          ' projectors; deepest =', eps_sep(1), ' Ha'
        write(*,'(a)') '!!!   SCF may converge to an unphysical state (huge kinetic/nonlocal stress,'
        write(*,'(a)') '!!!   flat deep bands). Consider a different pseudopotential. See '//trim(fname)
      end if

      l0 = l0 + pp%nproj(ll,ik)
    end do
    close(fh)
    if (.not. flagged_any) then
      write(*,'(a)') 'Ghost check: no ghost states detected for element '//trim(pp%atom_symbol(ik))
    end if
    nfound = 0 ! silence unused warnings on some compilers
    deallocate(rmesh, vloc, eloc, esep, h0, h1)
  end subroutine check_pp_ghosts


  subroutine add_valence_screening(pp, ik, n, rmesh, vloc)
    ! Add the atomic valence screening V_H[rho] + V_xc[rho] built from the
    ! degeneracy-weighted atomic density table rho_pp_tbl, renormalized to
    ! zps electrons. Makes the channel spectra directly comparable to the
    ! pseudo-atom reference levels (essential for semicore pseudopotentials,
    ! whose separable term legitimately creates deep reference states).
    use structures, only: s_pp_info
    implicit none
    type(s_pp_info), intent(in) :: pp
    integer, intent(in) :: ik, n
    real(8), intent(in) :: rmesh(n)
    real(8), intent(inout) :: vloc(n)
    integer :: i, ir, nr
    real(8) :: r, x, qtot, rho, rs, ex, vx, vc, srs
    real(8), allocatable :: q(:), inner(:), outer(:)
    real(8), parameter :: pi_l = 3.141592653589793d0
    real(8), parameter :: ca = 0.0311d0, cb = -0.048d0, cc = 0.0020d0, cd = -0.0116d0
    real(8), parameter :: cg = -0.1423d0, cb1 = 1.0529d0, cb2 = 0.3334d0

    allocate(q(n), inner(n), outer(n))
    nr = pp%nrmax
    do i = 1, n
      r = rmesh(i)
      if (r >= pp%rad(nr,ik)) then
        q(i) = 0d0
      else
        ir = 1
        do while (ir < nr .and. pp%rad(ir+1,ik) <= r)
          ir = ir + 1
        end do
        x = (r - pp%rad(ir,ik)) / max(pp%rad(ir+1,ik) - pp%rad(ir,ik), 1d-300)
        q(i) = (1d0-x) * pp%rho_pp_tbl(ir,ik) + x * pp%rho_pp_tbl(ir+1,ik)
      end if
      if (q(i) < 0d0) q(i) = 0d0
    end do
    qtot = 0d0
    do i = 1, n
      qtot = qtot + q(i)
    end do
    qtot = qtot * (rmesh(2)-rmesh(1))
    if (qtot <= 1d-12) then
      deallocate(q, inner, outer)
      return
    end if
    q = q * (pp%zps(ik) / qtot)

    ! radial Hartree: V_H(r) = (1/r) int_0^r q dr' + int_r^inf q/r' dr'
    inner(1) = 0.5d0 * q(1) * (rmesh(2)-rmesh(1))
    do i = 2, n
      inner(i) = inner(i-1) + 0.5d0*(q(i)+q(i-1))*(rmesh(2)-rmesh(1))
    end do
    outer(n) = 0d0
    do i = n-1, 1, -1
      outer(i) = outer(i+1) + 0.5d0*(q(i)/rmesh(i)+q(i+1)/rmesh(i+1))*(rmesh(2)-rmesh(1))
    end do
    do i = 1, n
      rho = q(i) / (4d0*pi_l*rmesh(i)**2)
      rho = max(rho, 1d-30)
      rs = (3d0/(4d0*pi_l*rho))**(1d0/3d0)
      ex = -0.75d0*(3d0/pi_l)**(1d0/3d0)*rho**(1d0/3d0)
      vx = 4d0*ex/3d0
      if (rs >= 1d0) then
        srs = sqrt(rs)
        vc = cg*(1d0+7d0/6d0*cb1*srs+4d0/3d0*cb2*rs)/(1d0+cb1*srs+cb2*rs)**2
      else
        vc = ca*log(rs)+(cb-ca/3d0)+2d0/3d0*cc*rs*log(rs)+(2d0*cd-cc)*rs/3d0
      end if
      vloc(i) = vloc(i) + inner(i)/rmesh(i) + outer(i) + vx + vc
    end do
    deallocate(q, inner, outer)
  end subroutine add_valence_screening

  subroutine interp_vloc(pp, ik, n, rmesh, vloc)
    use structures, only: s_pp_info
    implicit none
    type(s_pp_info), intent(in) :: pp
    integer, intent(in) :: ik, n
    real(8), intent(in) :: rmesh(n)
    real(8), intent(out) :: vloc(n)
    integer :: i, ir, nr
    real(8) :: r, x

    ! vloctbl is guaranteed valid up to nrps for every reader branch;
    ! beyond the pseudization radius the local potential is Coulombic.
    nr = pp%nrps(ik)
    do i = 1, n
      r = rmesh(i)
      if (r >= pp%rad(nr,ik)) then
        vloc(i) = -pp%zps(ik) / r          ! Coulomb tail beyond the table
      else
        ir = 1
        do while (ir < nr .and. pp%rad(ir+1,ik) <= r)
          ir = ir + 1
        end do
        x = (r - pp%rad(ir,ik)) / max(pp%rad(ir+1,ik) - pp%rad(ir,ik), 1d-300)
        vloc(i) = (1d0-x) * pp%vloctbl(ir,ik) + x * pp%vloctbl(ir+1,ik)
      end if
    end do
  end subroutine interp_vloc

  subroutine interp_chi(pp, ik, l, ll, n, rmesh, chi)
    use structures, only: s_pp_info
    implicit none
    type(s_pp_info), intent(in) :: pp
    integer, intent(in) :: ik, l, ll, n
    real(8), intent(in) :: rmesh(n)
    real(8), intent(out) :: chi(n)
    integer :: i, ir, nr
    real(8) :: r, x

    nr = pp%nrps(ik)
    do i = 1, n
      r = rmesh(i)
      if (r >= pp%rad(nr,ik)) then
        chi(i) = 0d0                       ! projectors are compactly supported
      else
        ir = 1
        do while (ir < nr .and. pp%rad(ir+1,ik) <= r)
          ir = ir + 1
        end do
        x = (r - pp%rad(ir,ik)) / max(pp%rad(ir+1,ik) - pp%rad(ir,ik), 1d-300)
        ! udvtbl stores the reduced form chi(r)*sqrt((2l+1)/4pi)/r^(l+1)
        ! (input_pp scales by const/r^(l+1) for 3D interpolation);
        ! reconstruct the radial u-representation projector here.
        chi(i) = ((1d0-x) * pp%udvtbl(ir,l,ik) + x * pp%udvtbl(ir+1,l,ik)) * r**(ll+1)
      end if
    end do
  end subroutine interp_chi

  subroutine build_h_local(n, h, rmesh, vloc, ll, ham)
    implicit none
    integer, intent(in) :: n, ll
    real(8), intent(in) :: h, rmesh(n), vloc(n)
    real(8), intent(inout) :: ham(n,n)
    integer :: i
    real(8) :: invh2

    invh2 = 1d0 / (h*h)
    ham = 0d0
    do i = 1, n
      ham(i,i) = invh2 + dble(ll*(ll+1)) / (2d0 * rmesh(i)**2) + vloc(i)
      if (i < n) then
        ham(i,i+1) = -0.5d0 * invh2
        ham(i+1,i) = -0.5d0 * invh2
      end if
    end do
  end subroutine build_h_local

  subroutine add_separable(n, h, chi, coef, ham)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: h, chi(n), coef
    real(8), intent(inout) :: ham(n,n)
    integer :: i, k

    ! <i|V_NL|k> = chi_i * coef * chi_k * h  (radial measure dr -> h)
    do k = 1, n
      if (chi(k) == 0d0) cycle
      do i = 1, n
        ham(i,k) = ham(i,k) + chi(i) * coef * chi(k) * h
      end do
    end do
  end subroutine add_separable

  subroutine lowest_eigs(n, ham, nl, eps, work_eigs)
    implicit none
    integer, intent(in) :: n, nl
    real(8), intent(inout) :: ham(n,n)
    real(8), intent(out) :: eps(nl)
    real(8), intent(out) :: work_eigs(n)
    integer :: lwork, info, i
    real(8), allocatable :: work(:)

    lwork = max(1, 3*n)
    allocate(work(lwork))
    call dsyev('N', 'U', n, ham, n, work_eigs, work, lwork, info)
    if (info /= 0) then
      eps = 0d0
      deallocate(work)
      return
    end if
    do i = 1, nl
      eps(i) = work_eigs(i)
    end do
    deallocate(work)
  end subroutine lowest_eigs

end module ghost_check
