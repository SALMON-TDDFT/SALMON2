module nstate_active_util
  implicit none
contains
  subroutine calc_nstate_active(rocc, nstate, nk, nspin, wtk, &
                                nstate_active_in, occ_threshold, nact, dropped)
    integer, intent(in)  :: nstate, nk, nspin, nstate_active_in
    real(8), intent(in)  :: rocc(nstate,nk,nspin), wtk(nk), occ_threshold
    integer, intent(out) :: nact
    real(8), intent(out) :: dropped
    integer :: ik, is, io, cnt, cmax
    if (nstate_active_in > 0) then
      nact = min(nstate_active_in, nstate)
    else if (occ_threshold > 0d0) then
      cmax = 0
      do is=1,nspin; do ik=1,nk
        cnt = count(rocc(:,ik,is) > occ_threshold)
        if (cnt > cmax) cmax = cnt
      end do; end do
      nact = max(1, min(cmax, nstate))
    else
      nact = nstate
    end if
    dropped = 0d0
    if (nact < nstate) then
      do is=1,nspin; do ik=1,nk; do io=nact+1,nstate
        dropped = dropped + rocc(io,ik,is)*wtk(ik)
      end do; end do; end do
    end if
  end subroutine
end module
