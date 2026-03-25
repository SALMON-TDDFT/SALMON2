  subroutine calculate_density_from_fragments(dg_frag, system, mg, rho, rho_s)
    use structures
    use salmon_global, only: nelec, nelec_spin
    use communication, only: comm_summation, comm_isend, comm_irecv, comm_wait_all
    use rt_dg_fragment_ops, only: fetch_remote_coef_pw_rows
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    type(s_dft_system),     intent(in)    :: system
    type(s_rgrid),          intent(in)    :: mg
    type(s_scalar),         intent(inout) :: rho
    type(s_scalar),         intent(inout) :: rho_s(system%nspin)

    integer, parameter :: rho_tag_base = 4101
    integer, parameter :: w_tag_base = 4201
    integer :: ifrag, io, i_local, ispin
    integer :: istate_frag
    integer :: ix, iy, iz, ixg, iyg, izg, owner_rank
    integer :: ig_i, nbf, ipw, n_pw
    integer :: nxyz(3), ixyz0(3)
    integer :: nocc_per_spin, nocc_spin
    integer :: irank, nreq_send, nreq_recv, ireq
    real(8) :: occ_factor
    complex(8) :: coef_i, psi_val, phase_pw, ci
    real(8) :: phi_i, rho_contrib
    real(8) :: total_charge, total_charge_local, scale_rho, elec_num_scaled_local
    real(8) :: rx, ry, rz, boxL(3), inv_sqrt_vol, theta
    real(8), allocatable :: w_local(:,:,:)
    integer, allocatable :: pw_row_ids(:), req_send(:), req_recv(:)
    logical, allocatable :: send_active(:), recv_active(:)
    complex(8), allocatable :: coef_pw_full(:,:,:)
    type(s_scalar), allocatable :: rho_send(:), w_send(:), rho_recv(:), w_recv(:)

    rho%f = 0.0d0
    do ispin = 1, system%nspin
      rho_s(ispin)%f = 0.0d0
    end do

    if (.not. allocated(dg_frag%phi_frag)) return

    allocate(w_local(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
    allocate(send_active(0:dg_frag%isize-1), recv_active(0:dg_frag%isize-1))
    allocate(rho_send(0:dg_frag%isize-1), w_send(0:dg_frag%isize-1))
    allocate(rho_recv(0:dg_frag%isize-1), w_recv(0:dg_frag%isize-1))

    rho%f = 0.0d0
    w_local = 0.0d0
    send_active = .false.
    recv_active = .false.

    ! Closed-shell fallback: nelec = 2 * nocc_per_spin.
    nocc_per_spin = min(dg_frag%nstate_tot, int(nelec / 2.0d0 + 1.0d-12))
    occ_factor = 2.0d0 / real(system%nspin, 8)
    n_pw = 0
    if (dg_frag%use_plane_wave_basis .and. allocated(dg_frag%coef_pw) .and. allocated(dg_frag%k_pw)) then
      n_pw = dg_frag%n_plane_waves
    end if
    boxL(1) = dg_frag%hgs(1) * real(mg%num(1), 8)
    boxL(2) = dg_frag%hgs(2) * real(mg%num(2), 8)
    boxL(3) = dg_frag%hgs(3) * real(mg%num(3), 8)
    inv_sqrt_vol = 1.0d0 / sqrt(max(1.0d-16, boxL(1) * boxL(2) * boxL(3)))

    do irank = 0, dg_frag%isize - 1
      if (irank == dg_frag%id) cycle
      if (dg_frag%is_frag_root) then
        if (local_fragments_overlap_rank_box(dg_frag, mg, irank)) then
          send_active(irank) = .true.
          allocate(rho_send(irank)%f(mg%is_all(1, irank):mg%ie_all(1, irank), &
                                     mg%is_all(2, irank):mg%ie_all(2, irank), &
                                     mg%is_all(3, irank):mg%ie_all(3, irank)))
          allocate(w_send(irank)%f(mg%is_all(1, irank):mg%ie_all(1, irank), &
                                   mg%is_all(2, irank):mg%ie_all(2, irank), &
                                   mg%is_all(3, irank):mg%ie_all(3, irank)))
          rho_send(irank)%f = 0.0d0
          w_send(irank)%f = 0.0d0
        end if
      end if
      if (rank_fragments_overlap_local_box(dg_frag, mg, irank)) then
        recv_active(irank) = .true.
        allocate(rho_recv(irank)%f(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
        allocate(w_recv(irank)%f(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3)))
        rho_recv(irank)%f = 0.0d0
        w_recv(irank)%f = 0.0d0
      end if
    end do

    if (dg_frag%is_frag_root) then
      if (n_pw > 0) then
        allocate(pw_row_ids(n_pw))
        do ipw = 1, n_pw
          pw_row_ids(ipw) = ipw
        end do
        allocate(coef_pw_full(n_pw, dg_frag%nstate_tot, dg_frag%nspin))
        coef_pw_full(:, :, :) = (0.0d0, 0.0d0)
        call fetch_remote_coef_pw_rows(dg_frag, pw_row_ids, coef_pw_full)
        deallocate(pw_row_ids)
      end if

      i_local = 0
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = i_local + 1

        nxyz(1:3) = dg_frag%nxyz_domain(1:3, ifrag)
        ixyz0(1:3) = dg_frag%ixyz_frag(1:3, ifrag)

        do iz = 1, nxyz(3)
          do iy = 1, nxyz(2)
            do ix = 1, nxyz(1)
              ixg = ixyz0(1) + ix - 1
              iyg = ixyz0(2) + iy - 1
              izg = ixyz0(3) + iz - 1

              ixg = mod(ixg - 1, mg%num(1)) + 1
              iyg = mod(iyg - 1, mg%num(2)) + 1
              izg = mod(izg - 1, mg%num(3)) + 1

              owner_rank = find_grid_owner(ixg, iyg, izg)
              if (owner_rank == dg_frag%id) then
                w_local(ixg, iyg, izg) = w_local(ixg, iyg, izg) + 1.0d0
              else if (send_active(owner_rank)) then
                w_send(owner_rank)%f(ixg, iyg, izg) = w_send(owner_rank)%f(ixg, iyg, izg) + 1.0d0
              end if
            end do
          end do
        end do

        do ispin = 1, system%nspin
          nocc_spin = nocc_per_spin
          if (system%nspin == 2 .and. sum(nelec_spin(:)) > 0) then
            nocc_spin = min(dg_frag%nstate_tot, nelec_spin(ispin))
          end if
          nbf = dg_frag%n_basis(ifrag, ispin)
          do io = 1, nocc_spin
            do iz = 1, nxyz(3)
              do iy = 1, nxyz(2)
                do ix = 1, nxyz(1)
                  ixg = ixyz0(1) + ix - 1
                  iyg = ixyz0(2) + iy - 1
                  izg = ixyz0(3) + iz - 1

                  ixg = mod(ixg - 1, mg%num(1)) + 1
                  iyg = mod(iyg - 1, mg%num(2)) + 1
                  izg = mod(izg - 1, mg%num(3)) + 1

                  psi_val = (0.0d0, 0.0d0)
                  do istate_frag = 1, nbf
                    ig_i = dg_frag%index_basis(istate_frag, ifrag, ispin)
                    if (ig_i < 1 .or. ig_i > dg_frag%n_mat_max) cycle
                    coef_i = dg_frag%coef(ig_i, io, ispin)
                    phi_i = dg_frag%phi_frag(ix, iy, iz, istate_frag, i_local)
                    psi_val = psi_val + coef_i * phi_i
                  end do
                  if (n_pw > 0) then
                    rx = real(ixg - 1, 8) * dg_frag%hgs(1)
                    ry = real(iyg - 1, 8) * dg_frag%hgs(2)
                    rz = real(izg - 1, 8) * dg_frag%hgs(3)
                    do ipw = 1, n_pw
                      theta = dg_frag%k_pw(1, ipw) * rx + dg_frag%k_pw(2, ipw) * ry + dg_frag%k_pw(3, ipw) * rz
                      phase_pw = cmplx(cos(theta), sin(theta), kind=8)
                      ci = coef_pw_full(ipw, io, ispin)
                      psi_val = psi_val + ci * phase_pw * inv_sqrt_vol
                    end do
                  end if

                  rho_contrib = occ_factor * real(conjg(psi_val) * psi_val, kind=8)

                  owner_rank = find_grid_owner(ixg, iyg, izg)
                  if (owner_rank == dg_frag%id) then
                    rho%f(ixg, iyg, izg) = rho%f(ixg, iyg, izg) + rho_contrib
                  else if (send_active(owner_rank)) then
                    rho_send(owner_rank)%f(ixg, iyg, izg) = rho_send(owner_rank)%f(ixg, iyg, izg) + rho_contrib
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
    end if

    nreq_recv = 0
    do irank = 0, dg_frag%isize - 1
      if (recv_active(irank)) nreq_recv = nreq_recv + 2
    end do
    if (nreq_recv > 0) then
      allocate(req_recv(nreq_recv))
      ireq = 0
      do irank = 0, dg_frag%isize - 1
        if (.not. recv_active(irank)) cycle
        ireq = ireq + 1
        req_recv(ireq) = comm_irecv(rho_recv(irank)%f, irank, rho_tag_base + dg_frag%id, dg_frag%icomm)
        ireq = ireq + 1
        req_recv(ireq) = comm_irecv(w_recv(irank)%f, irank, w_tag_base + dg_frag%id, dg_frag%icomm)
      end do
    end if

    nreq_send = 0
    do irank = 0, dg_frag%isize - 1
      if (send_active(irank)) nreq_send = nreq_send + 2
    end do
    if (nreq_send > 0) then
      allocate(req_send(nreq_send))
      ireq = 0
      do irank = 0, dg_frag%isize - 1
        if (.not. send_active(irank)) cycle
        ireq = ireq + 1
        req_send(ireq) = comm_isend(rho_send(irank)%f, irank, rho_tag_base + irank, dg_frag%icomm)
        ireq = ireq + 1
        req_send(ireq) = comm_isend(w_send(irank)%f, irank, w_tag_base + irank, dg_frag%icomm)
      end do
    end if

    if (nreq_recv > 0) then
      call comm_wait_all(req_recv)
      do irank = 0, dg_frag%isize - 1
        if (.not. recv_active(irank)) cycle
        rho%f = rho%f + rho_recv(irank)%f
        w_local = w_local + w_recv(irank)%f
        deallocate(rho_recv(irank)%f, w_recv(irank)%f)
      end do
      deallocate(req_recv)
    end if
    if (nreq_send > 0) then
      call comm_wait_all(req_send)
      do irank = 0, dg_frag%isize - 1
        if (.not. send_active(irank)) cycle
        deallocate(rho_send(irank)%f, w_send(irank)%f)
      end do
      deallocate(req_send)
    end if

    where (w_local > 0.5d0)
      rho%f = rho%f / w_local
    end where

    total_charge_local = sum(rho%f) * system%hvol
    call comm_summation(total_charge_local, total_charge, dg_frag%icomm)
    dg_frag%elec_num_raw = total_charge
    dg_frag%rho_scale_factor = 1.0d0
    if (total_charge > 1.0d-14 .and. total_charge == total_charge) then
      scale_rho = nelec / total_charge
      dg_frag%rho_scale_factor = scale_rho
      rho%f = rho%f * scale_rho
    end if
    elec_num_scaled_local = sum(rho%f) * system%hvol
    call comm_summation(elec_num_scaled_local, dg_frag%elec_num_scaled, dg_frag%icomm)

    do ispin = 1, system%nspin
      rho_s(ispin)%f = rho%f / real(system%nspin, 8)
    end do

    if (allocated(coef_pw_full)) deallocate(coef_pw_full)
    deallocate(w_local, send_active, recv_active, rho_send, w_send, rho_recv, w_recv)

  contains

  integer function find_grid_owner(ixg, iyg, izg) result(owner)
    implicit none
    integer, intent(in) :: ixg, iyg, izg
    integer :: jrank

    owner = dg_frag%id
    do jrank = 0, dg_frag%isize - 1
      if (ixg < mg%is_all(1, jrank) .or. ixg > mg%ie_all(1, jrank)) cycle
      if (iyg < mg%is_all(2, jrank) .or. iyg > mg%ie_all(2, jrank)) cycle
      if (izg < mg%is_all(3, jrank) .or. izg > mg%ie_all(3, jrank)) cycle
      owner = jrank
      return
    end do
  end function find_grid_owner

  logical function local_fragments_overlap_rank_box(dg_frag_local, mg_local, rank_box) result(overlap)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag_local
    type(s_rgrid), intent(in) :: mg_local
    integer, intent(in) :: rank_box
    integer :: jfrag

    overlap = .false.
    do jfrag = dg_frag_local%ifrag_start, dg_frag_local%ifrag_end
      if (fragment_overlaps_box(dg_frag_local%ixyz_frag(:, jfrag), dg_frag_local%nxyz_domain(:, jfrag), &
          mg_local%is_all(:, rank_box), mg_local%ie_all(:, rank_box), mg_local%num)) then
        overlap = .true.
        return
      end if
    end do
  end function local_fragments_overlap_rank_box

  logical function rank_fragments_overlap_local_box(dg_frag_local, mg_local, source_rank) result(overlap)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag_local
    type(s_rgrid), intent(in) :: mg_local
    integer, intent(in) :: source_rank
    integer :: jfrag

    overlap = .false.
    do jfrag = 1, dg_frag_local%n_frag
      if (dg_frag_local%id_array(jfrag) /= source_rank) cycle
      if (fragment_overlaps_box(dg_frag_local%ixyz_frag(:, jfrag), dg_frag_local%nxyz_domain(:, jfrag), &
          mg_local%is, mg_local%ie, mg_local%num)) then
        overlap = .true.
        return
      end if
    end do
  end function rank_fragments_overlap_local_box

  logical function fragment_overlaps_box(ixyz_frag, nxyz_domain, box_is, box_ie, num_grid) result(overlap)
    implicit none
    integer, intent(in) :: ixyz_frag(3), nxyz_domain(3), box_is(3), box_ie(3), num_grid(3)
    integer :: idir

    overlap = .true.
    do idir = 1, 3
      if (.not. periodic_range_overlaps(ixyz_frag(idir), nxyz_domain(idir), box_is(idir), box_ie(idir), num_grid(idir))) then
        overlap = .false.
        return
      end if
    end do
  end function fragment_overlaps_box

  logical function periodic_range_overlaps(start_idx, length_idx, box_is, box_ie, num_grid) result(overlap)
    implicit none
    integer, intent(in) :: start_idx, length_idx, box_is, box_ie, num_grid
    integer :: current_start, remaining_len, segment_len, seg_is, seg_ie

    overlap = .false.
    if (length_idx <= 0) return
    current_start = mod(start_idx - 1, num_grid) + 1
    remaining_len = length_idx
    do while (remaining_len > 0)
      segment_len = min(remaining_len, num_grid - current_start + 1)
      seg_is = current_start
      seg_ie = current_start + segment_len - 1
      if (max(seg_is, box_is) <= min(seg_ie, box_ie)) then
        overlap = .true.
        return
      end if
      remaining_len = remaining_len - segment_len
      current_start = 1
    end do
  end function periodic_range_overlaps

  end subroutine calculate_density_from_fragments
