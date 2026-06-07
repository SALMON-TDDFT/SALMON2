! DG fragment basis-file I/O for non-SOI RT.
#include "config.h"
module rt_dg_fragment_io
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use rt_dg_fragment_coefficients, only: rebuild_coef_owner_map
  use rt_dg_fragment_layout, only: get_fragment_group_root_rank
  implicit none
  private
  public :: read_fragment_basis_data

contains
  subroutine read_dg_occupation_seed(dg_frag)
    use filesystem, only: get_filehandle
    use communication, only: comm_is_root, comm_bcast
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer :: iunit, nspin_file, iostat_open
    integer, allocatable :: occ_file(:)
    logical :: has_file
    character(256) :: filename

    has_file = .false.
    nspin_file = dg_frag%nspin
    if (comm_is_root(dg_frag%id)) then
      filename = './data_dcdft/total/dg_occupation.bin'
      iunit = get_filehandle()
      open(iunit, file=filename, form='unformatted', access='stream', status='old', iostat=iostat_open)
      if (iostat_open == 0) then
        has_file = .true.
        read(iunit) nspin_file
        if (nspin_file >= 1) then
          allocate(occ_file(nspin_file))
          read(iunit) occ_file(1:nspin_file)
        end if
        close(iunit)
      end if
    end if

    call comm_bcast(has_file, dg_frag%icomm, 0)
    call comm_bcast(nspin_file, dg_frag%icomm, 0)
    if (.not. has_file) return
    if (nspin_file < 1) return

    if (.not. allocated(occ_file)) allocate(occ_file(nspin_file))
    call comm_bcast(occ_file, dg_frag%icomm, 0)

    if (nspin_file /= dg_frag%nspin) then
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,i0,a,i0,a)') "[WARN] dg_occupation.bin nspin mismatch: file=", nspin_file, &
          " runtime=", dg_frag%nspin, " (ignoring seed occupancy)"
      end if
      deallocate(occ_file)
      return
    end if

    if (.not. allocated(dg_frag%nocc_spin)) allocate(dg_frag%nocc_spin(dg_frag%nspin))
    dg_frag%nocc_spin(1:dg_frag%nspin) = min(dg_frag%nstate_tot, max(0, occ_file(1:dg_frag%nspin)))
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,10(1x,i0))') "[INFO] DG occupancy seed loaded:", dg_frag%nocc_spin(1:dg_frag%nspin)
    end if
    deallocate(occ_file)
  end subroutine read_dg_occupation_seed


  !=======================================================================
  ! Read fragment basis data from DC-LCFO calculation (MPI-parallelized)
  !=======================================================================
  subroutine read_fragment_basis_data(dg_frag, bdir_frag)
    use filesystem, only: get_filehandle
    use communication, only: comm_is_root, comm_bcast, comm_sync_all, comm_summation, comm_get_max
    use salmon_global, only: dg_nmat_cap_mode, dg_nmat_cap_fixed, dg_nmat_cap_multiple, nelec, nelec_spin, &
                             dg_subspace_extra_states
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    character(*), intent(in) :: bdir_frag

    character(32), parameter :: binfile_wf = "wavefunctions.bin"
    character(32), parameter :: binfile_bfb = "basis_functions_buffer.bin"
    character(32), parameter :: binfile_rg = "rgrid_index.bin"
    integer, parameter :: basis_buffer_magic = -22022212
    integer, parameter :: basis_buffer_version = 1
    character(256) :: filename
    integer :: iunit, ifrag, ispin, n_frag_file, nspin_file
    integer :: nstate_frag_file, nstate_tot_file, nstate_frag_keep
    integer, allocatable :: n_basis_tmp(:,:), index_basis_file(:,:,:)
    real(8), allocatable :: coef_state_file(:)
    integer :: n_mat_tmp(2)   ! nspin is expected to be 1 or 2
    integer :: ifrag_count, i_local, io, global_idx
    integer :: local_coef_max, local_idx
    integer :: nxyz_domain(3), nxyz_alloc(3), lgnum_frag(3), lgnum_total(3)
    integer :: nxyz_buffer_file(3), nxyz_box(3)
    integer :: magic_file, version_file
    integer, allocatable :: n_basis_frag(:)
    integer, allocatable :: jxyz_tot(:,:)
    integer :: ix, iy, iz, n
    integer :: ixg_store, iyg_store, izg_store
    integer :: ix_src, iy_src, iz_src
    integer :: ix_box, iy_box, iz_box
    integer :: nb  ! halo width
    integer :: nbasis_iter
    integer :: n_mat_cap, n_mat_cap_env, ienv
    integer :: nocc_max, nocc_eff, ifrag_best, occ_min, occ_max, cap_min, cap_max
    integer :: nocc_est, nocc_frag_est, n_basis_min, n_basis_max, n_basis_effective_min
    integer :: nvirt_est_min, min_virtual_required
    real(8) :: n_basis_avg, nvirt_est_avg
    integer :: state_col, occ_base, occ_extra, frag_occ_s, frag_occ_e, nocc_frag_seed
    integer :: env_status, env_len
    character(len=64) :: env_n_mat_cap
    logical :: warned_spin_discard, has_buffer_file, identity_seed_coefficients
    real(8) :: cap_avg, weight_best
    real(8) :: coef_diag_local(2), coef_diag_global(2)
    real(8), allocatable :: frag_weight_local(:,:,:), frag_weight_sum(:,:,:)
    integer, allocatable :: occ_count(:,:), cap_frag(:,:)

    ! Step 1: Root reads metadata from first fragment and broadcasts
    if (comm_is_root(dg_frag%id)) then
      ifrag = 1
      iunit = get_filehandle()
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), ifrag, '/', binfile_wf

      open(iunit, file=filename, form='unformatted', access='stream', status='old')
      read(iunit) n_frag_file, nspin_file, nstate_frag_file, nstate_tot_file

      if (n_frag_file /= dg_frag%n_frag .or. nspin_file /= dg_frag%nspin) then
        write(*,*) "Error: Fragment basis data mismatch"
        write(*,*) "  Expected n_frag=", dg_frag%n_frag, ", nspin=", dg_frag%nspin, &
                   ", nstate_frag=", dg_frag%nstate_frag
        write(*,*) "  Found    n_frag=", n_frag_file, ", nspin=", nspin_file, &
                   ", nstate_frag=", nstate_frag_file
        stop "DG-Fragment RT: Fragment basis data mismatch"
      end if

      close(iunit)
    end if

    ! Broadcast metadata to all ranks
    call comm_bcast(n_frag_file, dg_frag%icomm, 0)
    call comm_bcast(nspin_file, dg_frag%icomm, 0)
    call comm_bcast(nstate_frag_file, dg_frag%icomm, 0)
    call comm_bcast(nstate_tot_file, dg_frag%icomm, 0)

    identity_seed_coefficients = (nstate_tot_file < 0)
    if (identity_seed_coefficients) nstate_tot_file = -nstate_tot_file
    dg_frag%identity_seed_coefficients = identity_seed_coefficients
    dg_frag%fragment_basis_contracted = .false.
    dg_frag%defer_fragment_cap_to_local_eigen = .false.
    dg_frag%requested_fragment_basis_cap = 0
    dg_frag%seed_basis_runtime_capped = .false.
    if (dg_frag%id == 0) then
      if (identity_seed_coefficients) then
        write(*,'(1x,a)') "[WARN] DG identity seed detected: initial occupied states are fragment-local."
        write(*,'(1x,a)') "       This is the non-LCFO DG seed export, not a DGDFT/LCFO ground-state wavefunction."
      else
        write(*,'(1x,a)') "[INFO] DG LCFO coefficient seed detected: loading DC eigenvector coefficients."
      end if
    end if

    nstate_frag_keep = nstate_frag_file
    if (dg_frag%id == 0) then
      write(*,'(1x,a,a,a,i0,a,1pe12.4)') "[INFO] DG fragment cap input: mode='", trim(dg_nmat_cap_mode), &
        "' fixed=", dg_nmat_cap_fixed, " multiple=", dg_nmat_cap_multiple
    end if
    select case (trim(dg_nmat_cap_mode))
    case ('fragment_fixed', 'fixed_fragment')
      if (dg_nmat_cap_fixed < 1) then
        if (dg_frag%id == 0) then
          write(*,'(1x,a,a,a,i0)') "[FATAL] DG fragment cap mode='", trim(dg_nmat_cap_mode), &
            "' requires dg_nmat_cap_fixed >= 1; got ", dg_nmat_cap_fixed
        end if
        stop "DG-Fragment RT: invalid fragment basis cap"
      end if
      if (identity_seed_coefficients) then
        dg_frag%defer_fragment_cap_to_local_eigen = .true.
        dg_frag%requested_fragment_basis_cap = dg_nmat_cap_fixed
        if (dg_frag%id == 0) then
          write(*,'(1x,a,a,a,i0,a,i0,a,i0)') "[INFO] DG legacy identity fragment cap mode='", trim(dg_nmat_cap_mode), &
            "' file_nstate_frag=", nstate_frag_file, " keep_for_local_diag=", nstate_frag_keep, &
            " legacy_requested=", dg_nmat_cap_fixed
        end if
      else
        if (dg_frag%id == 0) then
          write(*,'(1x,a,a,a,i0,a)') "[WARN] DGDFT/LCFO seed ignores dg_nmat_cap_mode='", &
            trim(dg_nmat_cap_mode), "' fixed=", dg_nmat_cap_fixed, " and keeps the full DC-derived basis."
        end if
      end if
    case default
      if (nstate_frag_file /= dg_frag%nstate_frag) then
        if (dg_frag%id == 0) then
          write(*,'(1x,a,i0,a,i0,a)') "[INFO] nstate_frag differs: file=", nstate_frag_file, &
            " runtime=", dg_frag%nstate_frag, " (using fragment-state count from file)"
        end if
      end if
    end select
    dg_frag%nstate_frag = nstate_frag_keep

    ! Use the full state count stored in fragment files (disable occupied-state subset mode).
    if (nstate_tot_file /= dg_frag%nstate_tot) then
      if (dg_frag%id == 0) then
        write(*,'(1x,a,i0,a,i0,a)') "[INFO] nstate_tot differs: file=", nstate_tot_file, &
          " runtime=", dg_frag%nstate_tot, " (using full-state count from file)"
      end if
      dg_frag%nstate_tot = nstate_tot_file
    end if

    ! The identity seed encodes fragment-local basis states, not a dense global
    ! eigenvector matrix.  Load the occupied-column count before constructing
    ! coef so occupied columns can be distributed over every fragment.
    if (identity_seed_coefficients) call read_dg_occupation_seed(dg_frag)

    ! Allocate arrays
    allocate(dg_frag%n_basis(dg_frag%n_frag, dg_frag%nspin))
    if (.not. allocated(dg_frag%index_basis)) then
      allocate(dg_frag%index_basis(dg_frag%nstate_frag, dg_frag%n_frag, dg_frag%nspin))
    end if
    if (.not. allocated(dg_frag%n_mat)) then
      allocate(dg_frag%n_mat(dg_frag%nspin))
    end if
    allocate(n_basis_tmp(dg_frag%n_frag, dg_frag%nspin))
    allocate(index_basis_file(nstate_frag_file, dg_frag%n_frag, dg_frag%nspin))

    ! All ranks read global metadata from fragment 1.
    ! index_basis maps local->global indices across ALL fragments; every rank
    ! needs the full table or ranks > 0 will produce zero/NaN current.
    iunit = get_filehandle()
    write(filename, '(a, i6.6, a, a)') trim(bdir_frag), 1, '/', binfile_wf
    open(iunit, file=filename, form='unformatted', access='stream', status='old')
    read(iunit) n_frag_file, nspin_file, nstate_frag_file, nstate_tot_file
    if (nstate_tot_file < 0) nstate_tot_file = -nstate_tot_file
    read(iunit) n_mat_tmp(1:dg_frag%nspin)
    read(iunit) n_basis_tmp(1:dg_frag%n_frag, 1:dg_frag%nspin)
    read(iunit) index_basis_file(1:nstate_frag_file, 1:dg_frag%n_frag, 1:dg_frag%nspin)
    close(iunit)

    ! Step 3: Gather metadata (now consistent across all ranks)
    dg_frag%n_basis = min(n_basis_tmp, dg_frag%nstate_frag)
    dg_frag%index_basis(:, :, :) = 0
    dg_frag%index_basis(1:dg_frag%nstate_frag, 1:dg_frag%n_frag, 1:dg_frag%nspin) = &
      index_basis_file(1:dg_frag%nstate_frag, 1:dg_frag%n_frag, 1:dg_frag%nspin)
    dg_frag%n_mat(1:dg_frag%nspin) = n_mat_tmp(1:dg_frag%nspin)
    dg_frag%n_mat_max = maxval(n_mat_tmp(1:dg_frag%nspin))

    ! Optional legacy basis-size cap for fragment-comparison studies.
    ! DGDFT/LCFO coefficient seeds keep the full DC-derived basis; cap is a
    ! later optimization path, not part of the primary paper route.
    n_mat_cap = 0
    if (trim(dg_nmat_cap_mode) == 'fixed' .and. dg_nmat_cap_fixed >= 1) then
      if (identity_seed_coefficients) then
        n_mat_cap = dg_nmat_cap_fixed
        if (dg_frag%id == 0) then
          write(*,'(1x,a,a,a,i0)') "[INFO] DG fragment cap mode='", trim(dg_nmat_cap_mode), "' fixed=", n_mat_cap
        end if
      else
        if (dg_frag%id == 0) then
          write(*,'(1x,a,a,a,i0,a)') "[WARN] DGDFT/LCFO coefficient seed ignores dg_nmat_cap_mode='", &
            trim(dg_nmat_cap_mode), "' fixed=", dg_nmat_cap_fixed, " and keeps the full DC-derived basis."
        end if
      end if
    end if
    env_n_mat_cap = ""
    env_status = 1
    env_len = 0
    call get_environment_variable("SALMON_DG_NMAT_CAP", env_n_mat_cap, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      read(env_n_mat_cap(1:env_len), *, iostat=ienv) n_mat_cap_env
      if (ienv == 0 .and. n_mat_cap_env >= 1) then
        if (identity_seed_coefficients) then
          n_mat_cap = n_mat_cap_env
          if (dg_frag%id == 0) then
            write(*,'(1x,a,i0)') "[INFO] SALMON_DG_NMAT_CAP override applied: ", n_mat_cap
          end if
        else
          if (dg_frag%id == 0) then
            write(*,'(1x,a,i0,a)') "[WARN] DGDFT/LCFO coefficient seed ignores SALMON_DG_NMAT_CAP=", &
              n_mat_cap_env, " and keeps the full DC-derived basis."
          end if
        end if
      else
        if (dg_frag%id == 0) then
          write(*,'(1x,a,a,a)') "[WARN] Ignoring invalid SALMON_DG_NMAT_CAP='", &
                                trim(env_n_mat_cap(1:env_len)), "' (must be integer >= 1)"
        end if
      end if
    end if
    if (n_mat_cap >= 1) then
      if (any(dg_frag%n_mat(1:dg_frag%nspin) > n_mat_cap)) dg_frag%seed_basis_runtime_capped = .true.
      dg_frag%n_mat(1:dg_frag%nspin) = min(dg_frag%n_mat(1:dg_frag%nspin), n_mat_cap)
      dg_frag%n_mat_max = maxval(dg_frag%n_mat(1:dg_frag%nspin))
    end if

    ! Enforce index_basis consistency with the active matrix dimension:
    ! invalid indices are masked out to prevent OOB accesses.
    block
      integer :: ispin_cap, ifrag_cap, io_cap, idx_cap, max_keep
      do ispin_cap = 1, dg_frag%nspin
        max_keep = 0
        do ifrag_cap = 1, dg_frag%n_frag
          nbasis_iter = min(dg_frag%n_basis(ifrag_cap, ispin_cap), size(dg_frag%index_basis, 1))
          do io_cap = 1, nbasis_iter
            idx_cap = dg_frag%index_basis(io_cap, ifrag_cap, ispin_cap)
            if (idx_cap < 1 .or. idx_cap > dg_frag%n_mat(ispin_cap)) then
              dg_frag%index_basis(io_cap, ifrag_cap, ispin_cap) = 0
            else
              max_keep = max(max_keep, idx_cap)
            end if
          end do
        end do
        dg_frag%n_mat(ispin_cap) = max_keep
      end do
    end block
    dg_frag%n_mat_max = max(1, maxval(dg_frag%n_mat(1:dg_frag%nspin)))

    ! Compress fragmented/global basis indices to a dense contiguous range.
    ! The DC-LCFO metadata may retain large holes between fragment-local basis blocks,
    ! which inflates n_mat_max and all O(n_mat_max^2) operator matrices.
    block
      integer :: ispin_cmp, ifrag_cmp, io_cmp, idx_cmp, n_old, n_new
      integer, allocatable :: remap(:)
      do ispin_cmp = 1, dg_frag%nspin
        n_old = max(1, dg_frag%n_mat(ispin_cmp))
        allocate(remap(n_old))
        remap = 0
        n_new = 0
        do ifrag_cmp = 1, dg_frag%n_frag
          nbasis_iter = min(dg_frag%n_basis(ifrag_cmp, ispin_cmp), size(dg_frag%index_basis, 1))
          do io_cmp = 1, nbasis_iter
            idx_cmp = dg_frag%index_basis(io_cmp, ifrag_cmp, ispin_cmp)
            if (idx_cmp <= 0) cycle
            if (idx_cmp > n_old) then
              dg_frag%index_basis(io_cmp, ifrag_cmp, ispin_cmp) = 0
              cycle
            end if
            if (remap(idx_cmp) == 0) then
              n_new = n_new + 1
              remap(idx_cmp) = n_new
            end if
            dg_frag%index_basis(io_cmp, ifrag_cmp, ispin_cmp) = remap(idx_cmp)
          end do
        end do
        if (dg_frag%id == 0 .and. n_new < n_old) then
          write(*,'(1x,a,i0,a,i0,a,i0)') "[INFO] Compressed DG basis indices for ispin=", ispin_cmp, &
            ": old n_mat=", n_old, " new n_mat=", n_new
        end if
        if (n_new <= 0) then
          write(*,'(1x,a,i0,a)') "[FATAL] DG basis compression produced zero active basis for ispin=", ispin_cmp, "."
          stop "DG-Fragment RT: zero active basis after index compression"
        end if
        dg_frag%n_mat(ispin_cmp) = n_new
        deallocate(remap)
      end do
    end block
    dg_frag%n_mat_max = max(1, maxval(dg_frag%n_mat(1:dg_frag%nspin)))


    ! Validate index_basis uniqueness and coverage (root only)
    if (comm_is_root(dg_frag%id)) then
      block
        integer :: ispin_chk, ifrag_chk, io_chk, global_idx
        integer :: dup_count, out_count, miss_count
        integer, allocatable :: seen(:)
        do ispin_chk = 1, dg_frag%nspin
          allocate(seen(max(1, dg_frag%n_mat_max)))
          seen = 0
          dup_count = 0
          out_count = 0
          do ifrag_chk = 1, dg_frag%n_frag
            nbasis_iter = min(dg_frag%n_basis(ifrag_chk, ispin_chk), size(dg_frag%index_basis, 1))
            do io_chk = 1, nbasis_iter
              global_idx = dg_frag%index_basis(io_chk, ifrag_chk, ispin_chk)
              if (global_idx < 1 .or. global_idx > dg_frag%n_mat_max) then
                out_count = out_count + 1
              else
                if (seen(global_idx) == 1) dup_count = dup_count + 1
                seen(global_idx) = 1
              end if
            end do
          end do
          miss_count = count(seen == 0)
          if (dup_count > 0 .or. out_count > 0 .or. miss_count > 0) then
            write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[WARN] index_basis check (ispin=", ispin_chk, &
              "): dup=", dup_count, " out_of_range=", out_count, " missing=", miss_count
          end if
          deallocate(seen)
        end do
      end block
    end if

    dg_frag%owned_coef_start = 0
    dg_frag%owned_coef_end = -1

    ! Build the row-owner map before allocating coefficient storage.  In
    ! orbital mode this keeps only the subgroup-owned basis rows on each rank,
    ! avoiding a temporary full-fragment coefficient replica during seed load.
    call rebuild_coef_owner_map(dg_frag, "read-fragment-basis")
    local_coef_max = max(1, maxval(dg_frag%local_coef_count(1:dg_frag%nspin)))

    ! Step 4: nstate_tot was aligned to file metadata above (full-state mode).

    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    if (allocated(dg_frag%phi_frag_has_seed_buffer)) deallocate(dg_frag%phi_frag_has_seed_buffer)
    allocate(dg_frag%phi_frag_has_seed_buffer(ifrag_count))
    dg_frag%phi_frag_has_seed_buffer(:) = .false.

    if (n_mat_cap < 1 .and. trim(dg_nmat_cap_mode) == 'occ_multiple' .and. dg_nmat_cap_multiple > 0.0d0) then
      if (.not. identity_seed_coefficients) then
        if (dg_frag%id == 0) then
          write(*,'(1x,a)') "[WARN] DGDFT/LCFO coefficient seed ignores occ_multiple cap and keeps the full DC-derived basis."
        end if
      else
        if (dg_frag%nspin == 1) then
          nocc_max = max(1, min((nelec + 1) / 2, dg_frag%nstate_tot))
        else if (sum(nelec_spin(1:dg_frag%nspin)) > 0) then
          nocc_max = max(1, min(maxval(nelec_spin(1:dg_frag%nspin)), dg_frag%nstate_tot))
        else
          nocc_max = max(1, min(int(nelec / 2.0d0 + 1.0d-12), dg_frag%nstate_tot))
        end if

        allocate(frag_weight_local(dg_frag%n_frag, nocc_max, dg_frag%nspin))
        allocate(frag_weight_sum(dg_frag%n_frag, nocc_max, dg_frag%nspin))
        allocate(occ_count(dg_frag%n_frag, dg_frag%nspin))
        allocate(cap_frag(dg_frag%n_frag, dg_frag%nspin))
        frag_weight_local(:, :, :) = 0.0d0
        frag_weight_sum(:, :, :) = 0.0d0
        occ_count(:, :) = 0
        cap_frag(:, :) = 0

      do i_local = 1, ifrag_count
        ifrag = dg_frag%ifrag_start + i_local - 1
        do ispin = 1, dg_frag%nspin
          if (dg_frag%nspin == 1) then
            nocc_eff = max(1, min((nelec + 1) / 2, dg_frag%nstate_tot))
          else if (sum(nelec_spin(1:dg_frag%nspin)) > 0) then
            nocc_eff = max(1, min(nelec_spin(ispin), dg_frag%nstate_tot))
          else
            nocc_eff = max(1, min(int(nelec / 2.0d0 + 1.0d-12), dg_frag%nstate_tot))
          end if
          nbasis_iter = min(dg_frag%n_basis(ifrag, ispin), dg_frag%nstate_frag)
          do io = 1, nbasis_iter
            global_idx = dg_frag%index_basis(io, ifrag, ispin)
            if (global_idx >= 1 .and. global_idx <= nocc_eff) then
              frag_weight_local(ifrag, global_idx, ispin) = 1.0d0
            end if
          end do
        end do
      end do
      call comm_summation(frag_weight_local, frag_weight_sum, dg_frag%n_frag * nocc_max * dg_frag%nspin, dg_frag%icomm)

      do ispin = 1, dg_frag%nspin
        if (dg_frag%nspin == 1) then
          nocc_eff = max(1, min((nelec + 1) / 2, dg_frag%nstate_tot))
        else if (sum(nelec_spin(1:dg_frag%nspin)) > 0) then
          nocc_eff = max(1, min(nelec_spin(ispin), dg_frag%nstate_tot))
        else
          nocc_eff = max(1, min(int(nelec / 2.0d0 + 1.0d-12), dg_frag%nstate_tot))
        end if
        do io = 1, nocc_eff
          ifrag_best = 1
          weight_best = frag_weight_sum(1, io, ispin)
          do ifrag = 2, dg_frag%n_frag
            if (frag_weight_sum(ifrag, io, ispin) > weight_best) then
              ifrag_best = ifrag
              weight_best = frag_weight_sum(ifrag, io, ispin)
            end if
          end do
          occ_count(ifrag_best, ispin) = occ_count(ifrag_best, ispin) + 1
        end do
        do ifrag = 1, dg_frag%n_frag
          nbasis_iter = dg_frag%n_basis(ifrag, ispin)
          cap_frag(ifrag, ispin) = min(dg_frag%n_basis(ifrag, ispin), &
                                       int(floor(dg_nmat_cap_multiple * dble(occ_count(ifrag, ispin)))))
          cap_frag(ifrag, ispin) = max(1, cap_frag(ifrag, ispin))
          if (cap_frag(ifrag, ispin) < nbasis_iter) dg_frag%seed_basis_runtime_capped = .true.
          dg_frag%n_basis(ifrag, ispin) = cap_frag(ifrag, ispin)
          if (.false. .and. dg_frag%id == 0 .and. (cap_frag(ifrag, ispin) >= 60 .or. nbasis_iter >= 60)) then
            write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,f8.3,a,i0)') "        basis cap diag: rank=", dg_frag%id, &
              " ifrag=", ifrag, " ispin=", ispin, " before=", nbasis_iter, " occ_count=", occ_count(ifrag, ispin), &
              " multiple=", dg_nmat_cap_multiple, " after=", cap_frag(ifrag, ispin)
            flush(6)
          end if
          do io = cap_frag(ifrag, ispin) + 1, min(dg_frag%nstate_frag, size(dg_frag%index_basis, 1))
            dg_frag%index_basis(io, ifrag, ispin) = 0
          end do
        end do
      end do

      do ispin = 1, dg_frag%nspin
        do ifrag = 1, dg_frag%n_frag
          dg_frag%n_mat(ispin) = max(dg_frag%n_mat(ispin), 1)
        end do
      end do

      block
        integer :: ispin_cmp, ifrag_cmp, io_cmp, idx_cmp, n_old, n_new
        integer, allocatable :: remap(:)
        do ispin_cmp = 1, dg_frag%nspin
          n_old = max(1, dg_frag%n_mat(ispin_cmp))
          allocate(remap(n_old))
          remap = 0
          n_new = 0
          do ifrag_cmp = 1, dg_frag%n_frag
            nbasis_iter = min(dg_frag%n_basis(ifrag_cmp, ispin_cmp), size(dg_frag%index_basis, 1))
            do io_cmp = 1, nbasis_iter
              idx_cmp = dg_frag%index_basis(io_cmp, ifrag_cmp, ispin_cmp)
              if (idx_cmp <= 0) cycle
              if (idx_cmp > n_old) then
                dg_frag%index_basis(io_cmp, ifrag_cmp, ispin_cmp) = 0
                cycle
              end if
              if (remap(idx_cmp) == 0) then
                n_new = n_new + 1
                remap(idx_cmp) = n_new
              end if
              dg_frag%index_basis(io_cmp, ifrag_cmp, ispin_cmp) = remap(idx_cmp)
            end do
          end do
          dg_frag%n_mat(ispin_cmp) = max(1, n_new)
          deallocate(remap)
        end do
      end block
      dg_frag%n_mat_max = max(1, maxval(dg_frag%n_mat(1:dg_frag%nspin)))

      dg_frag%owned_coef_start = 0
      dg_frag%owned_coef_end = -1

      if (dg_frag%id == 0) then
        occ_min = minval(occ_count(:, 1:dg_frag%nspin))
        occ_max = maxval(occ_count(:, 1:dg_frag%nspin))
        cap_min = minval(cap_frag(:, 1:dg_frag%nspin))
        cap_max = maxval(cap_frag(:, 1:dg_frag%nspin))
        cap_avg = sum(dble(cap_frag(:, 1:dg_frag%nspin))) / dble(dg_frag%n_frag * dg_frag%nspin)
        write(*,'(1x,a,a,a,f8.3)') "[INFO] DG fragment cap mode='", trim(dg_nmat_cap_mode), &
          "' multiple=", dg_nmat_cap_multiple
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,f8.3,a,i0)') "[INFO] DG occ/cap summary: occ_min=", occ_min, &
          " occ_max=", occ_max, " cap_min=", cap_min, " cap_max=", cap_max, " cap_avg=", cap_avg, &
          " n_mat_max=", dg_frag%n_mat_max
      end if

        deallocate(frag_weight_local, frag_weight_sum, occ_count, cap_frag)
      end if
    end if

    if (.not. identity_seed_coefficients) then
      min_virtual_required = max(8, dg_subspace_extra_states)
      do ispin = 1, dg_frag%nspin
        if (dg_frag%nspin == 1) then
          nocc_est = max(0, min(int(nelec / 2.0d0 + 1.0d-12), dg_frag%nstate_tot))
        else if (sum(nelec_spin(1:dg_frag%nspin)) > 0) then
          nocc_est = max(0, min(nelec_spin(ispin), dg_frag%nstate_tot))
        else
          nocc_est = max(0, min(int(nelec / 2.0d0 + 1.0d-12), dg_frag%nstate_tot))
        end if
        nocc_frag_est = (nocc_est + max(1, dg_frag%n_frag) - 1) / max(1, dg_frag%n_frag)
        n_basis_min = minval(dg_frag%n_basis(1:dg_frag%n_frag, ispin))
        n_basis_max = maxval(dg_frag%n_basis(1:dg_frag%n_frag, ispin))
        n_basis_avg = sum(dble(dg_frag%n_basis(1:dg_frag%n_frag, ispin))) / dble(max(1, dg_frag%n_frag))
        n_basis_effective_min = n_basis_min
        nvirt_est_min = n_basis_effective_min - nocc_frag_est
        nvirt_est_avg = n_basis_avg - dble(nocc_frag_est)
        if (dg_frag%id == 0) then
          write(*,'(1x,a,i0,6(a,i0),2(a,f8.2))') "[INFO] DG LCFO seed basis capacity: ispin=", ispin, &
            " nocc=", nocc_est, " nocc_per_frag_est=", nocc_frag_est, &
            " n_basis_min=", n_basis_min, " n_basis_max=", n_basis_max, &
            " nvirt_est_min=", nvirt_est_min, " required_min=", min_virtual_required, &
            " n_basis_avg=", n_basis_avg, " nvirt_est_avg=", nvirt_est_avg
        end if
        if (nvirt_est_min < min_virtual_required) then
          if (dg_frag%id == 0) then
            write(*,'(1x,a)') "[FATAL] DGDFT/LCFO seed has too few unoccupied fragment-basis states for RT response."
            write(*,'(1x,a,i0,a,i0,a,i0)') "[FATAL] Estimated virtual states per fragment: min=", &
              nvirt_est_min, " required_min=", min_virtual_required, " ispin=", ispin
            write(*,'(1x,a)') "[FATAL] Increase the DC-LCFO exported fragment basis size or disable/reduce RT-side caps."
          end if
          stop "DG-Fragment RT: insufficient unoccupied fragment states in DGDFT seed"
        end if
      end do
    end if

    ! Reallocate coefficient arrays with correct n_mat_max dimension
    if (allocated(dg_frag%coef)) deallocate(dg_frag%coef)
    if (allocated(dg_frag%coef_new)) deallocate(dg_frag%coef_new)
    if (allocated(dg_frag%coef_work)) deallocate(dg_frag%coef_work)
    allocate(dg_frag%coef(local_coef_max, dg_frag%nstate_tot, dg_frag%nspin))
    if (identity_seed_coefficients) then
      allocate(dg_frag%coef_work(local_coef_max, dg_frag%nstate_tot, dg_frag%nspin))
      if (dg_frag%yn_adaptive_basis) then
        allocate(dg_frag%coef_new(local_coef_max, dg_frag%nstate_tot, dg_frag%nspin))
      end if
    end if
    dg_frag%coef = 0.0d0
    if (allocated(dg_frag%coef_new)) dg_frag%coef_new = 0.0d0
    if (allocated(dg_frag%coef_work)) dg_frag%coef_work = 0.0d0

    if (identity_seed_coefficients) then
      ! Step 5a: Build the initial occupied coefficient columns.  DC writes an
      ! identity seed to avoid a huge dense coefficient file; the occupied states
      ! must therefore be assigned fragment by fragment.  Using global basis row
      ! numbers directly would occupy only the first basis blocks and leave later
      ! fragments empty in large weak-scaling runs.
      if (dg_frag%id == 0) then
        do ispin = 1, dg_frag%nspin
          nocc_eff = dg_frag%nstate_tot
          if (allocated(dg_frag%nocc_spin)) then
            if (ispin <= size(dg_frag%nocc_spin)) nocc_eff = min(dg_frag%nstate_tot, max(0, dg_frag%nocc_spin(ispin)))
          end if
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') "[INFO] DG identity seed occupancy: ispin=", ispin, &
            " nocc=", nocc_eff, " nfrag=", dg_frag%n_frag, " base=", nocc_eff / max(1, dg_frag%n_frag), &
            " extra=", mod(nocc_eff, max(1, dg_frag%n_frag))
        end do
      end if
      do i_local = 1, ifrag_count
        ifrag = dg_frag%ifrag_start + i_local - 1
        do ispin = 1, dg_frag%nspin
          nocc_eff = dg_frag%nstate_tot
          if (allocated(dg_frag%nocc_spin)) then
            if (ispin <= size(dg_frag%nocc_spin)) nocc_eff = min(dg_frag%nstate_tot, max(0, dg_frag%nocc_spin(ispin)))
          end if
          occ_base = nocc_eff / max(1, dg_frag%n_frag)
          occ_extra = mod(nocc_eff, max(1, dg_frag%n_frag))
          frag_occ_s = (ifrag - 1) * occ_base + min(ifrag - 1, occ_extra) + 1
          frag_occ_e = ifrag * occ_base + min(ifrag, occ_extra)
          nocc_frag_seed = max(0, frag_occ_e - frag_occ_s + 1)
          nbasis_iter = min(dg_frag%n_basis(ifrag, ispin), dg_frag%nstate_frag)
          if (nocc_frag_seed > nbasis_iter) then
            write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') &
              "[FATAL] DG identity seed cannot represent occupied fragment states: rank=", dg_frag%id, &
              " ifrag=", ifrag, " ispin=", ispin, " nocc_frag=", nocc_frag_seed, " n_basis=", nbasis_iter
            stop "DG-Fragment RT: insufficient fragment basis for occupied seed"
          end if
          do io = 1, nbasis_iter
            global_idx = dg_frag%index_basis(io, ifrag, ispin)
            local_idx = 0
            if (global_idx > 0 .and. global_idx <= dg_frag%n_mat_max) local_idx = dg_frag%coef_global_to_local(global_idx, ispin)
            if (local_idx > 0 .and. local_idx <= size(dg_frag%coef, 1)) then
              state_col = 0
              if (io <= nocc_frag_seed) state_col = frag_occ_s + io - 1
              if (state_col >= 1 .and. state_col <= dg_frag%nstate_tot) then
                dg_frag%coef(local_idx, state_col, ispin) = (1.0d0, 0.0d0)
              end if
            end if
          end do
        end do
      end do
    else
      ! Step 5b: DC-LCFO coefficient seed.  wavefunctions.bin stores
      ! coef_wf(local_basis, state, spin), with each state column contiguous.
      ! Stream one column at a time so the RT rank keeps only its owned rows.
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        iunit = get_filehandle()
        write(filename, '(a, i6.6, a, a)') trim(bdir_frag), ifrag, '/', binfile_wf
        open(iunit, file=filename, form='unformatted', access='stream', status='old')
        read(iunit) n_frag_file, nspin_file, nstate_frag_file, nstate_tot_file
        read(iunit) n_mat_tmp(1:dg_frag%nspin)
        read(iunit) n_basis_tmp(1:dg_frag%n_frag, 1:dg_frag%nspin)
        read(iunit) index_basis_file(1:nstate_frag_file, 1:dg_frag%n_frag, 1:dg_frag%nspin)
        if (allocated(coef_state_file)) deallocate(coef_state_file)
        allocate(coef_state_file(nstate_frag_file))
        do ispin = 1, nspin_file
          do state_col = 1, nstate_tot_file
            read(iunit) coef_state_file(1:nstate_frag_file)
            if (ispin < 1 .or. ispin > dg_frag%nspin) cycle
            if (state_col > dg_frag%nstate_tot) cycle
            nbasis_iter = min(dg_frag%n_basis(ifrag, ispin), nstate_frag_file)
            do io = 1, nbasis_iter
              global_idx = dg_frag%index_basis(io, ifrag, ispin)
              local_idx = 0
              if (global_idx > 0 .and. global_idx <= dg_frag%n_mat_max) local_idx = dg_frag%coef_global_to_local(global_idx, ispin)
              if (local_idx > 0 .and. local_idx <= size(dg_frag%coef, 1)) then
                dg_frag%coef(local_idx, state_col, ispin) = dcmplx(coef_state_file(io), 0.0d0)
              end if
            end do
          end do
        end do
        close(iunit)
        deallocate(coef_state_file)
      end do
    end if

    coef_diag_local(:) = 0.0d0
    coef_diag_global(:) = 0.0d0
    if (allocated(dg_frag%coef)) then
      coef_diag_local(1) = sum(abs(dg_frag%coef)**2)
      coef_diag_local(2) = dble(count(abs(dg_frag%coef) > 0.0d0))
    end if
    call comm_summation(coef_diag_local, coef_diag_global, 2, dg_frag%icomm)
    if (comm_is_root(dg_frag%id)) then
      if (identity_seed_coefficients) then
        write(*,'(1x,a,2(a,1pe13.5))') "[INFO] DG identity seed coefficients loaded:", &
          " norm2=", coef_diag_global(1), " nonzero=", coef_diag_global(2)
      else
        write(*,'(1x,a,2(a,1pe13.5))') "[INFO] DG rank-distributed LCFO seed coefficients loaded:", &
          " norm2=", coef_diag_global(1), " nonzero=", coef_diag_global(2)
      end if
    end if
    if (coef_diag_global(1) <= 0.0d0 .or. coef_diag_global(2) <= 0.0d0) then
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a)') "[FATAL] DG seed coefficient matrix is empty after reading wavefunctions.bin."
      end if
      stop "DG-Fragment RT: empty seed coefficient matrix"
    end if

    ! Keep coefficients only on the owning fragment ranks.

    ! Step 6: Read real-space basis functions and grid mapping
    ! This enables density reconstruction: rho(r) = sum c*_i c_j phi_i(r) phi_j(r)
    ! Also extract fragment geometry information (ixyz_frag, nxyz_domain) for halo exchange
    nb = dg_frag%nxyz_buffer(1)  ! Assume uniform buffer width (4 for 4th-order stencil)
    nxyz_alloc = 0
    warned_spin_discard = .false.

    i_local = 0
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      i_local = i_local + 1

      ! Read grid index mapping
      iunit = get_filehandle()
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), ifrag, '/', binfile_rg

      open(iunit, file=filename, form='unformatted', access='stream', status='old')
      read(iunit) lgnum_frag(1:3), lgnum_total(1:3)

      if (.not. allocated(jxyz_tot)) then
        allocate(jxyz_tot(maxval(lgnum_frag), 3))
      end if
      do n = 1, 3
        read(iunit) jxyz_tot(1:lgnum_frag(n), n)
      end do
      close(iunit)

      ! Extract ixyz_frag (fragment origin in global coordinates, 1-based)
      ! jxyz_tot(1,:) gives the global index of the first grid point in this fragment
      dg_frag%ixyz_frag(1:3, ifrag) = jxyz_tot(1, 1:3)

      ! DG-RT requires the DC-exported buffer-aware basis.  The core-only
      ! basis_functions.bin cannot provide fragment-boundary stencil data.
      iunit = get_filehandle()
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), ifrag, '/', binfile_bfb
      inquire(file=filename, exist=has_buffer_file)

      if (.not. has_buffer_file) then
        write(*,'(1x,a,i0,a,a)') "[FATAL] missing DG buffer basis at ifrag=", ifrag, &
          " file=", trim(filename)
        write(*,'(1x,a)') "Regenerate the DC seed so basis_functions_buffer.bin is exported."
        stop "DG-Fragment RT: missing basis buffer file"
      end if
      open(iunit, file=filename, form='unformatted', access='stream', status='old')
      read(iunit) magic_file, version_file
      if (magic_file /= basis_buffer_magic .or. version_file /= basis_buffer_version) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "Error: invalid basis buffer header at ifrag=", ifrag, &
          " magic=", magic_file, " version=", version_file
        stop "DG-Fragment RT: invalid basis buffer file"
      end if
      read(iunit) nxyz_domain(1:3), nxyz_buffer_file(1:3), nspin_file, nstate_frag_file
      if (nspin_file < 1) then
        write(*,'(1x,a,i0,a,i0)') "Error: invalid nspin_file in basis buffer header at ifrag=", ifrag, &
                                   " nspin_file=", nspin_file
        stop "DG-Fragment RT: invalid nspin_file"
      end if
      do n = 1, 3
        if (dg_frag%num_fragment(n) > 1 .and. nxyz_buffer_file(n) < nb) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] DG seed buffer too small at ifrag=", ifrag, &
            " axis=", n, " seed_buffer=", nxyz_buffer_file(n), " required=", nb
          stop "DG-Fragment RT: insufficient basis buffer"
        end if
      end do
      nxyz_box(1:3) = nxyz_domain(1:3) + 2 * nxyz_buffer_file(1:3)
      if (allocated(n_basis_frag)) deallocate(n_basis_frag)
      allocate(n_basis_frag(nspin_file))
      read(iunit) n_basis_frag(1:nspin_file)

      ! Store domain size for this fragment
      dg_frag%nxyz_domain(1:3, ifrag) = nxyz_domain(1:3)

      ! Orbital fragment parallelism must not partition the fragment basis in
      ! real space.  Each subgroup rank keeps the same full fragment box, then
      ! matrix construction is split over basis rows/columns.
      if (.not. allocated(dg_frag%phi_frag)) then
        if (dg_frag%parallel_mode_orbital) then
          allocate(dg_frag%phi_frag(dg_frag%ixyz_frag(1, ifrag)-nb:dg_frag%ixyz_frag(1, ifrag)+nxyz_domain(1)-1+nb, &
                                     dg_frag%ixyz_frag(2, ifrag)-nb:dg_frag%ixyz_frag(2, ifrag)+nxyz_domain(2)-1+nb, &
                                     dg_frag%ixyz_frag(3, ifrag)-nb:dg_frag%ixyz_frag(3, ifrag)+nxyz_domain(3)-1+nb, &
                                     dg_frag%nstate_frag, ifrag_count))
        else
          allocate(dg_frag%phi_frag(dg_frag%mg%is(1)-nb:dg_frag%mg%ie(1)+nb, &
                                     dg_frag%mg%is(2)-nb:dg_frag%mg%ie(2)+nb, &
                                     dg_frag%mg%is(3)-nb:dg_frag%mg%ie(3)+nb, &
                                     dg_frag%nstate_frag, ifrag_count))
        end if
        dg_frag%phi_frag = 0.0d0  ! Initialize (including halo) to zero
      end if

      ! Read basis functions: f_basis(ix,iy,iz,ispin,istate)
      ! phi_frag has no spin dimension: keep spin-1 basis and discard extra spin channels
      ! while still consuming all records to keep stream alignment.
      block
        real(8), allocatable :: phi_box(:,:,:)
        if (nspin_file < 1 .or. nstate_frag_file < 1) then
          write(*,'(1x,a,i0,a,i0,a,i0)') "Error: invalid basis buffer header at ifrag=", ifrag, &
                                         " nspin_file=", nspin_file, " nstate_frag_file=", nstate_frag_file
          stop "DG-Fragment RT: invalid basis buffer header"
        end if
        allocate(phi_box(nxyz_box(1), nxyz_box(2), nxyz_box(3)))

        do ispin = 1, nspin_file
          do n = 1, nstate_frag_file
            read(iunit) phi_box(1:nxyz_box(1), 1:nxyz_box(2), 1:nxyz_box(3))

            if (ispin == 1 .and. n <= dg_frag%nstate_frag) then
              ! The buffer file is stored as an unwrapped box around the core.
              ! Unsplitted axes may still wrap within the fragment itself; split
              ! axes must be covered by the seed buffer checked above.
              do izg_store = lbound(dg_frag%phi_frag, 3), ubound(dg_frag%phi_frag, 3)
                iz_box = izg_store - (dg_frag%ixyz_frag(3, ifrag) - nxyz_buffer_file(3)) + 1
                if (iz_box < 1 .or. iz_box > nxyz_box(3)) then
                  iz_src = nxyz_buffer_file(3) + modulo(izg_store - dg_frag%ixyz_frag(3, ifrag), nxyz_domain(3)) + 1
                else
                  iz_src = iz_box
                end if
                do iyg_store = lbound(dg_frag%phi_frag, 2), ubound(dg_frag%phi_frag, 2)
                  iy_box = iyg_store - (dg_frag%ixyz_frag(2, ifrag) - nxyz_buffer_file(2)) + 1
                  if (iy_box < 1 .or. iy_box > nxyz_box(2)) then
                    iy_src = nxyz_buffer_file(2) + modulo(iyg_store - dg_frag%ixyz_frag(2, ifrag), nxyz_domain(2)) + 1
                  else
                    iy_src = iy_box
                  end if
                  do ixg_store = lbound(dg_frag%phi_frag, 1), ubound(dg_frag%phi_frag, 1)
                    ix_box = ixg_store - (dg_frag%ixyz_frag(1, ifrag) - nxyz_buffer_file(1)) + 1
                    if (ix_box < 1 .or. ix_box > nxyz_box(1)) then
                      ix_src = nxyz_buffer_file(1) + modulo(ixg_store - dg_frag%ixyz_frag(1, ifrag), nxyz_domain(1)) + 1
                    else
                      ix_src = ix_box
                    end if
                    dg_frag%phi_frag(ixg_store, iyg_store, izg_store, n, i_local) = &
                      phi_box(ix_src, iy_src, iz_src)
                  end do
                end do
              end do
              dg_frag%phi_frag_has_seed_buffer(i_local) = .true.
            end if
          end do
        end do

        if (nspin_file > 1 .and. .not. warned_spin_discard .and. comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,i0,a)') "[WARN] basis_functions_buffer.bin has nspin_file=", nspin_file, &
                                 "; using spin-1 basis only in phi_frag"
          warned_spin_discard = .true.
        end if

        if (allocated(phi_box)) deallocate(phi_box)
      end block

      close(iunit)
    end do

    ! Clean up temporary arrays
    if (allocated(jxyz_tot)) deallocate(jxyz_tot)
    if (allocated(n_basis_frag)) deallocate(n_basis_frag)

    ! Share fragment geometry metadata across all ranks.  The fragment root
    ! ranks are fixed by the stage-1 subgroup layout, so use the same
    ! deterministic mapping that initializes id_array before owner-map builds.
    block
      integer :: owner_rank
      do ifrag = 1, dg_frag%n_frag
        owner_rank = get_fragment_group_root_rank(ifrag, dg_frag%nproc_frag)
        call comm_bcast(dg_frag%ixyz_frag(1:3, ifrag), dg_frag%icomm, owner_rank)
        call comm_bcast(dg_frag%nxyz_domain(1:3, ifrag), dg_frag%icomm, owner_rank)
      end do
    end block

    ! Set flag indicating real-space basis functions are available
    dg_frag%has_real_space_basis = .true.

    ! Synchronize all ranks before proceeding
    call comm_sync_all(dg_frag%icomm)

    ! Clean up
    deallocate(n_basis_tmp, index_basis_file)

    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a)') "Fragment basis data loaded (coefficients + real-space basis)"
      write(*,'(1x,a,i0,a,i0,a,i0)') "  Domain size: ", nxyz_domain(1), " x ", nxyz_domain(2), " x ", nxyz_domain(3)
      write(*,'(1x,a,i0)') "  Number of basis functions per fragment: ", dg_frag%nstate_frag
      write(*,'(1x,a,i0)') "  Number of fragments loaded: ", ifrag_count
    end if

  end subroutine read_fragment_basis_data

end module rt_dg_fragment_io
