! DG fragment basis-file I/O for SOI RT.
#include "config.h"
module rt_dg_fragment_io_soi
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use rt_dg_fragment_coefficients, only: rebuild_coef_owner_map, get_fragment_coef_owner_rank
  use rt_dg_fragment_layout, only: get_fragment_group_root_rank
  implicit none
  private
  public :: read_fragment_basis_data

contains
  subroutine read_fragment_basis_data(dg_frag, bdir_frag)
    use filesystem, only: get_filehandle
    use communication, only: comm_is_root, comm_bcast, comm_sync_all, comm_summation
    use salmon_global, only: dg_nmat_cap_mode, dg_nmat_cap_fixed, dg_nmat_cap_multiple, nelec, nelec_spin, &
                             yn_adaptive_basis
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    character(*), intent(in) :: bdir_frag

    character(32), parameter :: binfile_wf = "wavefunctions_soi.bin"
    character(32), parameter :: binfile_bf = "basis_functions_soi.bin"
    character(32), parameter :: binfile_bfb = "basis_functions_buffer_soi.bin"
    character(32), parameter :: binfile_rg = "rgrid_index_soi.bin"
    integer, parameter :: basis_buffer_magic = -22022213
    integer, parameter :: basis_buffer_version = 1
    character(256) :: filename
    integer :: iunit, ifrag, ispin, n_frag_file, nspin_file
    integer :: nstate_frag_file, nstate_tot_file, nstate_frag_keep
    integer, allocatable :: n_basis_tmp(:,:), index_basis_file(:,:,:)
    complex(8), allocatable :: coef_state(:)
    integer :: n_mat_tmp(2)   ! nspin is expected to be 1 or 2
    integer :: ifrag_count, i_local, io, global_idx
    integer :: istate, ispin_file
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
    integer :: env_status, env_len
    character(len=64) :: env_n_mat_cap
    logical :: warned_spin_discard, warned_imag_discard, has_buffer_file
    real(8) :: cap_avg, weight_best
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
            "' requires dg_nmat_cap_fixed >= 1, got ", dg_nmat_cap_fixed
        end if
        stop "DG-Fragment RT: invalid fragment basis cap"
      end if
      nstate_frag_keep = min(nstate_frag_file, dg_nmat_cap_fixed)
      if (dg_frag%id == 0) then
        write(*,'(1x,a,a,a,i0,a,i0,a,i0)') "[INFO] DG fragment cap mode='", trim(dg_nmat_cap_mode), &
          "' file_nstate_frag=", nstate_frag_file, " keep=", nstate_frag_keep, &
          " requested=", dg_nmat_cap_fixed
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
    read(iunit) n_mat_tmp(1:dg_frag%nspin)
    read(iunit) n_basis_tmp(1:dg_frag%n_frag, 1:dg_frag%nspin)
    read(iunit) index_basis_file(1:nstate_frag_file, 1:dg_frag%n_frag, 1:dg_frag%nspin)
    close(iunit)

    ! Step 3: Gather metadata (now consistent across all ranks)
    dg_frag%n_basis = min(n_basis_tmp, dg_frag%nstate_frag)
    dg_frag%index_basis = 0
    dg_frag%index_basis(1:dg_frag%nstate_frag, 1:dg_frag%n_frag, 1:dg_frag%nspin) = &
      index_basis_file(1:dg_frag%nstate_frag, 1:dg_frag%n_frag, 1:dg_frag%nspin)
    dg_frag%n_mat(1:dg_frag%nspin) = n_mat_tmp(1:dg_frag%nspin)
    dg_frag%n_mat_max = maxval(n_mat_tmp(1:dg_frag%nspin))

    ! Optional basis-size cap for fragment-comparison studies.
    ! Normal production runs keep the full fragment basis unless
    ! SALMON_DG_NMAT_CAP is explicitly set.
    n_mat_cap = 0
    if (trim(dg_nmat_cap_mode) == 'fixed' .and. dg_nmat_cap_fixed >= 1) then
      n_mat_cap = dg_nmat_cap_fixed
      if (dg_frag%id == 0) then
        write(*,'(1x,a,a,a,i0)') "[INFO] DG fragment cap mode='", trim(dg_nmat_cap_mode), "' fixed=", n_mat_cap
      end if
    end if
    env_n_mat_cap = ""
    env_status = 1
    env_len = 0
    call get_environment_variable("SALMON_DG_NMAT_CAP", env_n_mat_cap, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      read(env_n_mat_cap(1:env_len), *, iostat=ienv) n_mat_cap_env
      if (ienv == 0 .and. n_mat_cap_env >= 1) then
        n_mat_cap = n_mat_cap_env
        if (dg_frag%id == 0) then
          write(*,'(1x,a,i0)') "[INFO] SALMON_DG_NMAT_CAP override applied: ", n_mat_cap
        end if
      else
        if (dg_frag%id == 0) then
          write(*,'(1x,a,a,a)') "[WARN] Ignoring invalid SALMON_DG_NMAT_CAP='", &
                                trim(env_n_mat_cap(1:env_len)), "' (must be integer >= 1)"
        end if
      end if
    end if
    if (n_mat_cap >= 1) then
      dg_frag%n_mat(1:dg_frag%nspin) = min(dg_frag%n_mat(1:dg_frag%nspin), n_mat_cap)
      dg_frag%n_mat_max = maxval(dg_frag%n_mat(1:dg_frag%nspin))
    end if

    ! Enforce cap consistency in index_basis:
    ! indices beyond capped n_mat are masked out to prevent OOB accesses.
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
    ! Keeping the DC-file numbering after a fragment-local cap leaves large
    ! holes between fragment blocks and inflates O(n_mat_max^2) work arrays.
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

    if (allocated(dg_frag%coef_owner)) deallocate(dg_frag%coef_owner)
    allocate(dg_frag%coef_owner(dg_frag%n_mat_max, dg_frag%nspin))
    dg_frag%coef_owner(:, :) = -1
    dg_frag%H_local_block_ids_valid = .false.
    do ispin = 1, dg_frag%nspin
      do ifrag = 1, dg_frag%n_frag
        nbasis_iter = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1))
        do io = 1, nbasis_iter
          global_idx = dg_frag%index_basis(io, ifrag, ispin)
          if (global_idx < 1 .or. global_idx > dg_frag%n_mat_max) cycle
          dg_frag%coef_owner(global_idx, ispin) = get_fragment_coef_owner_rank(dg_frag, ifrag, io, nbasis_iter)
        end do
      end do
    end do
    dg_frag%owned_coef_start = 0
    dg_frag%owned_coef_end = -1
    do global_idx = 1, dg_frag%n_mat_max
      if (any(dg_frag%coef_owner(global_idx, 1:dg_frag%nspin) == dg_frag%id)) then
        if (dg_frag%owned_coef_start == 0) dg_frag%owned_coef_start = global_idx
        dg_frag%owned_coef_end = global_idx
      end if
    end do

    ! Step 4: nstate_tot was aligned to file metadata above (full-state mode).

    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1

    if (n_mat_cap < 1 .and. trim(dg_nmat_cap_mode) == 'occ_multiple' .and. dg_nmat_cap_multiple > 0.0d0) then
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

      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        iunit = get_filehandle()
        write(filename, '(a, i6.6, a, a)') trim(bdir_frag), ifrag, '/', binfile_wf
        open(iunit, file=filename, form='unformatted', access='stream', status='old')
        read(iunit) n_frag_file, nspin_file, nstate_frag_file, nstate_tot_file
        read(iunit) n_mat_tmp(1:dg_frag%nspin)
        read(iunit) n_basis_tmp(1:dg_frag%n_frag, 1:dg_frag%nspin)
        read(iunit) index_basis_file(1:nstate_frag_file, 1:dg_frag%n_frag, 1:dg_frag%nspin)
        if (allocated(coef_state)) deallocate(coef_state)
        allocate(coef_state(nstate_frag_file))

        do ispin_file = 1, nspin_file
          do istate = 1, nstate_tot_file
            read(iunit) coef_state(1:nstate_frag_file)
            if (ispin_file < 1 .or. ispin_file > dg_frag%nspin) cycle
            if (dg_frag%nspin == 1) then
              nocc_eff = max(1, min((nelec + 1) / 2, dg_frag%nstate_tot))
            else if (sum(nelec_spin(1:dg_frag%nspin)) > 0) then
              nocc_eff = max(1, min(nelec_spin(ispin_file), dg_frag%nstate_tot))
            else
              nocc_eff = max(1, min(int(nelec / 2.0d0 + 1.0d-12), dg_frag%nstate_tot))
            end if
            if (istate > nocc_eff) cycle
            nbasis_iter = min(dg_frag%n_basis(ifrag, ispin_file), nstate_frag_file)
            frag_weight_local(ifrag, istate, ispin_file) = sum(abs(coef_state(1:nbasis_iter))**2)
          end do
        end do
        close(iunit)
        deallocate(coef_state)
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
          cap_frag(ifrag, ispin) = min(dg_frag%n_basis(ifrag, ispin), &
                                       int(floor(dg_nmat_cap_multiple * dble(occ_count(ifrag, ispin)))))
          cap_frag(ifrag, ispin) = max(1, cap_frag(ifrag, ispin))
          dg_frag%n_basis(ifrag, ispin) = cap_frag(ifrag, ispin)
          do io = cap_frag(ifrag, ispin) + 1, min(dg_frag%nstate_frag, size(dg_frag%index_basis, 1))
            dg_frag%index_basis(io, ifrag, ispin) = 0
          end do
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

      if (allocated(dg_frag%coef_owner)) deallocate(dg_frag%coef_owner)
      allocate(dg_frag%coef_owner(dg_frag%n_mat_max, dg_frag%nspin))
      dg_frag%coef_owner(:, :) = -1
      dg_frag%H_local_block_ids_valid = .false.
      do ispin = 1, dg_frag%nspin
        do ifrag = 1, dg_frag%n_frag
          nbasis_iter = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%index_basis, 1))
          do io = 1, nbasis_iter
            global_idx = dg_frag%index_basis(io, ifrag, ispin)
            if (global_idx < 1 .or. global_idx > dg_frag%n_mat_max) cycle
            dg_frag%coef_owner(global_idx, ispin) = get_fragment_coef_owner_rank(dg_frag, ifrag, io, nbasis_iter)
          end do
        end do
      end do
      dg_frag%owned_coef_start = 0
      dg_frag%owned_coef_end = -1
      do global_idx = 1, dg_frag%n_mat_max
        if (any(dg_frag%coef_owner(global_idx, 1:dg_frag%nspin) == dg_frag%id)) then
          if (dg_frag%owned_coef_start == 0) dg_frag%owned_coef_start = global_idx
          dg_frag%owned_coef_end = global_idx
        end if
      end do

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

    if (allocated(dg_frag%local_coef_count)) deallocate(dg_frag%local_coef_count)
    if (allocated(dg_frag%local_coef_global_ids)) deallocate(dg_frag%local_coef_global_ids)
    if (allocated(dg_frag%coef_global_to_local)) deallocate(dg_frag%coef_global_to_local)
    allocate(dg_frag%local_coef_count(dg_frag%nspin))
    dg_frag%local_coef_count(:) = 0
    do ispin = 1, dg_frag%nspin
      nbasis_iter = min(dg_frag%n_basis(dg_frag%ifrag_group, ispin), dg_frag%nstate_frag)
      do io = 1, nbasis_iter
        global_idx = dg_frag%index_basis(io, dg_frag%ifrag_group, ispin)
        if (global_idx < 1 .or. global_idx > dg_frag%n_mat_max) cycle
        if (dg_frag%coef_owner(global_idx, ispin) == dg_frag%id) then
          dg_frag%local_coef_count(ispin) = dg_frag%local_coef_count(ispin) + 1
        end if
      end do
    end do
    local_coef_max = max(1, maxval(dg_frag%local_coef_count(1:dg_frag%nspin)))
    allocate(dg_frag%local_coef_global_ids(local_coef_max, dg_frag%nspin))
    allocate(dg_frag%coef_global_to_local(dg_frag%n_mat_max, dg_frag%nspin))
    dg_frag%local_coef_global_ids(:, :) = 0
    dg_frag%coef_global_to_local(:, :) = 0
    do ispin = 1, dg_frag%nspin
      local_idx = 0
      nbasis_iter = min(dg_frag%n_basis(dg_frag%ifrag_group, ispin), dg_frag%nstate_frag)
      do io = 1, nbasis_iter
        global_idx = dg_frag%index_basis(io, dg_frag%ifrag_group, ispin)
        if (global_idx < 1 .or. global_idx > dg_frag%n_mat_max) cycle
        if (dg_frag%coef_owner(global_idx, ispin) /= dg_frag%id) cycle
        local_idx = local_idx + 1
        dg_frag%local_coef_global_ids(local_idx, ispin) = global_idx
        dg_frag%coef_global_to_local(global_idx, ispin) = local_idx
      end do
    end do

    ! Reallocate coefficient arrays with rank-local coefficient rows.
    if (allocated(dg_frag%coef)) deallocate(dg_frag%coef)
    if (allocated(dg_frag%coef_new)) deallocate(dg_frag%coef_new)
    if (allocated(dg_frag%coef_work)) deallocate(dg_frag%coef_work)
    allocate(dg_frag%coef(local_coef_max, dg_frag%nstate_tot, dg_frag%nspin))
    if (yn_adaptive_basis == 'y') then
      allocate(dg_frag%coef_new(local_coef_max, dg_frag%nstate_tot, dg_frag%nspin))
    end if
    allocate(dg_frag%coef_work(local_coef_max, dg_frag%nstate_tot, dg_frag%nspin))
    dg_frag%coef = 0.0d0
    if (allocated(dg_frag%coef_new)) dg_frag%coef_new = 0.0d0
    dg_frag%coef_work = 0.0d0

    ! Step 5: Reorganize coefficient data from fragment-local to rank-local rows.
    ! wavefunctions_soi.bin stores coef_wf(local_basis, state, spin); each state
    ! column is contiguous, so stream one local-basis vector at a time instead of
    ! keeping coef_local(nstate_frag,nstate_tot,nspin,ifrag_count) in memory.
    do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
      iunit = get_filehandle()
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), ifrag, '/', binfile_wf
      open(iunit, file=filename, form='unformatted', access='stream', status='old')
      read(iunit) n_frag_file, nspin_file, nstate_frag_file, nstate_tot_file
      read(iunit) n_mat_tmp(1:dg_frag%nspin)
      read(iunit) n_basis_tmp(1:dg_frag%n_frag, 1:dg_frag%nspin)
      read(iunit) index_basis_file(1:nstate_frag_file, 1:dg_frag%n_frag, 1:dg_frag%nspin)
      if (allocated(coef_state)) deallocate(coef_state)
      allocate(coef_state(nstate_frag_file))

      do ispin_file = 1, nspin_file
        do istate = 1, nstate_tot_file
          read(iunit) coef_state(1:nstate_frag_file)
          if (ispin_file < 1 .or. ispin_file > dg_frag%nspin) cycle
          if (istate > dg_frag%nstate_tot) cycle
          nbasis_iter = min(dg_frag%n_basis(ifrag, ispin_file), nstate_frag_file)
          do io = 1, nbasis_iter
            global_idx = dg_frag%index_basis(io, ifrag, ispin_file)
            local_idx = 0
            if (global_idx > 0 .and. global_idx <= dg_frag%n_mat_max) then
              local_idx = dg_frag%coef_global_to_local(global_idx, ispin_file)
            end if
            if (local_idx > 0 .and. local_idx <= size(dg_frag%coef, 1)) then
              dg_frag%coef(local_idx, istate, ispin_file) = coef_state(io)
            end if
          end do
        end do
      end do
      close(iunit)
      deallocate(coef_state)
    end do

    ! Keep coefficients only on the owning fragment ranks.

    ! Step 6: Read real-space basis functions and grid mapping
    ! This enables density reconstruction: rho(r) = sum c*_i c_j phi_i(r) phi_j(r)
    ! Also extract fragment geometry information (ixyz_frag, nxyz_domain) for halo exchange
    nb = dg_frag%nxyz_buffer(1)  ! Assume uniform buffer width (4 for 4th-order stencil)
    nxyz_alloc = 0
    warned_spin_discard = .false.
    warned_imag_discard = .false.

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

      ! Read buffer-aware basis functions.  DG-RT-SOI uses finite-difference
      ! traces at fragment faces, so the core-only basis file is insufficient.
      iunit = get_filehandle()
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), ifrag, '/', binfile_bfb
      inquire(file=filename, exist=has_buffer_file)
      if (.not. has_buffer_file) then
        write(*,'(1x,a,i0,a,a)') "[FATAL] missing DG-SOI buffer basis at ifrag=", ifrag, &
          " file=", trim(filename)
        write(*,'(1x,a)') "Regenerate the DC-SOI seed so basis_functions_buffer_soi.bin is exported."
        stop "DG-Fragment RT-SOI: missing basis buffer file"
      end if
      open(iunit, file=filename, form='unformatted', access='stream', status='old')
      read(iunit) magic_file, version_file
      if (magic_file /= basis_buffer_magic .or. version_file /= basis_buffer_version) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "Error: invalid SOI basis buffer header at ifrag=", ifrag, &
          " magic=", magic_file, " version=", version_file
        stop "DG-Fragment RT-SOI: invalid basis buffer file"
      end if
      read(iunit) nxyz_domain(1:3), nxyz_buffer_file(1:3), nspin_file, nstate_frag_file
      if (nspin_file < 1) then
        write(*,'(1x,a,i0,a,i0)') "Error: invalid nspin_file in basis buffer header at ifrag=", ifrag, &
                                   " nspin_file=", nspin_file
        stop "DG-Fragment RT: invalid nspin_file"
      end if
      do n = 1, 3
        if (dg_frag%num_fragment(n) > 1 .and. nxyz_buffer_file(n) < nb) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] DG-SOI seed buffer too small at ifrag=", ifrag, &
            " axis=", n, " seed_buffer=", nxyz_buffer_file(n), " required=", nb
          stop "DG-Fragment RT-SOI: insufficient basis buffer"
        end if
      end do
      nxyz_box(1:3) = nxyz_domain(1:3) + 2 * nxyz_buffer_file(1:3)
      if (allocated(n_basis_frag)) deallocate(n_basis_frag)
      allocate(n_basis_frag(nspin_file))
      read(iunit) n_basis_frag(1:nspin_file)

      ! Store domain size for this fragment
      dg_frag%nxyz_domain(1:3, ifrag) = nxyz_domain(1:3)

      if (.not. allocated(dg_frag%phi_frag)) then
        if (dg_frag%parallel_mode_orbital) then
          allocate(dg_frag%phi_frag(dg_frag%ixyz_frag(1, ifrag)-nb:dg_frag%ixyz_frag(1, ifrag)+nxyz_domain(1)-1+nb, &
                                     dg_frag%ixyz_frag(2, ifrag)-nb:dg_frag%ixyz_frag(2, ifrag)+nxyz_domain(2)-1+nb, &
                                     dg_frag%ixyz_frag(3, ifrag)-nb:dg_frag%ixyz_frag(3, ifrag)+nxyz_domain(3)-1+nb, &
                                     dg_frag%nstate_frag, ifrag_count))
          allocate(dg_frag%phi_frag_c(dg_frag%ixyz_frag(1, ifrag)-nb:dg_frag%ixyz_frag(1, ifrag)+nxyz_domain(1)-1+nb, &
                                      dg_frag%ixyz_frag(2, ifrag)-nb:dg_frag%ixyz_frag(2, ifrag)+nxyz_domain(2)-1+nb, &
                                      dg_frag%ixyz_frag(3, ifrag)-nb:dg_frag%ixyz_frag(3, ifrag)+nxyz_domain(3)-1+nb, &
                                      dg_frag%nstate_frag, ifrag_count))
          allocate(dg_frag%phi_frag_spinor_c(dg_frag%ixyz_frag(1, ifrag)-nb:dg_frag%ixyz_frag(1, ifrag)+nxyz_domain(1)-1+nb, &
                                             dg_frag%ixyz_frag(2, ifrag)-nb:dg_frag%ixyz_frag(2, ifrag)+nxyz_domain(2)-1+nb, &
                                             dg_frag%ixyz_frag(3, ifrag)-nb:dg_frag%ixyz_frag(3, ifrag)+nxyz_domain(3)-1+nb, &
                                             max(1, dg_frag%nspin), dg_frag%nstate_frag, ifrag_count))
        else
          allocate(dg_frag%phi_frag(dg_frag%mg%is(1)-nb:dg_frag%mg%ie(1)+nb, &
                                     dg_frag%mg%is(2)-nb:dg_frag%mg%ie(2)+nb, &
                                     dg_frag%mg%is(3)-nb:dg_frag%mg%ie(3)+nb, &
                                     dg_frag%nstate_frag, ifrag_count))
          allocate(dg_frag%phi_frag_c(dg_frag%mg%is(1)-nb:dg_frag%mg%ie(1)+nb, &
                                      dg_frag%mg%is(2)-nb:dg_frag%mg%ie(2)+nb, &
                                      dg_frag%mg%is(3)-nb:dg_frag%mg%ie(3)+nb, &
                                      dg_frag%nstate_frag, ifrag_count))
          allocate(dg_frag%phi_frag_spinor_c(dg_frag%mg%is(1)-nb:dg_frag%mg%ie(1)+nb, &
                                             dg_frag%mg%is(2)-nb:dg_frag%mg%ie(2)+nb, &
                                             dg_frag%mg%is(3)-nb:dg_frag%mg%ie(3)+nb, &
                                             max(1, dg_frag%nspin), dg_frag%nstate_frag, ifrag_count))
        end if
        dg_frag%phi_frag = 0.0d0  ! Initialize (including halo) to zero
        dg_frag%phi_frag_c = (0.0d0, 0.0d0)
        dg_frag%phi_frag_spinor_c = (0.0d0, 0.0d0)
      end if

      ! Read basis functions: f_basis(ix,iy,iz,ispin,istate)
      ! phi_frag has no spin dimension.
      ! Keep the full spinor in phi_frag_spinor_c; phi_frag_c/phi_frag retain
      ! component 1 as compatibility storage for legacy helper paths.
      block
        complex(8), allocatable :: phi_tmp(:,:,:)
        if (nspin_file < 1 .or. nstate_frag_file < 1) then
          write(*,'(1x,a,i0,a,i0,a,i0)') "Error: invalid basis buffer header at ifrag=", ifrag, &
                                         " nspin_file=", nspin_file, " nstate_frag_file=", nstate_frag_file
          stop "DG-Fragment RT-SOI: invalid basis buffer header"
        end if
        allocate(phi_tmp(nxyz_box(1), nxyz_box(2), nxyz_box(3)))

        do ispin = 1, nspin_file
          do n = 1, nstate_frag_file
            read(iunit) phi_tmp(1:nxyz_box(1), 1:nxyz_box(2), 1:nxyz_box(3))

            if (ispin <= size(dg_frag%phi_frag_spinor_c, 4) .and. n <= dg_frag%nstate_frag) then
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
                    dg_frag%phi_frag_spinor_c(ixg_store, iyg_store, izg_store, ispin, n, i_local) = &
                      phi_tmp(ix_src, iy_src, iz_src)
                    if (ispin == 1) then
                      dg_frag%phi_frag_c(ixg_store, iyg_store, izg_store, n, i_local) = phi_tmp(ix_src, iy_src, iz_src)
                      dg_frag%phi_frag(ixg_store, iyg_store, izg_store, n, i_local) = real(phi_tmp(ix_src, iy_src, iz_src))
                    end if
                  end do
                end do
              end do
            end if
          end do
        end do

        if (nspin_file > 1 .and. .not. warned_spin_discard .and. comm_is_root(dg_frag%id)) then
          write(*,'(1x,a,i0,a)') "[INFO] basis_functions_buffer_soi.bin has nspin_file=", nspin_file, &
                                 "; using full spinor storage for SOI operators"
          warned_spin_discard = .true.
        end if
        if (.not. warned_imag_discard .and. comm_is_root(dg_frag%id)) then
          write(*,'(1x,a)') "[INFO] DG-RT-SOI stores complex basis in phi_frag_c and real projection in phi_frag."
          warned_imag_discard = .true.
        end if

        deallocate(phi_tmp)
      end block

      close(iunit)
      if (.not. allocated(dg_frag%phi_frag_has_seed_buffer)) then
        allocate(dg_frag%phi_frag_has_seed_buffer(ifrag_count))
      end if
      dg_frag%phi_frag_has_seed_buffer(i_local) = .true.
    end do

    ! Clean up temporary arrays
    if (allocated(jxyz_tot)) deallocate(jxyz_tot)
    if (allocated(n_basis_frag)) deallocate(n_basis_frag)

    ! CRITICAL: Share fragment geometry metadata across all ranks for Halo initialization
    ! id_array is not initialized yet here, so reconstruct owner rank using
    ! the same block distribution rule as distribute_fragments().
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

end module rt_dg_fragment_io_soi
