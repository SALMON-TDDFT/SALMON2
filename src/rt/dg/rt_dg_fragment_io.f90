! DG fragment basis-file I/O for non-SOI RT.
#include "config.h"
module rt_dg_fragment_io
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use rt_dg_fragment_types, only: s_dg_fragment_rt, invalidate_coef_exchange_cache
  use rt_dg_fragment_coefficients, only: rebuild_coef_owner_map
  use rt_dg_fragment_layout, only: get_fragment_group_root_rank
  implicit none
  private
  public :: read_fragment_basis_data
  public :: initialize_coefficients_from_buffer_wannier_flux

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
    use salmon_global, only: nelec, nelec_spin, dg_subspace_extra_states, yn_dg_length_gauge
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    character(*), intent(in) :: bdir_frag

    character(32), parameter :: binfile_wf = "wavefunctions.bin"
    character(32), parameter :: binfile_bf = "basis_functions.bin"
    character(32), parameter :: binfile_bfb = "basis_functions_buffer.bin"
    character(32), parameter :: binfile_rg = "rgrid_index.bin"
    integer, parameter :: basis_buffer_magic = -22022212
    integer, parameter :: basis_buffer_version = 2
    character(256) :: filename
    integer :: iunit, ifrag, ispin, n_frag_file, nspin_file, iostat_read
    integer :: nstate_frag_file, nstate_tot_file, nstate_frag_keep
    integer, allocatable :: n_basis_tmp(:,:), index_basis_file(:,:,:)
    real(8), allocatable :: coef_state_file(:), esp_file(:,:)
    integer :: n_mat_tmp(2)   ! nspin is expected to be 1 or 2
    integer :: ifrag_count, i_local, io, global_idx
    integer :: local_coef_max, local_idx
    integer :: nxyz_domain(3), nxyz_alloc(3), lgnum_frag(3), lgnum_total(3)
    integer :: nxyz_buffer_file(3), nxyz_box(3)
    integer :: magic_file, version_file
    integer, allocatable :: n_basis_frag(:)
    integer, allocatable :: jxyz_tot(:,:)
    integer :: ix, iy, iz, n, iw
    integer :: ixg_store, iyg_store, izg_store
    integer :: ix_src, iy_src, iz_src
    integer :: ix_box, iy_box, iz_box
    integer :: nb  ! halo width
    integer :: nbasis_iter
    integer :: nocc_eff
    integer :: nocc_est, nocc_frag_est, n_basis_min, n_basis_max, n_basis_effective_min
    integer :: nvirt_est_min, min_virtual_required
    real(8) :: n_basis_avg, nvirt_est_avg
    integer :: state_col, occ_base, occ_extra, frag_occ_s, frag_occ_e, nocc_frag_seed
    logical :: warned_spin_discard, has_buffer_file, has_core_file, identity_seed_coefficients
    logical :: has_global_wannier_file
    logical :: esp_loaded
    real(8) :: coef_diag_local(2), coef_diag_global(2)

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
    if (identity_seed_coefficients) then
      if (dg_frag%id == 0) then
        write(*,'(1x,a)') "[FATAL] DG identity/local seed export is no longer supported."
        write(*,'(1x,a)') "[FATAL] Use DC-LCFO flux ground-state export with dense LCFO coefficients."
      end if
      stop "DG-Fragment RT: unsupported identity seed"
    end if
    dg_frag%identity_seed_coefficients = identity_seed_coefficients
    dg_frag%fragment_basis_contracted = .false.
    dg_frag%dc_lcfo_seed_basis_cleaned = .false.
    if (dg_frag%id == 0) then
      if (identity_seed_coefficients) then
        write(*,'(1x,a)') "[WARN] DG identity seed detected: initial occupied states are fragment-local."
        write(*,'(1x,a)') "       This is the non-LCFO DG seed export, not a DGDFT/LCFO ground-state wavefunction."
      else
        write(*,'(1x,a)') "[INFO] DG LCFO coefficient seed detected: loading DC eigenvector coefficients."
      end if
    end if

    nstate_frag_keep = nstate_frag_file
    if (nstate_frag_file /= dg_frag%nstate_frag) then
      if (dg_frag%id == 0) then
        write(*,'(1x,a,i0,a,i0,a)') "[INFO] nstate_frag differs: file=", nstate_frag_file, &
          " runtime=", dg_frag%nstate_frag, " (using fragment-state count from file)"
      end if
    end if
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
    call read_wannier_cluster_partition(dg_frag)

    call invalidate_coef_exchange_cache(dg_frag)

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
    if (.not. identity_seed_coefficients) then
      if (allocated(dg_frag%esp)) deallocate(dg_frag%esp)
      allocate(dg_frag%esp(dg_frag%nstate_tot, dg_frag%nspin))
      dg_frag%esp(:, :) = 0.0d0
      allocate(esp_file(nstate_tot_file, nspin_file))
    end if

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

    ! The global row numbering is part of the wavefunctions.bin contract:
    ! coef_wf(local_basis,state,spin), index_basis, and all DG operator blocks
    ! must refer to the same DC/EigenExa basis rows.  Sparse row ids are allowed
    ! because the DC basis keeps only the active n_basis rows in each fragment;
    ! reject only duplicate or out-of-range ids.
    block
      integer :: ispin_chk, ifrag_chk, io_chk, row_id
      integer :: dup_count, out_count, miss_count
      logical :: bad_index_basis
      integer, allocatable :: seen(:)
      bad_index_basis = .false.
      do ispin_chk = 1, dg_frag%nspin
        allocate(seen(max(1, dg_frag%n_mat(ispin_chk))))
        seen = 0
        dup_count = 0
        out_count = 0
        do ifrag_chk = 1, dg_frag%n_frag
          nbasis_iter = min(dg_frag%n_basis(ifrag_chk, ispin_chk), size(dg_frag%index_basis, 1))
          do io_chk = 1, nbasis_iter
            row_id = dg_frag%index_basis(io_chk, ifrag_chk, ispin_chk)
            if (row_id < 1 .or. row_id > dg_frag%n_mat(ispin_chk)) then
              out_count = out_count + 1
            else
              if (seen(row_id) == 1) dup_count = dup_count + 1
              seen(row_id) = 1
            end if
          end do
        end do
        miss_count = count(seen == 0)
        if (dup_count > 0 .or. out_count > 0) then
          bad_index_basis = .true.
          if (dg_frag%id == 0) then
            write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] invalid DG index_basis in wavefunctions.bin (ispin=", &
              ispin_chk, "): dup=", dup_count, " out_of_range=", out_count, " missing=", miss_count
          end if
        else if (miss_count > 0 .and. dg_frag%id == 0) then
          write(*,'(1x,a,i0,a,i0)') &
            "[WARN] sparse DG index_basis in wavefunctions.bin: ispin=", ispin_chk, &
            " missing_rows=", miss_count
        end if
        deallocate(seen)
      end do
      if (bad_index_basis) then
        if (dg_frag%id == 0) then
          write(*,'(1x,a)') "[FATAL] Regenerate the DC-LCFO/EigenExa seed; DG-RT no longer remaps coefficient rows."
        end if
        stop "DG-Fragment RT: invalid wavefunction index_basis"
      end if
    end block
    dg_frag%n_mat_max = max(1, maxval(dg_frag%n_mat(1:dg_frag%nspin)))
    if (yn_dg_length_gauge == 'y') then
      call read_buffer_periodic_wannier_basis_data(dg_frag, bdir_frag)
      inquire(file='./data_dcdft/total/wannier90_global_basis.bin', exist=has_global_wannier_file)
      if (.not. dg_frag%has_buffer_periodic_wannier_basis .and. .not. has_global_wannier_file) &
        call read_local_wannier_basis_data(dg_frag, bdir_frag)
      call read_wannier90_global_basis_metadata_if_requested(dg_frag)
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
            write(*,'(1x,a)') "[FATAL] Increase the DC-LCFO exported fragment basis size before DG-Fragment RT."
          end if
          stop "DG-Fragment RT: insufficient unoccupied fragment states in DGDFT seed"
        end if
      end do
    end if

    ! Reallocate coefficient arrays with correct n_mat_max dimension
    if (allocated(dg_frag%coef)) deallocate(dg_frag%coef)
    if (allocated(dg_frag%coef_new)) deallocate(dg_frag%coef_new)
    if (allocated(dg_frag%coef_work)) deallocate(dg_frag%coef_work)
    if (dg_frag%coef_state_block_mode) then
      allocate(dg_frag%coef(local_coef_max, max(1, dg_frag%coef_nstate_local), dg_frag%nspin))
    else
      allocate(dg_frag%coef(local_coef_max, dg_frag%nstate_tot, dg_frag%nspin))
    end if
    if (identity_seed_coefficients) then
      if (dg_frag%coef_state_block_mode) then
        allocate(dg_frag%coef_work(local_coef_max, max(1, dg_frag%coef_nstate_local), dg_frag%nspin))
      else
        allocate(dg_frag%coef_work(local_coef_max, dg_frag%nstate_tot, dg_frag%nspin))
      end if
      if (dg_frag%yn_adaptive_basis) then
        if (dg_frag%coef_state_block_mode) then
          allocate(dg_frag%coef_new(local_coef_max, max(1, dg_frag%coef_nstate_local), dg_frag%nspin))
        else
          allocate(dg_frag%coef_new(local_coef_max, dg_frag%nstate_tot, dg_frag%nspin))
        end if
      end if
    end if
    dg_frag%coef = 0.0d0
    if (allocated(dg_frag%coef_new)) dg_frag%coef_new = 0.0d0
    if (allocated(dg_frag%coef_work)) dg_frag%coef_work = 0.0d0
    if (yn_dg_length_gauge == 'y') call read_wannier90_global_basis_data(dg_frag, ifrag_count)

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
                if (dg_frag%coef_state_block_mode) then
                  if (state_col >= dg_frag%coef_state_start .and. state_col <= dg_frag%coef_state_end) then
                    dg_frag%coef(local_idx, state_col - dg_frag%coef_state_start + 1, ispin) = (1.0d0, 0.0d0)
                  end if
                else
                  dg_frag%coef(local_idx, state_col, ispin) = (1.0d0, 0.0d0)
                end if
              end if
            end if
          end do
        end do
      end do
    else
      ! Step 5b: DC-LCFO coefficient seed.  wavefunctions.bin stores
      ! coef_wf(local_basis, state, spin), with each state column contiguous.
      ! Stream one column at a time so the RT rank keeps only its owned rows.
      esp_loaded = .false.
      do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end
        i_local = ifrag - dg_frag%ifrag_start + 1
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
            nbasis_iter = min(dg_frag%n_basis(ifrag, ispin), nstate_frag_file)
            if (any(coef_state_file(1:nbasis_iter) /= coef_state_file(1:nbasis_iter)) .or. &
                any(abs(coef_state_file(1:nbasis_iter)) > huge(1.0d0))) then
              write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') &
                "[FATAL] non-finite coefficient in wavefunctions.bin: rank=", dg_frag%id, &
                " ifrag=", ifrag, " ispin=", ispin, " state=", state_col
              stop "DG-Fragment RT: non-finite wavefunction coefficient"
            end if
            if (ispin < 1 .or. ispin > dg_frag%nspin) cycle
            if (state_col > dg_frag%nstate_tot) cycle
            if (dg_frag%coef_state_block_mode) then
              if (state_col < dg_frag%coef_state_start .or. state_col > dg_frag%coef_state_end) cycle
            end if
            do io = 1, nbasis_iter
              global_idx = dg_frag%index_basis(io, ifrag, ispin)
              local_idx = 0
              if (global_idx > 0 .and. global_idx <= dg_frag%n_mat_max) local_idx = dg_frag%coef_global_to_local(global_idx, ispin)
              if (local_idx > 0 .and. local_idx <= size(dg_frag%coef, 1)) then
                if (dg_frag%coef_state_block_mode) then
                  dg_frag%coef(local_idx, state_col - dg_frag%coef_state_start + 1, ispin) = &
                    dcmplx(coef_state_file(io), 0.0d0)
                else
                  dg_frag%coef(local_idx, state_col, ispin) = dcmplx(coef_state_file(io), 0.0d0)
                end if
              end if
            end do
            if (dg_frag%has_global_wannier_basis .and. allocated(dg_frag%global_wannier_coef)) then
              if (state_col <= dg_frag%global_wannier_num_bands) then
                do iw = 1, dg_frag%global_wannier_num_wann
                  dg_frag%global_wannier_coef(1:nbasis_iter, iw, ispin, i_local) = &
                    dg_frag%global_wannier_coef(1:nbasis_iter, iw, ispin, i_local) + &
                    dcmplx(coef_state_file(1:nbasis_iter), 0.0d0) * &
                    dg_frag%global_wannier_transform(state_col, iw)
                end do
              end if
            end if
          end do
        end do
        if (allocated(esp_file) .and. .not. esp_loaded) then
          read(iunit, iostat=iostat_read) esp_file(1:nstate_tot_file, 1:nspin_file)
          if (iostat_read /= 0) then
            if (dg_frag%id == 0) then
              write(*,'(1x,a)') "[FATAL] DC-LCFO wavefunctions.bin has no EigenExa eigenvalue block."
              write(*,'(1x,a)') "[FATAL] Regenerate the DC-LCFO-Flux seed so RT can remove static eigenphases."
            end if
            stop "DG-Fragment RT: missing DC-LCFO eigenvalues"
          end if
          dg_frag%esp(1:min(dg_frag%nstate_tot,nstate_tot_file), 1:min(dg_frag%nspin,nspin_file)) = &
            esp_file(1:min(dg_frag%nstate_tot,nstate_tot_file), 1:min(dg_frag%nspin,nspin_file))
          esp_loaded = .true.
        end if
        close(iunit)
        deallocate(coef_state_file)
      end do
    end if

    if (yn_dg_length_gauge == 'y' .and. dg_frag%has_global_wannier_basis) then
      call apply_wannier_flux_eigen_seed_if_available(dg_frag, ifrag_count)
      call build_global_wannier_local_basis(dg_frag)
      call build_formal_dg_wannier_from_global_local(dg_frag)
    end if

    if (yn_dg_length_gauge == 'y' .and. dg_frag%has_buffer_periodic_wannier_basis) then
      call initialize_coefficients_from_buffer_wannier_flux(dg_frag)
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
    if (allocated(esp_file)) deallocate(esp_file)

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
        write(*,'(1x,a,i0,4(a,i0))') "Error: invalid basis buffer header at ifrag=", ifrag, &
          " magic=", magic_file, " expected_magic=", basis_buffer_magic, &
          " version=", version_file, " expected_version=", basis_buffer_version
        write(*,'(1x,a,i0,a)') "[FATAL] DG-Fragment RT requires basis_functions_buffer.bin version ", &
          basis_buffer_version, "."
        write(*,'(1x,a)') "[FATAL] Regenerate the DC-LCFO seed with the core-S-cleaned DG export path."
        stop "DG-Fragment RT: invalid basis buffer file"
      end if
      read(iunit) nxyz_domain(1:3), nxyz_buffer_file(1:3), nspin_file, nstate_frag_file
      if (nspin_file < 1) then
        write(*,'(1x,a,i0,a,i0)') "Error: invalid nspin_file in basis buffer header at ifrag=", ifrag, &
                                   " nspin_file=", nspin_file
        stop "DG-Fragment RT: invalid nspin_file"
      end if
      if (nspin_file /= dg_frag%nspin .or. nstate_frag_file /= dg_frag%nstate_frag) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') "[FATAL] DG basis buffer metadata mismatch at ifrag=", ifrag, &
          " nspin_file=", nspin_file, " expected=", dg_frag%nspin, &
          " nstate_frag_file=", nstate_frag_file, " expected=", dg_frag%nstate_frag
        stop "DG-Fragment RT: basis buffer metadata mismatch"
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
      if (any(n_basis_frag(1:dg_frag%nspin) /= dg_frag%n_basis(ifrag, 1:dg_frag%nspin))) then
        write(*,'(1x,a,i0)') "[FATAL] DG wavefunctions/basis buffer n_basis mismatch at ifrag=", ifrag
        write(*,'(1x,a,20(1x,i0))') "        wavefunctions:", dg_frag%n_basis(ifrag, 1:dg_frag%nspin)
        write(*,'(1x,a,20(1x,i0))') "        basis_buffer:", n_basis_frag(1:dg_frag%nspin)
        stop "DG-Fragment RT: inconsistent DC-LCFO seed files"
      end if

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

      ! The DC-LCFO eigenvector coefficients are defined in the core basis
      ! written to basis_functions.bin.  Keep the buffered file for halo/stencil
      ! values, but overwrite the core region from basis_functions.bin so the
      ! RT overlap/Hamiltonian integrals use exactly the DC diagonalization
      ! basis.
      iunit = get_filehandle()
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), ifrag, '/', binfile_bf
      inquire(file=filename, exist=has_core_file)
      if (.not. has_core_file) then
        write(*,'(1x,a,i0,a,a)') "[FATAL] missing DG core basis at ifrag=", ifrag, &
          " file=", trim(filename)
        write(*,'(1x,a)') "Regenerate the DC-LCFO seed so basis_functions.bin is exported."
        stop "DG-Fragment RT: missing core basis file"
      end if
      open(iunit, file=filename, form='unformatted', access='stream', status='old')
      block
        integer :: nxyz_core(3), nspin_core, nstate_core
        integer, allocatable :: n_basis_core(:)
        real(8), allocatable :: phi_core(:,:,:,:,:)
        read(iunit) nxyz_core(1:3), nspin_core, nstate_core
        if (any(nxyz_core(1:3) /= nxyz_domain(1:3)) .or. nspin_core /= nspin_file .or. &
            nstate_core /= nstate_frag_file) then
          write(*,'(1x,a,i0)') "[FATAL] DG core/buffer basis metadata mismatch at ifrag=", ifrag
          stop "DG-Fragment RT: inconsistent core and buffer basis files"
        end if
        allocate(n_basis_core(nspin_core))
        read(iunit) n_basis_core(1:nspin_core)
        if (any(n_basis_core(1:dg_frag%nspin) /= dg_frag%n_basis(ifrag, 1:dg_frag%nspin))) then
          write(*,'(1x,a,i0)') "[FATAL] DG wavefunctions/core basis n_basis mismatch at ifrag=", ifrag
          write(*,'(1x,a,20(1x,i0))') "        wavefunctions:", dg_frag%n_basis(ifrag, 1:dg_frag%nspin)
          write(*,'(1x,a,20(1x,i0))') "        basis_core:", n_basis_core(1:dg_frag%nspin)
          stop "DG-Fragment RT: inconsistent DC-LCFO core basis file"
        end if
        allocate(phi_core(nxyz_domain(1), nxyz_domain(2), nxyz_domain(3), nspin_core, nstate_core))
        read(iunit) phi_core(1:nxyz_domain(1), 1:nxyz_domain(2), 1:nxyz_domain(3), 1:nspin_core, 1:nstate_core)
        if (nspin_core >= 1) then
          do n = 1, min(dg_frag%nstate_frag, nstate_core)
            do iz = 1, nxyz_domain(3)
              izg_store = dg_frag%ixyz_frag(3, ifrag) + iz - 1
              do iy = 1, nxyz_domain(2)
                iyg_store = dg_frag%ixyz_frag(2, ifrag) + iy - 1
                do ix = 1, nxyz_domain(1)
                  ixg_store = dg_frag%ixyz_frag(1, ifrag) + ix - 1
                  dg_frag%phi_frag(ixg_store, iyg_store, izg_store, n, i_local) = phi_core(ix, iy, iz, 1, n)
                end do
              end do
            end do
          end do
        end if
        deallocate(phi_core, n_basis_core)
      end block
      close(iunit)
    end do

    if (.not. identity_seed_coefficients) then
      dg_frag%dc_lcfo_seed_basis_cleaned = .true.
      dg_frag%fragment_basis_contracted = .true.
    end if

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
    if (yn_dg_length_gauge == 'y' .and. dg_frag%has_local_wannier_basis) then
      call diagnose_local_wannier_tail(dg_frag)
    end if

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

  subroutine read_wannier_cluster_partition(dg_frag)
    use filesystem, only: get_filehandle
    use communication, only: comm_is_root, comm_bcast
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, parameter :: wannier_cluster_magic = -22022217
    integer, parameter :: wannier_cluster_version = 1
    integer :: iunit, iostat_open, io
    integer :: magic_file, version_file, n_frag_file, ncluster_file
    integer :: num_fragment_file(3), cluster_size_file(3), num_cluster_file(3)
    integer :: ifrag
    logical :: has_file
    character(256) :: filename

    has_file = .false.
    n_frag_file = dg_frag%n_frag
    ncluster_file = dg_frag%n_frag
    num_fragment_file(1:3) = dg_frag%num_fragment(1:3)
    cluster_size_file(1:3) = 1
    num_cluster_file(1:3) = dg_frag%num_fragment(1:3)

    if(allocated(dg_frag%wannier_cluster_id)) deallocate(dg_frag%wannier_cluster_id)
    if(allocated(dg_frag%wannier_cluster_range)) deallocate(dg_frag%wannier_cluster_range)

    if(comm_is_root(dg_frag%id)) then
      filename = './data_dcdft/total/wannier_cluster_partition.bin'
      iunit = get_filehandle()
      open(iunit,file=filename,form='unformatted',access='stream',status='old',iostat=iostat_open)
      if(iostat_open == 0) then
        read(iunit,iostat=io) magic_file, version_file
        if(io == 0 .and. magic_file == wannier_cluster_magic .and. version_file == wannier_cluster_version) then
          read(iunit,iostat=io) num_fragment_file(1:3), cluster_size_file(1:3), num_cluster_file(1:3)
          if(io == 0) read(iunit,iostat=io) n_frag_file, ncluster_file
          if(io == 0 .and. n_frag_file == dg_frag%n_frag .and. &
             all(num_fragment_file(1:3) == dg_frag%num_fragment(1:3)) .and. ncluster_file > 0) then
            allocate(dg_frag%wannier_cluster_id(dg_frag%n_frag))
            allocate(dg_frag%wannier_cluster_range(6,ncluster_file))
            read(iunit,iostat=io) dg_frag%wannier_cluster_id(1:dg_frag%n_frag)
            if(io == 0) read(iunit,iostat=io) dg_frag%wannier_cluster_range(1:6,1:ncluster_file)
            if(io == 0) has_file = .true.
          end if
        end if
        close(iunit)
      end if
      if(.not. has_file) then
        if(allocated(dg_frag%wannier_cluster_id)) deallocate(dg_frag%wannier_cluster_id)
        if(allocated(dg_frag%wannier_cluster_range)) deallocate(dg_frag%wannier_cluster_range)
        n_frag_file = dg_frag%n_frag
        ncluster_file = dg_frag%n_frag
        num_fragment_file(1:3) = dg_frag%num_fragment(1:3)
        cluster_size_file(1:3) = 1
        num_cluster_file(1:3) = dg_frag%num_fragment(1:3)
      end if
    end if

    call comm_bcast(has_file, dg_frag%icomm, 0)
    call comm_bcast(cluster_size_file, dg_frag%icomm, 0)
    call comm_bcast(num_cluster_file, dg_frag%icomm, 0)
    call comm_bcast(ncluster_file, dg_frag%icomm, 0)

    dg_frag%has_wannier_cluster_partition = has_file
    dg_frag%wannier_cluster_size(1:3) = cluster_size_file(1:3)
    dg_frag%num_wannier_cluster(1:3) = num_cluster_file(1:3)
    dg_frag%n_wannier_cluster = ncluster_file

    if(has_file) then
      if(.not. allocated(dg_frag%wannier_cluster_id)) allocate(dg_frag%wannier_cluster_id(dg_frag%n_frag))
      if(.not. allocated(dg_frag%wannier_cluster_range)) allocate(dg_frag%wannier_cluster_range(6,max(1,ncluster_file)))
      call comm_bcast(dg_frag%wannier_cluster_id, dg_frag%icomm, 0)
      call comm_bcast(dg_frag%wannier_cluster_range, dg_frag%icomm, 0)
    else
      dg_frag%wannier_cluster_size(1:3) = 1
      dg_frag%num_wannier_cluster(1:3) = dg_frag%num_fragment(1:3)
      dg_frag%n_wannier_cluster = dg_frag%n_frag
      if(allocated(dg_frag%wannier_cluster_id)) deallocate(dg_frag%wannier_cluster_id)
      if(allocated(dg_frag%wannier_cluster_range)) deallocate(dg_frag%wannier_cluster_range)
      allocate(dg_frag%wannier_cluster_id(dg_frag%n_frag))
      allocate(dg_frag%wannier_cluster_range(6,dg_frag%n_frag))
      do ifrag=1,dg_frag%n_frag
        dg_frag%wannier_cluster_id(ifrag) = ifrag
        dg_frag%wannier_cluster_range(1,ifrag) = fragment_coord_from_id(ifrag, dg_frag%num_fragment, 1)
        dg_frag%wannier_cluster_range(2,ifrag) = fragment_coord_from_id(ifrag, dg_frag%num_fragment, 1)
        dg_frag%wannier_cluster_range(3,ifrag) = fragment_coord_from_id(ifrag, dg_frag%num_fragment, 2)
        dg_frag%wannier_cluster_range(4,ifrag) = fragment_coord_from_id(ifrag, dg_frag%num_fragment, 2)
        dg_frag%wannier_cluster_range(5,ifrag) = fragment_coord_from_id(ifrag, dg_frag%num_fragment, 3)
        dg_frag%wannier_cluster_range(6,ifrag) = fragment_coord_from_id(ifrag, dg_frag%num_fragment, 3)
      end do
    end if

    if(comm_is_root(dg_frag%id)) then
      if(has_file) then
        write(*,'(1x,a,3(i0,1x),a,i0)') "[DG-WANNIER-CLUSTER] loaded cluster_size=", &
          dg_frag%wannier_cluster_size(1:3), " ncluster=", dg_frag%n_wannier_cluster
      else
        write(*,'(1x,a)') "[DG-WANNIER-CLUSTER] no cluster partition file; using one fragment per Wannier cluster."
      end if
    end if
  end subroutine read_wannier_cluster_partition

  integer function fragment_coord_from_id(ifrag, num_fragment, axis) result(coord)
    implicit none
    integer, intent(in) :: ifrag, num_fragment(3), axis
    integer :: rem
    coord = 1
    select case(axis)
    case(1)
      coord = (ifrag - 1) / max(1, num_fragment(2) * num_fragment(3)) + 1
    case(2)
      rem = modulo(ifrag - 1, max(1, num_fragment(2) * num_fragment(3)))
      coord = rem / max(1, num_fragment(3)) + 1
    case(3)
      rem = modulo(ifrag - 1, max(1, num_fragment(2) * num_fragment(3)))
      coord = modulo(rem, max(1, num_fragment(3))) + 1
    end select
  end function fragment_coord_from_id

  subroutine read_local_wannier_basis_data(dg_frag, bdir_frag)
    use filesystem, only: get_filehandle
    use communication, only: comm_is_root
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    character(*), intent(in) :: bdir_frag

    character(32), parameter :: binfile_lw = "local_wannier_basis.bin"
    integer, parameter :: local_wannier_magic = -22022214
    integer, parameter :: local_wannier_version = 2
    character(256) :: filename
    integer :: ifrag, i_local, ifrag_count, iunit, iostat_open
    integer :: magic_file, version_file, nxyz_domain(3), nspin_file
    integer :: nbasis_file, nproj_file, nkeep_file, nkeep_max, ispin, iw
    integer :: owned_count, total_count
    integer, allocatable :: proj_atom(:), proj_hybrid(:), keep_index(:)
    real(8), allocatable :: lambda_w(:), wcoef(:,:), r_wann(:,:,:), wcenter(:,:)

    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    if (ifrag_count <= 0) return

    if (allocated(dg_frag%local_wannier_nbasis)) deallocate(dg_frag%local_wannier_nbasis)
    if (allocated(dg_frag%local_wannier_nproj)) deallocate(dg_frag%local_wannier_nproj)
    if (allocated(dg_frag%local_wannier_nkeep)) deallocate(dg_frag%local_wannier_nkeep)
    if (allocated(dg_frag%local_wannier_coef)) deallocate(dg_frag%local_wannier_coef)
    if (allocated(dg_frag%local_wannier_r)) deallocate(dg_frag%local_wannier_r)
    if (allocated(dg_frag%local_wannier_center)) deallocate(dg_frag%local_wannier_center)
    if (allocated(dg_frag%local_wannier_owner_fragment)) deallocate(dg_frag%local_wannier_owner_fragment)
    if (allocated(dg_frag%local_wannier_owned)) deallocate(dg_frag%local_wannier_owned)
    allocate(dg_frag%local_wannier_nbasis(ifrag_count))
    allocate(dg_frag%local_wannier_nproj(ifrag_count))
    allocate(dg_frag%local_wannier_nkeep(ifrag_count))
    dg_frag%local_wannier_nbasis = 0
    dg_frag%local_wannier_nproj = 0
    dg_frag%local_wannier_nkeep = 0

    nkeep_max = 0
    do i_local = 1, ifrag_count
      ifrag = dg_frag%ifrag_start + i_local - 1
      iunit = get_filehandle()
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), ifrag, '/', binfile_lw
      open(iunit, file=filename, form='unformatted', access='stream', status='old', iostat=iostat_open)
      if (iostat_open /= 0) then
        write(*,'(1x,a,i0,a,a)') "[FATAL] missing DG local Wannier seed at ifrag=", ifrag, &
          " file=", trim(filename)
        stop "DG length gauge: missing local Wannier seed"
      end if
      read(iunit) magic_file, version_file
      if (magic_file /= local_wannier_magic .or. version_file /= local_wannier_version) then
        write(*,'(1x,a,i0,4(a,i0))') "[FATAL] invalid local Wannier seed header at ifrag=", ifrag, &
          " magic=", magic_file, " expected_magic=", local_wannier_magic, &
          " version=", version_file, " expected_version=", local_wannier_version
        stop "DG length gauge: invalid local Wannier seed"
      end if
      read(iunit) nxyz_domain(1:3), nspin_file, nbasis_file, nproj_file, nkeep_file
      close(iunit)
      if (nspin_file /= dg_frag%nspin) then
        write(*,'(1x,a,i0,a,i0,a,i0)') "[FATAL] local Wannier nspin mismatch at ifrag=", &
          ifrag, " file=", nspin_file, " runtime=", dg_frag%nspin
        stop "DG length gauge: local Wannier nspin mismatch"
      end if
      if (nbasis_file /= dg_frag%n_basis(ifrag, 1)) then
        write(*,'(1x,a,i0,a,i0,a,i0)') "[FATAL] local Wannier basis mismatch at ifrag=", &
          ifrag, " file=", nbasis_file, " wavefunctions=", dg_frag%n_basis(ifrag, 1)
        stop "DG length gauge: local Wannier basis mismatch"
      end if
      if (nproj_file <= 0 .or. nkeep_file <= 0 .or. nkeep_file > nproj_file) then
        write(*,'(1x,a,i0,a,i0,a,i0)') "[FATAL] invalid local Wannier rank at ifrag=", &
          ifrag, " nproj=", nproj_file, " nkeep=", nkeep_file
        stop "DG length gauge: invalid local Wannier rank"
      end if
      dg_frag%local_wannier_nbasis(i_local) = nbasis_file
      dg_frag%local_wannier_nproj(i_local) = nproj_file
      dg_frag%local_wannier_nkeep(i_local) = nkeep_file
      nkeep_max = max(nkeep_max, nkeep_file)
    end do

    allocate(dg_frag%local_wannier_coef(dg_frag%nstate_frag, nkeep_max, dg_frag%nspin, ifrag_count))
    allocate(dg_frag%local_wannier_r(3, nkeep_max, nkeep_max, dg_frag%nspin, ifrag_count))
    allocate(dg_frag%local_wannier_center(3, nkeep_max, dg_frag%nspin, ifrag_count))
    allocate(dg_frag%local_wannier_owner_fragment(nkeep_max, dg_frag%nspin, ifrag_count))
    allocate(dg_frag%local_wannier_owned(nkeep_max, dg_frag%nspin, ifrag_count))
    dg_frag%local_wannier_coef = 0.0d0
    dg_frag%local_wannier_r = 0.0d0
    dg_frag%local_wannier_center = 0.0d0
    dg_frag%local_wannier_owner_fragment = 0
    dg_frag%local_wannier_owned = .false.

    owned_count = 0
    total_count = 0
    do i_local = 1, ifrag_count
      ifrag = dg_frag%ifrag_start + i_local - 1
      nbasis_file = dg_frag%local_wannier_nbasis(i_local)
      nproj_file = dg_frag%local_wannier_nproj(i_local)
      nkeep_file = dg_frag%local_wannier_nkeep(i_local)
      allocate(proj_atom(nproj_file), proj_hybrid(nproj_file), keep_index(nkeep_file))
      allocate(lambda_w(nproj_file), wcoef(nbasis_file,nkeep_file), r_wann(3,nkeep_file,nkeep_file), &
               wcenter(3,nkeep_file))
      iunit = get_filehandle()
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), ifrag, '/', binfile_lw
      open(iunit, file=filename, form='unformatted', access='stream', status='old')
      read(iunit) magic_file, version_file
      read(iunit) nxyz_domain(1:3), nspin_file, nbasis_file, nproj_file, nkeep_file
      read(iunit) proj_atom(1:nproj_file), proj_hybrid(1:nproj_file)
      read(iunit) lambda_w(1:nproj_file), keep_index(1:nkeep_file)
      read(iunit) wcoef(1:nbasis_file,1:nkeep_file)
      read(iunit) r_wann(1:3,1:nkeep_file,1:nkeep_file)
      read(iunit) wcenter(1:3,1:nkeep_file)
      close(iunit)
      do ispin = 1, dg_frag%nspin
        dg_frag%local_wannier_coef(1:nbasis_file,1:nkeep_file,ispin,i_local) = &
          wcoef(1:nbasis_file,1:nkeep_file)
        dg_frag%local_wannier_r(1:3,1:nkeep_file,1:nkeep_file,ispin,i_local) = &
          r_wann(1:3,1:nkeep_file,1:nkeep_file)
        dg_frag%local_wannier_center(1:3,1:nkeep_file,ispin,i_local) = wcenter(1:3,1:nkeep_file)
        do iw = 1, nkeep_file
          dg_frag%local_wannier_owner_fragment(iw,ispin,i_local) = &
            find_wannier_owner_fragment(dg_frag, wcenter(1:3,iw))
          dg_frag%local_wannier_owned(iw,ispin,i_local) = &
            (dg_frag%local_wannier_owner_fragment(iw,ispin,i_local) == dg_frag%ifrag_group)
          total_count = total_count + 1
          if (dg_frag%local_wannier_owned(iw,ispin,i_local)) owned_count = owned_count + 1
        end do
      end do
      deallocate(proj_atom, proj_hybrid, keep_index, lambda_w, wcoef, r_wann, wcenter)
    end do

    dg_frag%has_local_wannier_basis = .true.
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0)') "[INFO] DG local Wannier seed loaded: local_fragments=", &
        ifrag_count, " keep_min=", minval(dg_frag%local_wannier_nkeep), &
        " keep_max=", maxval(dg_frag%local_wannier_nkeep), &
        " center_owned=", owned_count, " center_total=", total_count
    end if
  contains
    integer function find_wannier_owner_fragment(dg_frag_in, center) result(owner_frag)
      type(s_dg_fragment_rt), intent(in) :: dg_frag_in
      real(8), intent(in) :: center(3)
      integer :: ifrag_scan, idx(3), axis

      owner_frag = 0
      do axis = 1, 3
        if (dg_frag_in%hgs(axis) > 0.0d0) then
          idx(axis) = modulo(nint(center(axis) / dg_frag_in%hgs(axis)), dg_frag_in%lgnum_total(axis)) + 1
        else
          idx(axis) = 1
        end if
      end do
      do ifrag_scan = 1, dg_frag_in%n_frag
        if (.not. grid_in_fragment_core(dg_frag_in, idx, ifrag_scan)) cycle
        owner_frag = ifrag_scan
        return
      end do
      owner_frag = dg_frag_in%ifrag_group
    end function find_wannier_owner_fragment

    logical function grid_in_fragment_core(dg_frag_in, idx, ifrag_scan) result(is_inside)
      type(s_dg_fragment_rt), intent(in) :: dg_frag_in
      integer, intent(in) :: idx(3), ifrag_scan
      integer :: axis

      is_inside = .false.
      if (.not. allocated(dg_frag_in%frag_core_lo)) return
      if (ifrag_scan < 1 .or. ifrag_scan > size(dg_frag_in%frag_core_lo, 2)) return
      do axis = 1, 3
        if (idx(axis) < dg_frag_in%frag_core_lo(axis, ifrag_scan)) return
        if (idx(axis) > dg_frag_in%frag_core_hi(axis, ifrag_scan)) return
      end do
      is_inside = .true.
    end function grid_in_fragment_core
  end subroutine read_local_wannier_basis_data

  subroutine read_buffer_periodic_wannier_basis_data(dg_frag, bdir_frag)
    use eigen_subdiag_sub, only: eigen_dsyev
    use filesystem, only: get_filehandle
    use communication, only: comm_is_root, comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    character(*), intent(in) :: bdir_frag

    character(48), parameter :: binfile_bpw = "buffer_periodic_wannier_basis.bin"
    integer, parameter :: buffer_periodic_wannier_magic = -22022215
    integer, parameter :: buffer_periodic_wannier_version = 3
    character(256) :: filename
    integer :: ifrag, i_local, ifrag_count, iunit, iostat_open
    integer :: magic_file, version_file, ifrag_file
    integer :: aa_wann_source
    integer :: nxyz_domain(3), nxyz_buffer_file(3), nxyz_box(3), nspin_file
    integer :: nbasis_file, nproj_file, nkeep_file, nkeep_max, ispin_store
    integer :: iw_diag, jw_diag, irow, jcol, axis_center, idx_lo, idx_hi
    integer :: env_status, env_len
    integer, allocatable :: proj_atom(:), proj_hybrid(:), keep_index(:)
    real(8), allocatable :: lambda_w(:), spread_est(:), tail_est(:)
    real(8), allocatable :: wcoef(:,:), r_wann(:,:,:), wcenter(:,:), h_wann(:,:), v_wann(:,:,:)
    real(8), allocatable :: aa_wann(:,:,:)
    real(8), allocatable :: h_diag(:,:), evec(:,:), eval(:), v_eig(:,:), r_eig(:,:), r_eff(:,:)
    real(8) :: h_diag_abs_local, h_diag_abs_global, h_offdiag_abs_local, h_offdiag_abs_global
    real(8) :: r_diag_abs_local(3), r_diag_abs_global(3)
    real(8) :: r_offdiag_abs_local(3), r_offdiag_abs_global(3)
    real(8) :: position_scale
    character(32) :: env_position_mode
    character(32) :: env_position_scale
    logical :: use_direct_r_wann
    logical :: use_berry_position
    logical :: force_v_over_gap
    logical :: berry_position_used
    logical :: aa_available_axis
    logical :: has_first_seed

    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    if (ifrag_count <= 0) return
    env_position_mode = ''
    call get_environment_variable('SALMON_DG_BPW_POSITION_MODE', env_position_mode, length=env_len, status=env_status)
    use_direct_r_wann = .false.
    use_berry_position = .false.
    force_v_over_gap = .false.
    berry_position_used = .false.
    if (env_status == 0 .and. env_len > 0) then
      select case(trim(adjustl(env_position_mode(1:env_len))))
      case('rwann','RWANN','r_wann','R_WANN','direct','DIRECT','file','FILE')
        use_direct_r_wann = .true.
        use_berry_position = .false.
      case('berry','BERRY','aa','AA','aa_r','AA_R','a_wann','A_WANN')
        use_direct_r_wann = .false.
        use_berry_position = .true.
      case('v_over_gap','V_OVER_GAP','velocity','VELOCITY','legacy','LEGACY')
        use_direct_r_wann = .false.
        use_berry_position = .false.
        force_v_over_gap = .true.
      case default
        use_direct_r_wann = .false.
        use_berry_position = .false.
        force_v_over_gap = .true.
      end select
    end if
    position_scale = 1.0d0
    env_position_scale = ''
    call get_environment_variable('SALMON_DG_BPW_POSITION_SCALE', env_position_scale, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      read(env_position_scale(1:env_len), *, iostat=env_status) position_scale
      if (env_status /= 0) position_scale = 1.0d0
    end if

    if (allocated(dg_frag%buffer_wannier_nkeep)) deallocate(dg_frag%buffer_wannier_nkeep)
    if (allocated(dg_frag%buffer_wannier_coef)) deallocate(dg_frag%buffer_wannier_coef)
    if (allocated(dg_frag%buffer_wannier_spread)) deallocate(dg_frag%buffer_wannier_spread)
    if (allocated(dg_frag%buffer_wannier_tail)) deallocate(dg_frag%buffer_wannier_tail)
    if (allocated(dg_frag%buffer_wannier_h_flux)) deallocate(dg_frag%buffer_wannier_h_flux)
    if (allocated(dg_frag%buffer_wannier_v)) deallocate(dg_frag%buffer_wannier_v)
    if (allocated(dg_frag%buffer_wannier_frag_center)) deallocate(dg_frag%buffer_wannier_frag_center)
    if (allocated(dg_frag%buffer_wannier_xi_flux_blocks)) then
      do irow = 1, size(dg_frag%buffer_wannier_xi_flux_blocks)
        if (allocated(dg_frag%buffer_wannier_xi_flux_blocks(irow)%val)) &
          deallocate(dg_frag%buffer_wannier_xi_flux_blocks(irow)%val)
      end do
      deallocate(dg_frag%buffer_wannier_xi_flux_blocks)
    end if
    if (allocated(dg_frag%buffer_wannier_xi_flux_local_block_ids)) &
      deallocate(dg_frag%buffer_wannier_xi_flux_local_block_ids)
    dg_frag%n_buffer_wannier_xi_flux_blocks = 0
    dg_frag%buffer_wannier_xi_flux_available = .false.
    dg_frag%buffer_wannier_xi_flux_warned = .false.
    dg_frag%has_buffer_periodic_wannier_basis = .false.

    ifrag = dg_frag%ifrag_start
    write(filename, '(a, i6.6, a, a)') trim(bdir_frag), ifrag, '/', binfile_bpw
    inquire(file=filename, exist=has_first_seed)
    if (.not. has_first_seed) then
      if (comm_is_root(dg_frag%id)) &
        write(*,'(1x,a)') "[INFO] DG buffer-periodic Wannier seed not found; using legacy local Wannier data."
      return
    end if

    allocate(dg_frag%buffer_wannier_nkeep(ifrag_count))
    dg_frag%buffer_wannier_nkeep = 0
    nkeep_max = 0
    do i_local = 1, ifrag_count
      ifrag = dg_frag%ifrag_start + i_local - 1
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), ifrag, '/', binfile_bpw
      iunit = get_filehandle()
      open(iunit, file=filename, form='unformatted', access='stream', status='old', iostat=iostat_open)
      if (iostat_open /= 0) then
        write(*,'(1x,a,i0,a,a)') "[FATAL] missing DG buffer-periodic Wannier seed at ifrag=", &
          ifrag, " file=", trim(filename)
        stop "DG length gauge: incomplete buffer-periodic Wannier seed"
      end if
      read(iunit) magic_file, version_file
      if (magic_file /= buffer_periodic_wannier_magic .or. &
          version_file < 1 .or. version_file > buffer_periodic_wannier_version) then
        write(*,'(1x,a,i0,4(a,i0))') "[FATAL] invalid buffer-periodic Wannier seed at ifrag=", ifrag, &
          " magic=", magic_file, " expected_magic=", buffer_periodic_wannier_magic, &
          " version=", version_file, " max_supported_version=", buffer_periodic_wannier_version
        stop "DG length gauge: invalid buffer-periodic Wannier seed"
      end if
      read(iunit) ifrag_file, nxyz_domain(1:3), nxyz_buffer_file(1:3), nxyz_box(1:3)
      read(iunit) nspin_file, nbasis_file, nproj_file, nkeep_file
      close(iunit)
      if (ifrag_file /= ifrag .or. nspin_file /= dg_frag%nspin) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') "[FATAL] buffer-periodic Wannier metadata mismatch: ifrag=", &
          ifrag, " file_ifrag=", ifrag_file, " file_nspin=", nspin_file, " runtime_nspin=", dg_frag%nspin
        stop "DG length gauge: buffer-periodic Wannier metadata mismatch"
      end if
      if (nkeep_file <= 0 .or. nkeep_file > nproj_file) then
        write(*,'(1x,a,i0,a,i0,a,i0)') "[FATAL] invalid buffer-periodic Wannier rank at ifrag=", &
          ifrag, " nproj=", nproj_file, " nkeep=", nkeep_file
        stop "DG length gauge: invalid buffer-periodic Wannier rank"
      end if
      dg_frag%buffer_wannier_nkeep(i_local) = nkeep_file
      nkeep_max = max(nkeep_max, nkeep_file)
    end do

    allocate(dg_frag%buffer_wannier_spread(nkeep_max, ifrag_count))
    allocate(dg_frag%buffer_wannier_coef(dg_frag%nstate_frag, nkeep_max, dg_frag%nspin, ifrag_count))
    allocate(dg_frag%buffer_wannier_tail(nkeep_max, ifrag_count))
    allocate(dg_frag%buffer_wannier_h_flux(nkeep_max, nkeep_max, ifrag_count))
    allocate(dg_frag%buffer_wannier_v(3, nkeep_max, nkeep_max, ifrag_count))
    allocate(dg_frag%buffer_wannier_frag_center(3, ifrag_count))
    dg_frag%buffer_wannier_spread = 0.0d0
    dg_frag%buffer_wannier_coef = 0.0d0
    dg_frag%buffer_wannier_tail = 0.0d0
    dg_frag%buffer_wannier_h_flux = 0.0d0
    dg_frag%buffer_wannier_v = 0.0d0
    dg_frag%buffer_wannier_frag_center = 0.0d0
    h_diag_abs_local = 0.0d0
    h_offdiag_abs_local = 0.0d0
    r_diag_abs_local(:) = 0.0d0
    r_offdiag_abs_local(:) = 0.0d0

    do i_local = 1, ifrag_count
      ifrag = dg_frag%ifrag_start + i_local - 1
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), ifrag, '/', binfile_bpw
      iunit = get_filehandle()
      open(iunit, file=filename, form='unformatted', access='stream', status='old')
      read(iunit) magic_file, version_file
      read(iunit) ifrag_file, nxyz_domain(1:3), nxyz_buffer_file(1:3), nxyz_box(1:3)
      read(iunit) nspin_file, nbasis_file, nproj_file, nkeep_file
      allocate(proj_atom(nproj_file), proj_hybrid(nproj_file), keep_index(nkeep_file))
      allocate(lambda_w(nproj_file), spread_est(nkeep_file), tail_est(nkeep_file))
      allocate(wcoef(nbasis_file,nkeep_file), r_wann(3,nkeep_file,nkeep_file), &
               wcenter(3,nkeep_file), h_wann(nkeep_file,nkeep_file), v_wann(3,nkeep_file,nkeep_file))
      allocate(aa_wann(3,nkeep_file,nkeep_file))
      allocate(h_diag(nkeep_file,nkeep_file), evec(nkeep_file,nkeep_file), eval(nkeep_file))
      allocate(v_eig(nkeep_file,nkeep_file), r_eig(nkeep_file,nkeep_file), r_eff(nkeep_file,nkeep_file))
      read(iunit) proj_atom(1:nproj_file), proj_hybrid(1:nproj_file)
      read(iunit) lambda_w(1:nproj_file), keep_index(1:nkeep_file)
      read(iunit) spread_est(1:nkeep_file), tail_est(1:nkeep_file)
      read(iunit) wcoef(1:nbasis_file,1:nkeep_file)
      read(iunit) r_wann(1:3,1:nkeep_file,1:nkeep_file)
      read(iunit) wcenter(1:3,1:nkeep_file)
      read(iunit) h_wann(1:nkeep_file,1:nkeep_file)
      read(iunit) v_wann(1:3,1:nkeep_file,1:nkeep_file)
      aa_wann = 0.0d0
      if (version_file >= 2) read(iunit) aa_wann(1:3,1:nkeep_file,1:nkeep_file)
      aa_wann_source = 0
      if (version_file >= 3) read(iunit) aa_wann_source
      close(iunit)

      dg_frag%buffer_wannier_spread(1:nkeep_file,i_local) = spread_est(1:nkeep_file)
      do ispin_store = 1, dg_frag%nspin
        dg_frag%buffer_wannier_coef(1:nbasis_file,1:nkeep_file,ispin_store,i_local) = &
          wcoef(1:nbasis_file,1:nkeep_file)
      end do
      dg_frag%buffer_wannier_tail(1:nkeep_file,i_local) = tail_est(1:nkeep_file)
      dg_frag%buffer_wannier_h_flux(1:nkeep_file,1:nkeep_file,i_local) = h_wann(1:nkeep_file,1:nkeep_file)
      do axis_center = 1, 3
        if (allocated(dg_frag%ixyz_frag) .and. ifrag >= 1 .and. ifrag <= size(dg_frag%ixyz_frag, 2)) then
          idx_lo = dg_frag%ixyz_frag(axis_center,ifrag)
          idx_hi = idx_lo + nxyz_domain(axis_center) - 1
          if (associated(dg_frag%lg) .and. allocated(dg_frag%lg%coordinate) .and. &
              idx_lo >= lbound(dg_frag%lg%coordinate, 1) .and. &
              idx_hi <= ubound(dg_frag%lg%coordinate, 1)) then
            dg_frag%buffer_wannier_frag_center(axis_center,i_local) = &
              0.5d0 * (dg_frag%lg%coordinate(idx_lo, axis_center) + &
                       dg_frag%lg%coordinate(idx_hi, axis_center))
          else
            dg_frag%buffer_wannier_frag_center(axis_center,i_local) = dg_frag%hgs(axis_center) * &
              (dble(idx_lo - 1) + 0.5d0 * dble(nxyz_domain(axis_center) - 1))
          end if
        else
          dg_frag%buffer_wannier_frag_center(axis_center,i_local) = 0.0d0
        end if
      end do
      do iw_diag = 1, nkeep_file
        h_diag_abs_local = h_diag_abs_local + abs(h_wann(iw_diag,iw_diag))
        do jw_diag = 1, nkeep_file
          if (iw_diag == jw_diag) cycle
          h_offdiag_abs_local = h_offdiag_abs_local + abs(h_wann(iw_diag,jw_diag))
        end do
      end do
      h_diag(1:nkeep_file,1:nkeep_file) = h_wann(1:nkeep_file,1:nkeep_file)
      call eigen_dsyev(h_diag, eval, evec)
      do ispin_store = 1, 3
        aa_available_axis = version_file >= 3 .and. aa_wann_source == 1 .and. &
          maxval(abs(aa_wann(ispin_store,1:nkeep_file,1:nkeep_file))) > 0.0d0
        if (use_direct_r_wann) then
          r_eff(1:nkeep_file,1:nkeep_file) = r_wann(ispin_store,1:nkeep_file,1:nkeep_file)
        else if (use_berry_position .and. .not. aa_available_axis) then
          write(*,'(1x,a,i0,a,i0)') &
            "[FATAL] SALMON_DG_BPW_POSITION_MODE=AA_R was requested, but no Wannier90 AA_R block is available for ifrag=", &
            ifrag, " axis=", ispin_store
          stop "DG length gauge: requested BPW AA_R position matrix is unavailable"
        else if (use_berry_position .or. &
                 (.not. force_v_over_gap .and. .not. use_direct_r_wann .and. aa_available_axis)) then
          r_eff(1:nkeep_file,1:nkeep_file) = aa_wann(ispin_store,1:nkeep_file,1:nkeep_file)
          berry_position_used = .true.
        else
          v_eig(1:nkeep_file,1:nkeep_file) = matmul(transpose(evec(1:nkeep_file,1:nkeep_file)), &
            matmul(v_wann(ispin_store,1:nkeep_file,1:nkeep_file), evec(1:nkeep_file,1:nkeep_file)))
          r_eig(1:nkeep_file,1:nkeep_file) = 0.0d0
          do irow = 1, nkeep_file
            do jcol = 1, nkeep_file
              if (irow == jcol) cycle
              if (abs(eval(jcol) - eval(irow)) <= 1.0d-10) cycle
              r_eig(irow,jcol) = v_eig(irow,jcol) / &
                (eval(jcol) - eval(irow))
            end do
          end do
          r_eff(1:nkeep_file,1:nkeep_file) = matmul(evec(1:nkeep_file,1:nkeep_file), &
            matmul(r_eig(1:nkeep_file,1:nkeep_file), transpose(evec(1:nkeep_file,1:nkeep_file))))
        end if
        dg_frag%buffer_wannier_v(ispin_store,1:nkeep_file,1:nkeep_file,i_local) = &
          position_scale * r_eff(1:nkeep_file,1:nkeep_file)
        do iw_diag = 1, nkeep_file
          r_diag_abs_local(ispin_store) = r_diag_abs_local(ispin_store) + abs(position_scale * r_eff(iw_diag,iw_diag))
          do jw_diag = 1, nkeep_file
            if (iw_diag == jw_diag) cycle
            r_offdiag_abs_local(ispin_store) = r_offdiag_abs_local(ispin_store) + &
              abs(position_scale * r_eff(iw_diag,jw_diag))
          end do
        end do
      end do
      deallocate(proj_atom, proj_hybrid, keep_index, lambda_w, spread_est, tail_est)
      deallocate(wcoef, r_wann, wcenter, h_wann, v_wann, aa_wann)
      deallocate(h_diag, evec, eval, v_eig, r_eig, r_eff)
    end do

    dg_frag%has_buffer_periodic_wannier_basis = .true.
    call comm_summation(h_diag_abs_local, h_diag_abs_global, dg_frag%icomm)
    call comm_summation(h_offdiag_abs_local, h_offdiag_abs_global, dg_frag%icomm)
    call comm_summation(r_diag_abs_local, r_diag_abs_global, 3, dg_frag%icomm)
    call comm_summation(r_offdiag_abs_local, r_offdiag_abs_global, 3, dg_frag%icomm)
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,i0,a,i0,a,i0)') "[INFO] DG buffer-periodic Wannier seed loaded: local_fragments=", &
        ifrag_count, " keep_min=", minval(dg_frag%buffer_wannier_nkeep), &
        " keep_max=", maxval(dg_frag%buffer_wannier_nkeep)
      if (use_direct_r_wann) then
        write(*,'(1x,a)') "[DG-BPW-POSITION] mode=rwann"
      else if (berry_position_used .and. .not. force_v_over_gap) then
        write(*,'(1x,a)') "[DG-BPW-POSITION] mode=berry"
      else
        write(*,'(1x,a)') "[DG-BPW-POSITION] mode=v_over_gap"
      end if
      write(*,'(1x,a,1pe13.5)') "[DG-BPW-POSITION] scale=", position_scale
      write(*,'(1x,a,8(a,1pe13.5))') "[DG-BPW-MATRIX]", &
        " h_diag_abs=", h_diag_abs_global, " h_offdiag_abs=", h_offdiag_abs_global, &
        " r_diag_abs_x=", r_diag_abs_global(1), " r_offdiag_abs_x=", r_offdiag_abs_global(1), &
        " r_diag_abs_y=", r_diag_abs_global(2), " r_offdiag_abs_y=", r_offdiag_abs_global(2), &
        " r_diag_abs_z=", r_diag_abs_global(3), " r_offdiag_abs_z=", r_offdiag_abs_global(3)
    end if
  end subroutine read_buffer_periodic_wannier_basis_data

  subroutine read_wannier90_global_basis_metadata_if_requested(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    character(16) :: env_value
    integer :: env_status, env_len
    logical :: enabled

    enabled = .false.
    env_value = ''
    call get_environment_variable('SALMON_DG_W90_GLOBAL_METADATA', env_value, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      select case(trim(adjustl(env_value(1:env_len))))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        enabled = .true.
      end select
    end if
    if (.not. enabled) return
    call read_wannier90_global_basis_metadata(dg_frag)
  end subroutine read_wannier90_global_basis_metadata_if_requested

  subroutine read_wannier90_global_basis_metadata(dg_frag)
    use filesystem, only: get_filehandle
    use communication, only: comm_is_root
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    character(48), parameter :: binfile_w90g = "wannier90_global_basis.bin"
    integer, parameter :: wannier90_global_magic = -22022216
    integer, parameter :: wannier90_global_version = 2
    character(256) :: filename
    integer :: iunit, iostat_open, magic_file, version_file
    integer :: num_bands_file, num_wann_file, n_frag_file, iw, ifrag
    integer, allocatable :: owner_frag(:), owner_count(:)
    real(8), allocatable :: center_bohr(:,:), spread_aa2(:)
    real(8) :: spread_min, spread_max

    if (.not. comm_is_root(dg_frag%id)) return

    filename = './data_dcdft/total/'//binfile_w90g
    iunit = get_filehandle()
    open(iunit, file=filename, form='unformatted', access='stream', status='old', iostat=iostat_open)
    if (iostat_open /= 0) return

    read(iunit) magic_file, version_file
    if (magic_file /= wannier90_global_magic .or. version_file < 1 .or. version_file > wannier90_global_version) then
      close(iunit)
      write(*,'(1x,a,4(a,i0))') "[WARN] invalid Wannier90 global basis metadata:", &
        " magic=", magic_file, " expected_magic=", wannier90_global_magic, &
        " version=", version_file, " max_supported_version=", wannier90_global_version
      return
    end if
    read(iunit) num_bands_file, num_wann_file, n_frag_file
    allocate(owner_frag(num_wann_file), owner_count(max(1,n_frag_file)))
    allocate(center_bohr(3,num_wann_file), spread_aa2(num_wann_file))
    read(iunit) owner_frag(1:num_wann_file)
    read(iunit) center_bohr(1:3,1:num_wann_file)
    read(iunit) spread_aa2(1:num_wann_file)
    close(iunit)

    owner_count(:) = 0
    do iw=1,num_wann_file
      ifrag = owner_frag(iw)
      if (ifrag >= 1 .and. ifrag <= n_frag_file) owner_count(ifrag) = owner_count(ifrag) + 1
    end do
    spread_min = minval(spread_aa2(1:num_wann_file))
    spread_max = maxval(spread_aa2(1:num_wann_file))
    write(*,'(1x,a,i0,a,i0,a,i0,2(a,1pe12.4))') "[DG-W90-GLOBAL] metadata loaded: bands=", &
      num_bands_file, " wann=", num_wann_file, " fragments=", n_frag_file, &
      " spread_min=", spread_min, " spread_max=", spread_max
    write(*,'(1x,a,*(1x,i0))') "[DG-W90-GLOBAL] owner_count:", owner_count(1:n_frag_file)

    deallocate(owner_frag, owner_count, center_bohr, spread_aa2)
  end subroutine read_wannier90_global_basis_metadata

  subroutine apply_wannier_flux_eigen_seed_if_available(dg_frag, ifrag_count)
    use filesystem, only: get_filehandle
    use communication, only: comm_is_root, comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: ifrag_count

    character(48), parameter :: binfile_w90seed = "wannier_flux_eigen_seed.bin"
    integer, parameter :: wannier_flux_eigen_seed_magic = -22022219
    integer, parameter :: wannier_flux_eigen_seed_version = 1
    character(256) :: filename
    integer :: iunit, io
    integer :: magic_file, version_file
    integer :: num_bands_file, num_wann_file, nstate_seed, nspin_file, n_frag_file
    integer :: ifrag, i_local, ispin, io_basis, iw, state_col, local_idx, global_idx, local_state_col
    integer :: nbasis_iter
    real(8), allocatable :: eval_seed(:,:)
    complex(8), allocatable :: seed_wannier_to_eigen(:,:)
    complex(8) :: coeff
    real(8) :: norm_local(2), norm_global(2)
    logical :: has_file

    if (.not. dg_frag%has_global_wannier_basis) return
    if (.not. allocated(dg_frag%global_wannier_coef)) return
    if (ifrag_count <= 0) return

    filename = './data_dcdft/total/'//binfile_w90seed
    inquire(file=filename, exist=has_file)
    if (.not. has_file) return

    iunit = get_filehandle()
    open(iunit, file=filename, form='unformatted', access='stream', status='old', iostat=io)
    if (io /= 0) return
    read(iunit, iostat=io) magic_file, version_file
    if (io /= 0 .or. magic_file /= wannier_flux_eigen_seed_magic .or. &
        version_file /= wannier_flux_eigen_seed_version) then
      close(iunit)
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,4(a,i0))') "[WARN] invalid Wannier Flux eigen seed:", &
          " magic=", magic_file, " expected_magic=", wannier_flux_eigen_seed_magic, &
          " version=", version_file, " expected_version=", wannier_flux_eigen_seed_version
      end if
      return
    end if
    read(iunit, iostat=io) num_bands_file, num_wann_file, nstate_seed, nspin_file, n_frag_file
    if (io /= 0 .or. num_wann_file <= 0 .or. nstate_seed <= 0 .or. nspin_file <= 0) then
      close(iunit)
      if (comm_is_root(dg_frag%id)) &
        write(*,'(1x,a)') "[WARN] invalid Wannier Flux eigen seed dimensions; keeping LCFO seed."
      return
    end if
    allocate(eval_seed(nstate_seed,nspin_file))
    allocate(seed_wannier_to_eigen(num_wann_file,nstate_seed))
    read(iunit, iostat=io) eval_seed(1:nstate_seed,1:nspin_file)
    if (io == 0) read(iunit, iostat=io) seed_wannier_to_eigen(1:num_wann_file,1:nstate_seed)
    close(iunit)
    if (io /= 0) then
      if (comm_is_root(dg_frag%id)) &
        write(*,'(1x,a)') "[WARN] failed to read complete Wannier Flux eigen seed; keeping LCFO seed."
      deallocate(eval_seed, seed_wannier_to_eigen)
      return
    end if

    if (num_wann_file /= dg_frag%global_wannier_num_wann .or. &
        num_bands_file /= dg_frag%global_wannier_num_bands .or. &
        nstate_seed > dg_frag%nstate_tot .or. nspin_file /= dg_frag%nspin) then
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,8(a,i0))') "[FATAL] Wannier Flux eigen seed mismatch:", &
          " seed_bands=", num_bands_file, " rt_bands=", dg_frag%global_wannier_num_bands, &
          " seed_wann=", num_wann_file, " rt_wann=", dg_frag%global_wannier_num_wann, &
          " seed_states=", nstate_seed, " rt_states=", dg_frag%nstate_tot, &
          " seed_nspin=", nspin_file, " rt_nspin=", dg_frag%nspin
      end if
      stop "DG-Fragment RT: incompatible Wannier Flux eigen seed"
    end if
    if (nstate_seed /= dg_frag%nstate_tot .and. comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,i0,a,i0)') "[WARN] Wannier Flux eigen seed does not cover all RT states: seed=", &
        nstate_seed, " rt=", dg_frag%nstate_tot
    end if
    if (n_frag_file /= dg_frag%n_frag .and. comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,i0,a,i0)') "[WARN] Wannier Flux eigen seed fragment count differs: seed=", &
        n_frag_file, " rt=", dg_frag%n_frag
    end if

    dg_frag%coef = (0.0d0, 0.0d0)
    if (allocated(dg_frag%coef_work)) dg_frag%coef_work = (0.0d0, 0.0d0)
    if (allocated(dg_frag%coef_new)) dg_frag%coef_new = (0.0d0, 0.0d0)
    do i_local = 1, ifrag_count
      ifrag = dg_frag%ifrag_start + i_local - 1
      do ispin = 1, dg_frag%nspin
        nbasis_iter = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%global_wannier_coef, 1), &
                          dg_frag%nstate_frag)
        do state_col = 1, nstate_seed
          if (dg_frag%coef_state_block_mode) then
            if (state_col < dg_frag%coef_state_start .or. state_col > dg_frag%coef_state_end) cycle
            local_state_col = state_col - dg_frag%coef_state_start + 1
          else
            local_state_col = state_col
          end if
          if (local_state_col < 1 .or. local_state_col > size(dg_frag%coef, 2)) cycle
          do io_basis = 1, nbasis_iter
            global_idx = dg_frag%index_basis(io_basis, ifrag, ispin)
            local_idx = 0
            if (global_idx > 0 .and. global_idx <= dg_frag%n_mat_max) &
              local_idx = dg_frag%coef_global_to_local(global_idx, ispin)
            if (local_idx < 1 .or. local_idx > size(dg_frag%coef, 1)) cycle
            coeff = (0.0d0, 0.0d0)
            do iw = 1, num_wann_file
              coeff = coeff + dg_frag%global_wannier_coef(io_basis, iw, ispin, i_local) * &
                seed_wannier_to_eigen(iw, state_col)
            end do
            dg_frag%coef(local_idx, local_state_col, ispin) = coeff
          end do
        end do
      end do
    end do
    dg_frag%esp(1:nstate_seed,1:nspin_file) = eval_seed(1:nstate_seed,1:nspin_file)
    call invalidate_coef_exchange_cache(dg_frag)

    norm_local(:) = 0.0d0
    norm_local(1) = sum(abs(dg_frag%coef)**2)
    norm_local(2) = dble(count(abs(dg_frag%coef) > 0.0d0))
    call comm_summation(norm_local, norm_global, 2, dg_frag%icomm)
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,i0,a,i0,2(a,1pe13.5))') &
        "[DG-W90-SEED] applied prepared Flux eigen seed from global Wannier basis: states=", &
        nstate_seed, " wann=", num_wann_file, " norm2=", norm_global(1), " nonzero=", norm_global(2)
    end if
    if (norm_global(1) <= 0.0d0 .or. norm_global(2) <= 0.0d0) then
      if (comm_is_root(dg_frag%id)) &
        write(*,'(1x,a)') "[FATAL] Wannier Flux eigen seed produced an empty coefficient matrix."
      stop "DG-Fragment RT: empty Wannier Flux eigen seed"
    end if

    deallocate(eval_seed, seed_wannier_to_eigen)
  end subroutine apply_wannier_flux_eigen_seed_if_available

  subroutine read_wannier90_global_basis_data(dg_frag, ifrag_count)
    use filesystem, only: get_filehandle
    use communication, only: comm_is_root
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer, intent(in) :: ifrag_count
    character(48), parameter :: binfile_w90g = "wannier90_global_basis.bin"
    integer, parameter :: wannier90_global_magic = -22022216
    integer, parameter :: wannier90_global_version = 2
    character(256) :: filename
    integer :: iunit, iostat_open, io
    integer :: magic_file, version_file, position_available
    integer :: num_bands_file, num_wann_file, n_frag_file
    logical :: ok_position

    call clear_formal_dg_wannier_data(dg_frag)
    dg_frag%has_global_wannier_basis = .false.
    dg_frag%global_wannier_num_bands = 0
    dg_frag%global_wannier_num_wann = 0
    dg_frag%global_wannier_n_frag = 0
    if (allocated(dg_frag%global_wannier_owner_frag)) deallocate(dg_frag%global_wannier_owner_frag)
    if (allocated(dg_frag%global_wannier_center)) deallocate(dg_frag%global_wannier_center)
    if (allocated(dg_frag%global_wannier_spread)) deallocate(dg_frag%global_wannier_spread)
    if (allocated(dg_frag%global_wannier_transform)) deallocate(dg_frag%global_wannier_transform)
    if (allocated(dg_frag%global_wannier_position)) deallocate(dg_frag%global_wannier_position)
    if (allocated(dg_frag%global_wannier_coef)) deallocate(dg_frag%global_wannier_coef)
    if (allocated(dg_frag%global_wannier_local_nkeep)) deallocate(dg_frag%global_wannier_local_nkeep)
    if (allocated(dg_frag%global_wannier_local_ids)) deallocate(dg_frag%global_wannier_local_ids)
    if (allocated(dg_frag%global_wannier_local_owner_frag)) deallocate(dg_frag%global_wannier_local_owner_frag)
    if (allocated(dg_frag%global_wannier_local_center)) deallocate(dg_frag%global_wannier_local_center)
    if (allocated(dg_frag%global_wannier_local_coef)) deallocate(dg_frag%global_wannier_local_coef)
    if (allocated(dg_frag%global_wannier_local_position)) deallocate(dg_frag%global_wannier_local_position)
    dg_frag%has_global_wannier_local_basis = .false.
    dg_frag%has_global_wannier_position = .false.

    if (ifrag_count <= 0) return
    filename = './data_dcdft/total/'//binfile_w90g
    iunit = get_filehandle()
    open(iunit, file=filename, form='unformatted', access='stream', status='old', iostat=iostat_open)
    if (iostat_open /= 0) return

    read(iunit, iostat=io) magic_file, version_file
    if (io /= 0 .or. magic_file /= wannier90_global_magic .or. &
        version_file < 1 .or. version_file > wannier90_global_version) then
      close(iunit)
      if (comm_is_root(dg_frag%id)) then
        write(*,'(1x,a,4(a,i0))') "[WARN] invalid Wannier90 global basis data:", &
          " magic=", magic_file, " expected_magic=", wannier90_global_magic, &
          " version=", version_file, " max_supported_version=", wannier90_global_version
      end if
      return
    end if

    read(iunit, iostat=io) num_bands_file, num_wann_file, n_frag_file
    if (io /= 0 .or. num_bands_file <= 0 .or. num_wann_file <= 0 .or. n_frag_file <= 0) then
      close(iunit)
      if (comm_is_root(dg_frag%id)) &
        write(*,'(1x,a)') "[WARN] invalid Wannier90 global basis dimensions; ignoring global Wannier data."
      return
    end if
    if (num_bands_file > dg_frag%nstate_tot) then
      close(iunit)
      if (comm_is_root(dg_frag%id)) write(*,'(1x,a,i0,a,i0)') &
        "[WARN] Wannier90 global basis has more bands than RT states: bands=", &
        num_bands_file, " nstate_tot=", dg_frag%nstate_tot
      return
    end if
    if (n_frag_file /= dg_frag%n_frag .and. comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,i0,a,i0)') "[WARN] Wannier90 global basis fragment count differs: file=", &
        n_frag_file, " runtime=", dg_frag%n_frag
    end if

    allocate(dg_frag%global_wannier_owner_frag(num_wann_file))
    allocate(dg_frag%global_wannier_center(3,num_wann_file))
    allocate(dg_frag%global_wannier_spread(num_wann_file))
    allocate(dg_frag%global_wannier_transform(num_bands_file,num_wann_file))
    allocate(dg_frag%global_wannier_position(3,num_wann_file,num_wann_file))
    dg_frag%global_wannier_position = (0.0d0, 0.0d0)
    read(iunit, iostat=io) dg_frag%global_wannier_owner_frag(1:num_wann_file)
    if (io == 0) read(iunit, iostat=io) dg_frag%global_wannier_center(1:3,1:num_wann_file)
    if (io == 0) read(iunit, iostat=io) dg_frag%global_wannier_spread(1:num_wann_file)
    if (io == 0) read(iunit, iostat=io) dg_frag%global_wannier_transform(1:num_bands_file,1:num_wann_file)
    if (io == 0 .and. version_file >= 2) then
      position_available = 0
      read(iunit, iostat=io) position_available
      if (io == 0) read(iunit, iostat=io) dg_frag%global_wannier_position(1:3,1:num_wann_file,1:num_wann_file)
      dg_frag%has_global_wannier_position = (io == 0 .and. position_available /= 0)
    end if
    close(iunit)
    if (io == 0 .and. .not. dg_frag%has_global_wannier_position) then
      call read_wannier90_global_position_rdat(num_wann_file, dg_frag%global_wannier_position, ok_position)
      dg_frag%has_global_wannier_position = ok_position
    end if
    if (io /= 0) then
      if (comm_is_root(dg_frag%id)) &
        write(*,'(1x,a)') "[WARN] failed to read complete Wannier90 global basis data; ignoring it."
      if (allocated(dg_frag%global_wannier_owner_frag)) deallocate(dg_frag%global_wannier_owner_frag)
      if (allocated(dg_frag%global_wannier_center)) deallocate(dg_frag%global_wannier_center)
      if (allocated(dg_frag%global_wannier_spread)) deallocate(dg_frag%global_wannier_spread)
      if (allocated(dg_frag%global_wannier_transform)) deallocate(dg_frag%global_wannier_transform)
      if (allocated(dg_frag%global_wannier_position)) deallocate(dg_frag%global_wannier_position)
      return
    end if

    allocate(dg_frag%global_wannier_coef(dg_frag%nstate_frag, num_wann_file, dg_frag%nspin, ifrag_count))
    dg_frag%global_wannier_coef = (0.0d0, 0.0d0)
    dg_frag%global_wannier_num_bands = num_bands_file
    dg_frag%global_wannier_num_wann = num_wann_file
    dg_frag%global_wannier_n_frag = n_frag_file
    dg_frag%has_global_wannier_basis = .true.
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,i0,a,i0,a,i0)') "[DG-W90-GLOBAL] full basis loaded for RT: bands=", &
        num_bands_file, " wann=", num_wann_file, " fragments=", n_frag_file
      if (dg_frag%has_global_wannier_position) then
        write(*,'(1x,a)') "[DG-W90-GLOBAL] Wannier90 AA_R position matrix loaded."
      else
        write(*,'(1x,a)') "[DG-W90-GLOBAL] AA_R position matrix unavailable; center-only position is used."
      end if
    end if
  end subroutine read_wannier90_global_basis_data

  subroutine read_wannier90_global_position_rdat(num_wann_expected, position_aa_r, ok)
    use filesystem, only: get_filehandle
    use inputoutput, only: au_length_aa
    use salmon_global, only: sysname
    implicit none
    integer, intent(in) :: num_wann_expected
    complex(8), intent(inout) :: position_aa_r(:,:,:)
    logical, intent(out) :: ok
    character(256) :: filename, header, env_filename
    integer :: iunit, io, num_wann_file, nrpts_file, ir
    integer :: env_status, env_len
    integer :: rvec(3), n, m
    real(8) :: rx_re, rx_im, ry_re, ry_im, rz_re, rz_im
    logical :: exists

    ok = .false.
    if (num_wann_expected <= 0) return
    if (size(position_aa_r, 1) < 3 .or. size(position_aa_r, 2) < num_wann_expected .or. &
        size(position_aa_r, 3) < num_wann_expected) return

    env_filename = ''
    call get_environment_variable('SALMON_DG_W90_R_DAT', env_filename, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      filename = trim(env_filename(1:env_len))
    else
      filename = './data_dcdft/total/'//trim(sysname)//'_r.dat'
    end if
    inquire(file=filename, exist=exists)
    if (.not. exists) return

    iunit = get_filehandle()
    open(iunit, file=filename, status='old', action='read', iostat=io)
    if (io /= 0) return
    read(iunit, '(a)', iostat=io) header
    if (io == 0) read(iunit, *, iostat=io) num_wann_file
    if (io == 0) read(iunit, *, iostat=io) nrpts_file
    if (io /= 0 .or. num_wann_file /= num_wann_expected .or. nrpts_file <= 0) then
      close(iunit)
      return
    end if

    position_aa_r(1:3,1:num_wann_expected,1:num_wann_expected) = (0.0d0, 0.0d0)
    do ir = 1, nrpts_file * num_wann_file * num_wann_file
      read(iunit, *, iostat=io) rvec(1:3), n, m, rx_re, rx_im, ry_re, ry_im, rz_re, rz_im
      if (io /= 0) exit
      if (any(rvec(1:3) /= 0)) cycle
      if (n < 1 .or. n > num_wann_file .or. m < 1 .or. m > num_wann_file) cycle
      position_aa_r(1,n,m) = cmplx(rx_re, rx_im, kind=8) / au_length_aa
      position_aa_r(2,n,m) = cmplx(ry_re, ry_im, kind=8) / au_length_aa
      position_aa_r(3,n,m) = cmplx(rz_re, rz_im, kind=8) / au_length_aa
    end do
    close(iunit)
    ok = (io == 0)
  end subroutine read_wannier90_global_position_rdat

  subroutine clear_formal_dg_wannier_data(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag
    integer :: iblock

    if (allocated(dg_frag%dg_wannier_nkeep)) deallocate(dg_frag%dg_wannier_nkeep)
    if (allocated(dg_frag%dg_wannier_global_ids)) deallocate(dg_frag%dg_wannier_global_ids)
    if (allocated(dg_frag%dg_wannier_owner_frag)) deallocate(dg_frag%dg_wannier_owner_frag)
    if (allocated(dg_frag%dg_wannier_ref_center)) deallocate(dg_frag%dg_wannier_ref_center)
    if (allocated(dg_frag%dg_wannier_basis_coef)) deallocate(dg_frag%dg_wannier_basis_coef)
    if (allocated(dg_frag%dg_wannier_h0_local)) deallocate(dg_frag%dg_wannier_h0_local)
    if (allocated(dg_frag%dg_wannier_xi_local)) deallocate(dg_frag%dg_wannier_xi_local)
    if (allocated(dg_frag%dg_wannier_coef)) deallocate(dg_frag%dg_wannier_coef)
    if (allocated(dg_frag%dg_wannier_neighbor_blocks)) then
      do iblock = 1, size(dg_frag%dg_wannier_neighbor_blocks)
        if (allocated(dg_frag%dg_wannier_neighbor_blocks(iblock)%h_flux)) &
          deallocate(dg_frag%dg_wannier_neighbor_blocks(iblock)%h_flux)
        if (allocated(dg_frag%dg_wannier_neighbor_blocks(iblock)%xi_flux)) &
          deallocate(dg_frag%dg_wannier_neighbor_blocks(iblock)%xi_flux)
      end do
      deallocate(dg_frag%dg_wannier_neighbor_blocks)
    end if
    if (allocated(dg_frag%dg_wannier_local_neighbor_block_ids)) &
      deallocate(dg_frag%dg_wannier_local_neighbor_block_ids)
    dg_frag%has_formal_dg_wannier_basis = .false.
    dg_frag%n_dg_wannier_neighbor_blocks = 0
    dg_frag%dg_wannier_blocks_gauge_consistent = .false.
    dg_frag%dg_wannier_uses_full_global_position = .false.
  end subroutine clear_formal_dg_wannier_data

  subroutine build_formal_dg_wannier_from_global_local(dg_frag)
    use communication, only: comm_is_root
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag

    integer :: ifrag_count, nkeep_max, nstate_cols
    integer :: i_local, ifrag, ispin, iaxis, iw, jw, ib, jb
    integer :: ib_global, jb_global, nbf, nkeep
    complex(8) :: hsum

    call clear_formal_dg_wannier_data(dg_frag)
    if (.not. dg_frag%has_global_wannier_local_basis) return
    if (.not. allocated(dg_frag%global_wannier_local_nkeep)) return
    if (.not. allocated(dg_frag%global_wannier_local_ids)) return
    if (.not. allocated(dg_frag%global_wannier_local_coef)) return

    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    if (ifrag_count <= 0) return
    nkeep_max = max(1, maxval(dg_frag%global_wannier_local_nkeep))
    nstate_cols = 0
    if (allocated(dg_frag%coef)) nstate_cols = size(dg_frag%coef, 2)

    allocate(dg_frag%dg_wannier_nkeep(ifrag_count))
    allocate(dg_frag%dg_wannier_global_ids(nkeep_max, ifrag_count))
    allocate(dg_frag%dg_wannier_owner_frag(nkeep_max, ifrag_count))
    allocate(dg_frag%dg_wannier_ref_center(3, ifrag_count))
    allocate(dg_frag%dg_wannier_basis_coef(dg_frag%nstate_frag, nkeep_max, dg_frag%nspin, ifrag_count))
    allocate(dg_frag%dg_wannier_h0_local(nkeep_max, nkeep_max, dg_frag%nspin, ifrag_count))
    allocate(dg_frag%dg_wannier_xi_local(3, nkeep_max, nkeep_max, dg_frag%nspin, ifrag_count))
    if (nstate_cols > 0) allocate(dg_frag%dg_wannier_coef(nkeep_max, nstate_cols, dg_frag%nspin, ifrag_count))

    dg_frag%dg_wannier_nkeep = dg_frag%global_wannier_local_nkeep
    dg_frag%dg_wannier_global_ids = 0
    dg_frag%dg_wannier_owner_frag = 0
    dg_frag%dg_wannier_ref_center = 0.0d0
    dg_frag%dg_wannier_basis_coef = (0.0d0, 0.0d0)
    dg_frag%dg_wannier_h0_local = (0.0d0, 0.0d0)
    dg_frag%dg_wannier_xi_local = (0.0d0, 0.0d0)
    if (allocated(dg_frag%dg_wannier_coef)) dg_frag%dg_wannier_coef = (0.0d0, 0.0d0)

    dg_frag%dg_wannier_global_ids(1:size(dg_frag%global_wannier_local_ids,1),1:ifrag_count) = &
      dg_frag%global_wannier_local_ids(1:size(dg_frag%global_wannier_local_ids,1),1:ifrag_count)
    dg_frag%dg_wannier_owner_frag(1:size(dg_frag%global_wannier_local_owner_frag,1),1:ifrag_count) = &
      dg_frag%global_wannier_local_owner_frag(1:size(dg_frag%global_wannier_local_owner_frag,1),1:ifrag_count)
    dg_frag%dg_wannier_basis_coef(1:size(dg_frag%global_wannier_local_coef,1),1:size(dg_frag%global_wannier_local_coef,2), &
      1:dg_frag%nspin,1:ifrag_count) = dg_frag%global_wannier_local_coef(1:size(dg_frag%global_wannier_local_coef,1), &
      1:size(dg_frag%global_wannier_local_coef,2),1:dg_frag%nspin,1:ifrag_count)

    do i_local = 1, ifrag_count
      ifrag = dg_frag%ifrag_start + i_local - 1
      if (allocated(dg_frag%buffer_wannier_frag_center) .and. &
          i_local <= size(dg_frag%buffer_wannier_frag_center, 2)) then
        dg_frag%dg_wannier_ref_center(1:3, i_local) = dg_frag%buffer_wannier_frag_center(1:3, i_local)
      else if (dg_frag%dg_wannier_nkeep(i_local) > 0 .and. allocated(dg_frag%global_wannier_local_center)) then
        dg_frag%dg_wannier_ref_center(1:3, i_local) = &
          sum(dg_frag%global_wannier_local_center(1:3,1:dg_frag%dg_wannier_nkeep(i_local),i_local), dim=2) / &
          dble(dg_frag%dg_wannier_nkeep(i_local))
      end if

      nkeep = dg_frag%dg_wannier_nkeep(i_local)
      do ispin = 1, dg_frag%nspin
        nbf = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%dg_wannier_basis_coef, 1))
        if (nkeep <= 0 .or. nbf <= 0) cycle
        if (allocated(dg_frag%global_wannier_local_position) .and. dg_frag%has_global_wannier_position) then
          do jw = 1, nkeep
            do iw = 1, nkeep
              dg_frag%dg_wannier_xi_local(1:3, iw, jw, ispin, i_local) = &
                dg_frag%global_wannier_local_position(1:3, iw, jw, i_local)
            end do
          end do
        else if (allocated(dg_frag%global_wannier_local_center)) then
          do iw = 1, nkeep
            dg_frag%dg_wannier_xi_local(1:3, iw, iw, ispin, i_local) = &
              cmplx(dg_frag%global_wannier_local_center(1:3, iw, i_local), 0.0d0, kind=8)
          end do
        end if
        do iaxis = 1, 3
          do iw = 1, nkeep
            dg_frag%dg_wannier_xi_local(iaxis, iw, iw, ispin, i_local) = &
              dg_frag%dg_wannier_xi_local(iaxis, iw, iw, ispin, i_local) - &
              cmplx(dg_frag%dg_wannier_ref_center(iaxis, i_local), 0.0d0, kind=8)
          end do
        end do

        if (allocated(dg_frag%H_mat) .and. allocated(dg_frag%index_basis)) then
          do jw = 1, nkeep
            do iw = 1, nkeep
              hsum = (0.0d0, 0.0d0)
              do jb = 1, nbf
                jb_global = dg_frag%index_basis(jb, ifrag, ispin)
                if (jb_global < 1 .or. jb_global > size(dg_frag%H_mat, 2)) cycle
                do ib = 1, nbf
                  ib_global = dg_frag%index_basis(ib, ifrag, ispin)
                  if (ib_global < 1 .or. ib_global > size(dg_frag%H_mat, 1)) cycle
                  hsum = hsum + conjg(dg_frag%dg_wannier_basis_coef(ib, iw, ispin, i_local)) * &
                    cmplx(dg_frag%H_mat(ib_global, jb_global, ispin), 0.0d0, kind=8) * &
                    dg_frag%dg_wannier_basis_coef(jb, jw, ispin, i_local)
                end do
              end do
              dg_frag%dg_wannier_h0_local(iw, jw, ispin, i_local) = hsum
            end do
          end do
        end if
      end do
    end do

    dg_frag%has_formal_dg_wannier_basis = .true.
    dg_frag%dg_wannier_blocks_gauge_consistent = dg_frag%has_global_wannier_position
    dg_frag%dg_wannier_uses_full_global_position = dg_frag%has_global_wannier_position
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,i0,a,i0,a,l1)') &
        "[DG-WANNIER-FORMAL] local blocks initialized from global Wannier data: nkeep_max=", &
        nkeep_max, " local_fragments=", ifrag_count, " has_AA_R=", dg_frag%has_global_wannier_position
    end if
  end subroutine build_formal_dg_wannier_from_global_local

  subroutine build_global_wannier_local_basis(dg_frag)
    use communication, only: comm_is_root
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag

    integer :: ifrag_count, ifrag, i_local, iw, ilw, jlw, ispin, nbf
    integer :: nkeep_max, nkeep_min, nkeep_total
    integer :: center_range, env_status, env_len
    character(32) :: env_center_range

    dg_frag%has_global_wannier_local_basis = .false.
    if (.not. dg_frag%has_global_wannier_basis) return
    if (.not. allocated(dg_frag%global_wannier_coef)) return
    if (.not. allocated(dg_frag%global_wannier_owner_frag)) return
    if (.not. allocated(dg_frag%global_wannier_center)) return

    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    if (ifrag_count <= 0) return
    center_range = 1
    env_center_range = ''
    call get_environment_variable('SALMON_DG_GLOBAL_WANNIER_H_CENTER_RANGE', env_center_range, &
      length=env_len, status=env_status)
    if (env_status /= 0 .or. env_len <= 0) then
      call get_environment_variable('SALMON_DG_GLOBAL_WANNIER_CENTER_RANGE', env_center_range, &
        length=env_len, status=env_status)
    end if
    if (env_status == 0 .and. env_len > 0) then
      select case (adjustl(trim(env_center_range(1:env_len))))
      case ('all','ALL','global','GLOBAL','full','FULL')
        center_range = -1
      case default
        read(env_center_range(1:env_len), *, iostat=env_status) center_range
        if (env_status /= 0) center_range = 1
      end select
    end if

    if (allocated(dg_frag%global_wannier_local_nkeep)) deallocate(dg_frag%global_wannier_local_nkeep)
    if (allocated(dg_frag%global_wannier_local_ids)) deallocate(dg_frag%global_wannier_local_ids)
    if (allocated(dg_frag%global_wannier_local_owner_frag)) deallocate(dg_frag%global_wannier_local_owner_frag)
    if (allocated(dg_frag%global_wannier_local_center)) deallocate(dg_frag%global_wannier_local_center)
    if (allocated(dg_frag%global_wannier_local_coef)) deallocate(dg_frag%global_wannier_local_coef)
    if (allocated(dg_frag%global_wannier_local_position)) deallocate(dg_frag%global_wannier_local_position)

    allocate(dg_frag%global_wannier_local_nkeep(ifrag_count))
    dg_frag%global_wannier_local_nkeep(:) = 0
    do i_local = 1, ifrag_count
      ifrag = dg_frag%ifrag_start + i_local - 1
      do iw = 1, dg_frag%global_wannier_num_wann
        if (global_wannier_owner_in_h_range(dg_frag%global_wannier_owner_frag(iw), ifrag, &
            center_range, dg_frag%num_fragment)) then
          dg_frag%global_wannier_local_nkeep(i_local) = dg_frag%global_wannier_local_nkeep(i_local) + 1
        end if
      end do
    end do

    nkeep_max = max(1, maxval(dg_frag%global_wannier_local_nkeep))
    allocate(dg_frag%global_wannier_local_ids(nkeep_max, ifrag_count))
    allocate(dg_frag%global_wannier_local_owner_frag(nkeep_max, ifrag_count))
    allocate(dg_frag%global_wannier_local_center(3, nkeep_max, ifrag_count))
    allocate(dg_frag%global_wannier_local_coef(dg_frag%nstate_frag, nkeep_max, dg_frag%nspin, ifrag_count))
    allocate(dg_frag%global_wannier_local_position(3, nkeep_max, nkeep_max, ifrag_count))
    dg_frag%global_wannier_local_ids = 0
    dg_frag%global_wannier_local_owner_frag = 0
    dg_frag%global_wannier_local_center = 0.0d0
    dg_frag%global_wannier_local_coef = (0.0d0, 0.0d0)
    dg_frag%global_wannier_local_position = (0.0d0, 0.0d0)

    do i_local = 1, ifrag_count
      ifrag = dg_frag%ifrag_start + i_local - 1
      ilw = 0
      do iw = 1, dg_frag%global_wannier_num_wann
        if (.not. global_wannier_owner_in_h_range(dg_frag%global_wannier_owner_frag(iw), ifrag, &
            center_range, dg_frag%num_fragment)) cycle
        ilw = ilw + 1
        dg_frag%global_wannier_local_ids(ilw, i_local) = iw
        dg_frag%global_wannier_local_owner_frag(ilw, i_local) = dg_frag%global_wannier_owner_frag(iw)
        dg_frag%global_wannier_local_center(1:3, ilw, i_local) = dg_frag%global_wannier_center(1:3, iw)
        do ispin = 1, dg_frag%nspin
          nbf = min(dg_frag%n_basis(ifrag, ispin), size(dg_frag%global_wannier_coef, 1), &
                    size(dg_frag%global_wannier_local_coef, 1))
          if (nbf > 0) then
            dg_frag%global_wannier_local_coef(1:nbf, ilw, ispin, i_local) = &
              dg_frag%global_wannier_coef(1:nbf, iw, ispin, i_local)
          end if
        end do
      end do
      if (dg_frag%has_global_wannier_position .and. allocated(dg_frag%global_wannier_position)) then
        do ilw = 1, dg_frag%global_wannier_local_nkeep(i_local)
          iw = dg_frag%global_wannier_local_ids(ilw, i_local)
          do jlw = 1, dg_frag%global_wannier_local_nkeep(i_local)
            dg_frag%global_wannier_local_position(1:3, ilw, jlw, i_local) = &
              dg_frag%global_wannier_position(1:3, iw, dg_frag%global_wannier_local_ids(jlw, i_local))
          end do
        end do
      end if
    end do

    dg_frag%has_global_wannier_local_basis = .true.
    if (comm_is_root(dg_frag%id)) then
      nkeep_min = minval(dg_frag%global_wannier_local_nkeep)
      nkeep_total = sum(dg_frag%global_wannier_local_nkeep)
      write(*,'(1x,a,i0,a,i0,a,i0,a,i0)') &
        "[DG-W90-GLOBAL-LOCAL] built local-selected global Wannier basis: h_center_range=", &
        center_range, " keep_min=", nkeep_min, " keep_max=", maxval(dg_frag%global_wannier_local_nkeep), &
        " keep_local_total=", nkeep_total
    end if
  end subroutine build_global_wannier_local_basis

  logical function global_wannier_owner_in_h_range(owner_frag, center_frag, center_range, num_fragment) result(in_range)
    implicit none
    integer, intent(in) :: owner_frag, center_frag, center_range, num_fragment(3)

    if (owner_frag <= 0 .or. center_frag <= 0) then
      in_range = .false.
    else if (center_range < 0) then
      in_range = .true.
    else
      in_range = (fragment_periodic_manhattan_distance(owner_frag, center_frag, num_fragment) <= center_range)
    end if
  end function global_wannier_owner_in_h_range

  integer function fragment_periodic_manhattan_distance(ifrag_a, ifrag_b, num_fragment) result(dist)
    implicit none
    integer, intent(in) :: ifrag_a, ifrag_b, num_fragment(3)
    integer :: axis, da, ca(3), cb(3)

    dist = huge(1)
    if (any(num_fragment(1:3) <= 0)) return
    call fragment_coord3_from_id_io(ifrag_a, num_fragment, ca)
    call fragment_coord3_from_id_io(ifrag_b, num_fragment, cb)
    dist = 0
    do axis = 1, 3
      da = abs(ca(axis) - cb(axis))
      da = min(da, max(0, num_fragment(axis) - da))
      dist = dist + da
    end do
  end function fragment_periodic_manhattan_distance

  subroutine fragment_coord3_from_id_io(ifrag, num_fragment, coord)
    implicit none
    integer, intent(in) :: ifrag, num_fragment(3)
    integer, intent(out) :: coord(3)
    integer :: tmp

    coord(1:3) = 1
    if (ifrag <= 0 .or. any(num_fragment(1:3) <= 0)) return
    tmp = ifrag - 1
    coord(3) = mod(tmp, num_fragment(3)) + 1
    tmp = tmp / num_fragment(3)
    coord(2) = mod(tmp, num_fragment(2)) + 1
    tmp = tmp / num_fragment(2)
    coord(1) = mod(tmp, num_fragment(1)) + 1
  end subroutine fragment_coord3_from_id_io

  subroutine initialize_coefficients_from_buffer_wannier_flux(dg_frag)
    use eigen_subdiag_sub, only: eigen_dsyev
    use communication, only: comm_is_root, comm_summation
    implicit none
    type(s_dg_fragment_rt), intent(inout) :: dg_frag

    integer :: ifrag, i_local, ifrag_count, ispin, iw, io, k, state_col
    integer :: iaxis, occ, virt, env_status, env_len
    integer :: nocc_eff, occ_base, occ_extra, frag_occ_s, frag_occ_e, nocc_frag_seed
    integer :: nvirt_eff, virt_base, virt_extra, frag_virt_s, frag_virt_e, nvirt_frag_seed
    integer :: nseed_frag
    integer :: nw, nbf, global_idx, local_idx, local_state_col
    real(8), allocatable :: h_work(:,:), evec(:,:), eval(:), coef_vec(:), r_work(:,:), r_eig(:,:)
    real(8) :: diag_local(4), diag_global(4)
    real(8) :: trans_local(3,5), trans_global(3,5), amp, gap, strength_pair
    real(8) :: min_gap_local(3), max_pair_local(3), max_pair_gap_local(3)
    logical :: trace_bpw_transition
    character(32) :: env_bpw_transition

    if (.not. dg_frag%has_buffer_periodic_wannier_basis) return
    if (.not. allocated(dg_frag%buffer_wannier_coef)) return
    if (.not. allocated(dg_frag%buffer_wannier_h_flux)) return

    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    if (ifrag_count <= 0) return

    dg_frag%coef = (0.0d0, 0.0d0)
    dg_frag%buffer_wannier_flux_seed_applied = .false.
    diag_local(:) = 0.0d0
    diag_global(:) = 0.0d0
    trans_local(:, :) = 0.0d0
    trans_global(:, :) = 0.0d0
    trace_bpw_transition = .false.
    env_bpw_transition = ''
    call get_environment_variable('SALMON_DG_BPW_TRANSITION_TRACE', env_bpw_transition, length=env_len, status=env_status)
    if (env_status == 0 .and. env_len > 0) then
      select case(trim(adjustl(env_bpw_transition(1:env_len))))
      case('1','y','Y','yes','YES','true','TRUE','on','ON')
        trace_bpw_transition = .true.
      end select
    end if

    do i_local = 1, ifrag_count
      ifrag = dg_frag%ifrag_start + i_local - 1
      do ispin = 1, dg_frag%nspin
        nw = dg_frag%buffer_wannier_nkeep(i_local)
        nbf = min(dg_frag%n_basis(ifrag, ispin), dg_frag%nstate_frag, &
                  size(dg_frag%buffer_wannier_coef, 1))
        if (nw <= 0 .or. nbf <= 0) cycle

        nocc_eff = dg_frag%nstate_tot
        if (allocated(dg_frag%nocc_spin)) then
          if (ispin <= size(dg_frag%nocc_spin)) nocc_eff = min(dg_frag%nstate_tot, max(0, dg_frag%nocc_spin(ispin)))
        end if
        occ_base = nocc_eff / max(1, dg_frag%n_frag)
        occ_extra = mod(nocc_eff, max(1, dg_frag%n_frag))
        frag_occ_s = (ifrag - 1) * occ_base + min(ifrag - 1, occ_extra) + 1
        frag_occ_e = ifrag * occ_base + min(ifrag, occ_extra)
        nocc_frag_seed = max(0, frag_occ_e - frag_occ_s + 1)
        nvirt_eff = max(0, dg_frag%nstate_tot - nocc_eff)
        virt_base = nvirt_eff / max(1, dg_frag%n_frag)
        virt_extra = mod(nvirt_eff, max(1, dg_frag%n_frag))
        frag_virt_s = nocc_eff + (ifrag - 1) * virt_base + min(ifrag - 1, virt_extra) + 1
        frag_virt_e = nocc_eff + ifrag * virt_base + min(ifrag, virt_extra)
        nvirt_frag_seed = max(0, frag_virt_e - frag_virt_s + 1)
        nseed_frag = nocc_frag_seed + nvirt_frag_seed
        if (nseed_frag <= 0) cycle
        if (nseed_frag > nw) then
          write(*,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
            "[FATAL] buffer-periodic Wannier Flux seed has too few states: rank=", &
            dg_frag%id, " ifrag=", ifrag, " nocc_frag=", nocc_frag_seed, &
            " nvirt_frag=", nvirt_frag_seed, " needed=", nseed_frag, " nw=", nw
          stop "DG-Fragment RT: insufficient buffer-periodic Wannier states"
        end if

        allocate(h_work(nw,nw), evec(nw,nw), eval(nw), coef_vec(nbf))
        h_work(1:nw,1:nw) = dg_frag%buffer_wannier_h_flux(1:nw,1:nw,i_local)
        call eigen_dsyev(h_work, eval, evec)

        if (trace_bpw_transition .and. nocc_frag_seed > 0 .and. nvirt_frag_seed > 0) then
          min_gap_local(:) = 1.0d300
          max_pair_local(:) = 0.0d0
          max_pair_gap_local(:) = 0.0d0
          allocate(r_work(nw,nw), r_eig(nw,nw))
          do iaxis = 1, 3
            r_work(1:nw,1:nw) = matmul(dg_frag%buffer_wannier_v(iaxis,1:nw,1:nw,i_local), &
              evec(1:nw,1:nw))
            r_eig(1:nw,1:nw) = matmul(transpose(evec(1:nw,1:nw)), r_work(1:nw,1:nw))
            do virt = nocc_frag_seed + 1, nseed_frag
              do occ = 1, nocc_frag_seed
                gap = eval(virt) - eval(occ)
                amp = r_eig(occ,virt)
                strength_pair = amp * amp
                min_gap_local(iaxis) = min(min_gap_local(iaxis), gap)
                if (strength_pair > max_pair_local(iaxis)) then
                  max_pair_local(iaxis) = strength_pair
                  max_pair_gap_local(iaxis) = gap
                end if
                trans_local(iaxis,1) = trans_local(iaxis,1) + 1.0d0
                trans_local(iaxis,2) = trans_local(iaxis,2) + strength_pair
                if (gap > 1.0d-12) then
                  trans_local(iaxis,3) = trans_local(iaxis,3) + strength_pair / gap
                  trans_local(iaxis,4) = trans_local(iaxis,4) + 2.0d0 * gap * strength_pair
                  trans_local(iaxis,5) = trans_local(iaxis,5) + gap * strength_pair
                end if
              end do
            end do
          end do
          deallocate(r_work, r_eig)
          do iaxis = 1, 3
            write(*,'(1x,a,i0,a,i0,a,i0,a,i0,4(a,1pe13.5))') &
              "[DG-BPW-TRANSITION-LOCAL] rank=", dg_frag%id, " ifrag=", ifrag, " ispin=", ispin, &
              " axis=", iaxis, " min_gap_eV=", min_gap_local(iaxis) * 27.211386245988d0, &
              " max_pair_r2=", max_pair_local(iaxis), &
              " max_pair_gap_eV=", max_pair_gap_local(iaxis) * 27.211386245988d0, &
              " nvirt_frag=", dble(nvirt_frag_seed)
          end do
        end if

        do k = 1, nseed_frag
          if (k <= nocc_frag_seed) then
            state_col = frag_occ_s + k - 1
          else
            state_col = frag_virt_s + (k - nocc_frag_seed) - 1
          end if
          if (state_col < 1 .or. state_col > dg_frag%nstate_tot) cycle
          if (dg_frag%coef_state_block_mode) then
            if (state_col < dg_frag%coef_state_start .or. state_col > dg_frag%coef_state_end) cycle
            local_state_col = state_col - dg_frag%coef_state_start + 1
          else
            local_state_col = state_col
          end if
          if (local_state_col < 1 .or. local_state_col > size(dg_frag%coef, 2)) cycle

          coef_vec(1:nbf) = 0.0d0
          do iw = 1, nw
            coef_vec(1:nbf) = coef_vec(1:nbf) + &
              dg_frag%buffer_wannier_coef(1:nbf,iw,ispin,i_local) * evec(iw,k)
          end do

          do io = 1, nbf
            global_idx = dg_frag%index_basis(io, ifrag, ispin)
            local_idx = 0
            if (global_idx > 0 .and. global_idx <= dg_frag%n_mat_max) &
              local_idx = dg_frag%coef_global_to_local(global_idx, ispin)
            if (local_idx > 0 .and. local_idx <= size(dg_frag%coef, 1)) then
              dg_frag%coef(local_idx, local_state_col, ispin) = dcmplx(coef_vec(io), 0.0d0)
            end if
          end do
          if (allocated(dg_frag%esp)) then
            if (state_col <= size(dg_frag%esp, 1) .and. ispin <= size(dg_frag%esp, 2)) &
              dg_frag%esp(state_col, ispin) = eval(k)
          end if
          diag_local(1) = diag_local(1) + 1.0d0
          if (k > nocc_frag_seed) diag_local(4) = diag_local(4) + 1.0d0
        end do

        diag_local(2) = diag_local(2) + dble(nw)
        diag_local(3) = diag_local(3) + eval(1)
        deallocate(h_work, evec, eval, coef_vec)
      end do
    end do

    call comm_summation(diag_local, diag_global, 4, dg_frag%icomm)
    if (trace_bpw_transition) call comm_summation(trans_local, trans_global, 15, dg_frag%icomm)
    if (diag_global(1) > 0.0d0) dg_frag%buffer_wannier_flux_seed_applied = .true.
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,4(a,1pe13.5))') "[INFO] DG buffer-periodic Flux eigenstate seed applied:", &
        " seeded_states=", diag_global(1), " virtual_states=", diag_global(4), &
        " total_local_wannier=", diag_global(2), &
        " eval1_sum=", diag_global(3)
      if (trace_bpw_transition) then
        do iaxis = 1, 3
          write(*,'(1x,a,i0,6(a,1pe13.5))') "[DG-BPW-TRANSITION] axis=", iaxis, &
            " pairs=", trans_global(iaxis,1), &
            " sum_r2=", trans_global(iaxis,2), &
            " sum_r2_over_gap=", trans_global(iaxis,3), &
            " fsum_like=", trans_global(iaxis,4) / max(1.0d0, diag_global(1) - diag_global(4)), &
            " mean_gap_eV=", 27.211386245988d0 * trans_global(iaxis,5) / max(1.0d-300, trans_global(iaxis,2)), &
            " per_occ_r2=", trans_global(iaxis,2) / max(1.0d0, diag_global(1) - diag_global(4))
        end do
      end if
    end if
  end subroutine initialize_coefficients_from_buffer_wannier_flux

  subroutine diagnose_local_wannier_tail(dg_frag)
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag

    real(8), parameter :: tail_tol = 1.0d-8
    integer :: ifrag, i_local, ifrag_count, ispin, iw, ib
    integer :: ix, iy, iz, b, b_req, bmax, nbf, nkeep
    integer :: ibx, iby, ibz
    integer :: core_lo(3), core_hi(3)
    real(8) :: hvol, wval, rho_w, q_total, q_core, q_shell, max_shell_abs
    real(8) :: core_frac, shell_frac, outside_frac, required_tail
    real(8) :: min_core_frac, max_shell_frac, max_outside_frac, max_outside_at_buffer
    integer :: max_required_buffer, owned_count
    real(8), allocatable :: q_outside(:)
    character(16) :: tail_status

    if (.not. dg_frag%is_frag_root) return
    if (.not. allocated(dg_frag%phi_frag)) return
    if (.not. allocated(dg_frag%local_wannier_coef)) return
    if (.not. allocated(dg_frag%local_wannier_owned)) return
    ifrag_count = dg_frag%ifrag_end - dg_frag%ifrag_start + 1
    if (ifrag_count <= 0) return

    hvol = dg_frag%hgs(1) * dg_frag%hgs(2) * dg_frag%hgs(3)
    bmax = max(0, maxval(dg_frag%nxyz_buffer(1:3)))
    allocate(q_outside(0:bmax))

    do i_local = 1, ifrag_count
      ifrag = dg_frag%ifrag_start + i_local - 1
      nbf = min(dg_frag%n_basis(ifrag, 1), size(dg_frag%local_wannier_coef, 1))
      nkeep = min(dg_frag%local_wannier_nkeep(i_local), size(dg_frag%local_wannier_coef, 2))
      core_lo(1:3) = dg_frag%ixyz_frag(1:3, ifrag)
      core_hi(1:3) = dg_frag%ixyz_frag(1:3, ifrag) + dg_frag%nxyz_domain(1:3, ifrag) - 1
      min_core_frac = 1.0d0
      max_shell_frac = 0.0d0
      max_outside_frac = 0.0d0
      max_outside_at_buffer = 0.0d0
      max_required_buffer = 0
      owned_count = 0
      tail_status = 'OK'

      do ispin = 1, dg_frag%nspin
        do iw = 1, nkeep
          if (.not. dg_frag%local_wannier_owned(iw, ispin, i_local)) cycle
          owned_count = owned_count + 1
          q_total = 0.0d0
          q_core = 0.0d0
          q_shell = 0.0d0
          max_shell_abs = 0.0d0
          q_outside(0:bmax) = 0.0d0

          do iz = lbound(dg_frag%phi_frag, 3), ubound(dg_frag%phi_frag, 3)
            do iy = lbound(dg_frag%phi_frag, 2), ubound(dg_frag%phi_frag, 2)
              do ix = lbound(dg_frag%phi_frag, 1), ubound(dg_frag%phi_frag, 1)
                wval = 0.0d0
                do ib = 1, nbf
                  wval = wval + dg_frag%phi_frag(ix, iy, iz, ib, i_local) * &
                    dg_frag%local_wannier_coef(ib, iw, ispin, i_local)
                end do
                rho_w = wval * wval * hvol
                q_total = q_total + rho_w
                if (ix >= core_lo(1) .and. ix <= core_hi(1) .and. &
                    iy >= core_lo(2) .and. iy <= core_hi(2) .and. &
                    iz >= core_lo(3) .and. iz <= core_hi(3)) then
                  q_core = q_core + rho_w
                end if
                if (ix == lbound(dg_frag%phi_frag, 1) .or. ix == ubound(dg_frag%phi_frag, 1) .or. &
                    iy == lbound(dg_frag%phi_frag, 2) .or. iy == ubound(dg_frag%phi_frag, 2) .or. &
                    iz == lbound(dg_frag%phi_frag, 3) .or. iz == ubound(dg_frag%phi_frag, 3)) then
                  q_shell = q_shell + rho_w
                  max_shell_abs = max(max_shell_abs, abs(wval))
                end if
                do b = 0, bmax
                  ibx = merge(1, 0, ix < core_lo(1) - b .or. ix > core_hi(1) + b)
                  iby = merge(1, 0, iy < core_lo(2) - b .or. iy > core_hi(2) + b)
                  ibz = merge(1, 0, iz < core_lo(3) - b .or. iz > core_hi(3) + b)
                  if (ibx + iby + ibz > 0) q_outside(b) = q_outside(b) + rho_w
                end do
              end do
            end do
          end do

          if (q_total > 0.0d0) then
            core_frac = q_core / q_total
            shell_frac = q_shell / q_total
            outside_frac = max(0.0d0, 1.0d0 - core_frac)
            max_outside_at_buffer = max(max_outside_at_buffer, q_outside(bmax) / q_total)
            required_tail = tail_tol * q_total
            b_req = bmax
            do b = 0, bmax
              if (q_outside(b) <= required_tail) then
                b_req = b
                exit
              end if
            end do
            if (q_shell > required_tail) b_req = bmax + 1
            min_core_frac = min(min_core_frac, core_frac)
            max_shell_frac = max(max_shell_frac, shell_frac)
            max_outside_frac = max(max_outside_frac, outside_frac)
            max_required_buffer = max(max_required_buffer, b_req)
            if (q_outside(bmax) > required_tail .or. q_shell > required_tail) tail_status = 'INSUFFICIENT'
          end if
        end do
      end do

      if (owned_count > 0) then
        write(*,'(1x,a,i0,a,i0,a,i0,a,3(i0,1x),a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,i0,a,a,a,1pe12.4)') &
          "[DG-WANNIER-TAIL] rank=", dg_frag%id, " fragment=", ifrag, &
          " owned=", owned_count, " buffer=", dg_frag%nxyz_buffer(1), dg_frag%nxyz_buffer(2), &
          dg_frag%nxyz_buffer(3), " min_core_frac=", min_core_frac, &
          " max_outside_core=", max_outside_frac, " outer_shell_frac=", max_shell_frac, &
          " outside_at_buffer=", max_outside_at_buffer, &
          " recommended_buffer=", max_required_buffer, " status=", trim(tail_status), " tol=", tail_tol
      end if
    end do

    deallocate(q_outside)
  end subroutine diagnose_local_wannier_tail

end module rt_dg_fragment_io
