!
!  Copyright 2019-2020 SALMON developers
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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

#include "config.h"

! Ground state of the reference cell in dm_unfold_option='reference':
!   check_dm_unfold_reference : stop, before the SCF, if the run cannot produce the Hamiltonian file
!   write_dm_unfold_reference : write the physical arrays that define the Kohn-Sham Hamiltonian
!                               (local potential, non-local pseudopotential projectors, lattice, stencil)
!                               of the reference cell into <dir_out_restart>/dm_unfold_ref_hamiltonian.bin
!
! The file is a Fortran unformatted sequential file (native byte order, default integer, real(8)).
! It is meant to be kept together with wfn.bin etc. of the same data_for_restart directory.
! Record layout (format_version = 1); see also the header comment of each record in write_dm_unfold_reference:
!   1: file_magic(16), format_version
!   2: N(3), nspin, nion, nk, no, nd, nps, nlma
!   3: hgs(3), hvol, primitive_a(3,3)
!   4: if_orthogonal, coef_lap0, coef_lap(4,3), coef_nab(4,3), coef_f(6),
!      coef_lap0_nd1, coef_lap_nd1(1,3), coef_nab_nd1(1,3)
!   5: kion(nion), Rion(3,nion)
!   6: convergence(16), threshold, iter, sum1, E_tot      (informational only; see below)
!   7..6+nspin: V_local(N(1),N(2),N(3)) for each spin
!   next: mps(nion), jxyz(3,nps,nion), rxyz(3,nps,nion), ia_tbl(nlma), rinv_uvu(nlma), uv(nps,nlma)
! The file is meant to be read back by the same build: native byte order and default kinds are used, and the kinds and
! byte sizes are not recorded in the file. All the values are in atomic units. The unused tail (j > mps(ia)) of
! jxyz, rxyz and uv is zero.
! Hamiltonian terms that are not contained in the file (DFT+U, XC operator depending on the kinetic-energy density,
! HSE06, absorbing boundary) are rejected in check_dm_unfold_reference.
! Record 6 is not a convergence flag: "iter" is the number of SCF iterations done, "sum1" is the last value of the
! convergence measure selected by "convergence" (1d9 if no iteration was done), and the reason why the SCF stopped is
! not recorded. The accuracy of the Hamiltonian/orbitals has to be judged by the band residuals on the user side.
!
! Restore side (first stage): the Hamiltonian H(q) of the reference cell for an arbitrary Bloch vector q, rebuilt from
! the file and applied with the ordinary hpsi (no change in hpsi and the other existing routines):
!   read_dm_unfold_reference    : read the file (the root process of the communicator reads and broadcasts)
!   build_dm_unfold_reference   : check against the super cell and the current settings, set up the structures hpsi needs
!   apply_h_dm_unfold_reference : H(q)*psi for no_work vectors given on the real-space grid of the reference cell
!   release_dm_unfold_reference : free everything (read and build can be done again)
! The Hamiltonian is not divided in real space: every process which calls these routines has its own complete copy.
! build rejects the builds/settings which hpsi would handle with the tables of the grid of this run: OpenACC and the
! stencil with explicit vectorization. A build without MPI is rejected too. The other optional libraries (ScaLAPACK,
! FFTE, FFTW) are not used by the H action of these routines and are not rejected.
! Errors in the restore routines terminate the whole job (MPI_Abort), because they are not necessarily called by all the
! processes at the same time.
! The Bloch vector q is the Cartesian one (the same convention as system%vec_k), and no vector potential is added.
!
! Eigenstates of H(q) (second stage):
!   solve_dm_unfold_reference   : no_work states of H(q) from the bottom of the spectrum, found by the iteration which the
!                                 ground state calculation and the projection of the real-time calculation use
!                                 (subspace diagonalization -> conjugate gradient -> Gram-Schmidt), starting from the vectors
!                                 given on entry (for example, the states for a q nearby). The convergence is judged by
!                                 the residual |H phi - e phi| of the lowest no_conv of the no_work states, after a last subspace
!                                 diagonalization. Only these no_conv states are to be used as the eigenstates: the others are
!                                 buffer (work) states which help the next call and are not eigenstates of H(q).
!                                 The routine needs yn_scalapack='n' (the processes of this module are not set up for ScaLAPACK).
!
! Matrix elements (third stage):
!   calc_tm_dm_unfold_reference : -i <psi_i|d/dr|psi_j> and the commutator term of the non-local pseudopotential for given
!                                 vectors at an arbitrary Bloch vector q, with the definitions of write_tm_data (tm.bin of
!                                 the ground state), so that the matrices of the states of H(q) from solve_dm_unfold_reference
!                                 are consistent with the Hamiltonian they are eigenstates of.

module dm_unfold_reference_sub
  use structures, only: s_rgrid, s_parallel_info, s_dft_system, s_stencil, s_scalar, s_pp_grid, s_sendrecv_grid, s_orbital, s_cg
  implicit none
  private

  public :: s_dm_unfold_ref
  public :: check_dm_unfold_reference
  public :: write_dm_unfold_reference
  public :: read_dm_unfold_reference
  public :: build_dm_unfold_reference
  public :: apply_h_dm_unfold_reference
  public :: calc_tm_dm_unfold_reference
  public :: solve_dm_unfold_reference
  public :: release_dm_unfold_reference

  ! The Hamiltonian of the reference cell. The components are private: they are set only by the routines of this module.
  ! Immutable after build: lattice, V_local, ppg (except zekr_uV), stencil, grid.
  ! Changed by apply and solve: system%vec_k (= q), ppg%zekr_uV, and the work arrays/dimensions which depend on the number of
  ! vectors.
  type s_dm_unfold_ref
    private
    logical :: is_read  = .false.
    logical :: is_built = .false.
    integer :: no_work     = 0    ! number of vectors the work arrays are made for
    integer :: file_nk  = 0    ! number of k points and bands of the reference ground state (informational)
    integer :: file_no  = 0
    real(8) :: file_hvol = 0d0 ! hvol written in the file
    type(s_dft_system)    :: system
    type(s_stencil)       :: stencil
    type(s_rgrid)         :: lg, mg
    type(s_parallel_info) :: info
    type(s_scalar), allocatable :: V_local(:)
    type(s_pp_grid)       :: ppg
    type(s_orbital)       :: tpsi, htpsi
    type(s_orbital)       :: gpsi ! work arrays of solve_dm_unfold_reference (allocated when it is used)
    type(s_cg)            :: cg
  end type s_dm_unfold_ref

contains

  subroutine check_dm_unfold_reference(info,system)
    use structures, only: s_parallel_info, s_dft_system
    use salmon_global, only: theory, iperiodic, yn_opt, yn_dc, yn_jm, yn_spinorbit, xc, absorbing_boundary, &
                           & write_gs_restart_data, method_wf_distributor
    use plusU_global, only: PLUS_U_ON
    use communication, only: comm_is_root
    use parallelization, only: nproc_id_global, end_parallel
    implicit none
    type(s_parallel_info), intent(in) :: info
    type(s_dft_system),    intent(in) :: system
    character(256) :: msg

    msg = ''
    if (theory /= 'dft') then
      msg = "theory must be 'dft'."
    else if (iperiodic /= 3) then
      msg = 'a periodic system (iperiodic=3) is required.'
    else if (yn_opt /= 'n' .or. yn_dc /= 'n' .or. yn_jm /= 'n') then
      msg = "yn_opt, yn_dc and yn_jm must be 'n'."
    else if (yn_spinorbit /= 'n') then
      msg = "yn_spinorbit must be 'n'."
    else if (system%nspin /= 1) then
      msg = 'nspin /= 1 is not supported.'
    else if (PLUS_U_ON) then
      msg = 'DFT+U (hubbard_param.dat) is not supported.'
    else if (system%xc_payload%use_tau_operator) then
      msg = 'the XC operator depending on the kinetic-energy density (meta-GGA) is not supported.'
    else if (trim(xc) == 'hse06') then
      msg = "xc='hse06' is not supported."
    else if (absorbing_boundary /= 'none') then
      msg = "absorbing_boundary /= 'none' is not supported."
    else if (any(info%nprgrid /= 1)) then
      msg = 'real-space parallelization (nproc_rgrid) is not allowed; use only nproc_k and nproc_ob.'
    else if (write_gs_restart_data /= 'all' .and. write_gs_restart_data /= 'wfn') then
      msg = "write_gs_restart_data must be 'all' or 'wfn'."
    else if (method_wf_distributor /= 'single') then
      msg = "method_wf_distributor must be 'single'."
    end if
#ifndef USE_MPI
    msg = 'a build with MPI support is required.'
#endif

    if (len_trim(msg) > 0) then
      if (comm_is_root(nproc_id_global)) then
        write(*,'(a)') 'Error: dm_unfold_option=reference (ground state of the reference cell): '//trim(msg)
      end if
      call end_parallel
      stop
    end if
  end subroutine check_dm_unfold_reference

  subroutine write_dm_unfold_reference(odir,iter,sum1,system,energy,lg,mg,stencil,V_local,ppg)
    use structures, only: s_rgrid, s_dft_system, s_dft_energy, s_stencil, s_scalar, s_pp_grid
    use salmon_global, only: convergence, threshold
    use communication, only: comm_is_root, comm_bcast, comm_sync_all
    use parallelization, only: nproc_id_global, nproc_group_global, end_parallel
    implicit none
    character(*),       intent(in) :: odir
    integer,            intent(in) :: iter
    real(8),            intent(in) :: sum1
    type(s_dft_system), intent(in) :: system
    type(s_dft_energy), intent(in) :: energy
    type(s_rgrid),      intent(in) :: lg, mg
    type(s_stencil),    intent(in) :: stencil
    type(s_scalar),     intent(in) :: V_local(:)
    type(s_pp_grid),    intent(in) :: ppg
    integer :: fh, ios, ios_close, ispin, ilma, ia, nps, nlma, nion
    real(8), allocatable :: uv_w(:,:)
    character(256) :: fname, msg
    ! the name and the format of the file: the same values in read_dm_unfold_reference
    character(*),  parameter :: dm_unfold_ref_file_name = 'dm_unfold_ref_hamiltonian.bin'
    character(16), parameter :: file_magic = 'SALMON_DMUNF_REF'
    integer,       parameter :: format_version = 1

    fname = trim(odir)//dm_unfold_ref_file_name
    nion = system%nion
    ios  = 0
    msg  = ''

    ! only the global root checks and writes (ppg and V_local are identical on all the processes, and the grid is
    ! not divided), and its status is shared by all the processes below
    if (comm_is_root(nproc_id_global)) then
      nps  = size(ppg%jxyz,2)
      nlma = ppg%nlma
      if (any(mg%is /= lg%is) .or. any(mg%ie /= lg%ie)) then
        ios = -1
        msg = 'internal error: the real-space grid is divided.'
      else if (size(V_local) /= system%nspin .or. size(ppg%jxyz,3) /= nion .or. size(ppg%uv,1) /= nps &
             & .or. size(ppg%uv,2) < nlma .or. size(ppg%ia_tbl) < nlma .or. size(ppg%rinv_uvu) < nlma) then
        ios = -1
        msg = 'internal error: inconsistent array sizes.'
      else if (system%xc_payload%use_tau_operator) then
        ios = -1
        msg = 'the XC operator depending on the kinetic-energy density is not stored in the file.'
      end if

      if (ios == 0) then
        ! uv is set only for j <= mps(ia); the unused tail is zeroed in the copy to be written
        allocate(uv_w(nps,nlma))
        uv_w = ppg%uv(:,1:nlma)
        do ilma = 1, nlma
          ia = ppg%ia_tbl(ilma)
          uv_w(ppg%mps(ia)+1:nps,ilma) = 0d0
        end do

        open(newunit=fh, file=trim(fname), form='unformatted', access='sequential', status='replace', action='write', &
           & iostat=ios)
        if (ios == 0) then
          write_records: do
            write(fh,iostat=ios) file_magic, format_version
            if (ios /= 0) exit write_records
            write(fh,iostat=ios) lg%num(1:3), system%nspin, nion, system%nk, system%no, lg%nd, nps, nlma
            if (ios /= 0) exit write_records
            write(fh,iostat=ios) system%hgs, system%hvol, system%primitive_a
            if (ios /= 0) exit write_records
            write(fh,iostat=ios) stencil%if_orthogonal, stencil%coef_lap0, stencil%coef_lap, stencil%coef_nab, &
                               & stencil%coef_f, stencil%coef_lap0_nd1, stencil%coef_lap_nd1, stencil%coef_nab_nd1
            if (ios /= 0) exit write_records
            write(fh,iostat=ios) system%kion(1:nion), system%Rion(1:3,1:nion)
            if (ios /= 0) exit write_records
            write(fh,iostat=ios) convergence, threshold, iter, sum1, energy%E_tot
            if (ios /= 0) exit write_records
            do ispin = 1, system%nspin
              write(fh,iostat=ios) V_local(ispin)%f(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
              if (ios /= 0) exit write_records
            end do
            write(fh,iostat=ios) ppg%mps(1:nion), ppg%jxyz, ppg%rxyz, ppg%ia_tbl(1:nlma), ppg%rinv_uvu(1:nlma), uv_w
            exit write_records
          end do write_records
          ! the status of close (flush of the last buffer) is also shared; a failure of write is kept as it is
          close(fh, iostat=ios_close)
          if (ios == 0) ios = ios_close
        end if
        if (ios /= 0) msg = 'failed to write '//trim(fname)
      end if
    end if

    call comm_bcast(ios, nproc_group_global)
    if (ios /= 0) then
      if (comm_is_root(nproc_id_global)) then
        if (len_trim(msg) == 0) msg = 'failed to write '//trim(fname)
        write(*,'(a,i0,a)') 'Error: dm_unfold reference Hamiltonian file (iostat=', ios, '): '//trim(msg)
      end if
      call end_parallel
      stop
    end if

    call comm_sync_all
    if (comm_is_root(nproc_id_global)) write(*,'(a)') '  dm_unfold reference Hamiltonian written: '//trim(fname)

  end subroutine write_dm_unfold_reference

  !---------------------------------------------------------------------------------------------------------------------
  ! restore side
  !---------------------------------------------------------------------------------------------------------------------

  ! Terminate the whole job. The restore routines are not necessarily called by all the processes at the same time
  ! (apply is called with a different q and number of vectors on each process), so MPI_Finalize (end_parallel) is not
  ! used: a process which stopped alone would leave the other processes waiting. The message is written by the process
  ! which calls this routine (nothing is written if msg is empty).
  subroutine abort_job(msg)
    use, intrinsic :: iso_fortran_env, only: output_unit
    use parallelization, only: nproc_id_global
#ifdef USE_MPI
    use mpi, only: MPI_COMM_WORLD
#endif
    implicit none
    character(*), intent(in) :: msg
    integer :: ierr, ios
    if (len_trim(msg) > 0) then
      write(*,'(a,i0,a)') 'Error: dm_unfold reference Hamiltonian (process ', nproc_id_global, '): '//trim(msg)
    end if
    flush(output_unit, iostat=ios) ! a failure of flush must not prevent the abort
#ifdef USE_MPI
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
#endif
    stop 1
  end subroutine abort_job

  ! ihead = [N(3), nspin, nion, nk, no, nd, nps, nlma] (record 2 of the file)
  subroutine set_dimensions(ref,ihead)
    implicit none
    type(s_dm_unfold_ref), intent(inout) :: ref
    integer,               intent(in)    :: ihead(10)
    ref%lg%num(1:3) = ihead(1:3)
    ref%lg%is(1:3)  = 1
    ref%lg%ie(1:3)  = ihead(1:3)
    ref%system%nspin = ihead(4)
    ref%system%nion  = ihead(5)
    ref%file_nk      = ihead(6)
    ref%file_no      = ihead(7)
    ref%ppg%nps      = ihead(9)
    ref%ppg%nlma     = ihead(10)
  end subroutine set_dimensions

  subroutine allocate_arrays(ref)
    implicit none
    type(s_dm_unfold_ref), intent(inout) :: ref
    integer :: ispin, nps, nlma, nion
    nps  = ref%ppg%nps
    nlma = ref%ppg%nlma
    nion = ref%system%nion
    allocate(ref%V_local(ref%system%nspin))
    do ispin = 1, ref%system%nspin
      allocate(ref%V_local(ispin)%f(ref%lg%is(1):ref%lg%ie(1),ref%lg%is(2):ref%lg%ie(2),ref%lg%is(3):ref%lg%ie(3)))
    end do
    allocate(ref%system%kion(nion), ref%system%Rion(3,nion))
    allocate(ref%ppg%mps(nion), ref%ppg%jxyz(3,nps,nion), ref%ppg%rxyz(3,nps,nion))
    allocate(ref%ppg%ia_tbl(nlma), ref%ppg%rinv_uvu(nlma), ref%ppg%uv(nps,nlma))
  end subroutine allocate_arrays

  ! Read the file. The process with rank 0 of icomm reads it and broadcasts the contents to the processes of icomm:
  ! all the processes of icomm have to call this routine.
  subroutine read_dm_unfold_reference(idir,icomm,ref)
    use communication, only: comm_is_root, comm_bcast, comm_get_groupinfo
    implicit none
    character(*),          intent(in)    :: idir ! directory (with the last '/') which contains the file
    integer,               intent(in)    :: icomm
    type(s_dm_unfold_ref), intent(inout) :: ref
    type(s_dm_unfold_ref) :: blank
    integer :: myrank, nprocs, fh, ios, ios_close, ispin, ihead(10), version, iter
    character(16)  :: magic, conv
    real(8)        :: thr, sum1, etot
    character(256) :: fname, msg
    integer, allocatable :: ibuf(:)
    ! the name and the format of the file: the same values in write_dm_unfold_reference
    character(*),  parameter :: dm_unfold_ref_file_name = 'dm_unfold_ref_hamiltonian.bin'
    character(16), parameter :: file_magic = 'SALMON_DMUNF_REF'
    integer,       parameter :: format_version = 1
    integer,       parameter :: nd_ref = 4 ! half width of the finite-difference stencil (fixed in hpsi); the same in build

#ifndef USE_MPI
    call abort_job('read: a build with MPI support is required.')
#endif
    ref = blank
    fname = trim(idir)//dm_unfold_ref_file_name
    call comm_get_groupinfo(icomm, myrank, nprocs)
    ios  = 0
    msg  = ''
    ihead = 0

    ! the root reads all the records; the status is shared before the contents are broadcast
    if (comm_is_root(myrank)) then
      open(newunit=fh, file=trim(fname), form='unformatted', access='sequential', status='old', action='read', iostat=ios)
      if (ios == 0) then
        read_records: do
          read(fh,iostat=ios) magic, version
          if (ios /= 0) exit read_records
          if (magic /= file_magic .or. version /= format_version) then
            ios = 1
            msg = 'unknown file (magic or format_version does not match).'
            exit read_records
          end if
          read(fh,iostat=ios) ihead
          if (ios /= 0) exit read_records
          if (any(ihead(1:3) < 1) .or. ihead(4) /= 1 .or. any(ihead(5:7) < 1) .or. ihead(8) /= nd_ref &
              & .or. ihead(9) < 0 .or. ihead(10) < 0) then
            ios = 1
            msg = 'invalid dimensions in the header (nspin /= 1 or nd /= 4 is not supported).'
            exit read_records
          end if
          call set_dimensions(ref,ihead)
          call allocate_arrays(ref)
          read(fh,iostat=ios) ref%system%hgs, ref%file_hvol, ref%system%primitive_a
          if (ios /= 0) exit read_records
          read(fh,iostat=ios) ref%stencil%if_orthogonal, ref%stencil%coef_lap0, ref%stencil%coef_lap, &
                            & ref%stencil%coef_nab, ref%stencil%coef_f, ref%stencil%coef_lap0_nd1, &
                            & ref%stencil%coef_lap_nd1, ref%stencil%coef_nab_nd1
          if (ios /= 0) exit read_records
          read(fh,iostat=ios) ref%system%kion, ref%system%Rion
          if (ios /= 0) exit read_records
          read(fh,iostat=ios) conv, thr, iter, sum1, etot ! informational record, not used
          if (ios /= 0) exit read_records
          do ispin = 1, ref%system%nspin
            read(fh,iostat=ios) ref%V_local(ispin)%f
            if (ios /= 0) exit read_records
          end do
          ! the last record is one record: it has to be read by a single read statement
          read(fh,iostat=ios) ref%ppg%mps, ref%ppg%jxyz, ref%ppg%rxyz, ref%ppg%ia_tbl, ref%ppg%rinv_uvu, ref%ppg%uv
          exit read_records
        end do read_records
        close(fh, iostat=ios_close)
        if (ios == 0) ios = ios_close
      end if
      if (ios /= 0 .and. len_trim(msg) == 0) msg = 'failed to read '//trim(fname)
    end if

    call comm_bcast(ios, icomm)
    if (ios /= 0) then
      if (comm_is_root(myrank)) then
        write(*,'(a,i0,a)') 'Error: dm_unfold reference Hamiltonian file (iostat=', ios, '): '//trim(msg)
      end if
      call abort_job('') ! the message has been written by the root of icomm
    end if

    ! dimensions first, then the allocation on the other processes, then the contents
    call comm_bcast(ihead, icomm)
    if (.not. comm_is_root(myrank)) then
      call set_dimensions(ref,ihead)
      call allocate_arrays(ref)
    end if
    call comm_bcast(ref%system%hgs, icomm)
    call comm_bcast(ref%file_hvol, icomm)
    call comm_bcast(ref%system%primitive_a, icomm)
    call comm_bcast(ref%stencil%if_orthogonal, icomm)
    call comm_bcast(ref%stencil%coef_lap0, icomm)
    call comm_bcast(ref%stencil%coef_lap, icomm)
    call comm_bcast(ref%stencil%coef_nab, icomm)
    call comm_bcast(ref%stencil%coef_f, icomm)
    call comm_bcast(ref%stencil%coef_lap0_nd1, icomm)
    call comm_bcast(ref%stencil%coef_lap_nd1, icomm)
    call comm_bcast(ref%stencil%coef_nab_nd1, icomm)
    call comm_bcast(ref%system%kion, icomm)
    call comm_bcast(ref%system%Rion, icomm)
    do ispin = 1, ref%system%nspin
      call comm_bcast(ref%V_local(ispin)%f, icomm)
    end do
    call comm_bcast(ref%ppg%mps, icomm)
    allocate(ibuf(size(ref%ppg%jxyz))) ! there is no bcast for an integer array of rank 3
    if (comm_is_root(myrank)) ibuf = reshape(ref%ppg%jxyz, [size(ref%ppg%jxyz)])
    call comm_bcast(ibuf, icomm)
    ref%ppg%jxyz = reshape(ibuf, shape(ref%ppg%jxyz))
    deallocate(ibuf)
    call comm_bcast(ref%ppg%rxyz, icomm)
    call comm_bcast(ref%ppg%ia_tbl, icomm)
    call comm_bcast(ref%ppg%rinv_uvu, icomm)
    call comm_bcast(ref%ppg%uv, icomm)

    ref%is_read = .true.
    if (comm_is_root(myrank)) write(*,'(a)') '  dm_unfold reference Hamiltonian read: '//trim(fname)
  end subroutine read_dm_unfold_reference

  ! Set up the structures which hpsi needs. The super cell is the cell of the current run:
  !   a_super(:,i) = nmul(i)*a_ref(:,i),  num_super(i) = nmul(i)*N_ref(i),  hgs_super = hgs_ref
  ! (nmul = 1 and the arguments of the reference cell itself when the file is used in the reference run).
  ! All the processes which called read_dm_unfold_reference have to call this routine.
  subroutine build_dm_unfold_reference(ref,primitive_a_super,num_super,hgs_super,nmul)
    use salmon_global, only: iperiodic, yn_periodic, yn_jm, yn_spinorbit, xc, absorbing_boundary, quiet
    use plusU_global, only: PLUS_U_ON
    use lattice, only: init_lattice
    use initialization_sub, only: init_grid_parallel
    use code_optimization, only: optimized_stencil_is_callable
#ifdef USE_MPI
    use mpi, only: MPI_COMM_SELF
#endif
    implicit none
    type(s_dm_unfold_ref), intent(inout) :: ref
    real(8),               intent(in)    :: primitive_a_super(3,3), hgs_super(3)
    integer,               intent(in)    :: num_super(3), nmul(3)
    type(s_stencil) :: st_chk
    integer :: i, comm_self
    real(8) :: tol
    logical :: quiet_save
    character(256) :: msg
    integer, parameter :: nd_ref = 4 ! half width of the finite-difference stencil (fixed in hpsi); the same in read

    if (.not. ref%is_read) call abort_job('build: the file has not been read.')
    if (ref%is_built) call abort_job('build: already built (call release_dm_unfold_reference before building again).')

    ! supported builds and current settings (they are used by hpsi/zpseudo as they are; none of them is changed here)
    msg = ''
#ifdef USE_OPENACC
    msg = 'a build with OpenACC is not supported.'
#endif
    if (len_trim(msg) == 0) then
      if (optimized_stencil_is_callable) then
        ! zstencil then uses the tables of the periodic indices (modx, mody, modz) made for the grid of this run,
        ! not those of the reference cell
        msg = 'the stencil with explicit vectorization is in use (yn_want_stencil_hand_vectorization).'
      else if (iperiodic /= 3 .or. yn_periodic /= 'y') then
        msg = "a periodic system (iperiodic=3, yn_periodic='y') is required."
      else if (yn_jm /= 'n' .or. yn_spinorbit /= 'n') then
        msg = "yn_jm and yn_spinorbit must be 'n'."
      else if (PLUS_U_ON) then
        msg = 'DFT+U (hubbard_param.dat) is not supported.'
      else if (trim(xc) == 'hse06') then
        msg = "xc='hse06' is not supported."
      else if (absorbing_boundary /= 'none') then
        msg = "absorbing_boundary /= 'none' is not supported."
      end if
    end if
    if (len_trim(msg) > 0) call abort_job('build: '//trim(msg))

    ! the cell of the file and the cell of the run: a_super = nmul*a_ref, N_super = nmul*N_ref, hgs_super = hgs_ref
    tol = 1d-8*maxval(abs(ref%system%primitive_a))
    do i = 1, 3
      if (nmul(i) < 1) then
        msg = 'nmul must be a positive integer.'
        exit
      else if (num_super(i) /= nmul(i)*ref%lg%num(i)) then
        msg = 'the number of grid points of the super cell is not nmul times that of the reference cell in the file.'
        exit
      else if (any(abs(primitive_a_super(:,i) - dble(nmul(i))*ref%system%primitive_a(:,i)) > tol)) then
        msg = 'the lattice vectors of the super cell are not nmul times those of the reference cell in the file.'
        exit
      else if (abs(hgs_super(i) - ref%system%hgs(i)) > 1d-8*ref%system%hgs(i)) then
        msg = 'the grid spacing of the super cell differs from that of the reference cell in the file.'
        exit
      end if
    end do
    if (len_trim(msg) > 0) call abort_job('build: '//trim(msg))

    ! lattice: coef_F and the volume element are derived again and compared with the file; the other stencil
    ! coefficients are those in the file
    ref%system%ngrid = product(ref%lg%num)
    call init_lattice(ref%system,st_chk)
    if (any(abs(st_chk%coef_f - ref%stencil%coef_f) > 1d-10*max(1d0,maxval(abs(st_chk%coef_f))))) then
      call abort_job('build: coef_F in the file is inconsistent with the lattice vectors.')
    end if
    if (abs(ref%system%hvol - ref%file_hvol) > 1d-10*ref%file_hvol) then
      call abort_job('build: hvol in the file is inconsistent with the lattice vectors and the grid.')
    end if

    ! one Bloch vector and one spin at a time; no vector potential
    ref%system%nk = 1
    ref%system%no = 1
    if (allocated(ref%system%vec_k)) deallocate(ref%system%vec_k)
    if (allocated(ref%system%wtk))   deallocate(ref%system%wtk)
    allocate(ref%system%vec_k(3,1), ref%system%wtk(1))
    ref%system%vec_k = 0d0
    ref%system%wtk   = 1d0
    ref%system%vec_Ac = 0d0
    ref%system%if_real_orbital = .false.
    ref%system%xc_payload%use_tau_operator = .false.

    ! grid: not divided. lg is set again from the physical region, so that the padding below does not accumulate.
    ref%lg%is_array(1:3) = ref%lg%is(1:3)
    ref%lg%ie_array(1:3) = ref%lg%ie(1:3)
    ref%lg%nd   = nd_ref
    ref%lg%ndir = 3

    ref%info%npk       = 1
    ref%info%nporbital = 1
    ref%info%nprgrid   = 1
    ref%info%id_rko    = 0
    ref%info%isize_rko = 1
    if (allocated(ref%info%imap)) deallocate(ref%info%imap)
    allocate(ref%info%imap(0:0,0:0,0:0,0:0,0:0))
    ref%info%imap = 0
    ref%info%if_divide_rspace = .false.
    ref%info%if_divide_orbit  = .false.
#ifdef USE_MPI
    comm_self = MPI_COMM_SELF
#else
    comm_self = 0
#endif
    ref%info%icomm_r   = comm_self ; ref%info%id_r   = 0 ; ref%info%isize_r   = 1
    ref%info%icomm_k   = comm_self ; ref%info%id_k   = 0 ; ref%info%isize_k   = 1
    ref%info%icomm_o   = comm_self ; ref%info%id_o   = 0 ; ref%info%isize_o   = 1
    ref%info%icomm_ro  = comm_self ; ref%info%id_ro  = 0 ; ref%info%isize_ro  = 1
    ref%info%icomm_ko  = comm_self ; ref%info%id_ko  = 0 ; ref%info%isize_ko  = 1
    ref%info%icomm_rko = comm_self
    ref%info%im_s = 1 ; ref%info%im_e = 1 ; ref%info%numm = 1
    ref%info%ik_s = 1 ; ref%info%ik_e = 1 ; ref%info%numk = 1
    ref%info%io_s = 1 ; ref%info%io_e = 1 ; ref%info%numo = 1

    ! init_grid_parallel writes a line from every process whose id_rko is 0 unless quiet. This is the only place where a
    ! global variable of salmon_global is changed: it is set only during this call and restored (the routine is not
    ! meant to be called by several threads of a process at the same time).
    quiet_save = quiet
    quiet = .true.
    call init_grid_parallel(ref%info,ref%lg,ref%mg)
    quiet = quiet_save

    ! V_local has no padding: it is read on the physical region, which is the region hpsi expects
    do i = 1, ref%system%nspin
      if (any(lbound(ref%V_local(i)%f) /= ref%mg%is) .or. any(ubound(ref%V_local(i)%f) /= ref%mg%ie)) then
        call abort_job('build: internal error, the region of V_local differs from that of the grid.')
      end if
    end do

    ! non-local pseudopotential: zekr_uV(j,ilma,1) is made for each q in apply (j > mps(ia) is not used)
    if (allocated(ref%ppg%zekr_uV)) deallocate(ref%ppg%zekr_uV)
    allocate(ref%ppg%zekr_uV(ref%ppg%nps,ref%ppg%nlma,1:1))
    ref%ppg%zekr_uV = (0d0,0d0)

    ref%no_work = 0
    ref%is_built = .true.
  end subroutine build_dm_unfold_reference

  ! H(q)*psi for no_work vectors. psi and hpsi_out are given on the physical region of the grid of the reference cell:
  ! (N1,N2,N3,no_work). q is the Cartesian Bloch vector (system%vec_k convention). There is no communication.
  subroutine apply_h_dm_unfold_reference(ref,q,no_work,psi,hpsi_out)
    use hamiltonian, only: hpsi, update_kvector_nonlocalpt
    implicit none
    type(s_dm_unfold_ref), intent(inout) :: ref
    real(8),               intent(in)    :: q(3)
    integer,               intent(in)    :: no_work
    complex(8),            intent(in)    :: psi(:,:,:,:)
    complex(8),            intent(out)   :: hpsi_out(:,:,:,:)
    type(s_sendrecv_grid) :: srg ! not used (the grid is not divided); it is only an argument of hpsi
    integer :: io, is(3), ie(3)

    if (.not. ref%is_built) call abort_job('apply: the Hamiltonian has not been built.')
    if (no_work < 1) call abort_job('apply: no_work must be positive.')
    is = ref%mg%is
    ie = ref%mg%ie
    if (any(shape(psi) /= [ref%mg%num(1),ref%mg%num(2),ref%mg%num(3),no_work]) .or. &
        any(shape(hpsi_out) /= [ref%mg%num(1),ref%mg%num(2),ref%mg%num(3),no_work])) then
      call abort_job('apply: the shape of psi or hpsi_out is not (N1,N2,N3,no_work).')
    end if

    ! work arrays and the dimensions which depend on the number of vectors
    call set_work_arrays(ref,no_work)

    ! k = q, A = 0 (k + A is not added twice), and the phase of the non-local projectors for this q
    ref%system%vec_k(1:3,1) = q(1:3)
    ref%system%vec_Ac = 0d0
    call update_kvector_nonlocalpt(1,1,ref%system,ref%ppg)

    do io = 1, no_work
      ref%tpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),1,io,1,1) = psi(:,:,:,io)
    end do

    call hpsi(ref%tpsi,ref%htpsi,ref%info,ref%mg,ref%V_local,ref%system,ref%stencil,srg,ref%ppg)

    do io = 1, no_work
      hpsi_out(:,:,:,io) = ref%htpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),1,io,1,1)
    end do
  end subroutine apply_h_dm_unfold_reference

  ! Matrix elements of the current for no_mat vectors given on the real-space grid of the reference cell, (N1,N2,N3,no_mat),
  ! at the Bloch vector q (Cartesian, the same convention as apply and solve; no vector potential is added). The definitions
  ! and the order of the indices are those of write_tm_data (the transition moments tm.bin of the ground state), with the first
  ! index of the matrix for the bra: for i, j = 1 .. no_mat and c = x, y, z,
  !   upu(c,i,j)           = -i <psi_i| d/dr_c |psi_j>
  !   u_rVnl_Vnlr_u(c,i,j) = -i ( <psi_i| e^{-iqr} r_c V_nl e^{iqr} |psi_j> - <psi_i| e^{-iqr} V_nl r_c e^{iqr} |psi_j> )
  ! (the volume element is included in both). psi_j are the periodic parts u_j of the Bloch orbitals, as the states of the
  ! ground state calculation and of solve. The gradient is the finite-difference one with the Cartesian transformation of the
  ! oblique cells (calc_gradient_psi), and r in the non-local term is the position of the grid point of the projector of the
  ! file (rxyz), as in write_tm_data. The matrix of dH(q)/dq_c, the velocity, between the vectors is upu + u_rVnl_Vnlr_u
  ! + q_c <psi_i|psi_j>, where the last term is q_c delta_ij only if the vectors are orthonormal (it is not in the returned
  ! matrices). The vectors do not have to be eigenstates or orthonormal, and the order of the bands does not matter.
  ! No communication; every process which calls it uses its own copy of the Hamiltonian. This routine does not use the work
  ! arrays of apply and solve (their contents and their dimensions are not changed), but changes system%vec_k (= q) and
  ! ppg%zekr_uV as apply does.
  subroutine calc_tm_dm_unfold_reference(ref,q,no_mat,psi,upu,u_rVnl_Vnlr_u)
    use hamiltonian, only: update_kvector_nonlocalpt
    use stencil_sub, only: calc_gradient_psi
    use math_constants, only: zI
    implicit none
    type(s_dm_unfold_ref), intent(inout) :: ref
    real(8),               intent(in)    :: q(3)
    integer,               intent(in)    :: no_mat
    complex(8),            intent(in)    :: psi(:,:,:,:)
    complex(8),            intent(out)   :: upu(:,:,:), u_rVnl_Vnlr_u(:,:,:)
    complex(8), allocatable :: wpsi(:,:,:), gpsi_c(:,:,:,:), uVpsi(:), uVrpsi(:,:)
    complex(8) :: veik, wrk(3)
    integer :: is(3), ie(3), i, j, c, ilma, ia, jj, ix, iy, iz
    real(8) :: x, y, z

    if (.not. ref%is_built) call abort_job('calc_tm: the Hamiltonian has not been built.')
    if (ref%system%nspin /= 1) call abort_job('calc_tm: only nspin = 1 is supported.')
    if (no_mat < 1) call abort_job('calc_tm: no_mat must be positive.')
    is = ref%mg%is
    ie = ref%mg%ie
    if (any(shape(psi) /= [ref%mg%num(1),ref%mg%num(2),ref%mg%num(3),no_mat])) then
      call abort_job('calc_tm: the shape of psi is not (N1,N2,N3,no_mat).')
    end if
    if (any(shape(upu) /= [3,no_mat,no_mat]) .or. any(shape(u_rVnl_Vnlr_u) /= [3,no_mat,no_mat])) then
      call abort_job('calc_tm: the shape of upu or u_rVnl_Vnlr_u is not (3,no_mat,no_mat).')
    end if

    ! -i <psi_i|nabla|psi_j>: the gradient of one vector at a time (as write_tm_data does)
    allocate(wpsi(ref%mg%is_array(1):ref%mg%ie_array(1),ref%mg%is_array(2):ref%mg%ie_array(2), &
                  ref%mg%is_array(3):ref%mg%ie_array(3)))
    allocate(gpsi_c(3,ref%mg%is_array(1):ref%mg%ie_array(1),ref%mg%is_array(2):ref%mg%ie_array(2), &
                  ref%mg%is_array(3):ref%mg%ie_array(3)))
    wpsi = (0d0,0d0)
    gpsi_c = (0d0,0d0)
    do j = 1, no_mat
      wpsi(is(1):ie(1),is(2):ie(2),is(3):ie(3)) = psi(:,:,:,j)
      call calc_gradient_psi(wpsi,gpsi_c,ref%mg%is_array,ref%mg%ie_array,ref%mg%is,ref%mg%ie, &
                             ref%mg%idx,ref%mg%idy,ref%mg%idz,ref%stencil%coef_nab,ref%system%rmatrix_B)
      !$omp parallel do private(i,c,wrk)
      do i = 1, no_mat
        do c = 1, 3
          wrk(c) = sum(conjg(psi(:,:,:,i)) * gpsi_c(c,is(1):ie(1),is(2):ie(2),is(3):ie(3)))
        end do
        upu(1:3,i,j) = -zI * wrk(1:3) * ref%system%hvol
      end do
      !$omp end parallel do
    end do
    deallocate(wpsi, gpsi_c)

    ! the non-local term for this q: zekr_uV = e^{-iqr} uV (update_kvector_nonlocalpt), vec_Ac = 0
    ref%system%vec_k(1:3,1) = q(1:3)
    ref%system%vec_Ac = 0d0
    call update_kvector_nonlocalpt(1,1,ref%system,ref%ppg)

    allocate(uVpsi(no_mat), uVrpsi(3,no_mat))
    u_rVnl_Vnlr_u = (0d0,0d0)
    do ilma = 1, ref%ppg%nlma
      ia = ref%ppg%ia_tbl(ilma)
      uVpsi  = (0d0,0d0)
      uVrpsi = (0d0,0d0)
      !$omp parallel do private(jj,x,y,z,ix,iy,iz,veik)
      do i = 1, no_mat
        do jj = 1, ref%ppg%mps(ia)
          x  = ref%ppg%rxyz(1,jj,ia)
          y  = ref%ppg%rxyz(2,jj,ia)
          z  = ref%ppg%rxyz(3,jj,ia)
          ix = ref%ppg%jxyz(1,jj,ia) - is(1) + 1 ! the index of psi
          iy = ref%ppg%jxyz(2,jj,ia) - is(2) + 1 ! the index of psi
          iz = ref%ppg%jxyz(3,jj,ia) - is(3) + 1 ! the index of psi
          veik = conjg(ref%ppg%zekr_uV(jj,ilma,1))
          uVpsi(i)    = uVpsi(i)    + veik     * psi(ix,iy,iz,i)  ! <v|e^{iqr}|u>
          uVrpsi(1,i) = uVrpsi(1,i) + veik * x * psi(ix,iy,iz,i)  ! <v|e^{iqr} x|u>
          uVrpsi(2,i) = uVrpsi(2,i) + veik * y * psi(ix,iy,iz,i)
          uVrpsi(3,i) = uVrpsi(3,i) + veik * z * psi(ix,iy,iz,i)
        end do
      end do
      !$omp end parallel do
      uVpsi = uVpsi * ref%ppg%rinv_uvu(ilma)
      !$omp parallel do private(i,j) collapse(2)
      do j = 1, no_mat
        do i = 1, no_mat
          ! -i ( <u_i|e^{-iqr} r|v><v|e^{iqr}|u_j> - <u_i|e^{-iqr}|v><v|e^{iqr} r|u_j> ) * hvol
          u_rVnl_Vnlr_u(1:3,i,j) = u_rVnl_Vnlr_u(1:3,i,j) &
            - zI * ( conjg(uVrpsi(1:3,i)) * uVpsi(j) - conjg(uVpsi(i)) * uVrpsi(1:3,j) ) * ref%system%hvol
        end do
      end do
      !$omp end parallel do
    end do
    deallocate(uVpsi, uVrpsi)
  end subroutine calc_tm_dm_unfold_reference

  ! Work arrays and the dimensions which depend on the number of vectors. Nothing is done if no_work is the same as before.
  subroutine set_work_arrays(ref,no_work)
    implicit none
    type(s_dm_unfold_ref), intent(inout) :: ref
    integer,               intent(in)    :: no_work

    if (no_work == ref%no_work) return
    if (allocated(ref%tpsi%zwf))  deallocate(ref%tpsi%zwf)
    if (allocated(ref%htpsi%zwf)) deallocate(ref%htpsi%zwf)
    allocate(ref%tpsi%zwf (ref%mg%is_array(1):ref%mg%ie_array(1),ref%mg%is_array(2):ref%mg%ie_array(2), &
                           ref%mg%is_array(3):ref%mg%ie_array(3),1,1:no_work,1:1,1:1))
    allocate(ref%htpsi%zwf(ref%mg%is_array(1):ref%mg%ie_array(1),ref%mg%is_array(2):ref%mg%ie_array(2), &
                           ref%mg%is_array(3):ref%mg%ie_array(3),1,1:no_work,1:1,1:1))
    ref%tpsi%zwf  = (0d0,0d0)
    ref%htpsi%zwf = (0d0,0d0)
    ref%tpsi%update_zwf_overlap = .true.
    ref%htpsi%update_zwf_overlap = .true.
    ! the work arrays of solve_dm_unfold_reference are made again (gscg_zwf allocates those of cg when they are not allocated)
    if (allocated(ref%gpsi%zwf))      deallocate(ref%gpsi%zwf)
    if (allocated(ref%cg%pk%zwf))     deallocate(ref%cg%pk%zwf)
    if (allocated(ref%cg%pre_gk%zwf)) deallocate(ref%cg%pre_gk%zwf)
    if (allocated(ref%cg%hpk%zwf))    deallocate(ref%cg%hpk%zwf)

    ref%system%no = no_work
    ref%info%io_s = 1
    ref%info%io_e = no_work
    ref%info%numo = no_work
    ! the tables of the orbital index which the subspace diagonalization and the Gram-Schmidt routines refer to,
    ! for one orbital process
    ref%info%numo_max = no_work
    if (allocated(ref%info%io_s_all)) deallocate(ref%info%io_s_all)
    if (allocated(ref%info%io_e_all)) deallocate(ref%info%io_e_all)
    if (allocated(ref%info%numo_all)) deallocate(ref%info%numo_all)
    if (allocated(ref%info%irank_io)) deallocate(ref%info%irank_io)
    allocate(ref%info%io_s_all(0:0), ref%info%io_e_all(0:0), ref%info%numo_all(0:0), ref%info%irank_io(1:no_work))
    ref%info%io_s_all(0) = 1
    ref%info%io_e_all(0) = no_work
    ref%info%numo_all(0) = no_work
    ref%info%irank_io    = 0
    ref%no_work = no_work
  end subroutine set_work_arrays

  ! no_work states of H(q) from the bottom of the spectrum: H(q) phi_n = e_n phi_n, <phi_n|phi_m> = delta_nm (the volume element
  ! included) for the lowest no_conv states. The caller uses only the states 1..no_conv as the eigenstates.
  !   phi (N1,N2,N3,no_work): the starting vectors on entry (linearly independent; they are orthonormalized first), the
  !                           eigenstates in the order of increasing e_n on exit
  !   q                     : the Cartesian Bloch vector (the same convention as apply), no vector potential is added
  !   no_conv               : the number of the lowest states which have to converge (1 <= no_conv <= no_work). The states
  !                           no_conv+1, .., no_work are the buffer states (no_buffer = no_work - no_conv of them): they are
  !                           iterated and returned, but their residuals do not enter the judgement and they are not guaranteed
  !                           to be eigenstates (the upper bands may cross the bands which are not in the set, and then they
  !                           converge slowly). The residual of the no_conv states does not show that no lower state of H(q)
  !                           has been missed. If the boundary between no_conv and no_conv+1 cuts a degenerate group, the
  !                           subspace of the lowest no_conv states is not unique: the caller chooses no_conv at a gap of the
  !                           spectrum (or takes the whole group), or knows that the output depends on the cut.
  !   tol_res               : the states are accepted when max_{n<=no_conv} res_n <= tol_res and orth <= tol_orth (orth: all
  !                           no_work states). tol_res is in Hartree (the residual is the norm of an energy). There is no lower
  !                           limit of tol_res. A tol_res which cannot be reached (for example 1d-20) makes the routine iterate
  !                           nsweep_max sweeps (in the tests without preconditioning the residual then stayed at the level of
  !                           1d-14 .. 1d-13). This relies on the correction of the conjugate gradient in conjugate_gradient.f90
  !                           (guard of the step: a band rejected by the guard is not updated any more in the call): without it
  !                           the residual of the converged states suddenly grew to 1d-5 .. 1d-1 after 20 .. 50 sweeps when the
  !                           iteration was continued.
  !   nsweep_max            : the maximum number of the sweeps (a sweep = conjugate gradient with ncg steps + Gram-Schmidt).
  !                           Every judgement is made after a subspace diagonalization, so with nsweep_max = 0 the starting
  !                           vectors are rotated to the eigenbasis of H(q) in their span and only evaluated
  !   ncg                   : the number of the steps of the conjugate gradient in a sweep (the caller gives the input ncg)
  !   eps, res              : e_n = <phi_n|H|phi_n> and res_n = ||(H - e_n) phi_n|| (the norm with the volume element), of the
  !                           returned states
  !   orth                  : max| <phi_n|phi_m> - delta_nm | of the returned states
  !   nsweep_done           : the number of the sweeps done
  !   converged             : whether the condition on tol_res and tol_orth is satisfied. The routine does not stop the job when
  !                           it is not, and the last states are returned in any case (as the iteration of the projection in
  !                           write.f90 does): the caller checks converged and res, and decides. A state which has not converged
  !                           is not an eigenstate of H(q) to the accuracy of tol_res; how much it is off depends on res and on
  !                           the gaps to the other states, not on res alone. If eps, res or orth contain NaN or Inf, the job is
  !                           aborted: such a state cannot be used as a basis, whatever the tol_res.
  ! The routine follows solve_orbitals in the ground state calculation (scf_iteration.f90) and the iteration of the projection in
  ! write.f90, with the judgement by the residual of the lowest no_conv states instead of the kinetic energy of the occupied
  ! states. The number of the applications of H (each for all the no_work states) is 2 + 3 nsweep_done + (the sum of the numbers
  ! of the steps the conjugate gradient actually did), which is 2 + nsweep_done*(ncg + 3) if every conjugate gradient runs
  ! ncg steps (subspace diagonalization and evaluation: 1 each, conjugate gradient: 1 + the steps).
  ! There is no communication outside of the process.
  subroutine solve_dm_unfold_reference(ref,q,no_work,no_conv,phi,tol_res,nsweep_max,ncg,eps,res,orth,nsweep_done,converged)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use hamiltonian, only: update_kvector_nonlocalpt
    use salmon_global, only: yn_scalapack
    use structures, only: allocate_orbital_complex
    use subspace_diagonalization, only: ssdg
    use Conjugate_Gradient, only: gscg_zwf
    use gram_schmidt_orth, only: gram_schmidt
    implicit none
    type(s_dm_unfold_ref), intent(inout) :: ref
    real(8),               intent(in)    :: q(3), tol_res
    integer,               intent(in)    :: no_work, no_conv, nsweep_max, ncg
    complex(8),            intent(inout) :: phi(:,:,:,:)
    real(8),               intent(out)   :: eps(no_work), res(no_work), orth
    integer,               intent(out)   :: nsweep_done
    logical,               intent(out)   :: converged
    type(s_sendrecv_grid) :: srg ! not used (the grid is not divided); it is only an argument
    integer :: io, iter, is(3), ie(3)
    real(8) :: rmax
    real(8), parameter :: tol_orth = 1d-10 ! allowed max| <phi_n|phi_m> - delta_nm | of the states returned

    if (.not. ref%is_built) call abort_job('solve: the Hamiltonian has not been built.')
    if (yn_scalapack == 'y') call abort_job("solve: yn_scalapack must be 'n'.")
    if (no_work < 1) call abort_job('solve: no_work must be positive.')
    if (no_conv < 1 .or. no_conv > no_work) call abort_job('solve: no_conv must satisfy 1 <= no_conv <= no_work.')
    if (tol_res <= 0d0 .or. nsweep_max < 0 .or. ncg < 0) then
      call abort_job('solve: tol_res must be positive, nsweep_max and ncg must not be negative.')
    end if
    is = ref%mg%is
    ie = ref%mg%ie
    if (any(shape(phi) /= [ref%mg%num(1),ref%mg%num(2),ref%mg%num(3),no_work])) then
      call abort_job('solve: the shape of phi is not (N1,N2,N3,no_work).')
    end if

    call set_work_arrays(ref,no_work)
    if (.not. allocated(ref%gpsi%zwf)) then
      call allocate_orbital_complex(1,ref%mg,ref%info,ref%gpsi)
      ref%gpsi%update_zwf_overlap = .true.
    end if

    ! k = q, A = 0, and the phase of the non-local projectors for this q (the same as apply)
    ref%system%vec_k(1:3,1) = q(1:3)
    ref%system%vec_Ac = 0d0
    call update_kvector_nonlocalpt(1,1,ref%system,ref%ppg)

    do io = 1, no_work
      ref%tpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),1,io,1,1) = phi(:,:,:,io)
    end do
    call gram_schmidt(ref%system,ref%mg,ref%info,ref%tpsi)

    nsweep_done = 0
    converged = .false.
    do iter = 0, nsweep_max
      call ssdg(ref%mg,ref%system,ref%info,ref%stencil,ref%tpsi,ref%htpsi,ref%ppg,ref%V_local,srg)
      call evaluate_bands(ref,srg,eps,res,orth)
      if (.not. (all(ieee_is_finite(eps)) .and. all(ieee_is_finite(res)) .and. ieee_is_finite(orth))) then
        call abort_job('solve: eps, res or orth is NaN or Inf.')
      end if
      rmax = maxval(res(1:no_conv))
      if (rmax <= tol_res .and. orth <= tol_orth) then
        converged = .true.
        exit
      end if
      if (iter == nsweep_max) exit
      call gscg_zwf(ncg,ref%mg,ref%system,ref%info,ref%stencil,ref%ppg,ref%V_local,srg,ref%tpsi,ref%htpsi,ref%gpsi,ref%cg)
      call gram_schmidt(ref%system,ref%mg,ref%info,ref%tpsi)
      nsweep_done = nsweep_done + 1
    end do

    do io = 1, no_work
      phi(:,:,:,io) = ref%tpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),1,io,1,1)
    end do
  end subroutine solve_dm_unfold_reference

  ! e_n, the residual and the orthonormality of the states in ref%tpsi (ref%htpsi is used as a work array)
  subroutine evaluate_bands(ref,srg,eps,res,orth)
    use hamiltonian, only: hpsi
    implicit none
    type(s_dm_unfold_ref), intent(inout) :: ref
    type(s_sendrecv_grid), intent(inout) :: srg
    real(8),               intent(out)   :: eps(:), res(:), orth
    complex(8), allocatable :: w(:,:), hw(:,:), s(:,:)
    integer :: nv, nr, io, jo, is(3), ie(3)
    real(8) :: hvol, nrm
    complex(8) :: zc

    nv = ref%no_work
    nr = product(ref%mg%num)
    is = ref%mg%is
    ie = ref%mg%ie
    hvol = ref%system%hvol
    call hpsi(ref%tpsi,ref%htpsi,ref%info,ref%mg,ref%V_local,ref%system,ref%stencil,srg,ref%ppg)

    allocate(w(nr,nv), hw(nr,nv), s(nv,nv))
    do io = 1, nv
      w(:,io)  = reshape(ref%tpsi%zwf (is(1):ie(1),is(2):ie(2),is(3):ie(3),1,io,1,1), [nr])
      hw(:,io) = reshape(ref%htpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),1,io,1,1), [nr])
    end do
    do io = 1, nv
      nrm = hvol*sum(abs(w(:,io))**2)
      zc  = hvol*sum(conjg(w(:,io))*hw(:,io))/nrm
      eps(io) = real(zc,8)
      res(io) = sqrt(hvol*sum(abs(hw(:,io) - eps(io)*w(:,io))**2)/nrm)
    end do
    call zgemm('C','N',nv,nv,nr,cmplx(hvol,0d0,8),w,nr,w,nr,(0d0,0d0),s,nv)
    orth = 0d0
    do jo = 1, nv
    do io = 1, nv
      if (io == jo) then
        orth = max(orth, abs(s(io,jo) - (1d0,0d0)))
      else
        orth = max(orth, abs(s(io,jo)))
      end if
    end do
    end do
    deallocate(w, hw, s)
  end subroutine evaluate_bands

  ! Free everything. read and build can be done again after this.
  subroutine release_dm_unfold_reference(ref)
    implicit none
    type(s_dm_unfold_ref), intent(inout) :: ref
    type(s_dm_unfold_ref) :: blank
    ref = blank ! the allocatable components of ref are deallocated by the assignment
  end subroutine release_dm_unfold_reference

end module dm_unfold_reference_sub
