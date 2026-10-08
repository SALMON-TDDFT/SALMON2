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

module dm_unfold_reference_sub
  implicit none
  private

  public :: check_dm_unfold_reference
  public :: write_dm_unfold_reference

  character(*), parameter, public :: dm_unfold_ref_file_name = 'dm_unfold_ref_hamiltonian.bin'
  character(16), parameter :: file_magic = 'SALMON_DMUNF_REF'
  integer, parameter :: format_version = 1

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

end module dm_unfold_reference_sub
