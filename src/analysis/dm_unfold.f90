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

module dm_unfold_sub
  implicit none

contains

  subroutine init_dm_unfold(lg,system,info,ofl,unfold)

  use structures
  use communication, only: comm_is_root, comm_bcast, comm_summation
  use parallelization, only: nproc_id_global, nproc_group_global, end_parallel
  use salmon_global, only: dm_unfold_option, no_ref, base_directory, sysname, natom, izatom, kion, &
                         & yn_out_mom_distr_gs, dq_mom, nq_mom, num_kgrid, &
                         & al_pr, al_vec1_pr, al_vec2_pr, al_vec3_pr
  use filesystem, only: open_filehandle
  use inputoutput, only: t_unit_time, t_unit_ac, t_unit_current
  use math_constants, only: zI,pi
  use lattice, only: calc_inverse
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  type(s_rgrid),           intent(in)    :: lg
  type(s_dft_system),      intent(inout) :: system
  type(s_parallel_info),   intent(in)    :: info
  type(s_ofile),           intent(out)   :: ofl
  type(s_unfold),          intent(inout) :: unfold
#ifdef USE_MPI
  character(256) :: iofile
  integer :: icomm
  integer :: gsize(7), lsize(7), lstart(7)
  integer :: gsize4(4), lsize4(4), lstart4(4)
  integer :: iopen_flag, minfo, mfile
  integer :: source_type, local_type, global_type, file_type
  integer :: ierr, n_count
  integer(kind=MPI_OFFSET_KIND) :: disp_upu, disp_vnl
  integer :: nspin,ispin,nsk_se,isk_s,isk_e,ie_ref(3),ihk,ir1,ir2,ir3
  integer :: ih1,ih2,ih3,ir1_ref,ir2_ref,ir3_ref,ig1_ref,ig2_ref,ig3_ref,io_ref,isk,ilk
  integer :: ig1,ig2,ig3,iqx,iqy,iqz,fp,iatom,j,ik
  real(8) :: omega_ref,B_ref(3,3),rsum,rsum_l,rj_l(3),rj(3),gx,gy,gz
  real(8) :: qx,qy,qz,dqx,dqy,dqz,value,nq_sum
  complex(8) :: zsum
  real(8),allocatable :: reta_uu(:,:,:,:),nq_l(:,:,:),nq_l_private(:,:,:)
  logical :: exists, e_occupation, e_wfn, e_tm
  integer :: i
  real(8) :: a_pr(3,3), ainv_pr(3,3), detA_pr, pmat_r(3,3), pmat_resid
  real(8) :: norm_pr(3), A_ref(3,3)
  integer :: pmat_i(3,3)

  ! -- primitive-to-reference translation-phase labeling (Phase A; unfolding.tex sec.9.5) --
  integer,parameter :: nhprk_max = 4096
  real(8),parameter :: hprk_thresh = 1d-2
  integer :: nhprk, ibox_c, n1c, n2c, n3c, nfound_c, jshift, icand, ibest_c, nbad_l, nbad
  integer :: nvec_pr(3,nhprk_max), hvec_pr(3,nhprk_max)
  real(8) :: fcoset_n(3,nhprk_max), fcoset_h(3,nhprk_max), fc3(3)
  real(8) :: pinv(3,3), pinvT(3,3), detP_r, tol_coset, pmat_r8(3,3)
  real(8) :: tc_cart(3,nhprk_max)
  complex(8),allocatable :: phi_pred(:,:), phase_gj(:,:,:,:)
  complex(8) :: cscore, phi_meas(nhprk_max)
  real(8) :: score_abs, csum_g2
  integer,allocatable :: hprk_label_l(:,:)
  real(8),allocatable :: hprk_score_l(:,:)
  logical :: found_dup

  if( dm_unfold_option /= 'super' ) then
    if (comm_is_root(nproc_id_global)) then
      write(*,"(A)") "dm_unfold_option /= 'super' at init_dm_unfold"
    end if
    call end_parallel
    stop
  end if

  if( system%nspin /= 1 ) then
    if (comm_is_root(nproc_id_global)) then
      write(*,"(A)") 'nspin /= 1 not allowed in dm_unfold calculations'
    end if
    call end_parallel
    stop
  end if

  if(comm_is_root(nproc_id_global)) then
    write(*,*) 'Entering init_dm_unfold'
  end if

  ! reference-cell lattice vectors: the 'super'-stage run's own cell
  ! (system%primitive_a) divided by unfold%num_hkgrid along each direction
  ! -- the same integer ratio already used to build B_ref/omega_ref below.
  ! This keeps al_pr/al_vec1-3_pr expressed at the natural (small) size of
  ! the reference cell, matching the 'reference'-stage run, rather than
  ! forcing the user to specify them at the much larger supercell scale.
  A_ref(:,1) = system%primitive_a(:,1) / dble(unfold%num_hkgrid(1))
  A_ref(:,2) = system%primitive_a(:,2) / dble(unfold%num_hkgrid(2))
  A_ref(:,3) = system%primitive_a(:,3) / dble(unfold%num_hkgrid(3))

  ! primitive-to-reference correspondence: build a^P from al_pr/al_vec1-3_pr
  ! (same orthogonal-vs-general convention as al/al_vec1-3 in init_dft_system),
  ! then compute the integer matrix P from a^R = A_ref (the reference cell,
  ! see above) and a^P. The note stores lattice vectors as matrix ROWS, but
  ! this code stores them as matrix COLUMNS (a_pr(1:3,j) = a^P_j), so here
  ! A_R = A_P * P and P = (A_P)^{-1} A_R. A diagonal/orthogonal test cell
  ! cannot distinguish this from the note's own row-convention formula; a
  ! non-diagonal example is needed.
  ! al_pr/al_vec1-3_pr are optional: only when ALL of them are exactly
  ! zero is this cell treated as its own primitive cell (a^P = a^R, P = I)
  ! -- this is a "no R->P decomposition performed" default, not a claim
  ! that the cell is physically primitive. Any other partial/incomplete
  ! specification is rejected as an error, never silently defaulted.
  if( sum(abs(al_pr))==0d0 .and. &
      sum(abs(al_vec1_pr))+sum(abs(al_vec2_pr))+sum(abs(al_vec3_pr))==0d0 ) then

    a_pr = A_ref
    pmat_i = 0
    pmat_i(1,1) = 1 ; pmat_i(2,2) = 1 ; pmat_i(3,3) = 1

  else

    if(al_vec1_pr(2)==0d0 .and. al_vec1_pr(3)==0d0 .and. al_vec2_pr(1)==0d0 .and. &
       al_vec2_pr(3)==0d0 .and. al_vec3_pr(1)==0d0 .and. al_vec3_pr(2)==0d0) then
      ! orthogonal case: accept al_pr, or else a diagonal-only al_vec1-3_pr
      ! (mirrors the al/al_vec1(1),al_vec2(2),al_vec3(3) fallback in
      ! init_dft_system) -- anything less complete than that is an error.
      if( al_pr(1)*al_pr(2)*al_pr(3) /= 0d0 ) then
        a_pr = 0d0
        a_pr(1,1) = al_pr(1)
        a_pr(2,2) = al_pr(2)
        a_pr(3,3) = al_pr(3)
      else if( al_vec1_pr(1)*al_vec2_pr(2)*al_vec3_pr(3) /= 0d0 ) then
        a_pr = 0d0
        a_pr(1,1) = al_vec1_pr(1)
        a_pr(2,2) = al_vec2_pr(2)
        a_pr(3,3) = al_vec3_pr(3)
      else
        stop 'incomplete primitive-cell lattice vectors (al_pr/al_vec1-3_pr) in dm_unfold, super'
      end if
    else
      a_pr(1:3,1) = al_vec1_pr
      a_pr(1:3,2) = al_vec2_pr
      a_pr(1:3,3) = al_vec3_pr
    end if

    ! reject a zero-length vector, and a degenerate/near-singular cell,
    ! BEFORE calc_inverse (which divides by det unconditionally). The
    ! degeneracy test is the dimensionless |det(A_P)|/(|a1||a2||a3|),
    ! not a bare |det(A_P)| threshold, so it does not depend on the
    ! overall size of the primitive cell.
    norm_pr(1) = sqrt(sum(a_pr(1:3,1)**2))
    norm_pr(2) = sqrt(sum(a_pr(1:3,2)**2))
    norm_pr(3) = sqrt(sum(a_pr(1:3,3)**2))
    if( any(norm_pr(:) < 1d-12) ) then
      stop 'zero-length primitive-cell lattice vector (al_pr/al_vec1-3_pr) in dm_unfold, super'
    end if

    detA_pr = a_pr(1,1)*a_pr(2,2)*a_pr(3,3) + a_pr(2,1)*a_pr(3,2)*a_pr(1,3) + a_pr(3,1)*a_pr(1,2)*a_pr(2,3) &
            - a_pr(1,3)*a_pr(2,2)*a_pr(3,1) - a_pr(2,3)*a_pr(3,2)*a_pr(1,1) - a_pr(3,3)*a_pr(1,2)*a_pr(2,1)
    if( abs(detA_pr)/(norm_pr(1)*norm_pr(2)*norm_pr(3)) < 1d-8 ) then
      stop 'degenerate primitive-cell lattice vectors (al_pr/al_vec1-3_pr) in dm_unfold, super'
    end if

    call calc_inverse(a_pr, ainv_pr, detA_pr)
    pmat_r = matmul(ainv_pr, A_ref)
    pmat_i = nint(pmat_r)
    pmat_resid = maxval(abs(pmat_r - dble(pmat_i)))
    if( pmat_resid > 1d-6 ) then
      stop 'reference-cell lattice vectors (derived from al/al_vec1-3 and num_skgrid/num_kgrid) are not &
        &an integer multiple of al_pr/al_vec1-3_pr (P is not an integer matrix) in dm_unfold, super'
    end if

  end if

  unfold%a_pr = a_pr
  unfold%pmat = pmat_i

  if (comm_is_root(nproc_id_global)) then
    write(*,"(A)") 'dm_unfold_option=super: primitive-to-reference matrix P (A_R = A_P * P, a^R_i = sum_j P(j,i) a^P_j):'
    do i = 1,3
      write(*,"(3I6)") pmat_i(i,1:3)
    end do
  end if

  if (comm_is_root(nproc_id_global)) then
    inquire(file='reference/wfn.bin', exist = e_wfn)
    inquire(file='reference/occupation.bin', exist = e_occupation)
    inquire(file='reference/tm.bin', exist = e_tm)
  end if
  exists = e_wfn .and. e_occupation .and. e_tm
  call comm_bcast(exists, nproc_group_global)
  if( .not. exists ) then
    if (comm_is_root(nproc_id_global)) then
      write(*,"(A)") 'Error: file not found, reference/wfn.bin, occupation.bin, tm.bin'
    end if
    call end_parallel
    stop
  end if

  nspin = 1
  isk_s = (info%ik_s-1) * unfold%nhrsk + 1
  isk_e = info%ik_e * unfold%nhrsk
  nsk_se = isk_e - isk_s + 1

  if( any(mod(lg%ie(:), unfold%num_hkgrid(:)) /= 0) ) then
    if (comm_is_root(nproc_id_global)) then
      write(*,"(A)") 'ie_ref(1:3) error at init_dm_unfold'
    end if
    call end_parallel
    stop
  end if

  ie_ref(1:3)=lg%ie(1:3)/unfold%num_hkgrid(1:3) ! size of reference cell grid

  allocate( unfold%psi_ref(1:ie_ref(1),1:ie_ref(2),1:ie_ref(3),nspin,no_ref,isk_s:isk_e,1) )

  source_type = MPI_DOUBLE_COMPLEX
  minfo = MPI_INFO_NULL
  iopen_flag = MPI_MODE_RDONLY
  iofile = "reference/wfn.bin"
  icomm = info%icomm_k

! create MPI_Type (Window) of process-local wave function
  gsize  = [ie_ref(1:3), nspin, no_ref, nsk_se, 1]
  lsize  = [ie_ref(1:3), nspin, no_ref, nsk_se, 1]
  lstart = [1,1,1,      1,     1,     1,      1] - 1

  call MPI_Type_create_subarray(7, gsize, lsize, lstart, MPI_ORDER_FORTRAN, source_type, local_type, ierr)
  call MPI_Type_commit(local_type, ierr)

! create MPI_Type (Window) of global wave function
  gsize  = [ie_ref(1:3), nspin, no_ref, unfold%nsk,     1]
  lstart = [1,1,1,      1,     1,     isk_s,          1] - 1

  call MPI_Type_create_subarray(7, gsize, lsize, lstart, MPI_ORDER_FORTRAN, source_type, global_type, ierr)
  call MPI_Type_commit(global_type, ierr)

  call MPI_File_open(icomm, iofile, iopen_flag, minfo, mfile, ierr)
  call MPI_File_set_view(mfile, 0_MPI_OFFSET_KIND, local_type, global_type, 'native', MPI_INFO_NULL, ierr)

  call MPI_File_read_all(mfile, unfold%psi_ref, 1, local_type, MPI_STATUS_IGNORE, ierr)

  call MPI_File_close(mfile, ierr)
  call MPI_Type_free(global_type, ierr)
  call MPI_Type_free( local_type, ierr)

  if(comm_is_root(nproc_id_global)) then
    write(*,*) 'End reading reference/wfn.bin'
  end if

! read occupation

  allocate( unfold%rocc_ref(no_ref, unfold%nsk, nspin) )

  if(comm_is_root(nproc_id_global)) then
     iofile = "reference/occupation.bin"
     open(888,file=iofile,form='unformatted')
     read(888) unfold%rocc_ref(1:no_ref,1:unfold%nsk,1:nspin)
     close(888)
  end if
  call comm_bcast(unfold%rocc_ref,nproc_group_global)

  if(comm_is_root(nproc_id_global))then
    write(*,*) 'End reading reference/occupation.bin'
  end if

! read transition dipole matrix elements

    allocate( unfold%upu_ref(3, no_ref, no_ref, isk_s:isk_e), unfold%u_rVnl_Vnlr_u_ref(3, no_ref, no_ref, isk_s:isk_e) )

    iofile = "reference/tm.bin"
    icomm = info%icomm_k

    gsize4  = [3, no_ref , no_ref , unfold%nsk       ]
    lsize4  = [3, no_ref , no_ref , isk_e - isk_s + 1]
    lstart4 = [1, 1     , 1     , isk_s            ] - 1

    n_count = lsize4(1)*lsize4(2)*lsize4(3)*lsize4(4)

    call MPI_Type_create_subarray(4, gsize4, lsize4, lstart4, MPI_ORDER_FORTRAN, source_type, file_type, ierr)
    call MPI_Type_commit(file_type, ierr)
    call MPI_File_open(icomm, iofile, iopen_flag, minfo, mfile, ierr)

    disp_upu = 0_MPI_OFFSET_KIND
    disp_vnl = int(3,8)*no_ref*no_ref*unfold%nsk*int(16,8)

    call MPI_File_set_view(mfile, disp_upu, source_type, file_type, 'native', minfo, ierr)
    call MPI_File_read_all(mfile, unfold%upu_ref(1,1,1,isk_s), n_count, source_type, MPI_STATUS_IGNORE, ierr)

    call MPI_File_set_view(mfile, disp_vnl, source_type, file_type, 'native', minfo, ierr)
    call MPI_File_read_all(mfile, unfold%u_rVnl_Vnlr_u_ref(1,1,1,isk_s), n_count, source_type, MPI_STATUS_IGNORE, ierr)

    call MPI_File_close(mfile, ierr) 
    call MPI_Type_free(file_type, ierr)

  if(comm_is_root(nproc_id_global))then
    write(*,*) 'End reading reference/tm.bin'
  end if

! exp(i hat_k r) table
  allocate( unfold%eihkr_tbl(lg%ie(1),lg%ie(2),lg%ie(3),unfold%nhrsk) )
  !$omp parallel do private(ihk,ih1,ih2,ih3,ir1_ref,ir2_ref,ir3_ref,ir1,ir2,ir3) collapse(4)
  do ihk = 1, unfold%nhrsk
  do ih1 = 1, unfold%num_hkgrid(1)
  do ih2 = 1, unfold%num_hkgrid(2)
  do ih3 = 1, unfold%num_hkgrid(3)
  do ir1_ref = 1, ie_ref(1)
  do ir2_ref = 1, ie_ref(2)
  do ir3_ref = 1, ie_ref(3)
    ir1 = ir1_ref + (ih1-1) * ie_ref(1)
    ir2 = ir2_ref + (ih2-1) * ie_ref(2)
    ir3 = ir3_ref + (ih3-1) * ie_ref(3)
    unfold%eihkr_tbl(ir1,ir2,ir3,ihk) = exp( zI * (unfold%vec_hrsk(1,ihk)*(ir1-1)*system%hgs(1) &
    &  + unfold%vec_hrsk(2,ihk)*(ir2-1)*system%hgs(2) + unfold%vec_hrsk(3,ihk)*(ir3-1)*system%hgs(3) ) )
  end do
  end do
  end do
  end do
  end do
  end do
  end do

  if(comm_is_root(nproc_id_global)) then
    write(ofl%file_dm_unfold,"(2A,'_dm_unfold.data')") trim(base_directory),trim(SYSname)
    ofl%fh_dm_unfold = open_filehandle(ofl%file_dm_unfold)
    open(ofl%fh_dm_unfold,file=ofl%file_dm_unfold)
        write(ofl%fh_dm_unfold, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1,  "Time", trim(t_unit_time%name), &
        & 2,  "Ac_tot_x", trim(t_unit_ac%name), &
        & 3,  "Ac_tot_y", trim(t_unit_ac%name), &
        & 4,  "Ac_tot_z", trim(t_unit_ac%name), &
        & 5,  "Tr[rho(t)]", "none", &
        & 6,  "Jx:<upu>-d", trim(t_unit_current%name), &
        & 7,  "Jy:<upu>-d", trim(t_unit_current%name), &
        & 8,  "Jz:<upu>-d", trim(t_unit_current%name), &
        & 9,  "Jx:<u[r,V]u>-d", trim(t_unit_current%name), &
        & 10, "Jy:<u[r,V]u>-d", trim(t_unit_current%name), &
        & 11, "Jz:<u[r,V]u>-d", trim(t_unit_current%name), &
        & 12, "Jx:<upu>-nd", trim(t_unit_current%name), &
        & 13, "Jy:<upu>-nd", trim(t_unit_current%name), &
        & 14, "Jz:<upu>-nd", trim(t_unit_current%name), &
        & 15, "Jx:<u[r,V]u>-nd", trim(t_unit_current%name), &
        & 16, "Jy:<u[r,V]u>-nd", trim(t_unit_current%name), &
        & 17, "Jz:<u[r,V]u>-nd", trim(t_unit_current%name), &
        & 18, "Jx:k", trim(t_unit_current%name), &
        & 19, "Jy:k", trim(t_unit_current%name), &
        & 20, "Jz:k", trim(t_unit_current%name), &
        & 21, "Jx:A", trim(t_unit_current%name), &
        & 22, "Jy:A", trim(t_unit_current%name), &
        & 23, "Jz:A", trim(t_unit_current%name), &
        & 24, "Jx:rho woV", trim(t_unit_current%name), &
        & 25, "Jy:rho woV", trim(t_unit_current%name), &
        & 26, "Jz:rho woV", trim(t_unit_current%name), &
        & 27, "Jx:rho", trim(t_unit_current%name), &
        & 28, "Jy:rho", trim(t_unit_current%name), &
        & 29, "Jz:rho", trim(t_unit_current%name), &
        & 30, "rho(t)|u(G)|2", "none", &
        & 31, "Jx:uG2-d", trim(t_unit_current%name), &
        & 32, "Jy:uG2-d", trim(t_unit_current%name), &
        & 33, "Jz:uG2-d", trim(t_unit_current%name), &
        & 34, "Jx:uG2-nd", trim(t_unit_current%name), &
        & 35, "Jy:uG2-nd", trim(t_unit_current%name), &
        & 36, "Jz:uG2-nd", trim(t_unit_current%name), &
        & 37, "Jx:uG2k-d", trim(t_unit_current%name), &
        & 38, "Jy:uG2k-d", trim(t_unit_current%name), &
        & 39, "Jz:uG2k-d", trim(t_unit_current%name), &
        & 40, "Jx:uG2A-d", trim(t_unit_current%name), &
        & 41, "Jy:uG2A-d", trim(t_unit_current%name), &
        & 42, "Jz:uG2A-d", trim(t_unit_current%name), &
        & 43, "Jx:uG2", trim(t_unit_current%name), &
        & 44, "Jy:uG2", trim(t_unit_current%name), &
        & 45, "Jz:uG2", trim(t_unit_current%name), &
        & 46, "int nq(q)", "none", &
        & 47, "Jx:nq-d", trim(t_unit_current%name), &
        & 48, "Jy:nq-d", trim(t_unit_current%name), &
        & 49, "Jz:nq-d", trim(t_unit_current%name), &
        & 50, "Jx:nq-nd", trim(t_unit_current%name), &
        & 51, "Jy:nq-nd", trim(t_unit_current%name), &
        & 52, "Jz:nq-nd", trim(t_unit_current%name), &
        & 53, "Jx:nq", trim(t_unit_current%name), &
        & 54, "Jy:nq", trim(t_unit_current%name), &
        & 55, "Jz:nq", trim(t_unit_current%name)

  end if

! Fourier transform of reference cell Bloch orbital
  allocate( unfold%psi_refG(1:ie_ref(1),1:ie_ref(2),1:ie_ref(3),nspin,no_ref,isk_s:isk_e,1) )

  ispin = 1
  omega_ref = system%hvol * system%ngrid / dble(unfold%num_hkgrid(1)*unfold%num_hkgrid(2)*unfold%num_hkgrid(3))

  !$omp parallel do private(ilk,ihk,isk,io_ref,ig1_ref,ig2_ref,ig3_ref,zsum,ir1_ref,ir2_ref,ir3_ref) collapse(2)
  do ilk = info%ik_s, info%ik_e ! large k
  do ihk = 1, unfold%nhrsk   ! hat k
    isk = unfold%isk_tbl(ilk,ihk) !small k = large k + hat k
  do io_ref = 1, no_ref
  do ig1_ref = 1, ie_ref(1)
  do ig2_ref = 1, ie_ref(2)
  do ig3_ref = 1, ie_ref(3)
    zsum = 0d0
    do ir1_ref = 1, ie_ref(1)
    do ir2_ref = 1, ie_ref(2)
    do ir3_ref = 1, ie_ref(3)
      zsum = zsum + unfold%psi_ref(ir1_ref,ir2_ref,ir3_ref,ispin,io_ref,isk,1) &
      & * exp( -2*pi*zI*((ig1_ref-1)*(ir1_ref-1)/dble(ie_ref(1))+(ig2_ref-1)*(ir2_ref-1)/dble(ie_ref(2)) &
      &                 +(ig3_ref-1)*(ir3_ref-1)/dble(ie_ref(3))) )
    end do
    end do
    end do
    unfold%psi_refG(ig1_ref,ig2_ref,ig3_ref,ispin,io_ref,isk,1) = zsum *system%hvol / omega_ref
  end do
  end do
  end do
  end do

  end do
  end do

  allocate( reta_uu(1:ie_ref(1),1:ie_ref(2),1:ie_ref(3),isk_s:isk_e) )
!$omp parallel do private(ilk,ihk,isk,ig1_ref,ig2_ref,ig3_ref,io_ref,rsum) collapse(2)
  do ilk = info%ik_s, info%ik_e
  do ihk = 1, unfold%nhrsk
    isk = unfold%isk_tbl(ilk,ihk)
  do ig1_ref = 1, ie_ref(1)
  do ig2_ref = 1, ie_ref(2)
  do ig3_ref = 1, ie_ref(3)
    rsum = 0d0
    do io_ref = 1, no_ref
      rsum = rsum + unfold%rocc_ref(io_ref,isk,ispin) * abs(unfold%psi_refG(ig1_ref,ig2_ref,ig3_ref,ispin,io_ref,isk,1))**2
    end do
    reta_uu(ig1_ref,ig2_ref,ig3_ref,isk) = rsum * omega_ref
  enddo
  enddo
  enddo

  enddo
  enddo

  B_ref(:,:) = system%primitive_b(:,:)
  B_ref(:,1) = B_ref(:,1) * unfold%num_hkgrid(1)
  B_ref(:,2) = B_ref(:,2) * unfold%num_hkgrid(2)
  B_ref(:,3) = B_ref(:,3) * unfold%num_hkgrid(3)
  rsum_l = 0d0
  rj_l(1:3) = 0d0
!$omp parallel do private(ilk,ihk,isk,ig1_ref,ig2_ref,ig3_ref,ig1,ig2,ig3,gx,gy,gz) reduction(+:rsum_l,rj_l) collapse(2)
  do ilk = info%ik_s, info%ik_e
  do ihk = 1, unfold%nhrsk
    isk = unfold%isk_tbl(ilk,ihk)
  do ig1_ref = 1, ie_ref(1)
  do ig2_ref = 1, ie_ref(2)
  do ig3_ref = 1, ie_ref(3)
    ig1 = ig1_ref - 1
    ig2 = ig2_ref - 1
    ig3 = ig3_ref - 1
    if( ig1 > ie_ref(1)/2 ) ig1 = ig1 - ie_ref(1)
    if( ig2 > ie_ref(2)/2 ) ig2 = ig2 - ie_ref(2)
    if( ig3 > ie_ref(3)/2 ) ig3 = ig3 - ie_ref(3)
    gx = ig1*B_ref(1,1) + ig2*B_ref(1,2) + ig3*B_ref(1,3)
    gy = ig1*B_ref(2,1) + ig2*B_ref(2,2) + ig3*B_ref(2,3)
    gz = ig1*B_ref(3,1) + ig2*B_ref(3,2) + ig3*B_ref(3,3)
    rsum_l = rsum_l + unfold%wtk_ref(isk) * reta_uu(ig1_ref,ig2_ref,ig3_ref,isk)
    rj_l(1) = rj_l(1) + unfold%wtk_ref(isk) * reta_uu(ig1_ref,ig2_ref,ig3_ref,isk) &
    & * (gx + system%vec_k(1,ilk) + unfold%vec_hrsk(1,ihk))
    rj_l(2) = rj_l(2) + unfold%wtk_ref(isk) * reta_uu(ig1_ref,ig2_ref,ig3_ref,isk) &
    & * (gy + system%vec_k(2,ilk) + unfold%vec_hrsk(2,ihk))
    rj_l(3) = rj_l(3) + unfold%wtk_ref(isk) * reta_uu(ig1_ref,ig2_ref,ig3_ref,isk) &
    & * (gz + system%vec_k(3,ilk) + unfold%vec_hrsk(3,ihk))
  end do
  end do
  end do

  end do
  end do

  rsum = 0d0
  rj(:) = 0d0
  call comm_summation(rsum_l,rsum,icomm)
  call comm_summation(rj_l,rj,3,icomm)

  rj = rj / omega_ref

  if(comm_is_root(nproc_id_global))then
    write(*,'(A,7x,2f17.12)') 'N:sum rho_uu(k,G)       ',rsum
    write(*,'(A,7x,6f17.12)') 'J:sum rho_uu(k,G)(G+k)  ',rj(1:3)
  end if

! ===================================================================
! Primitive-to-reference translation-phase labeling (Phase A).
! See unfolding.tex sec.9.5 (eq:tc-overlap-G) for the underlying
! algorithm; only the labeling step is implemented here -- the
! energy-eigenbasis-recovery step for same-hat_k mixing (sec.9.5,
! "Recovering the energy eigenbasis within a shared-hat_k block") is
! a separate, later step and is NOT performed by this block.
!
! unfold%nhprk = |det(pmat)|, the number of candidate primitive-cell
! unfolding vectors hat_k per reference-cell band. Both the candidates
! vec_hprk(:,c) = B_ref*h_c (h_c integer, c=1..nhprk) and the
! coset-representative real-space shifts t_c(j) = a_pr*n_j (n_j
! integer, j=1..nhprk) enumerate the same finite abelian group
! Z^3/(pmat)Z^3 (order nhprk): n,n' give the same real-space coset
! (mod the reference lattice) iff pinv*n == pinv*n' (mod 1,
! componentwise, pinv = pmat^{-1}); h,h' give the same reciprocal
! coset (mod the primitive reciprocal lattice) iff
! transpose(pinv)*h == transpose(pinv)*h' (mod 1). With this,
! hat_k(h_c).t_c(n_j) = 2*pi * dot_product(transpose(pinv)*h_c, n_j)
! exactly (unfolding.tex sec.9.5 derivation), so the nhprk x nhprk
! predicted-phase table phi_pred(c,j) = exp(i*hat_k(h_c).t_c(n_j)) is
! built from integers and pinv alone -- no Cartesian hat_k/t_c dot
! products are needed for it (Cartesian vectors are only needed for
! vec_hprk itself and for the measured side, via the G-grid, below).
!
! Enumeration uses an adaptive search box (doubled until nhprk
! distinct cosets are found, via floor()-based wrapping to [0,1) --
! NOT round()-based, which is ambiguous exactly on the half-integer
! coset values that occur whenever |det(pmat)| is even) rather than a
! fixed a priori bound, since a safe bound in terms of pmat's entries
! is not simple to state for a general integer matrix; nhprk_max is a
! generous sanity cap, not a physical limit.
!
! Per band (io_ref,isk), the measured phase for shift j is the
! G-space overlap (eq:tc-overlap-G):
!   phi_meas(j) = sum_g |c_g|^2 exp(i G_g.t_c(j)) / sum_g |c_g|^2,
! c_g = psi_refG(...,io_ref,isk,1); this is exact for the finite
! trigonometric interpolant psi_refG defines (see unfolding.tex
! sec.9.5 for the scope of this exactness statement). The band's
! best-matching candidate and match quality come from
!   score(c) = (1/nhprk) * sum_j phi_meas(j) * conjg(phi_pred(c,j)),
! using character orthogonality of the nhprk phi_pred(c,:) rows
! (distinct candidates are exactly orthogonal over the full coset
! set): |score(c)| <= 1 always, with equality only when phi_meas
! matches phi_pred(c,:) exactly. hprk_label is the score-maximizing
! c; when 1-|score| exceeds hprk_thresh the job does NOT stop -- it
! warns, stores the sentinel label 0, and continues (per the
! 2026-09-18 design discussion: a runtime degeneracy/near-degeneracy
! must never halt the job).
! ===================================================================

  pmat_r8 = dble(pmat_i)
  call calc_inverse(pmat_r8, pinv, detP_r)
  nhprk = nint(abs(detP_r))
  if( nhprk < 1 ) then
    if (comm_is_root(nproc_id_global)) then
      write(*,"(A)") 'Error: det(pmat) rounds to < 1 in dm_unfold primitive-to-reference labeling'
    end if
    call end_parallel
    stop
  end if
  if( nhprk > nhprk_max ) then
    if (comm_is_root(nproc_id_global)) then
      write(*,"(A,I0)") 'Error: |det(pmat)| exceeds the coset-enumeration cap in dm_unfold, nhprk=', nhprk
    end if
    call end_parallel
    stop
  end if
  unfold%nhprk = nhprk
  pinvT = transpose(pinv)
  tol_coset = 1d-6

! -- enumerate nhprk distinct real-space cosets (n -> t_c(n) = a_pr*n) --
  ibox_c = nhprk
  do
    nfound_c = 0
    coset_n_search: do n1c = 0, ibox_c-1
    do n2c = 0, ibox_c-1
    do n3c = 0, ibox_c-1
      fc3(:) = pinv(:,1)*n1c + pinv(:,2)*n2c + pinv(:,3)*n3c
      fc3(:) = fc3(:) - dble(floor(fc3(:)))
      found_dup = .false.
      do j = 1, nfound_c
        if( sum(abs(fc3(:)-fcoset_n(:,j))) < tol_coset ) then
          found_dup = .true.
          exit
        end if
      end do
      if( .not. found_dup ) then
        nfound_c = nfound_c + 1
        nvec_pr(1,nfound_c) = n1c
        nvec_pr(2,nfound_c) = n2c
        nvec_pr(3,nfound_c) = n3c
        fcoset_n(:,nfound_c) = fc3(:)
        if( nfound_c >= nhprk ) exit coset_n_search
      end if
    end do
    end do
    end do coset_n_search
    if( nfound_c >= nhprk ) exit
    ibox_c = ibox_c * 2
    if( ibox_c > 64*nhprk+64 ) then
      if (comm_is_root(nproc_id_global)) then
        write(*,"(A)") 'Error: could not enumerate a full set of real-space cosets in dm_unfold primitive-to-reference labeling'
      end if
      call end_parallel
      stop
    end if
  end do

! -- enumerate nhprk distinct reciprocal cosets (h -> hat_k(h) = B_ref*h) --
  ibox_c = nhprk
  do
    nfound_c = 0
    coset_h_search: do n1c = 0, ibox_c-1
    do n2c = 0, ibox_c-1
    do n3c = 0, ibox_c-1
      fc3(:) = pinvT(:,1)*n1c + pinvT(:,2)*n2c + pinvT(:,3)*n3c
      fc3(:) = fc3(:) - dble(floor(fc3(:)))
      found_dup = .false.
      do j = 1, nfound_c
        if( sum(abs(fc3(:)-fcoset_h(:,j))) < tol_coset ) then
          found_dup = .true.
          exit
        end if
      end do
      if( .not. found_dup ) then
        nfound_c = nfound_c + 1
        hvec_pr(1,nfound_c) = n1c
        hvec_pr(2,nfound_c) = n2c
        hvec_pr(3,nfound_c) = n3c
        fcoset_h(:,nfound_c) = fc3(:)
        if( nfound_c >= nhprk ) exit coset_h_search
      end if
    end do
    end do
    end do coset_h_search
    if( nfound_c >= nhprk ) exit
    ibox_c = ibox_c * 2
    if( ibox_c > 64*nhprk+64 ) then
      if (comm_is_root(nproc_id_global)) then
        write(*,"(A)") 'Error: could not enumerate a full set of reciprocal cosets in dm_unfold primitive-to-reference labeling'
      end if
      call end_parallel
      stop
    end if
  end do

  allocate( unfold%vec_hprk(3,nhprk) )
  allocate( phi_pred(nhprk,nhprk) )
  do icand = 1, nhprk
    unfold%vec_hprk(1,icand) = B_ref(1,1)*hvec_pr(1,icand) + B_ref(1,2)*hvec_pr(2,icand) + B_ref(1,3)*hvec_pr(3,icand)
    unfold%vec_hprk(2,icand) = B_ref(2,1)*hvec_pr(1,icand) + B_ref(2,2)*hvec_pr(2,icand) + B_ref(2,3)*hvec_pr(3,icand)
    unfold%vec_hprk(3,icand) = B_ref(3,1)*hvec_pr(1,icand) + B_ref(3,2)*hvec_pr(2,icand) + B_ref(3,3)*hvec_pr(3,icand)
    do jshift = 1, nhprk
      phi_pred(icand,jshift) = exp( zI * (2d0*pi) * dot_product(fcoset_h(:,icand), dble(nvec_pr(:,jshift))) )
    end do
  end do

  do jshift = 1, nhprk
    tc_cart(1,jshift) = a_pr(1,1)*nvec_pr(1,jshift) + a_pr(1,2)*nvec_pr(2,jshift) + a_pr(1,3)*nvec_pr(3,jshift)
    tc_cart(2,jshift) = a_pr(2,1)*nvec_pr(1,jshift) + a_pr(2,2)*nvec_pr(2,jshift) + a_pr(2,3)*nvec_pr(3,jshift)
    tc_cart(3,jshift) = a_pr(3,1)*nvec_pr(1,jshift) + a_pr(3,2)*nvec_pr(2,jshift) + a_pr(3,3)*nvec_pr(3,jshift)
  end do

! -- precompute the G-grid phase table exp(i G_g.t_c(j)), shared by all bands --
  allocate( phase_gj(ie_ref(1),ie_ref(2),ie_ref(3),nhprk) )
  !$omp parallel do private(ig1_ref,ig2_ref,ig3_ref,ig1,ig2,ig3,gx,gy,gz,jshift) collapse(3)
  do ig1_ref = 1, ie_ref(1)
  do ig2_ref = 1, ie_ref(2)
  do ig3_ref = 1, ie_ref(3)
    ig1 = ig1_ref - 1
    ig2 = ig2_ref - 1
    ig3 = ig3_ref - 1
    if( ig1 > ie_ref(1)/2 ) ig1 = ig1 - ie_ref(1)
    if( ig2 > ie_ref(2)/2 ) ig2 = ig2 - ie_ref(2)
    if( ig3 > ie_ref(3)/2 ) ig3 = ig3 - ie_ref(3)
    gx = ig1*B_ref(1,1) + ig2*B_ref(1,2) + ig3*B_ref(1,3)
    gy = ig1*B_ref(2,1) + ig2*B_ref(2,2) + ig3*B_ref(2,3)
    gz = ig1*B_ref(3,1) + ig2*B_ref(3,2) + ig3*B_ref(3,3)
    do jshift = 1, nhprk
      phase_gj(ig1_ref,ig2_ref,ig3_ref,jshift) = &
        & exp( zI * (gx*tc_cart(1,jshift) + gy*tc_cart(2,jshift) + gz*tc_cart(3,jshift)) )
    end do
  end do
  end do
  end do

! -- per-band labeling --
  allocate( hprk_label_l(no_ref,unfold%nsk), hprk_score_l(no_ref,unfold%nsk) )
  hprk_label_l = 0
  hprk_score_l = 0d0
  nbad_l = 0

  ! isk (the reference-cell k-point index that indexes psi_refG) ranges over
  ! exactly the contiguous block isk_s:isk_e for this rank -- see isk_s/isk_e
  ! above, and unfold%isk_tbl's construction in lattice.f90 (isk is a plain
  ! running count over the ilk-outer/ihk-inner loop there). Phase A's own
  ! computation below only ever uses isk and io_ref; it does not depend on
  ! the reference-to-super-reference decomposition (ilk,ihk) at all, so we
  ! loop directly over isk rather than reconstructing it via isk_tbl(ilk,ihk).
  !$omp parallel do private(isk,io_ref,csum_g2,ig1_ref,ig2_ref,ig3_ref,jshift,icand,cscore, &
  !$omp   score_abs,ibest_c,phi_meas) reduction(+:nbad_l)
  do isk = isk_s, isk_e
  do io_ref = 1, no_ref
    csum_g2 = 0d0
    do ig1_ref = 1, ie_ref(1)
    do ig2_ref = 1, ie_ref(2)
    do ig3_ref = 1, ie_ref(3)
      csum_g2 = csum_g2 + abs( unfold%psi_refG(ig1_ref,ig2_ref,ig3_ref,1,io_ref,isk,1) )**2
    end do
    end do
    end do

    do jshift = 1, nhprk
      phi_meas(jshift) = 0d0
      do ig1_ref = 1, ie_ref(1)
      do ig2_ref = 1, ie_ref(2)
      do ig3_ref = 1, ie_ref(3)
        phi_meas(jshift) = phi_meas(jshift) + abs( unfold%psi_refG(ig1_ref,ig2_ref,ig3_ref,1,io_ref,isk,1) )**2 &
          & * phase_gj(ig1_ref,ig2_ref,ig3_ref,jshift)
      end do
      end do
      end do
      if( csum_g2 > 0d0 ) phi_meas(jshift) = phi_meas(jshift) / csum_g2
    end do

    ibest_c = 0
    score_abs = -1d0
    do icand = 1, nhprk
      cscore = sum( phi_meas(1:nhprk) * conjg(phi_pred(icand,1:nhprk)) ) / dble(nhprk)
      if( abs(cscore) > score_abs ) then
        score_abs = abs(cscore)
        ibest_c = icand
      end if
    end do

    if( 1d0 - score_abs > hprk_thresh ) then
      nbad_l = nbad_l + 1
      hprk_label_l(io_ref,isk) = 0   ! sentinel: no confident label; job continues
    else
      hprk_label_l(io_ref,isk) = ibest_c
    end if
    hprk_score_l(io_ref,isk) = score_abs
  end do
  end do

  allocate( unfold%hprk_label(no_ref,unfold%nsk), unfold%hprk_score(no_ref,unfold%nsk) )
  call comm_summation(hprk_label_l, unfold%hprk_label, no_ref*unfold%nsk, info%icomm_k)
  call comm_summation(hprk_score_l, unfold%hprk_score, no_ref*unfold%nsk, info%icomm_k)
  call comm_summation(nbad_l, nbad, info%icomm_k)

  deallocate( phase_gj, phi_pred, hprk_label_l, hprk_score_l )

  if(comm_is_root(nproc_id_global)) then
    write(*,"(A,I0)") 'primitive-to-reference labeling (Phase A): nhprk = ', nhprk
    if( nbad > 0 ) then
      write(*,"(A,I0,A,ES10.3,A)") 'Warning: ', nbad, ' reference-cell band(s) had 1-|score| above ', &
        & hprk_thresh, '; hat_k label set to the sentinel value 0 for these bands (job continues). Detail:'
      do isk = 1, unfold%nsk
      do io_ref = 1, no_ref
        if( unfold%hprk_label(io_ref,isk) == 0 ) then
          write(*,"(A,I0,A,I0,A,F10.6)") '  isk=', isk, '  io_ref=', io_ref, '  |score|=', unfold%hprk_score(io_ref,isk)
        end if
      end do
      end do
    else
      write(*,"(A)") 'primitive-to-reference labeling (Phase A): all reference-cell bands matched a single hat_k candidate cleanly.'
    end if
  end if

  if( yn_out_mom_distr_gs == 'y' ) then

!   grid for momentum distribution, -nq_mom < iq < nq_mom with dq spacing
    if (dq_mom < 1d-9) dq_mom = (((2*pi)**3/system%det_a)/(num_kgrid(1)*num_kgrid(2)*num_kgrid(3)))**(1d0/3d0)
    if (nq_mom <= 0) nq_mom = 2*int((real(num_kgrid(1),kind=8)*num_kgrid(2)*num_kgrid(3))**(1.0d0/3.0d0))
    allocate( nq_l(-nq_mom:nq_mom,-nq_mom:nq_mom,-nq_mom:nq_mom) )
    nq_l = 0.0d0

!$omp parallel default(shared) &
!$omp private(ilk,ihk,isk,ig1_ref,ig2_ref,ig3_ref,ig1,ig2,ig3,gx,gy,gz,qx,qy,qz,iqx,iqy,iqz, &
!$omp dqx,dqy,dqz,value,nq_l_private)
    allocate( nq_l_private(-nq_mom:nq_mom,-nq_mom:nq_mom,-nq_mom:nq_mom) )
    nq_l_private = 0.0d0

!$omp do collapse(2) 
    do ilk = info%ik_s, info%ik_e
    do ihk = 1, unfold%nhrsk
      isk = unfold%isk_tbl(ilk,ihk)
    do ig1_ref = 1, ie_ref(1)
    do ig2_ref = 1, ie_ref(2)
    do ig3_ref = 1, ie_ref(3)
      ig1 = ig1_ref - 1
      ig2 = ig2_ref - 1
      ig3 = ig3_ref - 1
      if( ig1 > ie_ref(1)/2 ) ig1 = ig1 - ie_ref(1)
      if( ig2 > ie_ref(2)/2 ) ig2 = ig2 - ie_ref(2)
      if( ig3 > ie_ref(3)/2 ) ig3 = ig3 - ie_ref(3)
      gx = ig1*B_ref(1,1) + ig2*B_ref(1,2) + ig3*B_ref(1,3)
      gy = ig1*B_ref(2,1) + ig2*B_ref(2,2) + ig3*B_ref(2,3)
      gz = ig1*B_ref(3,1) + ig2*B_ref(3,2) + ig3*B_ref(3,3)
      qx = system%vec_k(1,ilk) + unfold%vec_hrsk(1,ihk) + gx
      qy = system%vec_k(2,ilk) + unfold%vec_hrsk(2,ihk) + gy
      qz = system%vec_k(3,ilk) + unfold%vec_hrsk(3,ihk) + gz
      iqx = floor( qx/dq_mom )
      iqy = floor( qy/dq_mom )
      iqz = floor( qz/dq_mom )
      if( iqx < -nq_mom .or. iqx >= nq_mom ) cycle
      if( iqy < -nq_mom .or. iqy >= nq_mom ) cycle
      if( iqz < -nq_mom .or. iqz >= nq_mom ) cycle
      dqx = qx/dq_mom - iqx
      dqy = qy/dq_mom - iqy
      dqz = qz/dq_mom - iqz
      value = unfold%wtk_ref(isk) * reta_uu(ig1_ref,ig2_ref,ig3_ref,isk) / dq_mom**3
      nq_l_private(iqx+1,iqy+1,iqz+1) = nq_l_private(iqx+1,iqy+1,iqz+1) + dqx    *dqy    *    dqz * value
      nq_l_private(iqx+1,iqy+1,iqz  ) = nq_l_private(iqx+1,iqy+1,iqz  ) + dqx    *dqy    *(1-dqz) * value
      nq_l_private(iqx+1,iqy,  iqz+1) = nq_l_private(iqx+1,iqy,  iqz+1) + dqx    *(1-dqy)*    dqz * value
      nq_l_private(iqx  ,iqy+1,iqz+1) = nq_l_private(iqx  ,iqy+1,iqz+1) + (1-dqx)*dqy    *    dqz * value
      nq_l_private(iqx+1,iqy  ,iqz  ) = nq_l_private(iqx+1,iqy  ,iqz  ) + dqx    *(1-dqy)*(1-dqz) * value
      nq_l_private(iqx  ,iqy+1,iqz  ) = nq_l_private(iqx  ,iqy+1,iqz  ) + (1-dqx)*dqy    *(1-dqz) * value
      nq_l_private(iqx  ,iqy  ,iqz+1) = nq_l_private(iqx  ,iqy  ,iqz+1) + (1-dqx)*(1-dqy)*    dqz * value
      nq_l_private(iqx  ,iqy  ,iqz  ) = nq_l_private(iqx  ,iqy  ,iqz  ) + (1-dqx)*(1-dqy)*(1-dqz) * value
    end do
    end do
    end do
!
    end do
    end do
!$omp end do

!$omp critical
    nq_l = nq_l + nq_l_private
!$omp end critical
    deallocate(nq_l_private)
!$omp end parallel

    allocate( unfold%nq_gs(-nq_mom:nq_mom,-nq_mom:nq_mom,-nq_mom:nq_mom) )
    call comm_summation(nq_l,unfold%nq_gs,(2*nq_mom+1)**3,icomm)

    rj(:) = 0d0
    nq_sum = sum(unfold%nq_gs(:,:,:)) * dq_mom**3
    do iqx = -nq_mom,nq_mom
    do iqy = -nq_mom,nq_mom
    do iqz = -nq_mom,nq_mom
      qx = iqx * dq_mom
      qy = iqy * dq_mom
      qz = iqz * dq_mom
      rj(1) = rj(1) + qx * unfold%nq_gs(iqx,iqy,iqz)
      rj(2) = rj(2) + qy * unfold%nq_gs(iqx,iqy,iqz)
      rj(3) = rj(3) + qz * unfold%nq_gs(iqx,iqy,iqz)
    end do
    end do
    end do
    rj = rj * dq_mom**3/omega_ref

    if(comm_is_root(nproc_id_global))then
      write(*,'(A,7x,f17.12)')  'N:int n(q,t)            ',nq_sum
      write(*,'(A,7x,3f17.12)') 'J:int q n(q,t)          ',rj(1:3)
    end if

    if(comm_is_root(nproc_id_global))then
      iofile = trim(base_directory)//trim(sysname)//"_momgs.cube"
      fp = 271
      open(fp,file=iofile)
      write(fp,*) "momentum distribution: reference cell, ground state"
      write(fp,*) "All values here are in a.u."

      write(fp,'(i5,3f12.6)') natom, -nq_mom*dq_mom, -nq_mom*dq_mom, -nq_mom*dq_mom
      write(fp,'(i5,3f12.6)') 2*nq_mom+1, dq_mom, 0.0d0, 0.0d0
      write(fp,'(i5,3f12.6)') 2*nq_mom+1, 0.0d0, dq_mom, 0.0d0
      write(fp,'(i5,3f12.6)') 2*nq_mom+1, 0.0d0, 0.0d0, dq_mom
      do iatom=1,natom
        ik=Kion(iatom)
        write(fp,'(i5,4f12.6)') izatom(ik),dble(izatom(ik)),(system%Rion(j,iatom),j=1,3)
      end do
      do iqx = -nq_mom,nq_mom
      do iqy = -nq_mom,nq_mom
        write(fp,'(6(1X,E23.15E3))', advance="yes") (unfold%nq_gs(iqx,iqy,iqz),iqz = -nq_mom,nq_mom)
      end do
      end do
      close(fp)
    end if

  end if ! yn_out_mom_distr_gs
#else
  stop 'dm_unfold_option requires a build with MPI support.'
#endif

  end subroutine init_dm_unfold

!===================================================================================================================================
  subroutine dm_unfold(itt,system,info,lg,ofl,psi_t,unfold)

    use structures
    use communication, only: comm_is_root, comm_summation
    use parallelization, only: nproc_id_global
    use salmon_global, only: no_ref, dt, num_kgrid, num_skgrid, sysname, base_directory, natom,izatom,kion, &
                            & dq_mom, nq_mom, yn_out_mom_distr_rt, out_mom_distr_rt_step
    use inputoutput, only: t_unit_time, t_unit_ac, t_unit_current
    use math_constants, only: pi
    implicit none
    integer                 ,intent(in) :: itt
    type(s_ofile)           ,intent(in) :: ofl
    type(s_rgrid)           ,intent(in) :: lg
    type(s_dft_system)      ,intent(in) :: system
    type(s_parallel_info)   ,intent(in) :: info
    type(s_orbital)         ,intent(in) :: psi_t
    type(s_unfold)          ,intent(in) :: unfold

    character(10) :: filenum
    character(60) :: iofile
    integer fp,j,iatom,ik
    integer ie_ref(3),ilk,ihk,isk,io,io_ref,ih1,ih2,ih3,ir1_ref,ir2_ref,ir3_ref
    integer ir1,ir2,ir3,ig1,ig2,ig3,ig1_ref,ig2_ref,ig3_ref,ispin,ig_ref(3),iqx,iqy,iqz
    integer io_ref1,io_ref2, icomm, isk_s, isk_e, nsk_se
    real(8) :: B_ref(3,3),gx,gy,gz,omega_ref,dqx,dqy,dqz,qx,qy,qz,value
    complex(8),allocatable :: mat(:,:,:,:), eta(:,:,:), eta_l(:,:,:)
    complex(8),allocatable :: eta_uu_d(:,:,:,:),eta_uu_nd(:,:,:,:)
    complex(8),allocatable :: nq_d_l(:,:,:),nq_d(:,:,:),nq_nd_l(:,:,:),nq_nd(:,:,:),nq_d_l_private(:,:,:),nq_nd_l_private(:,:,:)
    complex(8) :: zj1(3),zj2(3),zj3(3),zj4(3),zj5(3),zj6(3),zj1_uu(3),zj2_uu(3),zj3_uu(3),zj4_uu(3),zj_d(3),zj_nd(3)
    complex(8) :: zj1_l(3),zj2_l(3),zj3_l(3),zj4_l(3),zj5_l(3),zj6_l(3),zsum,zsum_d,zsum_nd,zsum_l,zsum_uu
      
    allocate( mat(no_ref,unfold%nhrsk,info%io_s:info%io_e,info%ik_s:info%ik_e))
    ie_ref(1:3) = lg%ie(1:3)/unfold%num_hkgrid(1:3)
    omega_ref = system%hvol * system%ngrid / dble(unfold%num_hkgrid(1)*unfold%num_hkgrid(2)*unfold%num_hkgrid(3))

    ispin = 1
    isk_s = (info%ik_s-1) * unfold%nhrsk + 1
    isk_e = info%ik_e * unfold%nhrsk
    nsk_se = isk_e - isk_s + 1

  !$omp parallel do private(ilk,ihk,isk,io,io_ref,zsum,ih1,ih2,ih3,ir1_ref,ir2_ref,ir3_ref,ir1,ir2,ir3) collapse(2)
    do ilk = info%ik_s, info%ik_e ! large k
    do ihk = 1, unfold%nhrsk   ! hat k
      isk = unfold%isk_tbl(ilk,ihk) !small k = large k + hat k
    do io = info%io_s, info%io_e     ! m, supercell
    do io_ref = 1, no_ref   ! n, reference
      zsum = 0d0
      do ih1 = 1, unfold%num_hkgrid(1)
      do ih2 = 1, unfold%num_hkgrid(2)
      do ih3 = 1, unfold%num_hkgrid(3)
      do ir1_ref = 1, ie_ref(1)
      do ir2_ref = 1, ie_ref(2)
      do ir3_ref = 1, ie_ref(3)
        ir1 = ir1_ref + (ih1-1) * ie_ref(1)
        ir2 = ir2_ref + (ih2-1) * ie_ref(2)
        ir3 = ir3_ref + (ih3-1) * ie_ref(3)
        zsum = zsum + conjg( unfold%psi_ref( ir1_ref, ir2_ref, ir3_ref, 1, io_ref, isk, 1) ) &
         &        * conjg( unfold%eihkr_tbl(ir1,ir2,ir3,ihk) ) &
         &        * psi_t%zwf( ir1, ir2, ir3, 1, io, ilk, 1 )
      end do
      end do
      end do
      end do
      end do
      end do
      zsum = zsum * system%hvol
      mat(io_ref, ihk, io, ilk) = zsum / (unfold%num_hkgrid(1)*unfold%num_hkgrid(2)*unfold%num_hkgrid(3))
    end do
    end do
    end do
    end do

    allocate( eta(no_ref,no_ref,isk_s:isk_e),eta_l(no_ref,no_ref,isk_s:isk_e) )
    eta_l = 0.0d0
  !$omp parallel do private(ilk,ihk,isk,io_ref1,io_ref2,zsum,io) collapse(2)
    do ilk = info%ik_s, info%ik_e
    do ihk = 1, unfold%nhrsk
      isk = unfold%isk_tbl(ilk,ihk)
      do io_ref1 = 1, no_ref
      do io_ref2 = 1, no_ref
        zsum = 0d0
        do io = info%io_s, info%io_e
          zsum = zsum + system%rocc(io,ilk,1) * mat(io_ref1, ihk, io, ilk) * conjg( mat(io_ref2, ihk, io, ilk) )
        end do
        eta_l(io_ref1, io_ref2, isk) = zsum
      end do
      end do
    end do
    end do
    eta_l = eta_l * (unfold%num_hkgrid(1)*unfold%num_hkgrid(2)*unfold%num_hkgrid(3))
    eta = 0.0d0
    call comm_summation(eta_l,eta,no_ref*no_ref*nsk_se,info%icomm_o)

    zsum_l = 0d0
    zj1_l(1:3) = 0d0
    zj2_l(1:3) = 0d0
    zj3_l(1:3) = 0d0
    zj4_l(1:3) = 0d0
    zj5_l(1:3) = 0d0
    zj6_l(1:3) = 0d0
  !$omp parallel do private(ilk,ihk,isk,io_ref1,io_ref2) reduction(+:zsum_l,zj1_l,zj2_l,zj3_l,zj4_l,zj5_l,zj6_l) collapse(2)
    do ilk = info%ik_s, info%ik_e
    do ihk = 1, unfold%nhrsk
       isk = unfold%isk_tbl(ilk,ihk)
    do io_ref1 = 1, no_ref
       zsum_l = zsum_l + eta(io_ref1, io_ref1, isk) * unfold%wtk_ref(isk)
       zj1_l(:) = zj1_l(:) + eta(io_ref1, io_ref1, isk) * unfold%wtk_ref(isk) &
        & * unfold%upu_ref(:, io_ref1, io_ref1, isk)
       zj2_l(:) = zj2_l(:) + eta(io_ref1, io_ref1, isk) * unfold%wtk_ref(isk) &
        & * unfold%u_rVnl_Vnlr_u_ref(:, io_ref1, io_ref1, isk)
       zj5_l(:) = zj5_l(:) + eta(io_ref1, io_ref1, isk) * unfold%wtk_ref(isk) &
        & * (system%vec_k(:,ilk) + unfold%vec_hrsk(:,ihk))
       zj6_l(:) = zj6_l(:) + eta(io_ref1, io_ref1, isk) * unfold%wtk_ref(isk) &
        & * system%vec_Ac(:)
    do io_ref2 = 1, no_ref
       if( io_ref1 /= io_ref2) then
         zj3_l(:) = zj3_l(:) + eta(io_ref2, io_ref1, isk) * unfold%wtk_ref(isk) &
        & * unfold%upu_ref(:, io_ref1, io_ref2, isk)
         zj4_l(:) = zj4_l(:) + eta(io_ref2, io_ref1, isk) * unfold%wtk_ref(isk) &
        & * unfold%u_rVnl_Vnlr_u_ref(:, io_ref1, io_ref2, isk)
       end if
    end do
    end do
    end do
    end do
    icomm = info%icomm_k

    zsum = 0d0
    zj1(:) = 0d0
    zj2(:) = 0d0
    zj3(:) = 0d0
    zj4(:) = 0d0
    zj5(:) = 0d0
    zj6(:) = 0d0
    call comm_summation(zsum_l,zsum,icomm)
    call comm_summation(zj1_l,zj1,3,icomm)
    call comm_summation(zj2_l,zj2,3,icomm)
    call comm_summation(zj3_l,zj3,3,icomm)
    call comm_summation(zj4_l,zj4,3,icomm)
    call comm_summation(zj5_l,zj5,3,icomm)
    call comm_summation(zj6_l,zj6,3,icomm)

    zj1 = zj1 / omega_ref
    zj2 = zj2 / omega_ref
    zj3 = zj3 / omega_ref
    zj4 = zj4 / omega_ref
    zj5 = zj5 / omega_ref
    zj6 = zj6 / omega_ref

  if(comm_is_root(nproc_id_global))then
    write(*,'(A,2x,i7,2x,A,3f17.12)') 'dm_unfold  it=', itt, '     Ac(t)=',system%vec_Ac(:)
    write(*,'(A,7x,f17.12)')          'N:Tr[rho(t)]             ', real(zsum)
    write(*,'(A,7x,3f17.12)')         'J:rho(t)<unk|i[h,r]|unk> ', real(zj1(:)+zj2(:)+zj3(:)+zj4(:)+zj5(:)+zj6(:))
  end if

  allocate( eta_uu_d(1:ie_ref(1),1:ie_ref(2),1:ie_ref(3),isk_s:isk_e) )
  allocate( eta_uu_nd(1:ie_ref(1),1:ie_ref(2),1:ie_ref(3),isk_s:isk_e) )
!$omp parallel do private(ilk,ihk,isk,ig1_ref,ig2_ref,ig3_ref,io_ref1,io_ref2,zsum_d,zsum_nd) collapse(2)
  do ilk = info%ik_s, info%ik_e
  do ihk = 1, unfold%nhrsk
    isk = unfold%isk_tbl(ilk,ihk)
  do ig1_ref = 1, ie_ref(1)
  do ig2_ref = 1, ie_ref(2)
  do ig3_ref = 1, ie_ref(3)
    zsum_d = 0d0
    zsum_nd = 0d0
    do io_ref1 = 1, no_ref
    do io_ref2 = 1, no_ref
      if( io_ref1 == io_ref2 ) then
        zsum_d = zsum_d + eta(io_ref1, io_ref2, isk) * unfold%psi_refG(ig1_ref,ig2_ref,ig3_ref,ispin,io_ref1,isk,1) &
        &                               * conjg( unfold%psi_refG(ig1_ref,ig2_ref,ig3_ref,ispin,io_ref2,isk,1) )
      else
        zsum_nd = zsum_nd + eta(io_ref1, io_ref2, isk) * unfold%psi_refG(ig1_ref,ig2_ref,ig3_ref,ispin,io_ref1,isk,1) &
        &                               * conjg( unfold%psi_refG(ig1_ref,ig2_ref,ig3_ref,ispin,io_ref2,isk,1) )
      end if
    end do
    end do
    eta_uu_d(ig1_ref,ig2_ref,ig3_ref,isk) = zsum_d * omega_ref
    eta_uu_nd(ig1_ref,ig2_ref,ig3_ref,isk) = zsum_nd * omega_ref
  enddo
  enddo
  enddo

  enddo
  enddo

  B_ref(:,:) = system%primitive_b(:,:)
  B_ref(:,1) = B_ref(:,1) * unfold%num_hkgrid(1)
  B_ref(:,2) = B_ref(:,2) * unfold%num_hkgrid(2)
  B_ref(:,3) = B_ref(:,3) * unfold%num_hkgrid(3)
  zsum_l = 0d0
  zj1_l(1:3) = 0d0
  zj2_l(1:3) = 0d0
  zj3_l(1:3) = 0d0
  zj4_l(1:3) = 0d0
!$omp parallel do private(ilk,ihk,isk,ig1_ref,ig2_ref,ig3_ref,ig1,ig2,ig3,gx,gy,gz) &
!$omp reduction(+:zsum_l,zj1_l,zj2_l,zj3_l,zj4_l) collapse(2)
  do ilk = info%ik_s, info%ik_e
  do ihk = 1, unfold%nhrsk
    isk = unfold%isk_tbl(ilk,ihk)
  do ig1_ref = 1, ie_ref(1)
  do ig2_ref = 1, ie_ref(2)
  do ig3_ref = 1, ie_ref(3)
    ig1 = ig1_ref - 1
    ig2 = ig2_ref - 1
    ig3 = ig3_ref - 1
    if( ig1 > ie_ref(1)/2 ) ig1 = ig1 - ie_ref(1)
    if( ig2 > ie_ref(2)/2 ) ig2 = ig2 - ie_ref(2)
    if( ig3 > ie_ref(3)/2 ) ig3 = ig3 - ie_ref(3)
    gx = ig1*B_ref(1,1) + ig2*B_ref(1,2) + ig3*B_ref(1,3)
    gy = ig1*B_ref(2,1) + ig2*B_ref(2,2) + ig3*B_ref(2,3)
    gz = ig1*B_ref(3,1) + ig2*B_ref(3,2) + ig3*B_ref(3,3)
    zsum_l = zsum_l + unfold%wtk_ref(isk) * eta_uu_d(ig1_ref,ig2_ref,ig3_ref,isk)
    zj1_l(1) = zj1_l(1) + unfold%wtk_ref(isk) * eta_uu_d(ig1_ref,ig2_ref,ig3_ref,isk) * gx
    zj1_l(2) = zj1_l(2) + unfold%wtk_ref(isk) * eta_uu_d(ig1_ref,ig2_ref,ig3_ref,isk) * gy
    zj1_l(3) = zj1_l(3) + unfold%wtk_ref(isk) * eta_uu_d(ig1_ref,ig2_ref,ig3_ref,isk) * gz
    zj2_l(1) = zj2_l(1) + unfold%wtk_ref(isk) * eta_uu_nd(ig1_ref,ig2_ref,ig3_ref,isk) * gx
    zj2_l(2) = zj2_l(2) + unfold%wtk_ref(isk) * eta_uu_nd(ig1_ref,ig2_ref,ig3_ref,isk) * gy
    zj2_l(3) = zj2_l(3) + unfold%wtk_ref(isk) * eta_uu_nd(ig1_ref,ig2_ref,ig3_ref,isk) * gz
    zj3_l(:) = zj3_l(:) + unfold%wtk_ref(isk) * eta_uu_d(ig1_ref,ig2_ref,ig3_ref,isk) &
        & * (system%vec_k(:,ilk) + unfold%vec_hrsk(:,ihk))
    zj4_l(:) = zj4_l(:) + unfold%wtk_ref(isk) * eta_uu_d(ig1_ref,ig2_ref,ig3_ref,isk) * system%vec_Ac(:)
  end do
  end do
  end do

  end do
  end do

  zsum_uu = 0d0
  zj1_uu(:) = 0d0
  zj2_uu(:) = 0d0
  zj3_uu(:) = 0d0
  zj4_uu(:) = 0d0
  call comm_summation(zsum_l,zsum_uu,icomm)
  call comm_summation(zj1_l,zj1_uu,3,icomm)
  call comm_summation(zj2_l,zj2_uu,3,icomm)
  call comm_summation(zj3_l,zj3_uu,3,icomm)
  call comm_summation(zj4_l,zj4_uu,3,icomm)

  zj1_uu = zj1_uu / omega_ref
  zj2_uu = zj2_uu / omega_ref
  zj3_uu = zj3_uu / omega_ref
  zj4_uu = zj4_uu / omega_ref

  if(comm_is_root(nproc_id_global))then
    write(*,'(A,7x,3f17.12)')         'J:rho_uu(k,G)(G+k+A)     ', real(zj1_uu(:)+zj2_uu(:)+zj3_uu(:)+zj4_uu(:))
  end if

  if( yn_out_mom_distr_rt == 'y' .and. (itt==1 .or. mod(itt,out_mom_distr_rt_step)==0)) then

!   grid for momentum distribution, -nq_mom < iq < nq_mom with dq spacing
  if (dq_mom < 1d-9) dq_mom = (((2*pi)**3/system%det_a)/(num_kgrid(1)*num_kgrid(2)*num_kgrid(3)))**(1d0/3d0)
  if (nq_mom <= 0) nq_mom = 2*int((real(num_kgrid(1),kind=8)*num_kgrid(2)*num_kgrid(3))**(1.0d0/3.0d0))
  allocate( nq_d_l( -nq_mom:nq_mom,-nq_mom:nq_mom,-nq_mom:nq_mom) )
  allocate( nq_nd_l(-nq_mom:nq_mom,-nq_mom:nq_mom,-nq_mom:nq_mom) )
  nq_d_l = 0.0d0
  nq_nd_l = 0.0d0

!$omp parallel default(shared) &
!$omp private(ilk,ihk,isk,ig1_ref,ig2_ref,ig3_ref,ig1,ig2,ig3,gx,gy,gz,qx,qy,qz,iqx,iqy,iqz, &
!$omp dqx,dqy,dqz,value,nq_d_l_private,nq_nd_l_private)
    allocate( nq_d_l_private(-nq_mom:nq_mom,-nq_mom:nq_mom,-nq_mom:nq_mom) )
    allocate( nq_nd_l_private(-nq_mom:nq_mom,-nq_mom:nq_mom,-nq_mom:nq_mom) )
    nq_d_l_private = 0.0d0
    nq_nd_l_private = 0.0d0

!$omp do collapse(2) 
    do ilk = info%ik_s, info%ik_e
    do ihk = 1, unfold%nhrsk
      isk = unfold%isk_tbl(ilk,ihk)
    do ig1_ref = 1, ie_ref(1)
    do ig2_ref = 1, ie_ref(2)
    do ig3_ref = 1, ie_ref(3)
      ig1 = ig1_ref - 1
      ig2 = ig2_ref - 1
      ig3 = ig3_ref - 1
      if( ig1 > ie_ref(1)/2 ) ig1 = ig1 - ie_ref(1)
      if( ig2 > ie_ref(2)/2 ) ig2 = ig2 - ie_ref(2)
      if( ig3 > ie_ref(3)/2 ) ig3 = ig3 - ie_ref(3)
      gx = ig1*B_ref(1,1) + ig2*B_ref(1,2) + ig3*B_ref(1,3)
      gy = ig1*B_ref(2,1) + ig2*B_ref(2,2) + ig3*B_ref(2,3)
      gz = ig1*B_ref(3,1) + ig2*B_ref(3,2) + ig3*B_ref(3,3)
      qx = system%vec_Ac(1) + system%vec_k(1,ilk) + unfold%vec_hrsk(1,ihk) + gx
      qy = system%vec_Ac(2) + system%vec_k(2,ilk) + unfold%vec_hrsk(2,ihk) + gy
      qz = system%vec_Ac(3) + system%vec_k(3,ilk) + unfold%vec_hrsk(3,ihk) + gz
      iqx = floor( qx/dq_mom )
      iqy = floor( qy/dq_mom )
      iqz = floor( qz/dq_mom )
      if( iqx < -nq_mom .or. iqx >= nq_mom ) cycle
      if( iqy < -nq_mom .or. iqy >= nq_mom ) cycle
      if( iqz < -nq_mom .or. iqz >= nq_mom ) cycle
      dqx = qx/dq_mom - iqx
      dqy = qy/dq_mom - iqy
      dqz = qz/dq_mom - iqz

      value = unfold%wtk_ref(isk) * eta_uu_d(ig1_ref,ig2_ref,ig3_ref,isk) / dq_mom**3
      nq_d_l_private(iqx+1,iqy+1,iqz+1) = nq_d_l_private(iqx+1,iqy+1,iqz+1) + dqx    *dqy    *    dqz * value
      nq_d_l_private(iqx+1,iqy+1,iqz  ) = nq_d_l_private(iqx+1,iqy+1,iqz  ) + dqx    *dqy    *(1-dqz) * value
      nq_d_l_private(iqx+1,iqy,  iqz+1) = nq_d_l_private(iqx+1,iqy,  iqz+1) + dqx    *(1-dqy)*    dqz * value
      nq_d_l_private(iqx  ,iqy+1,iqz+1) = nq_d_l_private(iqx  ,iqy+1,iqz+1) + (1-dqx)*dqy    *    dqz * value
      nq_d_l_private(iqx+1,iqy  ,iqz  ) = nq_d_l_private(iqx+1,iqy  ,iqz  ) + dqx    *(1-dqy)*(1-dqz) * value
      nq_d_l_private(iqx  ,iqy+1,iqz  ) = nq_d_l_private(iqx  ,iqy+1,  iqz) + (1-dqx)*dqy    *(1-dqz) * value
      nq_d_l_private(iqx  ,iqy  ,iqz+1) = nq_d_l_private(iqx  ,iqy  ,iqz+1) + (1-dqx)*(1-dqy)*    dqz * value
      nq_d_l_private(iqx  ,iqy  ,iqz  ) = nq_d_l_private(iqx  ,iqy  ,iqz  ) + (1-dqx)*(1-dqy)*(1-dqz) * value

      value = unfold%wtk_ref(isk) * eta_uu_nd(ig1_ref,ig2_ref,ig3_ref,isk) / dq_mom**3
      nq_nd_l_private(iqx+1,iqy+1,iqz+1) = nq_nd_l_private(iqx+1,iqy+1,iqz+1) + dqx    *dqy    *    dqz * value
      nq_nd_l_private(iqx+1,iqy+1,iqz  ) = nq_nd_l_private(iqx+1,iqy+1,iqz  ) + dqx    *dqy    *(1-dqz) * value
      nq_nd_l_private(iqx+1,iqy,  iqz+1) = nq_nd_l_private(iqx+1,iqy,  iqz+1) + dqx    *(1-dqy)*    dqz * value
      nq_nd_l_private(iqx  ,iqy+1,iqz+1) = nq_nd_l_private(iqx  ,iqy+1,iqz+1) + (1-dqx)*dqy    *    dqz * value
      nq_nd_l_private(iqx+1,iqy  ,iqz  ) = nq_nd_l_private(iqx+1,iqy  ,iqz  ) + dqx    *(1-dqy)*(1-dqz) * value
      nq_nd_l_private(iqx  ,iqy+1,iqz  ) = nq_nd_l_private(iqx  ,iqy+1,  iqz) + (1-dqx)*dqy    *(1-dqz) * value
      nq_nd_l_private(iqx  ,iqy  ,iqz+1) = nq_nd_l_private(iqx  ,iqy  ,iqz+1) + (1-dqx)*(1-dqy)*    dqz * value
      nq_nd_l_private(iqx  ,iqy  ,iqz  ) = nq_nd_l_private(iqx  ,iqy  ,iqz  ) + (1-dqx)*(1-dqy)*(1-dqz) * value
    end do
    end do
    end do
!
    end do
    end do
!$omp end do

!$omp critical
    nq_d_l = nq_d_l + nq_d_l_private
    nq_nd_l = nq_nd_l + nq_nd_l_private
!$omp end critical
    deallocate(nq_d_l_private, nq_nd_l_private)
!$omp end parallel

    allocate( nq_d(-nq_mom:nq_mom,-nq_mom:nq_mom,-nq_mom:nq_mom), nq_nd(-nq_mom:nq_mom,-nq_mom:nq_mom,-nq_mom:nq_mom) )
    nq_d = 0.0d0
    nq_nd = 0.0d0
    call comm_summation(nq_d_l,nq_d,(2*nq_mom+1)**3,icomm)
    call comm_summation(nq_nd_l,nq_nd,(2*nq_mom+1)**3,icomm)

    zsum_d = sum(nq_d(:,:,:))*dq_mom**3
    zsum_nd = sum(nq_nd(:,:,:))*dq_mom**3

    zj_d(:) = 0d0
    zj_nd(:) = 0d0
    do iqx = -nq_mom,nq_mom
    do iqy = -nq_mom,nq_mom
    do iqz = -nq_mom,nq_mom
      qx = iqx * dq_mom
      qy = iqy * dq_mom
      qz = iqz * dq_mom
      zj_d(1) = zj_d(1) + qx * nq_d(iqx,iqy,iqz)
      zj_d(2) = zj_d(2) + qy * nq_d(iqx,iqy,iqz)
      zj_d(3) = zj_d(3) + qz * nq_d(iqx,iqy,iqz)
      zj_nd(1) = zj_nd(1) + qx * nq_nd(iqx,iqy,iqz)
      zj_nd(2) = zj_nd(2) + qy * nq_nd(iqx,iqy,iqz)
      zj_nd(3) = zj_nd(3) + qz * nq_nd(iqx,iqy,iqz)
    end do
    end do
    end do
    zj_d = zj_d * dq_mom**3/omega_ref
    zj_nd = zj_nd * dq_mom**3/omega_ref

  if(comm_is_root(nproc_id_global))then
    write(*,'(A,7x,f17.12)')          'J:int n(q,t)             ', real(zsum_d)
    write(*,'(A,7x,3f17.12)')         'J:int q n(q,t)           ', real(zj_d(:)+zj_nd(:))
  end if

    if(comm_is_root(nproc_id_global))then
      write(filenum, '(i6.6)') itt
      iofile = trim(base_directory)//trim(sysname)//"_momdiff_"//trim(adjustl(filenum))//".cube"
      fp = 271
      open(fp,file=iofile)
      write(fp,*) "momentum distribution: reference cell, difference from ground state"
      write(fp,*) "All values here are in a.u."

      write(fp,'(i5,3f12.6)') natom, -nq_mom*dq_mom, -nq_mom*dq_mom, -nq_mom*dq_mom
      write(fp,'(i5,3f12.6)') 2*nq_mom+1, dq_mom, 0.0d0, 0.0d0
      write(fp,'(i5,3f12.6)') 2*nq_mom+1, 0.0d0, dq_mom, 0.0d0
      write(fp,'(i5,3f12.6)') 2*nq_mom+1, 0.0d0, 0.0d0, dq_mom
      do iatom=1,natom
        ik=Kion(iatom)
        write(fp,'(i5,4f12.6)') izatom(ik),dble(izatom(ik)),(system%Rion(j,iatom),j=1,3)
      end do
      do iqx = -nq_mom,nq_mom
      do iqy = -nq_mom,nq_mom
        write(fp,'(6(1X,E23.15E3))', advance="yes") &
        & (real(nq_d(iqx,iqy,iqz)+nq_nd(iqx,iqy,iqz)) - unfold%nq_gs(iqx,iqy,iqz),iqz = -nq_mom,nq_mom)
      end do
      end do
      close(fp)
    end if

  end if ! out_dm_unfold_mom_step

  if(comm_is_root(nproc_id_global))then

    write(ofl%fh_dm_unfold,'(55f17.12)') itt*dt, system%vec_Ac(1:3)*t_unit_ac%conv, &
    & real(zsum),real(zj1(1:3))*t_unit_current%conv,real(zj2(1:3))*t_unit_current%conv,real(zj3(1:3))*t_unit_current%conv, &
    & real(zj4(1:3))*t_unit_current%conv,real(zj5(1:3))*t_unit_current%conv,real(zj6(1:3))*t_unit_current%conv, &
    & real(zj1(1:3)+zj3(1:3)+zj5(1:3)+zj6(1:3))*t_unit_current%conv, &
    & real(zj1(1:3)+zj2(1:3)+zj3(1:3)+zj4(1:3)+zj5(1:3)+zj6(1:3))*t_unit_current%conv, &
    & real(zsum_uu),real(zj1_uu(1:3))*t_unit_current%conv,real(zj2_uu(1:3))*t_unit_current%conv, &
    & real(zj3_uu(1:3))*t_unit_current%conv,real(zj4_uu(1:3))*t_unit_current%conv, &
    & real(zj1_uu(1:3)+zj2_uu(1:3)+zj3_uu(1:3)+zj4_uu(1:3))*t_unit_current%conv, &
    & real(zsum_d+zsum_nd),real(zj_d(1:3))*t_unit_current%conv,real(zj_nd(1:3))*t_unit_current%conv, &
    & real(zj_d(1:3)+zj_nd(1:3))*t_unit_current%conv

  end if

  return
           
  end subroutine dm_unfold

end module dm_unfold_sub
