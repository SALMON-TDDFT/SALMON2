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
module dm_unfold_sub
  implicit none

contains

  subroutine init_dm_unfold(lg,system,info,ofl,unfold)

  use structures
  use communication, only: comm_is_root, comm_bcast, comm_summation
  use parallelization, only: nproc_id_global, nproc_group_global, end_parallel
  use salmon_global, only: dm_unfold_option, no_pr, base_directory, sysname, natom, izatom, kion, &
                         & yn_out_mom_distr_gs, dq_mom, nq_mom, num_kgrid
  use filesystem, only: open_filehandle
  use inputoutput, only: t_unit_time, t_unit_ac, t_unit_current
  use math_constants, only: zI,pi
  use mpi
  implicit none
  type(s_rgrid),           intent(in)    :: lg
  type(s_dft_system),      intent(inout) :: system
  type(s_parallel_info),   intent(in)    :: info
  type(s_ofile),           intent(out)   :: ofl
  type(s_unfold),          intent(inout) :: unfold
  character(256) :: iofile
  integer :: icomm
  integer :: gsize(7), lsize(7), lstart(7)
  integer :: gsize4(4), lsize4(4), lstart4(4)
  integer :: iopen_flag, minfo, mfile
  integer :: source_type, local_type, global_type, file_type
  integer :: ierr, n_count
  integer(kind=MPI_OFFSET_KIND) :: disp_upu, disp_vnl
  integer :: nspin,ispin,nsk_se,isk_s,isk_e,ie_pr(3),ihk,ir1,ir2,ir3
  integer :: ih1,ih2,ih3,ir1_pr,ir2_pr,ir3_pr,ig1_pr,ig2_pr,ig3_pr,io_pr,isk,ilk
  integer :: ig1,ig2,ig3,iqx,iqy,iqz,fp,iatom,j,ik
  real(8) :: omega_pr,B_pr(3,3),rsum,rsum_l,rj_l(3),rj(3),gx,gy,gz
  real(8) :: qx,qy,qz,dqx,dqy,dqz,value,nq_sum
  complex(8) :: zsum
  real(8),allocatable :: reta_uu(:,:,:,:),nq_l(:,:,:),nq_l_private(:,:,:)
  logical :: exists, e_occupation, e_wfn, e_tm

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

  if (comm_is_root(nproc_id_global)) then
    inquire(file='primitive/wfn.bin', exist = e_wfn)
    inquire(file='primitive/occupation.bin', exist = e_occupation)
    inquire(file='primitive/tm.bin', exist = e_tm)
  end if
  exists = e_wfn .and. e_occupation .and. e_tm
  call comm_bcast(exists, nproc_group_global)
  if( .not. exists ) then
    if (comm_is_root(nproc_id_global)) then
      write(*,"(A)") 'Error: file not found, primitive/wfn.bin, occupation.bin, tm.bin'
    end if
    call end_parallel
    stop
  end if

  nspin = 1
  isk_s = (info%ik_s-1) * unfold%nhk + 1
  isk_e = info%ik_e * unfold%nhk
  nsk_se = isk_e - isk_s + 1

  if( any(mod(lg%ie(:), unfold%num_hkgrid(:)) /= 0) ) then
    if (comm_is_root(nproc_id_global)) then
      write(*,"(A)") 'ie_pr(1:3) error at init_dm_unfold'
    end if
    call end_parallel
    stop
  end if

  ie_pr(1:3)=lg%ie(1:3)/unfold%num_hkgrid(1:3) ! size of primitive cell grid

  allocate( unfold%psi_pr(1:ie_pr(1),1:ie_pr(2),1:ie_pr(3),nspin,no_pr,isk_s:isk_e,1) )

  source_type = MPI_DOUBLE_COMPLEX
  minfo = MPI_INFO_NULL
  iopen_flag = MPI_MODE_RDONLY
  iofile = "primitive/wfn.bin"
  icomm = info%icomm_k

! create MPI_Type (Window) of process-local wave function
  gsize  = [ie_pr(1:3), nspin, no_pr, nsk_se, 1]
  lsize  = [ie_pr(1:3), nspin, no_pr, nsk_se, 1]
  lstart = [1,1,1,      1,     1,     1,      1] - 1

  call MPI_Type_create_subarray(7, gsize, lsize, lstart, MPI_ORDER_FORTRAN, source_type, local_type, ierr)
  call MPI_Type_commit(local_type, ierr)

! create MPI_Type (Window) of global wave function
  gsize  = [ie_pr(1:3), nspin, no_pr, unfold%nsk,     1]
  lstart = [1,1,1,      1,     1,     isk_s,          1] - 1

  call MPI_Type_create_subarray(7, gsize, lsize, lstart, MPI_ORDER_FORTRAN, source_type, global_type, ierr)
  call MPI_Type_commit(global_type, ierr)

  call MPI_File_open(icomm, iofile, iopen_flag, minfo, mfile, ierr)
  call MPI_File_set_view(mfile, 0_MPI_OFFSET_KIND, local_type, global_type, 'native', MPI_INFO_NULL, ierr)

  call MPI_File_read_all(mfile, unfold%psi_pr, 1, local_type, MPI_STATUS_IGNORE, ierr)

  call MPI_File_close(mfile, ierr)
  call MPI_Type_free(global_type, ierr)
  call MPI_Type_free( local_type, ierr)

  if(comm_is_root(nproc_id_global)) then
    write(*,*) 'End reading primitive/wfn.bin'
  end if

! read occupation

  allocate( unfold%rocc_pr(no_pr, unfold%nsk, nspin) )

  if(comm_is_root(nproc_id_global)) then
     iofile = "primitive/occupation.bin"
     open(888,file=iofile,form='unformatted')
     read(888) unfold%rocc_pr(1:no_pr,1:unfold%nsk,1:nspin)
     close(888)
  end if
  call comm_bcast(unfold%rocc_pr,nproc_group_global)

  if(comm_is_root(nproc_id_global))then
    write(*,*) 'End reading primitive/occupation.bin'
  end if

! read transition dipole matrix elements

    allocate( unfold%upu_pr(3, no_pr, no_pr, isk_s:isk_e), unfold%u_rVnl_Vnlr_u_pr(3, no_pr, no_pr, isk_s:isk_e) )

    iofile = "primitive/tm.bin"
    icomm = info%icomm_k

    gsize4  = [3, no_pr , no_pr , unfold%nsk       ]
    lsize4  = [3, no_pr , no_pr , isk_e - isk_s + 1]
    lstart4 = [1, 1     , 1     , isk_s            ] - 1

    n_count = lsize4(1)*lsize4(2)*lsize4(3)*lsize4(4)

    call MPI_Type_create_subarray(4, gsize4, lsize4, lstart4, MPI_ORDER_FORTRAN, source_type, file_type, ierr)
    call MPI_Type_commit(file_type, ierr)
    call MPI_File_open(icomm, iofile, iopen_flag, minfo, mfile, ierr)

    disp_upu = 0_MPI_OFFSET_KIND
    disp_vnl = int(3,8)*no_pr*no_pr*unfold%nsk*int(16,8)

    call MPI_File_set_view(mfile, disp_upu, source_type, file_type, 'native', minfo, ierr)
    call MPI_File_read_all(mfile, unfold%upu_pr(1,1,1,isk_s), n_count, source_type, MPI_STATUS_IGNORE, ierr)

    call MPI_File_set_view(mfile, disp_vnl, source_type, file_type, 'native', minfo, ierr)
    call MPI_File_read_all(mfile, unfold%u_rVnl_Vnlr_u_pr(1,1,1,isk_s), n_count, source_type, MPI_STATUS_IGNORE, ierr)

    call MPI_File_close(mfile, ierr) 
    call MPI_Type_free(file_type, ierr)

  if(comm_is_root(nproc_id_global))then
    write(*,*) 'End reading primitive/tm.bin'
  end if

! exp(i hat_k r) table
  allocate( unfold%eihkr_tbl(lg%ie(1),lg%ie(2),lg%ie(3),unfold%nhk) )
  !$omp parallel do private(ihk,ih1,ih2,ih3,ir1_pr,ir2_pr,ir3_pr,ir1,ir2,ir3) collapse(4)
  do ihk = 1, unfold%nhk
  do ih1 = 1, unfold%num_hkgrid(1)
  do ih2 = 1, unfold%num_hkgrid(2)
  do ih3 = 1, unfold%num_hkgrid(3)
  do ir1_pr = 1, ie_pr(1)
  do ir2_pr = 1, ie_pr(2)
  do ir3_pr = 1, ie_pr(3)
    ir1 = ir1_pr + (ih1-1) * ie_pr(1)
    ir2 = ir2_pr + (ih2-1) * ie_pr(2)
    ir3 = ir3_pr + (ih3-1) * ie_pr(3)
    unfold%eihkr_tbl(ir1,ir2,ir3,ihk) = exp( zI * (unfold%vec_hk(1,ihk)*(ir1-1)*system%hgs(1) &
    &  + unfold%vec_hk(2,ihk)*(ir2-1)*system%hgs(2) + unfold%vec_hk(3,ihk)*(ir3-1)*system%hgs(3) ) )
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
        & 11, "Jz:,u[r,V]u>-d", trim(t_unit_current%name), &
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

! Fourier transform of primitive cell Bloch orbital
  allocate( unfold%psi_prG(1:ie_pr(1),1:ie_pr(2),1:ie_pr(3),nspin,no_pr,isk_s:isk_e,1) )

  ispin = 1
  omega_pr = system%hvol * system%ngrid / dble(unfold%num_hkgrid(1)*unfold%num_hkgrid(2)*unfold%num_hkgrid(3))

  !$omp parallel do private(ilk,ihk,isk,io_pr,ig1_pr,ig2_pr,ig3_pr,zsum,ir1_pr,ir2_pr,ir3_pr) collapse(2)
  do ilk = info%ik_s, info%ik_e ! large k
  do ihk = 1, unfold%nhk   ! hat k
    isk = unfold%isk_tbl(ilk,ihk) !small k = large k + hat k
  do io_pr = 1, no_pr
  do ig1_pr = 1, ie_pr(1)
  do ig2_pr = 1, ie_pr(2)
  do ig3_pr = 1, ie_pr(3)
    zsum = 0d0
    do ir1_pr = 1, ie_pr(1)
    do ir2_pr = 1, ie_pr(2)
    do ir3_pr = 1, ie_pr(3)
      zsum = zsum + unfold%psi_pr(ir1_pr,ir2_pr,ir3_pr,ispin,io_pr,isk,1) &
      & * exp( -2*pi*zI*((ig1_pr-1)*(ir1_pr-1)/dble(ie_pr(1))+(ig2_pr-1)*(ir2_pr-1)/dble(ie_pr(2)) &
      &                 +(ig3_pr-1)*(ir3_pr-1)/dble(ie_pr(3))) )
    end do
    end do
    end do
    unfold%psi_prG(ig1_pr,ig2_pr,ig3_pr,ispin,io_pr,isk,1) = zsum *system%hvol / omega_pr
  end do
  end do
  end do
  end do

  end do
  end do

  allocate( reta_uu(1:ie_pr(1),1:ie_pr(2),1:ie_pr(3),isk_s:isk_e) )
!$omp parallel do private(ilk,ihk,isk,ig1_pr,ig2_pr,ig3_pr,io_pr,rsum) collapse(2)
  do ilk = info%ik_s, info%ik_e
  do ihk = 1, unfold%nhk
    isk = unfold%isk_tbl(ilk,ihk)
  do ig1_pr = 1, ie_pr(1)
  do ig2_pr = 1, ie_pr(2)
  do ig3_pr = 1, ie_pr(3)
    rsum = 0d0
    do io_pr = 1, no_pr
      rsum = rsum + unfold%rocc_pr(io_pr,isk,ispin) * abs(unfold%psi_prG(ig1_pr,ig2_pr,ig3_pr,ispin,io_pr,isk,1))**2
    end do
    reta_uu(ig1_pr,ig2_pr,ig3_pr,isk) = rsum * omega_pr
  enddo
  enddo
  enddo

  enddo
  enddo

  B_pr(:,:) = system%primitive_b(:,:)
  B_pr(:,1) = B_pr(:,1) * unfold%num_hkgrid(1)
  B_pr(:,2) = B_pr(:,2) * unfold%num_hkgrid(2)
  B_pr(:,3) = B_pr(:,3) * unfold%num_hkgrid(3)
  rsum_l = 0d0
  rj_l(1:3) = 0d0
!$omp parallel do private(ilk,ihk,isk,ig1_pr,ig2_pr,ig3_pr,ig1,ig2,ig3,gx,gy,gz) reduction(+:rsum_l,rj_l) collapse(2)
  do ilk = info%ik_s, info%ik_e
  do ihk = 1, unfold%nhk
    isk = unfold%isk_tbl(ilk,ihk)
  do ig1_pr = 1, ie_pr(1)
  do ig2_pr = 1, ie_pr(2)
  do ig3_pr = 1, ie_pr(3)
    ig1 = ig1_pr - 1
    ig2 = ig2_pr - 1
    ig3 = ig3_pr - 1
    if( ig1 > ie_pr(1)/2 ) ig1 = ig1 - ie_pr(1)
    if( ig2 > ie_pr(2)/2 ) ig2 = ig2 - ie_pr(2)
    if( ig3 > ie_pr(3)/2 ) ig3 = ig3 - ie_pr(3)
    gx = ig1*B_pr(1,1) + ig2*B_pr(1,2) + ig3*B_pr(1,3)
    gy = ig1*B_pr(2,1) + ig2*B_pr(2,2) + ig3*B_pr(2,3)
    gz = ig1*B_pr(3,1) + ig2*B_pr(3,2) + ig3*B_pr(3,3)
    rsum_l = rsum_l + unfold%wtk_pr(isk) * reta_uu(ig1_pr,ig2_pr,ig3_pr,isk)
    rj_l(1) = rj_l(1) + unfold%wtk_pr(isk) * reta_uu(ig1_pr,ig2_pr,ig3_pr,isk) &
    & * (gx + system%vec_k(1,ilk) + unfold%vec_hk(1,ihk))
    rj_l(2) = rj_l(2) + unfold%wtk_pr(isk) * reta_uu(ig1_pr,ig2_pr,ig3_pr,isk) &
    & * (gy + system%vec_k(2,ilk) + unfold%vec_hk(2,ihk))
    rj_l(3) = rj_l(3) + unfold%wtk_pr(isk) * reta_uu(ig1_pr,ig2_pr,ig3_pr,isk) &
    & * (gz + system%vec_k(3,ilk) + unfold%vec_hk(3,ihk))
  end do
  end do
  end do

  end do
  end do

  rsum = 0d0
  rj(:) = 0d0
  call comm_summation(rsum_l,rsum,icomm)
  call comm_summation(rj_l,rj,3,icomm)

  rj = rj / omega_pr

  if(comm_is_root(nproc_id_global))then
    write(*,'(A,7x,2f17.12)') 'N:sum rho_uu(k,G)       ',rsum
    write(*,'(A,7x,6f17.12)') 'J:sum rho_uu(k,G)(G+k)  ',rj(1:3)
  end if

  if( yn_out_mom_distr_gs == 'y' ) then

!   grid for momentum distribution, -nq_mom < iq < nq_mom with dq spacing
    if (dq_mom < 1d-9) dq_mom = (((2*pi)**3/system%det_a)/(num_kgrid(1)*num_kgrid(2)*num_kgrid(3)))**(1d0/3d0)
    if (nq_mom <= 0) nq_mom = 2*int((real(num_kgrid(1),kind=8)*num_kgrid(2)*num_kgrid(3))**(1.0d0/3.0d0))
    allocate( nq_l(-nq_mom:nq_mom,-nq_mom:nq_mom,-nq_mom:nq_mom) )
    nq_l = 0.0d0

!$omp parallel default(shared) &
!$omp private(ilk,ihk,isk,ig1_pr,ig2_pr,ig3_pr,ig1,ig2,ig3,gx,gy,gz,qx,qy,qz,iqx,iqy,iqz, &
!$omp dqx,dqy,dqz,value,nq_l_private)
    allocate( nq_l_private(-nq_mom:nq_mom,-nq_mom:nq_mom,-nq_mom:nq_mom) )
    nq_l_private = 0.0d0

!$omp do collapse(2) 
    do ilk = info%ik_s, info%ik_e
    do ihk = 1, unfold%nhk
      isk = unfold%isk_tbl(ilk,ihk)
    do ig1_pr = 1, ie_pr(1)
    do ig2_pr = 1, ie_pr(2)
    do ig3_pr = 1, ie_pr(3)
      ig1 = ig1_pr - 1
      ig2 = ig2_pr - 1
      ig3 = ig3_pr - 1
      if( ig1 > ie_pr(1)/2 ) ig1 = ig1 - ie_pr(1)
      if( ig2 > ie_pr(2)/2 ) ig2 = ig2 - ie_pr(2)
      if( ig3 > ie_pr(3)/2 ) ig3 = ig3 - ie_pr(3)
      gx = ig1*B_pr(1,1) + ig2*B_pr(1,2) + ig3*B_pr(1,3)
      gy = ig1*B_pr(2,1) + ig2*B_pr(2,2) + ig3*B_pr(2,3)
      gz = ig1*B_pr(3,1) + ig2*B_pr(3,2) + ig3*B_pr(3,3)
      qx = system%vec_k(1,ilk) + unfold%vec_hk(1,ihk) + gx
      qy = system%vec_k(2,ilk) + unfold%vec_hk(2,ihk) + gy
      qz = system%vec_k(3,ilk) + unfold%vec_hk(3,ihk) + gz
      iqx = floor( qx/dq_mom )
      iqy = floor( qy/dq_mom )
      iqz = floor( qz/dq_mom )
      if( iqx < -nq_mom .or. iqx >= nq_mom ) cycle
      if( iqy < -nq_mom .or. iqy >= nq_mom ) cycle
      if( iqz < -nq_mom .or. iqz >= nq_mom ) cycle
      dqx = qx/dq_mom - iqx
      dqy = qy/dq_mom - iqy
      dqz = qz/dq_mom - iqz
      value = unfold%wtk_pr(isk) * reta_uu(ig1_pr,ig2_pr,ig3_pr,isk) / dq_mom**3
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
    rj = rj * dq_mom**3/omega_pr

    if(comm_is_root(nproc_id_global))then
      write(*,'(A,7x,f17.12)')  'N:int n(q,t)            ',nq_sum
      write(*,'(A,7x,3f17.12)') 'J:int q n(q,t)          ',rj(1:3)
    end if

    if(comm_is_root(nproc_id_global))then
      iofile = trim(base_directory)//trim(sysname)//"_momgs.cube"
      fp = 271
      open(fp,file=iofile)
      write(fp,*) "momentum distribution: primitive cell, ground state"
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

  end subroutine init_dm_unfold

!===================================================================================================================================
  subroutine dm_unfold(itt,system,info,lg,ofl,psi_t,unfold)

    use structures
    use communication, only: comm_is_root, comm_summation
    use parallelization, only: nproc_id_global
    use salmon_global, only: no_pr, dt, num_kgrid, num_skgrid, sysname, base_directory, natom,izatom,kion, &
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
    integer ie_pr(3),ilk,ihk,isk,io,io_pr,ih1,ih2,ih3,ir1_pr,ir2_pr,ir3_pr
    integer ir1,ir2,ir3,ig1,ig2,ig3,ig1_pr,ig2_pr,ig3_pr,ispin,ig_pr(3),iqx,iqy,iqz
    integer io_pr1,io_pr2, icomm, isk_s, isk_e, nsk_se
    real(8) :: B_pr(3,3),gx,gy,gz,omega_pr,dqx,dqy,dqz,qx,qy,qz,value
    complex(8),allocatable :: mat(:,:,:,:), eta(:,:,:), eta_l(:,:,:)
    complex(8),allocatable :: eta_uu_d(:,:,:,:),eta_uu_nd(:,:,:,:)
    complex(8),allocatable :: nq_d_l(:,:,:),nq_d(:,:,:),nq_nd_l(:,:,:),nq_nd(:,:,:),nq_d_l_private(:,:,:),nq_nd_l_private(:,:,:)
    complex(8) :: zj1(3),zj2(3),zj3(3),zj4(3),zj5(3),zj6(3),zj1_uu(3),zj2_uu(3),zj3_uu(3),zj4_uu(3),zj_d(3),zj_nd(3)
    complex(8) :: zj1_l(3),zj2_l(3),zj3_l(3),zj4_l(3),zj5_l(3),zj6_l(3),zsum,zsum_d,zsum_nd,zsum_l,zsum_uu
      
    allocate( mat(no_pr,unfold%nhk,info%io_s:info%io_e,info%ik_s:info%ik_e))
    ie_pr(1:3) = lg%ie(1:3)/unfold%num_hkgrid(1:3)
    omega_pr = system%hvol * system%ngrid / dble(unfold%num_hkgrid(1)*unfold%num_hkgrid(2)*unfold%num_hkgrid(3))

    ispin = 1
    isk_s = (info%ik_s-1) * unfold%nhk + 1
    isk_e = info%ik_e * unfold%nhk
    nsk_se = isk_e - isk_s + 1

  !$omp parallel do private(ilk,ihk,isk,io,io_pr,zsum,ih1,ih2,ih3,ir1_pr,ir2_pr,ir3_pr,ir1,ir2,ir3) collapse(2)
    do ilk = info%ik_s, info%ik_e ! large k
    do ihk = 1, unfold%nhk   ! hat k
      isk = unfold%isk_tbl(ilk,ihk) !small k = large k + hat k
    do io = info%io_s, info%io_e     ! m, supercell
    do io_pr = 1, no_pr   ! n, primitive
      zsum = 0d0
      do ih1 = 1, unfold%num_hkgrid(1)
      do ih2 = 1, unfold%num_hkgrid(2)
      do ih3 = 1, unfold%num_hkgrid(3)
      do ir1_pr = 1, ie_pr(1)
      do ir2_pr = 1, ie_pr(2)
      do ir3_pr = 1, ie_pr(3)
        ir1 = ir1_pr + (ih1-1) * ie_pr(1)
        ir2 = ir2_pr + (ih2-1) * ie_pr(2)
        ir3 = ir3_pr + (ih3-1) * ie_pr(3)
        zsum = zsum + conjg( unfold%psi_pr( ir1_pr, ir2_pr, ir3_pr, 1, io_pr, isk, 1) ) &
         &        * conjg( unfold%eihkr_tbl(ir1,ir2,ir3,ihk) ) &
         &        * psi_t%zwf( ir1, ir2, ir3, 1, io, ilk, 1 )
      end do
      end do
      end do
      end do
      end do
      end do
      zsum = zsum * system%hvol
      mat(io_pr, ihk, io, ilk) = zsum / (unfold%num_hkgrid(1)*unfold%num_hkgrid(2)*unfold%num_hkgrid(3))
    end do
    end do
    end do
    end do

    allocate( eta(no_pr,no_pr,isk_s:isk_e),eta_l(no_pr,no_pr,isk_s:isk_e) )
    eta_l = 0.0d0
  !$omp parallel do private(ilk,ihk,isk,io_pr1,io_pr2,zsum,io) collapse(2)
    do ilk = info%ik_s, info%ik_e
    do ihk = 1, unfold%nhk
      isk = unfold%isk_tbl(ilk,ihk)
      do io_pr1 = 1, no_pr
      do io_pr2 = 1, no_pr
        zsum = 0d0
        do io = info%io_s, info%io_e
          zsum = zsum + system%rocc(io,ilk,1) * mat(io_pr1, ihk, io, ilk) * conjg( mat(io_pr2, ihk, io, ilk) )
        end do
        eta_l(io_pr1, io_pr2, isk) = zsum
      end do
      end do
    end do
    end do
    eta_l = eta_l * (unfold%num_hkgrid(1)*unfold%num_hkgrid(2)*unfold%num_hkgrid(3))
    eta = 0.0d0
    call comm_summation(eta_l,eta,no_pr*no_pr*nsk_se,info%icomm_o)

    zsum_l = 0d0
    zj1_l(1:3) = 0d0
    zj2_l(1:3) = 0d0
    zj3_l(1:3) = 0d0
    zj4_l(1:3) = 0d0
    zj5_l(1:3) = 0d0
    zj6_l(1:3) = 0d0
  !$omp parallel do private(ilk,ihk,isk,io_pr1,io_pr2) reduction(+:zsum_l,zj1_l,zj2_l,zj3_l,zj4_l,zj5_l,zj6_l) collapse(2)
    do ilk = info%ik_s, info%ik_e
    do ihk = 1, unfold%nhk
       isk = unfold%isk_tbl(ilk,ihk)
    do io_pr1 = 1, no_pr
       zsum_l = zsum_l + eta(io_pr1, io_pr1, isk) * unfold%wtk_pr(isk)
       zj1_l(:) = zj1_l(:) + eta(io_pr1, io_pr1, isk) * unfold%wtk_pr(isk) &
        & * unfold%upu_pr(:, io_pr1, io_pr1, isk)
       zj2_l(:) = zj2_l(:) + eta(io_pr1, io_pr1, isk) * unfold%wtk_pr(isk) &
        & * unfold%u_rVnl_Vnlr_u_pr(:, io_pr1, io_pr1, isk)
       zj5_l(:) = zj5_l(:) + eta(io_pr1, io_pr1, isk) * unfold%wtk_pr(isk) &
        & * (system%vec_k(:,ilk) + unfold%vec_hk(:,ihk))
       zj6_l(:) = zj6_l(:) + eta(io_pr1, io_pr1, isk) * unfold%wtk_pr(isk) &
        & * system%vec_Ac(:)
    do io_pr2 = 1, no_pr
       if( io_pr1 /= io_pr2) then
         zj3_l(:) = zj3_l(:) + eta(io_pr2, io_pr1, isk) * unfold%wtk_pr(isk) &
        & * unfold%upu_pr(:, io_pr1, io_pr2, isk)
         zj4_l(:) = zj4_l(:) + eta(io_pr2, io_pr1, isk) * unfold%wtk_pr(isk) &
        & * unfold%u_rVnl_Vnlr_u_pr(:, io_pr1, io_pr2, isk)
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

    zj1 = zj1 / omega_pr
    zj2 = zj2 / omega_pr
    zj3 = zj3 / omega_pr
    zj4 = zj4 / omega_pr
    zj5 = zj5 / omega_pr
    zj6 = zj6 / omega_pr

  if(comm_is_root(nproc_id_global))then
    write(*,'(A,2x,i7,2x,A,3f17.12)') 'dm_unfold  it=', itt, '     Ac(t)=',system%vec_Ac(:)
    write(*,'(A,7x,f17.12)')          'N:Tr[rho(t)]             ', real(zsum)
    write(*,'(A,7x,3f17.12)')         'J:rho(t)<unk|i[h,r]|unk> ', real(zj1(:)+zj2(:)+zj3(:)+zj4(:)+zj5(:)+zj6(:))
  end if

  allocate( eta_uu_d(1:ie_pr(1),1:ie_pr(2),1:ie_pr(3),isk_s:isk_e) )
  allocate( eta_uu_nd(1:ie_pr(1),1:ie_pr(2),1:ie_pr(3),isk_s:isk_e) )
!$omp parallel do private(ilk,ihk,isk,ig1_pr,ig2_pr,ig3_pr,io_pr1,io_pr2,zsum_d,zsum_nd) collapse(2)
  do ilk = info%ik_s, info%ik_e
  do ihk = 1, unfold%nhk
    isk = unfold%isk_tbl(ilk,ihk)
  do ig1_pr = 1, ie_pr(1)
  do ig2_pr = 1, ie_pr(2)
  do ig3_pr = 1, ie_pr(3)
    zsum_d = 0d0
    zsum_nd = 0d0
    do io_pr1 = 1, no_pr
    do io_pr2 = 1, no_pr
      if( io_pr1 == io_pr2 ) then
        zsum_d = zsum_d + eta(io_pr1, io_pr2, isk) * unfold%psi_prG(ig1_pr,ig2_pr,ig3_pr,ispin,io_pr1,isk,1) &
        &                               * conjg( unfold%psi_prG(ig1_pr,ig2_pr,ig3_pr,ispin,io_pr2,isk,1) )
      else
        zsum_nd = zsum_nd + eta(io_pr1, io_pr2, isk) * unfold%psi_prG(ig1_pr,ig2_pr,ig3_pr,ispin,io_pr1,isk,1) &
        &                               * conjg( unfold%psi_prG(ig1_pr,ig2_pr,ig3_pr,ispin,io_pr2,isk,1) )
      end if
    end do
    end do
    eta_uu_d(ig1_pr,ig2_pr,ig3_pr,isk) = zsum_d * omega_pr
    eta_uu_nd(ig1_pr,ig2_pr,ig3_pr,isk) = zsum_nd * omega_pr
  enddo
  enddo
  enddo

  enddo
  enddo

  B_pr(:,:) = system%primitive_b(:,:)
  B_pr(:,1) = B_pr(:,1) * unfold%num_hkgrid(1)
  B_pr(:,2) = B_pr(:,2) * unfold%num_hkgrid(2)
  B_pr(:,3) = B_pr(:,3) * unfold%num_hkgrid(3)
  zsum_l = 0d0
  zj1_l(1:3) = 0d0
  zj2_l(1:3) = 0d0
  zj3_l(1:3) = 0d0
  zj4_l(1:3) = 0d0
!$omp parallel do private(ilk,ihk,isk,ig1_pr,ig2_pr,ig3_pr,ig1,ig2,ig3,gx,gy,gz) &
!$omp reduction(+:zsum_l,zj1_l,zj2_l,zj3_l,zj4_l) collapse(2)
  do ilk = info%ik_s, info%ik_e
  do ihk = 1, unfold%nhk
    isk = unfold%isk_tbl(ilk,ihk)
  do ig1_pr = 1, ie_pr(1)
  do ig2_pr = 1, ie_pr(2)
  do ig3_pr = 1, ie_pr(3)
    ig1 = ig1_pr - 1
    ig2 = ig2_pr - 1
    ig3 = ig3_pr - 1
    if( ig1 > ie_pr(1)/2 ) ig1 = ig1 - ie_pr(1)
    if( ig2 > ie_pr(2)/2 ) ig2 = ig2 - ie_pr(2)
    if( ig3 > ie_pr(3)/2 ) ig3 = ig3 - ie_pr(3)
    gx = ig1*B_pr(1,1) + ig2*B_pr(1,2) + ig3*B_pr(1,3)
    gy = ig1*B_pr(2,1) + ig2*B_pr(2,2) + ig3*B_pr(2,3)
    gz = ig1*B_pr(3,1) + ig2*B_pr(3,2) + ig3*B_pr(3,3)
    zsum_l = zsum_l + unfold%wtk_pr(isk) * eta_uu_d(ig1_pr,ig2_pr,ig3_pr,isk)
    zj1_l(1) = zj1_l(1) + unfold%wtk_pr(isk) * eta_uu_d(ig1_pr,ig2_pr,ig3_pr,isk) * gx
    zj1_l(2) = zj1_l(2) + unfold%wtk_pr(isk) * eta_uu_d(ig1_pr,ig2_pr,ig3_pr,isk) * gy
    zj1_l(3) = zj1_l(3) + unfold%wtk_pr(isk) * eta_uu_d(ig1_pr,ig2_pr,ig3_pr,isk) * gz
    zj2_l(1) = zj2_l(1) + unfold%wtk_pr(isk) * eta_uu_nd(ig1_pr,ig2_pr,ig3_pr,isk) * gx
    zj2_l(2) = zj2_l(2) + unfold%wtk_pr(isk) * eta_uu_nd(ig1_pr,ig2_pr,ig3_pr,isk) * gy
    zj2_l(3) = zj2_l(3) + unfold%wtk_pr(isk) * eta_uu_nd(ig1_pr,ig2_pr,ig3_pr,isk) * gz
    zj3_l(:) = zj3_l(:) + unfold%wtk_pr(isk) * eta_uu_d(ig1_pr,ig2_pr,ig3_pr,isk) &
        & * (system%vec_k(:,ilk) + unfold%vec_hk(:,ihk))
    zj4_l(:) = zj4_l(:) + unfold%wtk_pr(isk) * eta_uu_d(ig1_pr,ig2_pr,ig3_pr,isk) * system%vec_Ac(:)
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

  zj1_uu = zj1_uu / omega_pr
  zj2_uu = zj2_uu / omega_pr
  zj3_uu = zj3_uu / omega_pr
  zj4_uu = zj4_uu / omega_pr

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
!$omp private(ilk,ihk,isk,ig1_pr,ig2_pr,ig3_pr,ig1,ig2,ig3,gx,gy,gz,qx,qy,qz,iqx,iqy,iqz, &
!$omp dqx,dqy,dqz,value,nq_d_l_private,nq_nd_l_private)
    allocate( nq_d_l_private(-nq_mom:nq_mom,-nq_mom:nq_mom,-nq_mom:nq_mom) )
    allocate( nq_nd_l_private(-nq_mom:nq_mom,-nq_mom:nq_mom,-nq_mom:nq_mom) )
    nq_d_l_private = 0.0d0
    nq_nd_l_private = 0.0d0

!$omp do collapse(2) 
    do ilk = info%ik_s, info%ik_e
    do ihk = 1, unfold%nhk
      isk = unfold%isk_tbl(ilk,ihk)
    do ig1_pr = 1, ie_pr(1)
    do ig2_pr = 1, ie_pr(2)
    do ig3_pr = 1, ie_pr(3)
      ig1 = ig1_pr - 1
      ig2 = ig2_pr - 1
      ig3 = ig3_pr - 1
      if( ig1 > ie_pr(1)/2 ) ig1 = ig1 - ie_pr(1)
      if( ig2 > ie_pr(2)/2 ) ig2 = ig2 - ie_pr(2)
      if( ig3 > ie_pr(3)/2 ) ig3 = ig3 - ie_pr(3)
      gx = ig1*B_pr(1,1) + ig2*B_pr(1,2) + ig3*B_pr(1,3)
      gy = ig1*B_pr(2,1) + ig2*B_pr(2,2) + ig3*B_pr(2,3)
      gz = ig1*B_pr(3,1) + ig2*B_pr(3,2) + ig3*B_pr(3,3)
      qx = system%vec_Ac(1) + system%vec_k(1,ilk) + unfold%vec_hk(1,ihk) + gx
      qy = system%vec_Ac(2) + system%vec_k(2,ilk) + unfold%vec_hk(2,ihk) + gy
      qz = system%vec_Ac(3) + system%vec_k(3,ilk) + unfold%vec_hk(3,ihk) + gz
      iqx = floor( qx/dq_mom )
      iqy = floor( qy/dq_mom )
      iqz = floor( qz/dq_mom )
      if( iqx < -nq_mom .or. iqx >= nq_mom ) cycle
      if( iqy < -nq_mom .or. iqy >= nq_mom ) cycle
      if( iqz < -nq_mom .or. iqz >= nq_mom ) cycle
      dqx = qx/dq_mom - iqx
      dqy = qy/dq_mom - iqy
      dqz = qz/dq_mom - iqz

      value = unfold%wtk_pr(isk) * eta_uu_d(ig1_pr,ig2_pr,ig3_pr,isk) / dq_mom**3
      nq_d_l_private(iqx+1,iqy+1,iqz+1) = nq_d_l_private(iqx+1,iqy+1,iqz+1) + dqx    *dqy    *    dqz * value
      nq_d_l_private(iqx+1,iqy+1,iqz  ) = nq_d_l_private(iqx+1,iqy+1,iqz  ) + dqx    *dqy    *(1-dqz) * value
      nq_d_l_private(iqx+1,iqy,  iqz+1) = nq_d_l_private(iqx+1,iqy,  iqz+1) + dqx    *(1-dqy)*    dqz * value
      nq_d_l_private(iqx  ,iqy+1,iqz+1) = nq_d_l_private(iqx  ,iqy+1,iqz+1) + (1-dqx)*dqy    *    dqz * value
      nq_d_l_private(iqx+1,iqy  ,iqz  ) = nq_d_l_private(iqx+1,iqy  ,iqz  ) + dqx    *(1-dqy)*(1-dqz) * value
      nq_d_l_private(iqx  ,iqy+1,iqz  ) = nq_d_l_private(iqx  ,iqy+1,  iqz) + (1-dqx)*dqy    *(1-dqz) * value
      nq_d_l_private(iqx  ,iqy  ,iqz+1) = nq_d_l_private(iqx  ,iqy  ,iqz+1) + (1-dqx)*(1-dqy)*    dqz * value
      nq_d_l_private(iqx  ,iqy  ,iqz  ) = nq_d_l_private(iqx  ,iqy  ,iqz  ) + (1-dqx)*(1-dqy)*(1-dqz) * value

      value = unfold%wtk_pr(isk) * eta_uu_nd(ig1_pr,ig2_pr,ig3_pr,isk) / dq_mom**3
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
    zj_d = zj_d * dq_mom**3/omega_pr
    zj_nd = zj_nd * dq_mom**3/omega_pr

  if(comm_is_root(nproc_id_global))then
    write(*,'(A,7x,f17.12)')          'J:int n(q,t)             ', real(zsum_d)
    write(*,'(A,7x,3f17.12)')         'J:int q n(q,t)           ', real(zj_d(:)+zj_nd(:))
  end if

    if(comm_is_root(nproc_id_global))then
      write(filenum, '(i6.6)') itt
      iofile = trim(base_directory)//trim(sysname)//"_momdiff_"//trim(adjustl(filenum))//".cube"
      fp = 271
      open(fp,file=iofile)
      write(fp,*) "momentum distribution: primitive cell, difference from ground state"
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
    & real(zj4(1:3))*t_unit_current%conv,real(zj5(1:3))*t_unit_current%conv,real(zj6(1:3)*t_unit_current%conv), &
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
