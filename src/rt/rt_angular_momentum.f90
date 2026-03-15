module rt_angular_momentum
  implicit none
  logical, save :: lz_state_initialized = .false.
  integer, save :: lz_num_x = 0, lz_num_y = 0
  integer, save :: lz_step_prev2 = -1, lz_step_prev1 = -1
  real(8), save :: lz_time_prev2 = 0d0, lz_time_prev1 = 0d0
  real(8), save :: lz_total_prev2 = 0d0, lz_total_prev1 = 0d0
  real(8), allocatable, save :: lz_xy_prev2(:,:), lz_xy_prev1(:,:)

contains

subroutine write_local_angular_momentum_xy(itt, lg, mg, system, info, singlescale)
  use structures, only: s_rgrid, s_dft_system, s_parallel_info, s_singlescale
  use communication, only: comm_summation
  use salmon_global, only: dt, optical_vortex_center_x, optical_vortex_center_y
  implicit none
  integer, intent(in) :: itt
  type(s_rgrid), intent(in) :: lg, mg
  type(s_dft_system), intent(in) :: system
  type(s_parallel_info), intent(in) :: info
  type(s_singlescale), intent(in) :: singlescale
  real(8), allocatable :: lz_xy_local(:,:), lz_xy_global(:,:)
  real(8) :: center_x, center_y, x, y, lz_local, lz_total_local, lz_total, time_now
  integer :: ix, iy, iz

  allocate(lz_xy_local(1:lg%num(1),1:lg%num(2)))
  allocate(lz_xy_global(1:lg%num(1),1:lg%num(2)))
  lz_xy_local = 0.0d0
  lz_xy_global = 0.0d0
  lz_total_local = 0.0d0
  time_now = dble(itt) * dt

  center_x = optical_vortex_center_x
  center_y = optical_vortex_center_y
  if (center_x < -1d20) center_x = 0.5d0 * lg%num(1) * system%hgs(1)
  if (center_y < -1d20) center_y = 0.5d0 * lg%num(2) * system%hgs(2)

  do iz = mg%is(3), mg%ie(3)
    do iy = mg%is(2), mg%ie(2)
      y = (dble(iy) - 0.5d0) * system%hgs(2) - center_y
      do ix = mg%is(1), mg%ie(1)
        x = (dble(ix) - 0.5d0) * system%hgs(1) - center_x
        lz_local = x * singlescale%curr(ix,iy,iz,2) - y * singlescale%curr(ix,iy,iz,1)
        lz_xy_local(ix,iy) = lz_xy_local(ix,iy) + lz_local * system%hgs(3)
        lz_total_local = lz_total_local + lz_local * system%Hvol
      end do
    end do
  end do

  call comm_summation(lz_xy_local, lz_xy_global, lg%num(1) * lg%num(2), info%icomm_r)
  call comm_summation(lz_total_local, lz_total, 1, info%icomm_r)

  if (.not. lz_state_initialized) then
    lz_num_x = lg%num(1)
    lz_num_y = lg%num(2)
    allocate(lz_xy_prev2(1:lz_num_x,1:lz_num_y))
    allocate(lz_xy_prev1(1:lz_num_x,1:lz_num_y))
    lz_xy_prev2 = lz_xy_global
    lz_step_prev2 = itt
    lz_time_prev2 = time_now
    lz_total_prev2 = lz_total
    lz_state_initialized = .true.
  else if (lz_step_prev1 < 0) then
    call dump_lz_record(lz_step_prev2, lz_time_prev2, lz_xy_prev2, &
      & (lz_xy_global - lz_xy_prev2) / max(time_now - lz_time_prev2, tiny(1d0)), &
      & lz_total_prev2, (lz_total - lz_total_prev2) / max(time_now - lz_time_prev2, tiny(1d0)), &
      & center_x, center_y, system)
    lz_xy_prev1 = lz_xy_global
    lz_step_prev1 = itt
    lz_time_prev1 = time_now
    lz_total_prev1 = lz_total
  else
    call dump_lz_record(lz_step_prev1, lz_time_prev1, lz_xy_prev1, &
      & (lz_xy_global - lz_xy_prev2) / max(time_now - lz_time_prev2, tiny(1d0)), &
      & lz_total_prev1, (lz_total - lz_total_prev2) / max(time_now - lz_time_prev2, tiny(1d0)), &
      & center_x, center_y, system)
    lz_xy_prev2 = lz_xy_prev1
    lz_step_prev2 = lz_step_prev1
    lz_time_prev2 = lz_time_prev1
    lz_total_prev2 = lz_total_prev1
    lz_xy_prev1 = lz_xy_global
    lz_step_prev1 = itt
    lz_time_prev1 = time_now
    lz_total_prev1 = lz_total
  end if

  deallocate(lz_xy_local, lz_xy_global)
end subroutine write_local_angular_momentum_xy

subroutine flush_local_angular_momentum_xy(system)
  use structures, only: s_dft_system
  use salmon_global, only: optical_vortex_center_x, optical_vortex_center_y
  implicit none
  type(s_dft_system), intent(in) :: system
  real(8) :: center_x, center_y

  if (.not. lz_state_initialized) return

  center_x = optical_vortex_center_x
  center_y = optical_vortex_center_y
  if (center_x < -1d20) center_x = 0.5d0 * lz_num_x * system%hgs(1)
  if (center_y < -1d20) center_y = 0.5d0 * lz_num_y * system%hgs(2)

  if (lz_step_prev1 >= 0) then
    call dump_lz_record(lz_step_prev1, lz_time_prev1, lz_xy_prev1, &
      & (lz_xy_prev1 - lz_xy_prev2) / max(lz_time_prev1 - lz_time_prev2, tiny(1d0)), &
      & lz_total_prev1, (lz_total_prev1 - lz_total_prev2) / max(lz_time_prev1 - lz_time_prev2, tiny(1d0)), &
      & center_x, center_y, system)
  else
    call dump_lz_record(lz_step_prev2, lz_time_prev2, lz_xy_prev2, 0d0 * lz_xy_prev2, &
      & lz_total_prev2, 0d0, center_x, center_y, system)
  end if
end subroutine flush_local_angular_momentum_xy

subroutine dump_lz_record(itt, time_now, lz_xy, dlzdt_xy, lz_total, dlzdt_total, center_x, center_y, system)
  use structures, only: s_dft_system
  use communication, only: comm_is_root
  use parallelization, only: nproc_id_global
  use salmon_global, only: base_directory, sysname
  implicit none
  integer, intent(in) :: itt
  real(8), intent(in) :: time_now, lz_xy(:,:), dlzdt_xy(:,:), lz_total, dlzdt_total, center_x, center_y
  type(s_dft_system), intent(in) :: system
  character(256) :: file_map, file_total, filenum
  integer :: ix, iy, iunit
  real(8) :: x, y

  if (.not. comm_is_root(nproc_id_global)) return

  write(filenum, '(i6.6)') itt
  file_map = trim(base_directory)//trim(sysname)//'_lz_xy_'//trim(adjustl(filenum))//'.data'
  open(newunit=iunit, file=trim(file_map), status='replace', action='write')
  write(iunit,'(a)') '# x y local_lz_zint dlocal_lz_dt_zint'
  do iy = 1, size(lz_xy, 2)
    y = (dble(iy) - 0.5d0) * system%hgs(2) - center_y
    do ix = 1, size(lz_xy, 1)
      x = (dble(ix) - 0.5d0) * system%hgs(1) - center_x
      write(iunit,'(4(1x,es24.16))') x, y, lz_xy(ix,iy), dlzdt_xy(ix,iy)
    end do
    write(iunit,*)
  end do
  close(iunit)

  file_total = trim(base_directory)//trim(sysname)//'_lz_total.data'
  if (itt == 0) then
    open(newunit=iunit, file=trim(file_total), status='replace', action='write')
    write(iunit,'(a)') '# step time local_lz_total dlocal_lz_total_dt'
  else
    open(newunit=iunit, file=trim(file_total), status='old', position='append', action='write')
  end if
  write(iunit,'(i10,1x,3(es24.16,1x))') itt, time_now, lz_total, dlzdt_total
  close(iunit)
end subroutine dump_lz_record

end module rt_angular_momentum
