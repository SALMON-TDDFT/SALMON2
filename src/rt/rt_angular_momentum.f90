module rt_angular_momentum
  implicit none
  logical, save :: angmom_state_initialized = .false.
  integer, save :: angmom_num_x = 0, angmom_num_y = 0
  integer, save :: angmom_step_prev2 = -1, angmom_step_prev1 = -1
  real(8), save :: angmom_time_prev2 = 0d0, angmom_time_prev1 = 0d0
  real(8), save :: lz_total_prev2 = 0d0, lz_total_prev1 = 0d0
  real(8), save :: sz_total_prev2 = 0d0, sz_total_prev1 = 0d0
  real(8), save :: jz_total_prev2 = 0d0, jz_total_prev1 = 0d0
  real(8), allocatable, save :: lz_xy_prev2(:,:), lz_xy_prev1(:,:)
  real(8), allocatable, save :: sz_xy_prev2(:,:), sz_xy_prev1(:,:)
  real(8), allocatable, save :: jz_xy_prev2(:,:), jz_xy_prev1(:,:)

contains

subroutine write_local_angular_momentum_xy(itt, lg, mg, system, info, singlescale, psi)
  use structures, only: s_rgrid, s_dft_system, s_parallel_info, s_singlescale, s_orbital
  use communication, only: comm_summation
  use salmon_global, only: dt, optical_vortex_center_x, optical_vortex_center_y, yn_spinorbit
  use noncollinear_module, only: calc_dm_noncollinear, calc_magnetization_micro
  implicit none
  integer, intent(in) :: itt
  type(s_rgrid), intent(in) :: lg, mg
  type(s_dft_system), intent(in) :: system
  type(s_parallel_info), intent(in) :: info
  type(s_singlescale), intent(in) :: singlescale
  type(s_orbital), intent(in) :: psi
  real(8), allocatable :: lz_xy_local(:,:), lz_xy_global(:,:)
  real(8), allocatable :: sz_xy_local(:,:), sz_xy_global(:,:)
  real(8), allocatable :: jz_xy_global(:,:), m_micro(:,:,:,:)
  real(8) :: center_x, center_y, x, y, lz_local, sz_local, time_now
  real(8) :: lz_total_local, lz_total, sz_total_local, sz_total, jz_total
  integer :: ix, iy, iz

  allocate(lz_xy_local(1:lg%num(1),1:lg%num(2)))
  allocate(lz_xy_global(1:lg%num(1),1:lg%num(2)))
  allocate(sz_xy_local(1:lg%num(1),1:lg%num(2)))
  allocate(sz_xy_global(1:lg%num(1),1:lg%num(2)))
  allocate(jz_xy_global(1:lg%num(1),1:lg%num(2)))
  lz_xy_local = 0.0d0
  lz_xy_global = 0.0d0
  sz_xy_local = 0.0d0
  sz_xy_global = 0.0d0
  jz_xy_global = 0.0d0
  lz_total_local = 0.0d0
  sz_total_local = 0.0d0
  time_now = dble(itt) * dt

  center_x = optical_vortex_center_x
  center_y = optical_vortex_center_y
  if (center_x < -1d20) center_x = 0.5d0 * lg%num(1) * system%hgs(1)
  if (center_y < -1d20) center_y = 0.5d0 * lg%num(2) * system%hgs(2)

  if (yn_spinorbit == 'y') then
    allocate(m_micro(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), 3))
    call calc_dm_noncollinear(psi, system, info, mg)
    call calc_magnetization_micro(mg, m_micro)
  end if

  do iz = mg%is(3), mg%ie(3)
    do iy = mg%is(2), mg%ie(2)
      y = (dble(iy) - 0.5d0) * system%hgs(2) - center_y
      do ix = mg%is(1), mg%ie(1)
        x = (dble(ix) - 0.5d0) * system%hgs(1) - center_x
        lz_local = x * singlescale%curr(ix,iy,iz,2) - y * singlescale%curr(ix,iy,iz,1)
        sz_local = 0.0d0
        if (yn_spinorbit == 'y') sz_local = m_micro(ix,iy,iz,3)
        lz_xy_local(ix,iy) = lz_xy_local(ix,iy) + lz_local * system%hgs(3)
        sz_xy_local(ix,iy) = sz_xy_local(ix,iy) + sz_local * system%hgs(3)
        lz_total_local = lz_total_local + lz_local * system%Hvol
        sz_total_local = sz_total_local + sz_local * system%Hvol
      end do
    end do
  end do

  call comm_summation(lz_xy_local, lz_xy_global, lg%num(1) * lg%num(2), info%icomm_r)
  call comm_summation(lz_total_local, lz_total, 1, info%icomm_r)
  call comm_summation(sz_xy_local, sz_xy_global, lg%num(1) * lg%num(2), info%icomm_r)
  call comm_summation(sz_total_local, sz_total, 1, info%icomm_r)
  jz_xy_global = lz_xy_global + sz_xy_global
  jz_total = lz_total + sz_total

  if (.not. angmom_state_initialized) then
    angmom_num_x = lg%num(1)
    angmom_num_y = lg%num(2)
    allocate(lz_xy_prev2(1:angmom_num_x,1:angmom_num_y))
    allocate(lz_xy_prev1(1:angmom_num_x,1:angmom_num_y))
    allocate(sz_xy_prev2(1:angmom_num_x,1:angmom_num_y))
    allocate(sz_xy_prev1(1:angmom_num_x,1:angmom_num_y))
    allocate(jz_xy_prev2(1:angmom_num_x,1:angmom_num_y))
    allocate(jz_xy_prev1(1:angmom_num_x,1:angmom_num_y))
    lz_xy_prev2 = lz_xy_global
    sz_xy_prev2 = sz_xy_global
    jz_xy_prev2 = jz_xy_global
    angmom_step_prev2 = itt
    angmom_time_prev2 = time_now
    lz_total_prev2 = lz_total
    sz_total_prev2 = sz_total
    jz_total_prev2 = jz_total
    angmom_state_initialized = .true.
  else if (angmom_step_prev1 < 0) then
    call dump_lz_record(angmom_step_prev2, angmom_time_prev2, lz_xy_prev2, sz_xy_prev2, jz_xy_prev2, &
      & (jz_xy_global - jz_xy_prev2) / max(time_now - angmom_time_prev2, tiny(1d0)), &
      & lz_total_prev2, sz_total_prev2, jz_total_prev2, &
      & (jz_total - jz_total_prev2) / max(time_now - angmom_time_prev2, tiny(1d0)), &
      & center_x, center_y, system)
    lz_xy_prev1 = lz_xy_global
    sz_xy_prev1 = sz_xy_global
    jz_xy_prev1 = jz_xy_global
    angmom_step_prev1 = itt
    angmom_time_prev1 = time_now
    lz_total_prev1 = lz_total
    sz_total_prev1 = sz_total
    jz_total_prev1 = jz_total
  else
    call dump_lz_record(angmom_step_prev1, angmom_time_prev1, lz_xy_prev1, sz_xy_prev1, jz_xy_prev1, &
      & (jz_xy_global - jz_xy_prev2) / max(time_now - angmom_time_prev2, tiny(1d0)), &
      & lz_total_prev1, sz_total_prev1, jz_total_prev1, &
      & (jz_total - jz_total_prev2) / max(time_now - angmom_time_prev2, tiny(1d0)), &
      & center_x, center_y, system)
    lz_xy_prev2 = lz_xy_prev1
    sz_xy_prev2 = sz_xy_prev1
    jz_xy_prev2 = jz_xy_prev1
    angmom_step_prev2 = angmom_step_prev1
    angmom_time_prev2 = angmom_time_prev1
    lz_total_prev2 = lz_total_prev1
    sz_total_prev2 = sz_total_prev1
    jz_total_prev2 = jz_total_prev1
    lz_xy_prev1 = lz_xy_global
    sz_xy_prev1 = sz_xy_global
    jz_xy_prev1 = jz_xy_global
    angmom_step_prev1 = itt
    angmom_time_prev1 = time_now
    lz_total_prev1 = lz_total
    sz_total_prev1 = sz_total
    jz_total_prev1 = jz_total
  end if

  if (allocated(m_micro)) deallocate(m_micro)
  deallocate(lz_xy_local, lz_xy_global, sz_xy_local, sz_xy_global, jz_xy_global)
end subroutine write_local_angular_momentum_xy

subroutine flush_local_angular_momentum_xy(system)
  use structures, only: s_dft_system
  use salmon_global, only: optical_vortex_center_x, optical_vortex_center_y
  implicit none
  type(s_dft_system), intent(in) :: system
  real(8) :: center_x, center_y

  if (.not. angmom_state_initialized) return

  center_x = optical_vortex_center_x
  center_y = optical_vortex_center_y
  if (center_x < -1d20) center_x = 0.5d0 * angmom_num_x * system%hgs(1)
  if (center_y < -1d20) center_y = 0.5d0 * angmom_num_y * system%hgs(2)

  if (angmom_step_prev1 >= 0) then
    call dump_lz_record(angmom_step_prev1, angmom_time_prev1, lz_xy_prev1, sz_xy_prev1, jz_xy_prev1, &
      & (jz_xy_prev1 - jz_xy_prev2) / max(angmom_time_prev1 - angmom_time_prev2, tiny(1d0)), &
      & lz_total_prev1, sz_total_prev1, jz_total_prev1, &
      & (jz_total_prev1 - jz_total_prev2) / max(angmom_time_prev1 - angmom_time_prev2, tiny(1d0)), &
      & center_x, center_y, system)
  else
    call dump_lz_record(angmom_step_prev2, angmom_time_prev2, lz_xy_prev2, sz_xy_prev2, jz_xy_prev2, &
      & 0d0 * jz_xy_prev2, lz_total_prev2, sz_total_prev2, jz_total_prev2, 0d0, center_x, center_y, system)
  end if
end subroutine flush_local_angular_momentum_xy

subroutine dump_lz_record(itt, time_now, lz_xy, sz_xy, jz_xy, djzdt_xy, lz_total, sz_total, jz_total, djzdt_total, center_x, center_y, system)
  use structures, only: s_dft_system
  use communication, only: comm_is_root
  use parallelization, only: nproc_id_global
  use salmon_global, only: base_directory, sysname
  use inputoutput, only: t_unit_length, t_unit_time, t_unit_time_inv
  implicit none
  integer, intent(in) :: itt
  real(8), intent(in) :: time_now, lz_xy(:,:), sz_xy(:,:), jz_xy(:,:), djzdt_xy(:,:)
  real(8), intent(in) :: lz_total, sz_total, jz_total, djzdt_total, center_x, center_y
  type(s_dft_system), intent(in) :: system
  character(256) :: file_map, file_total, filenum
  integer :: ix, iy, iunit
  real(8) :: x, y
  real(8) :: conv_ang_xy, conv_djzdt_xy, conv_ang_total, conv_djzdt_total
  character(64) :: unit_ang_xy, unit_djzdt_xy, unit_ang_total, unit_djzdt_total

  if (.not. comm_is_root(nproc_id_global)) return

  conv_ang_xy = t_unit_time_inv%conv
  conv_djzdt_xy = t_unit_time_inv%conv * t_unit_time_inv%conv
  conv_ang_total = t_unit_length%conv * t_unit_length%conv * t_unit_time_inv%conv
  conv_djzdt_total = t_unit_length%conv * t_unit_length%conv * t_unit_time_inv%conv * t_unit_time_inv%conv

  if (trim(t_unit_time_inv%name) == 'a.u.') then
    unit_ang_xy = 'a.u.'
    unit_djzdt_xy = 'a.u.'
  else
    unit_ang_xy = trim(t_unit_time_inv%name)
    unit_djzdt_xy = trim(t_unit_time_inv%name)//'^2'
  end if
  if (trim(t_unit_length%name) == 'a.u.' .or. trim(t_unit_time_inv%name) == 'a.u.') then
    unit_ang_total = 'a.u.'
    unit_djzdt_total = 'a.u.'
  else
    unit_ang_total = trim(t_unit_length%name)//'^2/'//trim(t_unit_time%name)
    unit_djzdt_total = trim(t_unit_length%name)//'^2/'//trim(t_unit_time%name)//'^2'
  end if

  write(filenum, '(i6.6)') itt
  file_map = trim(base_directory)//trim(sysname)//'_lz_xy_'//trim(adjustl(filenum))//'.data'
  open(newunit=iunit, file=trim(file_map), status='replace', action='write')
  write(iunit,'(a)') '# Local electronic angular momentum components integrated along z'
  write(iunit,'(a)') '# x: x coordinate measured from the chosen vortex center'
  write(iunit,'(a)') '# y: y coordinate measured from the chosen vortex center'
  write(iunit,'(a)') '# local_lz_zint: local Lz integrated over z'
  write(iunit,'(a)') '# local_sz_zint: local Sz integrated over z'
  write(iunit,'(a)') '# local_jz_zint: local Jz=Lz+Sz integrated over z'
  write(iunit,'(a)') '# dlocal_jz_dt_zint: time derivative of local_jz_zint'
  write(iunit,'(13a)') '# 1:x[', trim(t_unit_length%name), '] 2:y[', trim(t_unit_length%name), &
                     '] 3:local_lz_zint[', trim(unit_ang_xy), '] 4:local_sz_zint[', trim(unit_ang_xy), &
                     '] 5:local_jz_zint[', trim(unit_ang_xy), '] 6:dlocal_jz_dt_zint[', trim(unit_djzdt_xy), ']'
  do iy = 1, size(lz_xy, 2)
    y = ((dble(iy) - 0.5d0) * system%hgs(2) - center_y) * t_unit_length%conv
    do ix = 1, size(lz_xy, 1)
      x = ((dble(ix) - 0.5d0) * system%hgs(1) - center_x) * t_unit_length%conv
      write(iunit,'(6(1x,es24.16))') x, y, lz_xy(ix,iy) * conv_ang_xy, sz_xy(ix,iy) * conv_ang_xy, &
                                     jz_xy(ix,iy) * conv_ang_xy, djzdt_xy(ix,iy) * conv_djzdt_xy
    end do
    write(iunit,*)
  end do
  close(iunit)

  file_total = trim(base_directory)//trim(sysname)//'_lz_total.data'
  if (itt == 0) then
    open(newunit=iunit, file=trim(file_total), status='replace', action='write')
    write(iunit,'(a)') '# Total electronic angular momentum components'
    write(iunit,'(a)') '# step: RT step index'
    write(iunit,'(a)') '# time: RT time'
    write(iunit,'(a)') '# local_lz_total: total electronic Lz'
    write(iunit,'(a)') '# local_sz_total: total electronic Sz'
    write(iunit,'(a)') '# local_jz_total: total electronic Jz=Lz+Sz'
    write(iunit,'(a)') '# dlocal_jz_total_dt: time derivative of local_jz_total'
    write(iunit,'(11a)') '# 1:step[none] 2:time[', trim(t_unit_time%name), '] 3:local_lz_total[', trim(unit_ang_total), &
                        '] 4:local_sz_total[', trim(unit_ang_total), '] 5:local_jz_total[', trim(unit_ang_total), &
                        '] 6:dlocal_jz_total_dt[', trim(unit_djzdt_total), ']'
  else
    open(newunit=iunit, file=trim(file_total), status='old', position='append', action='write')
  end if
  write(iunit,'(i10,1x,5(es24.16,1x))') itt, time_now * t_unit_time%conv, lz_total * conv_ang_total, &
                                        sz_total * conv_ang_total, jz_total * conv_ang_total, djzdt_total * conv_djzdt_total
  close(iunit)
end subroutine dump_lz_record

end module rt_angular_momentum
