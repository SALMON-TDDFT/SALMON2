!
!  Copyright 2020 SALMON developers
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
!=======================================================================
module read_rtdata_file
  use filesystem, only: open_filehandle
  use inputoutput, only: t_unit_ac,t_unit_time
  implicit none

contains

!===================================================================================================================================

  function count_rows_from_rtdata_file(filename) result(icount)
    ! This function returns the number of data rows 
    ! in the `?????_rt.data` formatted file.
    ! Blank lines and comment lines are skipped.
    implicit none
    character(256), intent(in) :: filename
    character(1024) :: buf
    integer :: icount
    integer ::  flag, fh

    fh = open_filehandle(trim(filename), status='old')
    icount = 0
    do while (.true.)
        read(fh, '(a)', iostat=flag)  buf
        if (flag == 0) then
            buf = adjustl(buf)
            if (len_trim(buf) < 1) cycle
            if (buf(1:1) == '#' .or. buf(1:1) == '!') cycle
            icount = icount + 1
        else
            exit
        end if
    end do
    close(fh)
    return
  end function count_rows_from_rtdata_file

!===================================================================================================================================

  subroutine load_Ac_from_rtdata_file(filename, n_dat, t_dat, Ac_dat)
    ! This subroutine Reads the data of the given number of samples 
    ! from `???_rt.data` file. The empty and comment lines are skipped,
    ! and the values of the first four columns:
    ! time(1) and Ac_ext (2-4) are stored in t_dat, Ac_dat.
    implicit none
    character(256), intent(in) :: filename
    integer, intent(in) :: n_dat
    real(8), intent(out) :: t_dat(n_dat)
    real(8), intent(out) :: Ac_dat(1:3, n_dat)

    integer :: i, fh
    character(1024) :: buf

    fh = open_filehandle(trim(filename), status='old')
    i = 1
    do while (i <= n_dat)
        read(fh, '(a)')  buf
        buf = adjustl(buf)
        if (len_trim(buf) < 1) cycle
        if (buf(1:1) == '#' .or. buf(1:1) == '!') cycle
        read(buf, *) t_dat(i), Ac_dat(1, i), Ac_dat(2, i), Ac_dat(3, i)
        i = i + 1
    end do
    close(fh)

    ! Transform unit of Ac field
    Ac_dat(:, :) = Ac_dat(:, :) / t_unit_ac%conv
    t_dat(:) = t_dat(:) / t_unit_time%conv
    
    return
  end subroutine load_Ac_from_rtdata_file

  ! if you want use this subroutine, please switch subroutine name in
  ! rt/em_filed.f90
  subroutine load_Ac_from_rtdata_file_Ac_tot(filename, n_dat, t_dat, Ac_dat)
    ! This subroutine Reads the data of the given number of samples 
    ! from `???_rt.data` file. The empty and comment lines are skipped,
    ! and the values of the four columns:
    ! time(1) and Ac_tot (8-10) are stored in t_dat, Ac_dat.
    implicit none
    character(256), intent(in) :: filename
    integer, intent(in) :: n_dat
    real(8), intent(out) :: t_dat(n_dat)
    real(8), intent(out) :: Ac_dat(1:3, n_dat)
    real(8) :: tmp(10)

    integer :: i, fh
    character(1024) :: buf

    fh = open_filehandle(trim(filename), status='old')
    i = 1
    do while (i <= n_dat)
        read(fh, '(a)')  buf
        buf = adjustl(buf)
        if (len_trim(buf) < 1) cycle
        if (buf(1:1) == '#' .or. buf(1:1) == '!') cycle
        read(buf, *) t_dat(i), tmp(1:6), Ac_dat(1:3, i)  !Ac_tot
        i = i + 1
    end do
    close(fh)

    ! Transform unit of Ac field
    Ac_dat(:, :) = Ac_dat(:, :) / t_unit_ac%conv
    t_dat(:) = t_dat(:) / t_unit_time%conv
    
    return
  end subroutine load_Ac_from_rtdata_file_Ac_tot

!===================================================================================================================================

end module read_rtdata_file
