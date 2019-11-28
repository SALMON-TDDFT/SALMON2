!
!  Copyright 2019 SALMON developers
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
module misc_routines
  implicit none

  public :: string_lowercase
  public :: floor_pow2, ceiling_pow2
  public :: gen_logfilename
  public :: get_wtime

private
contains
  subroutine string_lowercase(str)
    ! Original code was written by David Frank  dave_frank@hotmail.com
    ! http://home.earthlink.net/~dave_gemini/strings.f90
    implicit none
    character(*),intent(inout) :: str

    integer,parameter :: dc = ichar('A') - ichar('a')
    character :: ch
    integer   :: i

    do i = 1,len(str)
      ch = str(i:i)
      if ('A' <= ch .and. ch <= 'Z') then
        ch = char(ichar(ch) - dc)
      end if
      str(i:i) = ch
    end do
  end subroutine

  function floor_pow2(n)
    implicit none
    integer            :: floor_pow2
    integer,intent(in) :: n
    integer :: k

    k = 1
    do while(k < n)
      k = k * 2
    end do
    if (k /= n) k = ishft(k, -1)
    floor_pow2 = k
  end function

  function ceiling_pow2(n)
    implicit none
    integer            :: ceiling_pow2
    integer,intent(in) :: n
    integer :: k

    k = 1
    do while(k < n)
      k = k * 2
    end do
    ceiling_pow2 = k
  end function ceiling_pow2

  ! input  : filename, extension
  ! output : <filename>_YYYYMMDD_hhmmss.<extension>
  function gen_logfilename(filename,extension)
    implicit none
    character(100)          :: gen_logfilename
    character(*),intent(in) :: filename
    character(*),intent(in),optional :: extension
    character(8)  :: d
    character(10) :: t
    call date_and_time(date=d,time=t)

    if (present(extension)) then
      write (gen_logfilename,'(A)') trim(filename)//'_'//d//'_'//t(1:6)//'.'//trim(extension)
    else
      write (gen_logfilename,'(A)') trim(filename)//'_'//d//'_'//t(1:6)//'.log'
    end if
  end function

  function get_wtime()
    implicit none
    real(8) :: get_wtime
    real(8) :: omp_get_wtime
    get_wtime = omp_get_wtime()
  end function
end module
