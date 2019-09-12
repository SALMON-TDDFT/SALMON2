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
module filesystem
  use, intrinsic :: iso_c_binding
  implicit none

  public :: file_exists, directory_exists
  public :: create_directory, remove_directory

  public :: atomic_create_directory ! for multi process

private
interface 
  subroutine posix_file_exists(filepath, retcode) bind(C,name='posix_file_exists')
    use, intrinsic :: iso_c_binding
    character(kind=c_char),intent(in) :: filepath
    integer(c_int),intent(out)        :: retcode
  end subroutine

  subroutine posix_directory_exists(dirpath, retcode) bind(C,name='posix_directory_exists')
    use, intrinsic :: iso_c_binding
    character(kind=c_char),intent(in) :: dirpath
    integer(c_int),intent(out)        :: retcode
  end subroutine

  subroutine posix_mkdir(dirpath, retcode) bind(C,name='posix_mkdir')
    use, intrinsic :: iso_c_binding
    character(kind=c_char),intent(in) :: dirpath
    integer(c_int),intent(out)        :: retcode
  end subroutine

  subroutine posix_rmdir(dirpath, retcode) bind(C,name='posix_mkdir')
    use, intrinsic :: iso_c_binding
    character(kind=c_char),intent(in) :: dirpath
    integer(c_int),intent(out)        :: retcode
  end subroutine
end interface

contains
  ! filepath: file path (relative or absolute)
  ! result: .true.  = file exists
  !         .false. = file not exists
  function file_exists(filepath) result(ret)
    implicit none
    character(*), intent(in) :: filepath
    logical :: ret
    integer :: retcode
    call posix_file_exists(adjustl(trim(filepath))//c_null_char, retcode)
    ret = (retcode == 1) ! success
  end function

  ! dirpath: directory path (relative or absolute)
  ! result: .true.  = directory exists
  !         .false. = directory not exists
  function directory_exists(dirpath) result(ret)
    implicit none
    character(*), intent(in) :: dirpath
    logical :: ret
    integer :: retcode
    call posix_directory_exists(adjustl(trim(dirpath))//c_null_char, retcode)
    ret = (retcode == 1) ! success
  end function

  ! dirpath: directory path (relative or absolute)
  subroutine create_directory(dirpath)
    implicit none
    character(*), intent(in) :: dirpath
    integer :: retcode
    call posix_mkdir(adjustl(trim(dirpath))//c_null_char, retcode)
    if (retcode /= 0) then
      stop 'fail: create_directory'
    end if
  end subroutine

  ! dirpath: directory path (relative or absolute)
  subroutine remove_directory(dirpath)
    implicit none
    character(*), intent(in) :: dirpath
    integer :: retcode
    call posix_rmdir(adjustl(trim(dirpath))//c_null_char, retcode)
    if (retcode /= 0) then
      stop 'fail: remove_directory'
    end if
  end subroutine

  ! dirpath:   directory path (relative or absolute)
  ! igroup:    communicator
  ! idelegate: A delegate process of creating directory
  subroutine atomic_create_directory(dirpath, igroup, idelegate)
    use salmon_communication, only: comm_is_root,comm_sync_all
    implicit none
    character(*), intent(in) :: dirpath
    integer, intent(in)      :: igroup, idelegate

    if (comm_is_root(idelegate)) then
      call create_directory(dirpath)
    end if
    call comm_sync_all(igroup) ! sync until directory created
  end subroutine
end module filesystem
