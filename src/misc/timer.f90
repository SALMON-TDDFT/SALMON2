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
module timer
  use misc_routines, only: get_wtime
  implicit none

  integer,public,parameter :: LOG_TOTAL                 = 0
  integer,public,parameter :: LOG_INIT                  = 1 ! general init.

  integer,public,parameter :: LOG_READ_LDA_DATA         = 10
  integer,public,parameter :: LOG_READ_RT_DATA          = 11
  integer,public,parameter :: LOG_WRITE_LDA_DATA        = 12
  integer,public,parameter :: LOG_WRITE_LDA_INFOS       = 13
  integer,public,parameter :: LOG_WRITE_RT_DATA         = 14
  integer,public,parameter :: LOG_WRITE_GS_RESULTS      = 15
  integer,public,parameter :: LOG_WRITE_RT_RESULTS      = 16

  integer,public,parameter :: LOG_INIT_GS               = 20
  integer,public,parameter :: LOG_INIT_GS_RESTART       = 21
  integer,public,parameter :: LOG_INIT_GS_ITERATION     = 22
  integer,public,parameter :: LOG_DEINIT_GS_ITERATION   = 23

  integer,public,parameter :: LOG_INIT_RT               = 25
  integer,public,parameter :: LOG_INIT_TIME_PROPAGATION = 26

  integer,public,parameter :: LOG_GS_ITERATION          = 30
  integer,public,parameter :: LOG_RT_ITERATION          = 31

  integer,public,parameter :: LOG_CALC_RHO              = 32
  integer,public,parameter :: LOG_CALC_HARTREE          = 33
  integer,public,parameter :: LOG_CALC_EXC_COR          = 34
  integer,public,parameter :: LOG_CALC_TOTAL_ENERGY     = 35
  integer,public,parameter :: LOG_CALC_ION_FORCE        = 36

  ! for GS
  integer,public,parameter :: LOG_CALC_GRAM_SCHMIDT     = 40
  integer,public,parameter :: LOG_CALC_SUBSPACE_DIAG    = 41
  integer,public,parameter :: LOG_CALC_MINIMIZATION     = 42
  integer,public,parameter :: LOG_CALC_CHANGE_ORDER     = 43
  integer,public,parameter :: LOG_CALC_ESP              = 44

  ! for RT
  integer,public,parameter :: LOG_CALC_EMFIELD          = 49
  integer,public,parameter :: LOG_CALC_VBOX             = 50
  integer,public,parameter :: LOG_CALC_TIME_PROPAGATION = 51
  integer,public,parameter :: LOG_HPSI                  = 52
  integer,public,parameter :: LOG_CALC_DP               = 55
  integer,public,parameter :: LOG_CALC_CURRENT          = 56
  integer,public,parameter :: LOG_CALC_VLOCAL           = 58 ! FIXME: wrong name
  integer,public,parameter :: LOG_CALC_PROJECTION       = 59
  integer,public,parameter :: LOG_CALC_QUADRUPOLE       = 60 ! FIXME: wrong name
  integer,public,parameter :: LOG_WRITE_ENERGIES        = 61
  integer,public,parameter :: LOG_WRITE_RT_INFOS        = 62
  integer,public,parameter :: LOG_RT_MISC               = 63
  integer,public,parameter :: LOG_RT_ANALYSIS           = 64

  integer,public,parameter :: LOG_ALLREDUCE_RHO             = 101
  integer,public,parameter :: LOG_ALLREDUCE_HARTREE         = 102
  integer,public,parameter :: LOG_ALLREDUCE_DIPOLE          = 103
  integer,public,parameter :: LOG_ALLREDUCE_TOTAL_ENERGY    = 104
  integer,public,parameter :: LOG_ALLREDUCE_CURRENT         = 105
  integer,public,parameter :: LOG_ALLREDUCE_INNER_PRODUCT3  = 106
  integer,public,parameter :: LOG_ALLREDUCE_INNER_PRODUCT5  = 107
  integer,public,parameter :: LOG_ALLREDUCE_INNER_PRODUCT7  = 108
  integer,public,parameter :: LOG_ALLREDUCE_ESP             = 109
  integer,public,parameter :: LOG_ALLREDUCE_ION_FORCE       = 110

  integer,public,parameter :: LOG_ALLGATHERV_TOTAL          = 120

  integer,public,parameter :: LOG_SENDRECV_TIME_PROPAGATION = 130
  integer,public,parameter :: LOG_SENDRECV_GRID             = 131

  ! for specific routines
  ! total_energy_periodic (GCEED part)
  integer,public,parameter :: LOG_TEP_TOTAL          = 200
  integer,public,parameter :: LOG_TEP_SENDRECV       = 201
  integer,public,parameter :: LOG_TEP_ORBITAL_ENERGY = 202
  integer,public,parameter :: LOG_TEP_ION_ION        = 203
  integer,public,parameter :: LOG_TEP_ION_ELECTRON   = 204
  integer,public,parameter :: LOG_TEP_NONLOCAL_1     = 205
  integer,public,parameter :: LOG_TEP_NONLOCAL_2     = 206

  ! current (GCEED part)
  integer,public,parameter :: LOG_CUR_TOTAL               = 210
  integer,public,parameter :: LOG_CUR_SENDRECV            = 211
  integer,public,parameter :: LOG_CUR_LOCAL               = 212
  integer,public,parameter :: LOG_CUR_NONLOCAL1           = 213
  integer,public,parameter :: LOG_CUR_NONLOCAL1_ALLREDUCE = 214
  integer,public,parameter :: LOG_CUR_NONLOCAL2           = 215
  integer,public,parameter :: LOG_CUR_NONLOCAL2_ALLREDUCE = 216

  ! conjugate gradient (gscg)
  integer,public,parameter :: LOG_GSCG_TOTAL          = 220
  integer,public,parameter :: LOG_GSCG_INIT           = 221
  integer,public,parameter :: LOG_GSCG_INIT_ITERATION = 222
  integer,public,parameter :: LOG_GSCG_ITERATION      = 223
  integer,public,parameter :: LOG_GSCG_DEINIT         = 224
  integer,public,parameter :: LOG_GSCG_ALLREDUCE      = 225

  ! subspace diag
  integer,public,parameter :: LOG_DIAG_TOTAL       = 230
  integer,public,parameter :: LOG_DIAG_INIT        = 231
  integer,public,parameter :: LOG_DIAG_VLOCAL      = 232
  integer,public,parameter :: LOG_DIAG_AMAT        = 233
  integer,public,parameter :: LOG_DIAG_ALLREDUCE   = 234
  integer,public,parameter :: LOG_DIAG_EIGEN       = 235
  integer,public,parameter :: LOG_DIAG_SET_ORBITAL = 236
  integer,public,parameter :: LOG_DIAG_UPDATE      = 237

  ! hpsi (ARTED part)
  integer,public,parameter :: LOG_HPSI_INIT    = 245
  integer,public,parameter :: LOG_HPSI_STENCIL = 246
  integer,public,parameter :: LOG_HPSI_PSEUDO  = 247
  integer,public,parameter :: LOG_HPSI_UPDATE  = 248

  ! ===============================================================


  public :: timer_initialize
  public :: timer_set, timer_reset
  public :: timer_get, timer_thread_get
  public :: timer_begin, timer_end
  public :: timer_thread_begin, timer_thread_end
  public :: timer_show_hour, timer_show_min
  public :: timer_show_current_hour, timer_show_current_min

  public :: timer_reentrance_read, timer_reentrance_write
  public :: timer_write, timer_thread_write


  integer,private,parameter   :: LOG_SIZE = 300
  real(8),private,allocatable :: log_time(:)
  real(8),private,allocatable :: log_temp(:)

  real(8),private,allocatable :: log_time_t(:,:)
  real(8),private,allocatable :: log_temp_t(:,:)

  character(*),private,parameter :: SHOW_FORMAT = '(a,f10.2,a,f10.2,a)'

private
contains
  subroutine timer_initialize
    use omp_lib, only: omp_get_max_threads
    implicit none
    allocate(log_time(0:LOG_SIZE - 1))
    allocate(log_temp(0:LOG_SIZE - 1))
    allocate(log_time_t(0:LOG_SIZE - 1, 0:omp_get_max_threads()-1))
    allocate(log_temp_t(0:LOG_SIZE - 1, 0:omp_get_max_threads()-1))
    call timer_reset
  end subroutine

  subroutine timer_set(e,t)
    implicit none
    integer,intent(in) :: e
    real(8),intent(in) :: t
    log_time(e) = t
    log_temp(e) = 0.d0
  end subroutine

  subroutine timer_reset(e)
    implicit none
    integer,intent(in),optional :: e
    if(present(e)) then
      log_time  (e)   = 0.d0
      log_temp  (e)   = 0.d0
      log_time_t(e,:) = 0.d0
      log_temp_t(e,:) = 0.d0
    else
      log_time  (:)   = 0.d0
      log_temp  (:)   = 0.d0
      log_time_t(:,:) = 0.d0
      log_time_t(:,:) = 0.d0
    end if
  end subroutine

  subroutine timer_reentrance_read(fd)
    use backup_routines, only: load_value
    implicit none
    integer,intent(in) :: fd
    deallocate(log_time);   call load_value(fd, log_time)
    deallocate(log_time_t); call load_value(fd, log_time_t)
  end subroutine

  subroutine timer_reentrance_write(fd)
    use backup_routines, only: save_value
    implicit none
    integer,intent(in) :: fd
    call save_value(fd, log_time)
    call save_value(fd, log_time_t)
  end subroutine

  subroutine timer_begin(id)
    implicit none
    integer,intent(in) :: id
    log_temp(id) = get_wtime()
  end subroutine

  subroutine timer_end(id)
    implicit none
    integer,intent(in) :: id
    log_time(id) = log_time(id) + get_wtime() - log_temp(id)
  end subroutine

  subroutine timer_thread_begin(id)
    use omp_lib
    implicit none
    integer,intent(in) :: id
    integer :: tid
    tid = omp_get_thread_num()
    if (tid == 0) then
      call timer_begin(id)
    end if
    log_temp_t(id,tid) = get_wtime()
  end subroutine

  subroutine timer_thread_end(id)
    use omp_lib
    implicit none
    integer,intent(in) :: id
    integer :: tid
    tid = omp_get_thread_num()
    if (tid == 0) then
      call timer_end(id)
    end if
    log_time_t(id,tid) = log_time_t(id,tid) + get_wtime() - log_temp_t(id,tid)
  end subroutine

  subroutine timer_show_hour(str, id)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: id
    real(8) :: time,hour
    time = log_time(id)
    hour = time / 3600
    write(*,SHOW_FORMAT) str,time,'sec =',hour,'hour'
  end subroutine

  subroutine timer_show_current_hour(str, id)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: id
    real(8) :: time,hour
    time = get_wtime() - log_temp(id) + log_time(id)
    hour = time / 3600
    write(*,SHOW_FORMAT) str,time,'sec =',hour,'hour'
  end subroutine

  subroutine timer_show_min(str, id)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: id
    real(8) :: time,mini
    time = log_time(id)
    mini = time / 60
    write(*,SHOW_FORMAT) str,time,'sec =',mini,'min'
  end subroutine

  subroutine timer_show_current_min(str, id)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: id
    real(8) :: time,mini
    time = get_wtime() - log_temp(id) + log_time(id)
    mini = time / 60
    write(*,SHOW_FORMAT) str,time,'sec =',mini,'min'
  end subroutine

  subroutine timer_write(fd,str,id)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: fd,id
    real(8) :: time
    time = log_time(id)
    write(fd,'(a,f16.8,a)') str,time,' [s]'
  end subroutine

  subroutine timer_thread_write(fd,str,id)
    use omp_lib
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: fd,id
    real(8) :: time
    integer :: i
    write(fd,*) str
    do i=0,omp_get_max_threads()-1
      time = log_time_t(id,i)
      write(fd,'(a,i4,a,f16.8,a)') 'tid =',i,': ',time,' [s]'
    end do
  end subroutine

  function timer_get(id)
    implicit none
    integer,intent(in) :: id
    real(8)            :: timer_get
    timer_get = log_time(id)
  end function

  function timer_thread_get(id,tid)
    implicit none
    integer,intent(in) :: id,tid
    real(8)            :: timer_thread_get
    timer_thread_get = log_time_t(id,tid)
  end function
end module
