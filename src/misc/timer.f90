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

  ! Calculation
  integer,public,parameter :: LOG_ALL          = 0
  integer,public,parameter :: LOG_STATIC       = 1
  integer,public,parameter :: LOG_GROUND_STATE = 2
  integer,public,parameter :: LOG_DYNAMICS     = 3
  integer,public,parameter :: LOG_IO           = 7

  ! GS routines
  integer,public,parameter :: LOG_CG           = 10
  integer,public,parameter :: LOG_DIAG         = 11
  integer,public,parameter :: LOG_SP_ENERGY    = 12
  integer,public,parameter :: LOG_GRAM_SCHMIDT = 13

  ! GS and RT routines
  integer,public,parameter :: LOG_DT_EVOLVE    = 20
  integer,public,parameter :: LOG_HPSI         = 21
  integer,public,parameter :: LOG_PSI_RHO      = 22
  integer,public,parameter :: LOG_HARTREE      = 23
  integer,public,parameter :: LOG_EXC_COR      = 24
  integer,public,parameter :: LOG_CURRENT      = 25
  integer,public,parameter :: LOG_TOTAL_ENERGY = 26
  integer,public,parameter :: LOG_ION_FORCE    = 27
  integer,public,parameter :: LOG_DT_EVOLVE_AC = 28
  integer,public,parameter :: LOG_ANA_RT_USEGS = 29
!  integer,public,parameter :: LOG_K_SHIFT_WF   = 29  !old name
  integer,public,parameter :: LOG_OTHER        = 30

  ! Hamiltonian
  integer,public,parameter :: LOG_HPSI_INIT    = 35
  integer,public,parameter :: LOG_HPSI_STENCIL = 36
  integer,public,parameter :: LOG_HPSI_PSEUDO  = 37
  integer,public,parameter :: LOG_HPSI_UPDATE  = 38

  ! Communication
  integer,public,parameter :: LOG_ALLREDUCE    = 40
  integer,public,parameter :: LOG_SENDRECV_GRID= 41


  ! for unified version
  ! ===============================================================
  integer,public,parameter :: LOG_TOTAL                 = 100

  integer,public,parameter :: LOG_READ_LDA_DATA         = 110
  integer,public,parameter :: LOG_READ_RT_DATA          = 111
  integer,public,parameter :: LOG_WRITE_LDA_DATA        = 112
  integer,public,parameter :: LOG_WRITE_LDA_INFOS       = 113
  integer,public,parameter :: LOG_WRITE_RT_DATA         = 114
  integer,public,parameter :: LOG_WRITE_GS_RESULTS      = 115
  integer,public,parameter :: LOG_WRITE_RT_RESULTS      = 116

  integer,public,parameter :: LOG_INIT_GS               = 120
  integer,public,parameter :: LOG_INIT_GS_RESTART       = 121
  integer,public,parameter :: LOG_INIT_GS_ITERATION     = 122
  integer,public,parameter :: LOG_DEINIT_GS_ITERATION   = 123

  integer,public,parameter :: LOG_INIT_RT               = 125
  integer,public,parameter :: LOG_INIT_TIME_PROPAGATION = 126

  integer,public,parameter :: LOG_GS_ITERATION          = 130
  integer,public,parameter :: LOG_RT_ITERATION          = 131

  integer,public,parameter :: LOG_CALC_RHO              = 132
  integer,public,parameter :: LOG_CALC_HARTREE          = 133
  integer,public,parameter :: LOG_CALC_EXC_COR          = 134
  integer,public,parameter :: LOG_CALC_TOTAL_ENERGY     = 135

  ! for GS
  integer,public,parameter :: LOG_CALC_GRAM_SCHMIDT     = 140
  integer,public,parameter :: LOG_CALC_SUBSPACE_DIAG    = 141
  integer,public,parameter :: LOG_CALC_MINIMIZATION     = 142
  integer,public,parameter :: LOG_CALC_CHANGE_ORDER     = 143

  ! for RT
  integer,public,parameter :: LOG_CALC_VBOX             = 150
  integer,public,parameter :: LOG_CALC_TIME_PROPAGATION = 151
  integer,public,parameter :: LOG_CALC_DP               = 155
  integer,public,parameter :: LOG_CALC_CURRENT          = 156
  integer,public,parameter :: LOG_CALC_VLOCAL           = 158 ! FIXME: wrong name
  integer,public,parameter :: LOG_CALC_PROJECTION       = 159
  integer,public,parameter :: LOG_CALC_QUADRUPOLE       = 160 ! FIXME: wrong name
  integer,public,parameter :: LOG_WRITE_ENERGIES        = 161
  integer,public,parameter :: LOG_WRITE_RT_INFOS        = 162

!  integer,public,parameter :: LOG_SENDRECV_TOTAL            = 200
  integer,public,parameter :: LOG_SENDRECV_TIME_PROPAGATION = 201

!  integer,public,parameter :: LOG_ALLREDUCE_TOTAL       = 300
  integer,public,parameter :: LOG_ALLREDUCE_RHO         = 301
  integer,public,parameter :: LOG_ALLREDUCE_HARTREE     = 302
  integer,public,parameter :: LOG_ALLREDUCE_DIPOLE      = 303
  integer,public,parameter :: LOG_ALLREDUCE_TOTAL_ENERGY= 304
  integer,public,parameter :: LOG_ALLREDUCE_CURRENT     = 305
  integer,public,parameter :: LOG_ALLREDUCE_INNER_PRODUCT3 = 306
  integer,public,parameter :: LOG_ALLREDUCE_INNER_PRODUCT5 = 307
  integer,public,parameter :: LOG_ALLREDUCE_INNER_PRODUCT7 = 308

  integer,public,parameter :: LOG_ALLGATHERV_TOTAL      = 400

  ! for specific routines
  ! total_energy_periodic (GCEED part)
  integer,public,parameter :: LOG_TEP_TOTAL             = 1000
  integer,public,parameter :: LOG_TEP_SENDRECV          = 1001
  integer,public,parameter :: LOG_TEP_ORBITAL_ENERGY    = 1002
  integer,public,parameter :: LOG_TEP_ION_ION           = 1003
  integer,public,parameter :: LOG_TEP_ION_ELECTRON      = 1004
  integer,public,parameter :: LOG_TEP_NONLOCAL_1        = 1005
  integer,public,parameter :: LOG_TEP_NONLOCAL_2        = 1006

  ! current (GCEED part)
  integer,public,parameter :: LOG_CUR_TOTAL               = 1010
  integer,public,parameter :: LOG_CUR_SENDRECV            = 1011
  integer,public,parameter :: LOG_CUR_LOCAL               = 1012
  integer,public,parameter :: LOG_CUR_NONLOCAL1           = 1013
  integer,public,parameter :: LOG_CUR_NONLOCAL1_ALLREDUCE = 1014
  integer,public,parameter :: LOG_CUR_NONLOCAL2           = 1015
  integer,public,parameter :: LOG_CUR_NONLOCAL2_ALLREDUCE = 1016

  ! subspace diag
  integer,public,parameter :: LOG_DIAG_TOTAL            = 1030
  integer,public,parameter :: LOG_DIAG_INIT             = 1031
  integer,public,parameter :: LOG_DIAG_VLOCAL           = 1032
  integer,public,parameter :: LOG_DIAG_AMAT             = 1033
  integer,public,parameter :: LOG_DIAG_ALLREDUCE        = 1034
  integer,public,parameter :: LOG_DIAG_EIGEN            = 1035
  integer,public,parameter :: LOG_DIAG_SET_ORBITAL      = 1036
  integer,public,parameter :: LOG_DIAG_UPDATE           = 1037

  ! conjugate gradient (gscg)
  integer,public,parameter :: LOG_GSCG_TOTAL            = 1020
  integer,public,parameter :: LOG_GSCG_INIT             = 1021
  integer,public,parameter :: LOG_GSCG_INIT_ITERATION   = 1022
  integer,public,parameter :: LOG_GSCG_ITERATION        = 1023
  integer,public,parameter :: LOG_GSCG_DEINIT           = 1024
  integer,public,parameter :: LOG_GSCG_ALLREDUCE        = 1025
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


  integer,private,parameter   :: LOG_SIZE = 2000
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
