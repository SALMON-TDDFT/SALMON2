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

  integer,public,parameter :: LOG_READ_GS_DATA         = 10
  integer,public,parameter :: LOG_READ_RT_DATA          = 11
  integer,public,parameter :: LOG_WRITE_GS_DATA        = 12
  integer,public,parameter :: LOG_WRITE_GS_INFO       = 13
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
  integer,public,parameter :: LOG_CALC_ALLGATHERV_VLOCAL = 58
  integer,public,parameter :: LOG_CALC_PROJECTION       = 59
  integer,public,parameter :: LOG_WRITE_ENERGIES        = 61
  integer,public,parameter :: LOG_CALC_CURRENT_ION      = 65
  integer,public,parameter :: LOG_CALC_EIGEN_ENERGY     = 66
  integer,public,parameter :: LOG_CALC_TOTAL_ENERGY_PERIODIC = 67
  integer,public,parameter :: LOG_CALC_SINGLESCALE      = 69
  integer,public,parameter :: LOG_CALC_DENSITY_MATRIX   = 68
  integer,public,parameter :: LOG_WRITE_RT_INFOS        = 62
  integer,public,parameter :: LOG_RT_MISC               = 63
  integer,public,parameter :: LOG_RT_ANALYSIS           = 64

  ! for MD
  integer,public,parameter :: LOG_MD_TEVOL_PART1        = 70
  integer,public,parameter :: LOG_MD_TEVOL_PART2        = 71
  integer,public,parameter :: LOG_MD_UPDATE_PSEUDO_PT   = 72

  integer,public,parameter :: LOG_CHECKPOINT_SELF       = 80
  integer,public,parameter :: LOG_CHECKPOINT_SYNC       = 81
  integer,public,parameter :: LOG_RESTART_SELF          = 82
  integer,public,parameter :: LOG_RESTART_SYNC          = 83

  ! =====================
  ! for specific routines

  ! singlescale maxwell-tddft
  integer,public,parameter :: LOG_SS_FDTD_CALC = 150
  integer,public,parameter :: LOG_SS_FDTD_COMM = 151
  integer,public,parameter :: LOG_SS_FDTD_COMM_COLL = 152
  integer,public,parameter :: LOG_SS_UPDATE_NONLOCALPT_MICROAC = 153

  ! total_energy module
  integer,public,parameter :: LOG_TE_ISOLATED_CALC       = 200
  integer,public,parameter :: LOG_TE_ISOLATED_COMM_COLL  = 201
  integer,public,parameter :: LOG_TE_PERIODIC_CALC       = 202
  integer,public,parameter :: LOG_TE_PERIODIC_COMM_COLL  = 203
  integer,public,parameter :: LOG_EIGEN_ENERGY_CALC      = 204
  integer,public,parameter :: LOG_EIGEN_ENERGY_HPSI      = 205
  integer,public,parameter :: LOG_EIGEN_ENERGY_COMM_COLL = 206

  ! density_matrix module
  integer,public,parameter :: LOG_DENSITY_MATRIX_CALC         = 210
  integer,public,parameter :: LOG_DENSITY_MATRIX_COMM_COLL    = 211
  integer,public,parameter :: LOG_DENSITY_MATRIX_COMM_HALO    = 212
  integer,public,parameter :: LOG_DENSITY_CALC                = 213
  integer,public,parameter :: LOG_DENSITY_COMM_COLL           = 214
  integer,public,parameter :: LOG_CURRENT_CALC                = 215
  integer,public,parameter :: LOG_CURRENT_CALC_UVPSI_RDIVIDED = 216
  integer,public,parameter :: LOG_CURRENT_COMM_COLL           = 217
  integer,public,parameter :: LOG_CURRENT_COMM_HALO           = 218
  integer,public,parameter :: LOG_MCURRENT_CALC               = 219
  integer,public,parameter :: LOG_MCURRENT_COMM_COLL          = 220

  ! conjugate_gradient module
  integer,public,parameter :: LOG_GSCG_ISOLATED_CALC      = 221
  integer,public,parameter :: LOG_GSCG_ISOLATED_HPSI      = 222
  integer,public,parameter :: LOG_GSCG_ISOLATED_COMM_COLL = 223
  integer,public,parameter :: LOG_GSCG_PERIODIC_CALC      = 224
  integer,public,parameter :: LOG_GSCG_PERIODIC_HPSI      = 225
  integer,public,parameter :: LOG_GSCG_PERIODIC_COMM_COLL = 226

  ! subspace_diagonalization module
  integer,public,parameter :: LOG_SSDG_ISOLATED_CALC      = 227
  integer,public,parameter :: LOG_SSDG_ISOLATED_HPSI      = 228
  integer,public,parameter :: LOG_SSDG_ISOLATED_COMM_COLL = 229
  integer,public,parameter :: LOG_SSDG_PERIODIC_CALC      = 230
  integer,public,parameter :: LOG_SSDG_PERIODIC_HPSI      = 231
  integer,public,parameter :: LOG_SSDG_PERIODIC_COMM_COLL = 232
  integer,public,parameter :: LOG_SSDG_PERIODIC_EIGEN     = 239

  ! subspace_diagonalization_so module
  integer,public,parameter :: LOG_SSDG_SO_ISOLATED_CALC      = 233
  integer,public,parameter :: LOG_SSDG_SO_ISOLATED_HPSI      = 234
  integer,public,parameter :: LOG_SSDG_SO_ISOLATED_COMM_COLL = 235
  integer,public,parameter :: LOG_SSDG_SO_PERIODIC_CALC      = 236
  integer,public,parameter :: LOG_SSDG_SO_PERIODIC_HPSI      = 237
  integer,public,parameter :: LOG_SSDG_SO_PERIODIC_COMM_COLL = 238

  ! force module
  integer,public,parameter :: LOG_CALC_FORCE_ION_ION  = 261
  integer,public,parameter :: LOG_CALC_FORCE_FOURIER  = 273
  integer,public,parameter :: LOG_CALC_FORCE_ELEC_ION = 262
  integer,public,parameter :: LOG_CALC_FORCE_GTPSI    = 263
  integer,public,parameter :: LOG_CALC_FORCE_DDEN     = 264
  integer,public,parameter :: LOG_CALC_FORCE_NONLOCAL = 265
  integer,public,parameter :: LOG_CALC_FORCE_LOCAL    = 266

  ! prep_pp module
  integer,public,parameter :: LOG_INIT_PS_TOTAL       = 267
  integer,public,parameter :: LOG_INIT_PS_CALC_NPS    = 268
  integer,public,parameter :: LOG_INIT_PS_CALC_JXYZ   = 269
  integer,public,parameter :: LOG_INIT_PS_LMA_UV      = 270
  integer,public,parameter :: LOG_INIT_PS_CALC_VPSL   = 271
  integer,public,parameter :: LOG_INIT_PS_UVPSI       = 272

  ! FIXME: modify later
  ! hamiltonian module
  integer,public,parameter :: LOG_ALLGATHERV_VLOCAL_CALC      = 240
  integer,public,parameter :: LOG_ALLGATHERV_VLOCAL_COMM_COLL = 241

  integer,public,parameter :: LOG_UHPSI_ALL            = 250
  integer,public,parameter :: LOG_UHPSI_UPDATE_OVERLAP = 251
  integer,public,parameter :: LOG_UHPSI_STENCIL        = 252
  integer,public,parameter :: LOG_UHPSI_SUBTRACTION    = 253
  integer,public,parameter :: LOG_UHPSI_PSEUDO         = 254
  integer,public,parameter :: LOG_UHPSI_PSEUDO_COMM    = 255

  integer,public,parameter :: LOG_UHPSI_OVL_PHASE1     = 256
  integer,public,parameter :: LOG_UHPSI_OVL_PHASE2     = 257
  integer,public,parameter :: LOG_UHPSI_OVL_PHASE2_COMM= 258
  integer,public,parameter :: LOG_UHPSI_OVL_PHASE3     = 259
  integer,public,parameter :: LOG_UHPSI_OVL_PHASE4     = 260


  ! FIXME: remove later
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

  public :: timer_enable_sub, timer_disable_sub
  public :: timer_get_sub


  integer,private,parameter   :: LOG_SIZE = 300
  integer,private             :: enable_range
  real(8),private,allocatable :: log_time(:,:)
  real(8),private,allocatable :: log_temp(:,:)
  logical,private,allocatable :: ticked(:)

  real(8),private,allocatable :: log_time_t(:,:)
  real(8),private,allocatable :: log_temp_t(:,:)
  logical,private,allocatable :: ticked_t(:,:)

  character(*),private,parameter :: SHOW_FORMAT = '(a,f10.2,a,f10.2,a)'

#define CHECK_TICKED(ID)        call now_ticked((ID),__LINE__,__FILE__)
#define CHECK_TICKED_T(ID,TID)  call now_ticked_t((ID),(TID),__LINE__,__FILE__)
#define CHECK_STOPPED(ID)       call now_stopped((ID),__LINE__,__FILE__)
#define CHECK_STOPPED_T(ID,TID) call now_stopped_t((ID),(TID),__LINE__,__FILE__)

  private :: now_ticked, now_ticked_t
  private :: now_stopped, now_stopped_t

private
contains
  subroutine timer_initialize
    use omp_lib, only: omp_get_max_threads
    implicit none
    allocate(log_time(0:LOG_SIZE - 1,2))
    allocate(log_temp(0:LOG_SIZE - 1,2))
    allocate(log_time_t(0:LOG_SIZE - 1, 0:omp_get_max_threads()-1))
    allocate(log_temp_t(0:LOG_SIZE - 1, 0:omp_get_max_threads()-1))
    allocate(ticked(0:LOG_SIZE - 1))
    allocate(ticked_t(0:LOG_SIZE - 1, 0:omp_get_max_threads()-1))
    ticked(:)     = .false.
    ticked_t(:,:) = .false.
    call timer_disable_sub
    call timer_reset
  end subroutine

  subroutine timer_set(e,t)
    implicit none
    integer,intent(in) :: e
    real(8),intent(in) :: t
    log_time(e,1:enable_range) = t
    log_temp(e,1:enable_range) = 0.d0
    ticked(e)   = .false.
  end subroutine

  subroutine timer_reset(e)
    implicit none
    integer,intent(in),optional :: e
    integer :: i,j
    if(present(e)) then
      CHECK_TICKED(e)
      do j=0,size(ticked_t,2)-1
        CHECK_TICKED_T(e,j)
      end do
      log_time  (e,1:enable_range) = 0.d0
      log_temp  (e,1:enable_range) = 0.d0
      log_time_t(e,:) = 0.d0
      log_temp_t(e,:) = 0.d0
      ticked(e)       = .false.
    else
      do i=0,LOG_SIZE-1
        CHECK_TICKED(i)
        do j=0,size(ticked_t,2)-1
          CHECK_TICKED_T(i,j)
        end do
      end do
      log_time  (:,1:enable_range) = 0.d0
      log_temp  (:,1:enable_range) = 0.d0
      log_time_t(:,:) = 0.d0
      log_time_t(:,:) = 0.d0
      ticked(:)       = .false.
      ticked_t(:,:)   = .false.
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
    CHECK_TICKED(id)
    log_temp(id,1:enable_range) = get_wtime()
    ticked(id)   = .true.
  end subroutine

  subroutine timer_end(id)
    implicit none
    integer,intent(in) :: id
    CHECK_STOPPED(id)
    log_time(id,1:enable_range) = log_time(id,1:enable_range) + get_wtime() - log_temp(id,1:enable_range)
    ticked(id)   = .false.
  end subroutine

  subroutine timer_thread_begin(id)
    use omp_lib
    implicit none
    integer,intent(in) :: id
    integer :: tid
    tid = omp_get_thread_num()
    CHECK_TICKED_T(id,tid)
    if (tid == 0) then
      call timer_begin(id)
    end if
    log_temp_t(id,tid) = get_wtime()
    ticked_t(id,tid)   = .true.
  end subroutine

  subroutine timer_thread_end(id)
    use omp_lib
    implicit none
    integer,intent(in) :: id
    integer :: tid
    tid = omp_get_thread_num()
    CHECK_STOPPED_T(id,tid)
    if (tid == 0) then
      call timer_end(id)
    end if
    log_time_t(id,tid) = log_time_t(id,tid) + get_wtime() - log_temp_t(id,tid)
    ticked_t(id,tid)   = .false.
  end subroutine

  subroutine timer_show_hour(str, id)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: id
    real(8) :: time,hour
    time = log_time(id,1)
    hour = time / 3600
    write(*,SHOW_FORMAT) str,time,'sec =',hour,'hour'
  end subroutine

  subroutine timer_show_current_hour(str, id)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: id
    real(8) :: time,hour
    time = get_wtime() - log_temp(id,1) + log_time(id,1)
    hour = time / 3600
    write(*,SHOW_FORMAT) str,time,'sec =',hour,'hour'
  end subroutine

  subroutine timer_show_min(str, id)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: id
    real(8) :: time,mini
    time = log_time(id,1)
    mini = time / 60
    write(*,SHOW_FORMAT) str,time,'sec =',mini,'min'
  end subroutine

  subroutine timer_show_current_min(str, id)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: id
    real(8) :: time,mini
    time = get_wtime() - log_temp(id,1) + log_time(id,1)
    mini = time / 60
    write(*,SHOW_FORMAT) str,time,'sec =',mini,'min'
  end subroutine

  subroutine timer_write(fd,str,id)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: fd,id
    real(8) :: time
    time = log_time(id,1)
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
    timer_get = log_time(id,1)
  end function

  subroutine timer_enable_sub()
    implicit none
    enable_range = 2
  end subroutine

  subroutine timer_disable_sub()
    implicit none
    enable_range = 1
  end subroutine

  function timer_get_sub(id)
    implicit none
    integer,intent(in) :: id
    real(8)            :: timer_get_sub
    timer_get_sub = log_time(id,2)
  end function

  function timer_thread_get(id,tid)
    implicit none
    integer,intent(in) :: id,tid
    real(8)            :: timer_thread_get
    timer_thread_get = log_time_t(id,tid)
  end function

  subroutine now_ticked(id,ln,sn)
    implicit none
    integer, intent(in)      :: id,ln
    character(*), intent(in) :: sn
    if (ticked(id)) then
      print '(A,I3,A,"(",A," line.",I4,")")', 'timer ',id,' now ticked... please check it.',sn,ln
#ifdef __INTEL_COMPILER
      call tracebackqq
#endif
      stop 'error'
    end if
  end subroutine

  subroutine now_ticked_t(id,tid,ln,sn)
    implicit none
    integer, intent(in)      :: id,tid,ln
    character(*), intent(in) :: sn
    if (ticked_t(id,tid)) then
      print '(A,I3,A,I3,A,"(",A," line.",I4,")")', 'thread ',tid,': timer ',id,' now ticked... please check it.',sn,ln
#ifdef __INTEL_COMPILER
!$omp critical
      call tracebackqq
!$omp end critical
#endif
      stop 'error'
    end if
  end subroutine

  subroutine now_stopped(id,ln,sn)
    implicit none
    integer, intent(in)      :: id,ln
    character(*), intent(in) :: sn
    if (.not.ticked(id)) then
      print '(A,I3,A,"(",A," line.",I4,")")', 'timer ',id,' now stopped... please check it.',sn,ln
#ifdef __INTEL_COMPILER
      call tracebackqq
#endif
      stop 'error'
    end if
  end subroutine

  subroutine now_stopped_t(id,tid,ln,sn)
    implicit none
    integer, intent(in)      :: id,tid,ln
    character(*), intent(in) :: sn
    if (.not.ticked_t(id,tid)) then
      print '(A,I3,A,I3,A,"(",A," line.",I4,")")', 'thread ',tid,': timer ',id,' now stopped... please check it.',sn,ln
#ifdef __INTEL_COMPILER
!$omp critical
      call tracebackqq
!$omp end critical
#endif
      stop 'error'
    end if
  end subroutine
end module
