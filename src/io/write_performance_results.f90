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
module write_performance_results
  implicit none

  private :: write_loadbalance
  public  :: write_gs_performance
  public  :: write_rt_performance

  integer,parameter,private :: mode_stdout = 1
  integer,parameter,private :: mode_csv    = 2

contains
  subroutine write_loadbalance(fd,nsize,tsrc,headers,write_mode)
    use salmon_parallel
    use salmon_communication
    use timer
    use math_constants
    implicit none
    integer,intent(in)       :: fd,nsize
    real(8),intent(in)       :: tsrc(nsize)
    character(30),intent(in) :: headers(nsize)
    integer,intent(in)       :: write_mode

    character(*), parameter :: time_format     = '(a30,6(f12.2))'
    character(*), parameter :: time_format_csv = '(a,",",5(f0.6,","),f0.6)'

    real(8) :: tmin(nsize),tmax(nsize),tdif(nsize),trel(nsize)
    real(8) :: tavg(nsize),tnorm(nsize),tstdev(nsize)
    integer :: i

    call comm_get_min(tsrc,tmin,nsize,nproc_group_global)
    call comm_get_max(tsrc,tmax,nsize,nproc_group_global)
    call comm_summation(tsrc,tavg,nsize,nproc_group_global)

    do i=1,nsize
      if (is_zero(tmax(i))) then
        tdif(i) = 0.d0
        trel(i) = 0.d0
        tavg(i) = 0.d0
      else
        tdif(i) = tmax(i) - tmin(i)
        trel(i) = tmax(i) / tmin(i)
        tavg(i) = tavg(i) / nproc_size_global
      end if
    end do

    do i=1,nsize
      tnorm(i) = tsrc(i)**2 - tavg(i)**2
    end do

    call comm_summation(tnorm,tstdev,nsize,nproc_group_global)

    do i=1,nsize
      tstdev(i) = sqrt(tstdev(i) / nproc_size_global)
    end do

    if (comm_is_root(nproc_id_global)) then
      if (write_mode == mode_stdout) then
        write (fd,'(a8,22x,6(a12))') 'Function','min','max','average','std. dev.','difference','relative'
      else
        write (fd,'(a,",",5(a,","),a)') 'Function','min','max','average','std. dev.','difference','relative'
      end if
      do i=1,nsize
        if (is_nonzero(trel(i))) then
          if (write_mode == mode_stdout) then
            write (fd,time_format) headers(i),tmin(i),tmax(i),tavg(i),tstdev(i),tdif(i),trel(i)
          else
            write (fd,time_format_csv) trim(headers(i)),tmin(i),tmax(i),tavg(i),tstdev(i),tdif(i),trel(i)
          end if
        end if
      end do
    end if
  end subroutine

  subroutine write_root(fd, str)
    use salmon_parallel
    use salmon_communication
    implicit none
    integer, intent(in)      :: fd
    character(*), intent(in) :: str
    if (comm_is_root(nproc_id_global)) then
      write (fd,'(a)') str
    end if
  end subroutine

  subroutine write_gs_performance(fd)
    use salmon_parallel
    use salmon_communication
    use timer
    implicit none
    integer, intent(in) :: fd

    integer,parameter :: LOG_SIZE = 10
    real(8)       :: tsrc(LOG_SIZE)
    character(30) :: headers(LOG_SIZE)
    integer :: mode

    mode = mode_csv

    call write_root(fd, '==================== elapsed time [s] ====================')
    if (comm_is_root(nproc_id_global)) then
      call timer_write(fd, 'elapsed time initializing          = ', LOG_INIT_GS)
      call timer_write(fd, 'elapsed time for reading data      = ', LOG_INIT_GS_RESTART)
      call timer_write(fd, 'elapsed time for init. scf iter.   = ', LOG_INIT_GS_ITERATION)
      call timer_write(fd, 'elapsed time for scf iterations    = ', LOG_GS_ITERATION)
      call timer_write(fd, 'elapsed time for deinit. scf iter. = ', LOG_DEINIT_GS_ITERATION)
      call timer_write(fd, 'elapsed time for writing data      = ', LOG_WRITE_RESULTS)
      call timer_write(fd, 'elapsed time for writing LDA data  = ', LOG_WRITE_LDA_DATA)
      call timer_write(fd, 'elapsed time for writing infos     = ', LOG_WRITE_INFOS)
      call timer_write(fd, 'total time                         = ', LOG_TOTAL)
    end if

    call write_root(fd, 'Load balance check [sec]')

    call set(1, LOG_CALC_GRAM_SCHMIDT , 'Gram Schmidt')
    call set(2, LOG_CALC_SUBSPACE_DIAG, 'subspace-diag.')
    call set(3, LOG_CALC_RHO          , 'calculating rho')
    call set(4, LOG_CALC_HARTREE      , 'Hartree routine')
    call set(5, LOG_CALC_EXC_COR      , 'Exc_Cor routine')
    call set(6, LOG_CALC_TOTAL_ENERGY , 'calculating Etot')
    call set(7, LOG_WRITE_RESULTS     , 'writing info.')
    call write_root(fd, '================== in scf iterations [s] =================')
    call write_loadbalance(fd, 7, tsrc, headers, mode)

    call set(1, LOG_DIAG_INIT       , 'initialization')
    call set(2, LOG_DIAG_VLOCAL     , 'Vlocal')
    call set(3, LOG_DIAG_AMAT       , 'Amat')
    call set(4, LOG_DIAG_ALLREDUCE  , 'allreduce')
    call set(5, LOG_DIAG_EIGEN      , 'eigen')
    call set(6, LOG_DIAG_SET_ORBITAL, 'set orbital')
    call set(7, LOG_DIAG_UPDATE     , 'set orbital')
    call write_root(fd, '================== in subspace-diag. [s] =================')
    call write_loadbalance(fd, 7, tsrc, headers, mode)

    call set(1, LOG_GSCG_TOTAL              , 'total')
    call set(2, LOG_GSCG_INIT               , 'init.')
    call set(3, LOG_GSCG_INIT_ITERATION     , 'init. iter.')
    call set(4, LOG_GSCG_ITERATION          , 'iterations')
    call set(5, LOG_GSCG_DEINIT             , 'deinit.')
    call set(6, LOG_ALLREDUCE_INNER_PRODUCT5, 'comm. for inner product(5)')
    call set(7, LOG_ALLREDUCE_INNER_PRODUCT7, 'comm. for inner product(7)')
    call write_root(fd, '================== in CG. [s] ============================')
    call write_loadbalance(fd, 7, tsrc, headers, mode)

  contains
    subroutine set(nid, tid, header)
      implicit none
      integer, intent(in)      :: nid, tid
      character(*), intent(in) :: header
      tsrc(nid) = timer_get(tid)
      write (headers(nid),'(a)') header
    end subroutine
  end subroutine

  subroutine write_rt_performance(fd)
    use salmon_parallel
    use salmon_communication
    use timer
    implicit none
    integer, intent(in) :: fd

    integer,parameter :: LOG_SIZE = 20
    real(8)       :: tsrc(LOG_SIZE)
    character(30) :: headers(LOG_SIZE)
    integer :: mode

    mode = mode_csv

    call write_root(fd, '==================== elapsed time [s] ====================')
    if (comm_is_root(nproc_id_global)) then
      call timer_write(fd, 'elapsed time for rt initialization  = ', LOG_INIT_RT)
      call timer_write(fd, 'elapsed time for reading lda data   = ', LOG_READ_LDA_DATA)
      call timer_write(fd, 'elapsed time for reading rt data    = ', LOG_READ_RT_DATA)
      call timer_write(fd, 'elapsed time for prep. time prop.   = ', LOG_INIT_TIME_PROPAGATION)
      call timer_write(fd, 'elapsed time for rt iterations      = ', LOG_RT_ITERATION)
      call timer_write(fd, 'elapsed time for writing rt data    = ', LOG_WRITE_RT_DATA)
      call timer_write(fd, 'elapsed time aft. writing rt data   = ', LOG_WRITE_RESULTS)
      call timer_write(fd, 'total time                          = ', LOG_TOTAL)
    end if

    call write_root(fd, '======================================================')
    call write_root(fd, '')

    call set( 1, LOG_CALC_VBOX            , 'Vbox')
    call set( 2, LOG_CALC_TIME_PROPAGATION, 'time propagation')
    call set( 3, LOG_CALC_RHO             , 'calculating rho')
    call set( 4, LOG_CALC_HARTREE         , 'Hartree routine')
    call set( 5, LOG_CALC_EXC_COR         , 'Exc_Cor routine')
    call set( 6, LOG_CALC_VLOCAL          , 'Vhxc')              ! FIXME: wrong name    
    call set( 7, LOG_CALC_DP              , 'calculating Dp')
    call set( 8, LOG_CALC_CURRENT         , 'calculating curr')
    call set( 9, LOG_CALC_TOTAL_ENERGY    , 'calculating Etot')
    call set(10, LOG_CALC_PROJECTION      , 'calc. projection')
    call set(11, LOG_CALC_QUADRUPOLE      , 'calc. quadrupole')  ! FIXME: wrong name
    call set(12, LOG_WRITE_ENERGIES       , 'writing energies')
    call set(13, LOG_WRITE_INFOS          , 'writing info etc.')
    call write_root(fd, '=========== elapsed time for rt iterations [s] ===========')
    call write_loadbalance(fd, 13, tsrc, headers, mode)

    call set(1, LOG_ALLREDUCE_RHO    , 'Allreduce in rho')
    call set(2, LOG_ALLREDUCE_HARTREE, 'Allreduce in Hartree')
    call set(3, LOG_ALLREDUCE_DIPOLE , 'Allreduce in dipole calc.')
    call set(4, LOG_ALLGATHERV_TOTAL , 'Allgatherv')
    call write_root(fd, '=========== communication time =======================')
    call write_loadbalance(fd, 4, tsrc, headers, mode)

    call set(1, LOG_TEP_SENDRECV      , 'sendrecv')
    call set(2, LOG_TEP_ORBITAL_ENERGY, 'orbital energy')
    call set(3, LOG_TEP_ION_ION       , 'ion-ion')
    call set(4, LOG_TEP_ION_ELECTRON  , 'ion-electron')
    call set(5, LOG_TEP_NONLOCAL_1    , 'nonlocal 1')
    call set(6, LOG_TEP_NONLOCAL_2    , 'nonlocal 2')
    call write_root(fd, '=========== total_energy_periodic ====================')
    call write_loadbalance(fd, 6, tsrc, headers, mode)

    call set(1, LOG_CUR_SENDRECV           , 'sendrecv')
    call set(2, LOG_CUR_LOCAL              , 'current (except nonlocal)')
    call set(3, LOG_CUR_NONLOCAL1          , 'current nonlocal (1)')
    call set(4, LOG_CUR_NONLOCAL1_ALLREDUCE, 'Allreduce nonlocal (1)')
    call set(5, LOG_CUR_NONLOCAL2          , 'current nonlocal (2)')
    call set(6, LOG_CUR_NONLOCAL2_ALLREDUCE, 'Allreduce nonlocal (2)')
    call write_root(fd, '=========== current ==================================')
    call write_loadbalance(fd, 6, tsrc, headers, mode)

  contains
    subroutine set(nid, tid, header)
      implicit none
      integer, intent(in)      :: nid, tid
      character(*), intent(in) :: header
      tsrc(nid) = timer_get(tid)
      write (headers(nid),'(a)') header
    end subroutine
  end subroutine write_rt_performance
end module write_performance_results
