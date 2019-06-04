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
module perflog
  implicit none

  public  :: write_performance

  integer,parameter,public :: write_mode_readable = 1
  integer,parameter,public :: write_mode_csv      = 2

private
contains
  subroutine write_loadbalance(fd,nsize,tsrc,headers,write_mode)
    use salmon_parallel
    use salmon_communication
    use timer
    use math_constants
    implicit none
    integer,intent(in)       :: fd,nsize
    real(8),intent(in)       :: tsrc(nsize)
    character(30),intent(in) :: headers(0:nsize)
    integer,intent(in)       :: write_mode

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
      if (write_mode == write_mode_readable) then
        write (fd,'(a30,6(a12))') headers(0),'min','max','average','std. dev.','difference','relative'
      else
        write (fd,'(a,",",5(a,","),a)') trim(headers(0)),'min','max','average','std. dev.','difference','relative'
      end if
      do i=1,nsize
        if (.not. is_zero(trel(i))) then
          if (write_mode == write_mode_readable) then
            write (fd,'(a30,6(f12.2))') headers(i),tmin(i),tmax(i),tavg(i),tstdev(i),tdif(i),trel(i)
          else
            write (fd,'(a,",",5(f0.6,","),f0.6)') trim(headers(i)),tmin(i),tmax(i),tavg(i),tstdev(i),tdif(i),trel(i)
          end if
        end if
      end do
    end if
  end subroutine

  subroutine write_flops(fd,write_mode)
    use salmon_parallel
    use salmon_communication
    use flops
    use math_constants
    use salmon_global, only: theory, iperiodic
    implicit none
    integer, intent(in) :: fd, write_mode

    real(8) :: lg(4),pg(4),mg(4),sg(4)

    select case(theory)
      case('TDDFT')
        select case(iperiodic)
          ! TODO: calc flops on unified hpsi
          ! 0: GCEED
          ! 3: ARTED
          case(3)
            call get_hamiltonian_flops(lg,pg,mg,sg)
            if (comm_is_root(nproc_id_global)) then
              if (write_mode == write_mode_readable) then
                write (fd,'(a30,4(a12))') 'hamiltonian','all','stencil','pseudo-pt','update'
                write (fd,'(a30,4(f12.2))') 'processor'       ,lg(4),lg(1),lg(2),lg(3)
                write (fd,'(a30,4(f12.2))') 'processor (best)',pg(4),pg(1),pg(2),pg(3)
                if (.not. is_zero(mg(4))) write (fd,'(a30,4(f12.2))') 'macro-grid',mg(4),mg(1),mg(2),mg(3)
                write (fd,'(a30,4(f12.2))') 'system'          ,sg(4),sg(1),sg(2),sg(3)
              else
                write (fd,'(a,",",3(a,","),a)') 'hamiltonian','all','stencil','pseudo-pt','update'
                write (fd,'(a,",",3(f0.6,","),f0.6)') 'processor'       ,lg(4),lg(1),lg(2),lg(3)
                write (fd,'(a,",",3(f0.6,","),f0.6)') 'processor (best)',pg(4),pg(1),pg(2),pg(3)
                if (.not. is_zero(mg(4))) write (fd,'(a,",",3(f0.6,","),f0.6)') 'macro-grid'      ,mg(4),mg(1),mg(2),mg(3)
                write (fd,'(a,",",3(f0.6,","),f0.6)') 'system'          ,sg(4),sg(1),sg(2),sg(3)
              end if
            end if
        end select
    end select
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

  subroutine write_performance(fd,mode)
    use salmon_parallel
    use salmon_communication
    use timer
    implicit none
    integer, intent(in) :: fd, mode

    integer,parameter :: LOG_SIZE = 20
    real(8)       :: tsrc(LOG_SIZE)
    character(30) :: headers(0:LOG_SIZE)

    call write_root(fd, '=== application breakdown [s] ===')
    if (comm_is_root(nproc_id_global)) then
      write (fd,'(a,f0.8)') 'total calculation time,',timer_get(LOG_TOTAL)
    end if

    call set(0, 0, 'scf calculation')
    call set(1, LOG_INIT_GS            , 'gs initialization')
    call set(2, LOG_INIT_GS_RESTART    , 'reading data')
    call set(3, LOG_INIT_GS_ITERATION  , 'init. scf iter.')
    call set(4, LOG_GS_ITERATION       , 'scf iterations')
    call set(5, LOG_DEINIT_GS_ITERATION, 'deinit. scf iter.')
    call set(6, LOG_WRITE_GS_RESULTS   , 'writing scf data')
    call set(7, LOG_WRITE_LDA_DATA     , 'writing LDA data')
    call set(8, LOG_WRITE_LDA_INFOS    , 'writing LDA infos')
    call write_loadbalance(fd, 8, tsrc, headers, mode)

    call set(0, 0, 'rt calculation')
    call set(1, LOG_INIT_RT              , 'rt initialization')
    call set(2, LOG_READ_LDA_DATA        , 'reading lda data')
    call set(3, LOG_READ_RT_DATA         , 'reading rt data')
    call set(4, LOG_INIT_TIME_PROPAGATION, 'prep. time prop.')
    call set(5, LOG_RT_ITERATION         , 'rt iterations')
    call set(6, LOG_WRITE_RT_DATA        , 'writing rt data')
    call set(7, LOG_WRITE_RT_RESULTS     , 'after writing rt data')
    call write_loadbalance(fd, 7, tsrc, headers, mode)

    call set(0, 0, 'in scf iterations')
    call set(1, LOG_CALC_GRAM_SCHMIDT , 'Gram Schmidt')
    call set(2, LOG_CALC_SUBSPACE_DIAG, 'subspace-diag.')
    call set(3, LOG_CALC_RHO          , 'calculating rho')
    call set(4, LOG_CALC_HARTREE      , 'Hartree routine')
    call set(5, LOG_CALC_EXC_COR      , 'Exc_Cor routine')
    call set(6, LOG_CALC_TOTAL_ENERGY , 'calculating Etot')
    call set(7, LOG_CALC_ESP          , 'calculating esp')
    call write_loadbalance(fd, 7, tsrc, headers, mode)

    call set(0, 0, 'in rt iterations')
    call set( 1, LOG_CALC_EMFIELD         , 'Emfield')
    call set( 2, LOG_CALC_VBOX            , 'Vbox')
    call set( 3, LOG_CALC_TIME_PROPAGATION, 'time propagation')
    call set( 4, LOG_CALC_RHO             , 'calculating rho')
    call set( 5, LOG_CALC_HARTREE         , 'Hartree routine')
    call set( 6, LOG_CALC_EXC_COR         , 'Exc_Cor routine')
    call set( 7, LOG_CALC_VLOCAL          , 'Vhxc')              ! FIXME: wrong name
    call set( 8, LOG_CALC_DP              , 'calculating Dp')
    call set( 9, LOG_CALC_CURRENT         , 'calculating curr')
    call set(10, LOG_CALC_TOTAL_ENERGY    , 'calculating Etot')
    call set(11, LOG_CALC_PROJECTION      , 'calc. projection')
    call set(12, LOG_CALC_QUADRUPOLE      , 'calc. quadrupole')  ! FIXME: wrong name
    call set(13, LOG_WRITE_ENERGIES       , 'writing energies')
    call set(14, LOG_WRITE_RT_INFOS       , 'writing info etc.')
    call set(15, LOG_RT_ANALYSIS          , 'analysis calc.')
    call set(16, LOG_RT_MISC              , 'misc.')
    call write_loadbalance(fd, 16, tsrc, headers, mode)

    call set(0, 0, 'in subspace-diag')
    call set(1, LOG_DIAG_INIT       , 'initialization')
    call set(2, LOG_DIAG_VLOCAL     , 'Vlocal')
    call set(3, LOG_DIAG_AMAT       , 'Amat')
    call set(4, LOG_DIAG_ALLREDUCE  , 'allreduce')
    call set(5, LOG_DIAG_EIGEN      , 'eigen')
    call set(6, LOG_DIAG_SET_ORBITAL, 'set orbital')
    call set(7, LOG_DIAG_UPDATE     , 'set orbital')
    call write_loadbalance(fd, 7, tsrc, headers, mode)

    call set(0, 0, 'in CG')
    call set(1, LOG_GSCG_TOTAL              , 'total')
    call set(2, LOG_GSCG_INIT               , 'init.')
    call set(3, LOG_GSCG_INIT_ITERATION     , 'init. iter.')
    call set(4, LOG_GSCG_ITERATION          , 'iterations')
    call set(5, LOG_GSCG_DEINIT             , 'deinit.')
    call set(6, LOG_ALLREDUCE_INNER_PRODUCT5, 'comm. for inner product(5)')
    call set(7, LOG_ALLREDUCE_INNER_PRODUCT7, 'comm. for inner product(7)')
    call set(8, LOG_GSCG_ALLREDUCE          , 'comm. for GSCG')
    call write_loadbalance(fd, 8, tsrc, headers, mode)

    call set(0, 0, 'in total_energy_periodic')
    call set(1, LOG_TEP_SENDRECV      , 'sendrecv')
    call set(2, LOG_TEP_ORBITAL_ENERGY, 'orbital energy')
    call set(3, LOG_TEP_ION_ION       , 'ion-ion')
    call set(4, LOG_TEP_ION_ELECTRON  , 'ion-electron')
    call set(5, LOG_TEP_NONLOCAL_1    , 'nonlocal 1')
    call set(6, LOG_TEP_NONLOCAL_2    , 'nonlocal 2')
    call write_loadbalance(fd, 6, tsrc, headers, mode)

    call set(0, 0, 'in current')
    call set(1, LOG_CUR_SENDRECV           , 'sendrecv')
    call set(2, LOG_CUR_LOCAL              , 'current (except nonlocal)')
    call set(3, LOG_CUR_NONLOCAL1          , 'current nonlocal (1)')
    call set(4, LOG_CUR_NONLOCAL1_ALLREDUCE, 'Allreduce nonlocal (1)')
    call set(5, LOG_CUR_NONLOCAL2          , 'current nonlocal (2)')
    call set(6, LOG_CUR_NONLOCAL2_ALLREDUCE, 'Allreduce nonlocal (2)')
    call write_loadbalance(fd, 6, tsrc, headers, mode)

    call set(0, 0, 'communication time')
    call set(1, LOG_ALLREDUCE_RHO    , 'Allreduce in rho')
    call set(2, LOG_ALLREDUCE_HARTREE, 'Allreduce in Hartree')
    call set(3, LOG_ALLREDUCE_DIPOLE , 'Allreduce in dipole calc.')
    call set(4, LOG_ALLREDUCE_ESP    , 'Allreduce in esp')
    call set(5, LOG_ALLGATHERV_TOTAL , 'Allgatherv')
    call write_loadbalance(fd, 5, tsrc, headers, mode)

    call write_root(fd, '=== performance [GFLOPS] ===')
    call write_flops(fd, mode)
  contains
    subroutine set(nid, tid, header)
      implicit none
      integer, intent(in)      :: nid, tid
      character(*), intent(in) :: header
      if (nid > 0) then
        tsrc(nid) = timer_get(tid)
      end if
      write (headers(nid),'(a)') header
    end subroutine
  end subroutine write_performance
end module perflog
