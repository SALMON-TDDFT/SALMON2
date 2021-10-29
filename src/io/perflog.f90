!
!  Copyright 2019-2020 SALMON developers
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
    use parallelization
    use communication
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
            write (fd,'(a30,6(e12.5e2))') headers(i),tmin(i),tmax(i),tavg(i),tstdev(i),tdif(i),trel(i)
          else
            write (fd,'(a,",",5(e12.5e2,","),e12.5e2)') trim(headers(i)),tmin(i),tmax(i),tavg(i),tstdev(i),tdif(i),trel(i)
          end if
        end if
      end do
    end if
  end subroutine

  subroutine write_flops(fd,write_mode)
    use parallelization
    use communication
    use math_constants
!    use flops
    use salmon_global, only: theory, iperiodic
    implicit none
    integer, intent(in) :: fd, write_mode

    real(8) :: lg(4),pg(4),mg(4),sg(4)

    select case(theory)
      case('tddft')
        select case(iperiodic)
          ! TODO: calc flops on unified hpsi
          ! 0: GCEED
          ! 3: ARTED
          case(3)
            !  if (yn_domain_parallel == 'y') then
            !  call get_hamiltonian_flops(lg,pg,mg,sg)
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
            !end if
        end select
    end select
  end subroutine

  subroutine write_root(fd, str)
    use parallelization
    use communication
    implicit none
    integer, intent(in)      :: fd
    character(*), intent(in) :: str
    if (comm_is_root(nproc_id_global)) then
      write (fd,'(a)') str
    end if
  end subroutine

  subroutine write_performance(fd,mode)
    use parallelization
    use communication
    use salmon_global, only: theory,method_singlescale,yn_ffte
    use timer
    implicit none
    integer, intent(in) :: fd, mode

    integer,parameter :: LOG_SIZE = 30
    real(8)       :: tsrc(LOG_SIZE)
    character(30) :: headers(0:LOG_SIZE)

    call write_root(fd, '=== application breakdown [s] ===')
    if (comm_is_root(nproc_id_global)) then
      write (fd,'(a,f0.8)') 'total calculation time,',timer_get(LOG_TOTAL)
    end if

    select case(theory)
    case('dft','dft_band','dft_md')
      call set(0, 0, 'scf calculation')
      call set(1, LOG_INIT_GS            , 'gs initialization')
      call set(2, LOG_INIT_GS_RESTART    , 'reading data')
      call set(3, LOG_INIT_GS_ITERATION  , 'init. scf iter.')
      call set(4, LOG_GS_ITERATION       , 'scf iterations')
      call set(5, LOG_DEINIT_GS_ITERATION, 'deinit. scf iter.')
      call set(6, LOG_WRITE_GS_RESULTS   , 'writing scf data')
      call set(7, LOG_WRITE_GS_DATA     , 'writing LDA data')
      call set(8, LOG_WRITE_GS_INFO     , 'writing LDA infos')
      call write_loadbalance(fd, 8, tsrc, headers, mode)

      call set(0, 0, 'in scf iterations')
      call set_sub( 1, LOG_CALC_MINIMIZATION  , 'Minimization')
      call set_sub( 2, LOG_CALC_GRAM_SCHMIDT  , 'Gram Schmidt')
      call set_sub( 3, LOG_CALC_SUBSPACE_DIAG , 'subspace-diag.')
      call set_sub( 4, LOG_CALC_RHO           , 'calculating rho')
      call set_sub( 5, LOG_CALC_HARTREE       , 'Hartree routine')
      call set_sub( 6, LOG_CALC_EXC_COR       , 'Exc_Cor routine')
      call set_sub( 7, LOG_CALC_TOTAL_ENERGY  , 'calculating Etot')
      call set_sub( 8, LOG_CALC_ESP           , 'calculating esp')
      call set_sub( 9, LOG_CALC_ION_FORCE     , 'calc force')
      call set_sub(10, LOG_MD_TEVOL_PART1     , 'MD-opt time-evol. part1')
      call set_sub(11, LOG_MD_TEVOL_PART2     , 'MD-opt time-evol. part2')
      call set_sub(12, LOG_MD_UPDATE_PSEUDO_PT, 'update pseudo-pt')
      call write_loadbalance(fd, 12, tsrc, headers, mode)
    case('dft2tddft')
      call set(0, 0, 'DFT data redistribution')
      call set(1, LOG_INIT_GS         , 'gs initialization')
      call set(2, LOG_INIT_GS_RESTART , 'reading data')
      call set(3, LOG_INIT_RT         , 'rt initialization')
      call set(4, LOG_WRITE_RT_DATA   , 'data redistribution')
      call set(5, LOG_WRITE_RT_RESULTS, 'writing data')
      call write_loadbalance(fd, 5, tsrc, headers, mode)
    case('tddft_response','tddft_pulse','single_scale_maxwell_tddft',&
         'multi_scale_maxwell_tddft','maxwell')
      call set(0, 0, 'rt calculation')
      call set(1, LOG_INIT_RT              , 'rt initialization')
      call set(2, LOG_READ_GS_DATA        , 'reading lda data')
      call set(3, LOG_READ_RT_DATA         , 'reading rt data')
      call set(4, LOG_INIT_TIME_PROPAGATION, 'prep. time prop.')
      call set(5, LOG_RT_ITERATION         , 'rt iterations')
      call set(6, LOG_WRITE_RT_DATA        , 'writing rt data')
      call set(7, LOG_WRITE_RT_RESULTS     , 'after writing rt data')
      call write_loadbalance(fd, 7, tsrc, headers, mode)

      call set(0, 0, 'in rt iterations')
      call set_sub( 1, LOG_CALC_EMFIELD         , 'Emfield')
      call set_sub( 2, LOG_CALC_VBOX            , 'Vbox')
      call set_sub( 3, LOG_CALC_TIME_PROPAGATION, 'time propagation')
      call set_sub( 4, LOG_CALC_RHO             , 'calculating rho')
      if (theory=='single_scale_maxwell_tddft' .and. method_singlescale=='1d_fourier' .and. yn_ffte=='y') then
        call set_sub( 5, LOG_CALC_HARTREE       , 'Hartree+FDTD by FFTE')
      else
        call set_sub( 5, LOG_CALC_HARTREE       , 'Hartree routine')
      end if
      call set_sub( 6, LOG_CALC_EXC_COR         , 'Exc_Cor routine')
      call set_sub( 7, LOG_CALC_ALLGATHERV_VLOCAL, 'allgatherv_vlocal')
      call set_sub( 8, LOG_CALC_DP              , 'calculating Dp')
      call set_sub( 9, LOG_CALC_DENSITY_MATRIX  , 'calculating density matrix')
      call set_sub(10, LOG_CALC_CURRENT         , 'calculating curr')
      call set_sub(11, LOG_CALC_TOTAL_ENERGY    , 'calculating Etot')
      call set_sub(12, LOG_CALC_PROJECTION      , 'calc. projection')
      call set_sub(13, LOG_WRITE_ENERGIES       , 'writing energies')
      call set_sub(14, LOG_CALC_EIGEN_ENERGY    , 'calc_eigen_energy')
      call set_sub(15, LOG_CALC_CURRENT_ION     , 'calc_current_ion')
      call set_sub(16, LOG_CALC_TOTAL_ENERGY_PERIODIC, 'calc_total_energy_periodic')
      call set_sub(17, LOG_CALC_SINGLESCALE     , 'calc singlescale')
      call set_sub(18, LOG_CALC_ION_FORCE       , 'calc force')
      call set_sub(19, LOG_MD_TEVOL_PART1       , 'MD-opt time-evol. part1')
      call set_sub(20, LOG_MD_TEVOL_PART2       , 'MD-opt time-evol. part2')
      call set_sub(21, LOG_MD_UPDATE_PSEUDO_PT  , 'update pseudo-pt')
      call set_sub(22, LOG_WRITE_RT_INFOS       , 'writing info etc.')
      call set_sub(23, LOG_RT_ANALYSIS          , 'analysis calc.')
      call set_sub(24, LOG_RT_MISC              , 'misc.')
      call write_loadbalance(fd, 24, tsrc, headers, mode)
    case default
      stop 'invalid theory @ perflog'
    end select

    if (theory == 'single_scale_maxwell_tddft') then
      call set(0, 0, 'in singlescale maxwell-tddft')
      call set_sub(1, LOG_SS_FDTD_CALC,      'FDTD calc')
      call set_sub(2, LOG_SS_FDTD_COMM,      'FDTD halo comm')
      call set_sub(3, LOG_SS_FDTD_COMM_COLL, 'FDTD coll comm')
      call set_sub(4, LOG_SS_UPDATE_NONLOCALPT_MICROAC, 'update nonlocal microAc')
      call write_loadbalance(fd, 4, tsrc, headers, mode)
    end if

    call set(0, 0, 'total_energy module')
    call set_sub(1, LOG_TE_ISOLATED_CALC,       'isolated calc.')
    call set_sub(2, LOG_TE_ISOLATED_COMM_COLL,  'isolated comm. coll')
    call set_sub(3, LOG_TE_PERIODIC_CALC,       'periodic calc.')
    call set_sub(4, LOG_TE_PERIODIC_COMM_COLL,  'periodic comm. coll.')
    call set_sub(5, LOG_EIGEN_ENERGY_CALC,      'eigen_energy calc.')
    call set_sub(6, LOG_EIGEN_ENERGY_COMM_COLL, 'eigen_energy comm. coll.')
    call set_sub(7, LOG_EIGEN_ENERGY_HPSI,      'eigen_energy hpsi')
    call write_loadbalance(fd, 7, tsrc, headers, mode)

    call set(0, 0, 'density_matrix module')
    call set_sub( 1, LOG_DENSITY_MATRIX_CALC,          'density_matrix calc.')
    call set_sub( 2, LOG_DENSITY_MATRIX_COMM_COLL,     'density_matrix comm. coll.')
    call set_sub( 3, LOG_DENSITY_MATRIX_COMM_HALO,     'density_matrix comm. halo')
    call set_sub( 4, LOG_DENSITY_CALC,                 'density calc.')
    call set_sub( 5, LOG_DENSITY_COMM_COLL,            'density comm. coll.')
    call set_sub( 6, LOG_CURRENT_CALC,                 'current calc.')
    call set_sub( 7, LOG_CURRENT_CALC_UVPSI_RDIVIDED,  'current rDivided uVpsi calc.')
    call set_sub( 8, LOG_CURRENT_COMM_COLL,            'current comm. coll.')
    call set_sub( 9, LOG_CURRENT_COMM_HALO,            'current comm. halo.')
    call set_sub(10, LOG_MCURRENT_CALC,                'micro-current calc.')
    call set_sub(11, LOG_MCURRENT_COMM_COLL,           'micro-current comm. coll.')
    call set_sub(12, LOG_CALC_STENCIL_CURRENT,           'calc stencil current')
    call write_loadbalance(fd, 12, tsrc, headers, mode)

    call set(0, 0, 'conjugate_gradient module')
    call set_sub(1, LOG_GSCG_ISOLATED_CALC,      'isolated calc.')
    call set_sub(2, LOG_GSCG_ISOLATED_COMM_COLL, 'isolated comm. coll.')
    call set_sub(3, LOG_GSCG_ISOLATED_HPSI,      'isolated hpsi')
    call set_sub(4, LOG_GSCG_PERIODIC_CALC,      'periodic calc.')
    call set_sub(5, LOG_GSCG_PERIODIC_COMM_COLL, 'periodic comm. coll.')
    call set_sub(6, LOG_GSCG_PERIODIC_HPSI,      'periodic hpsi')
    call write_loadbalance(fd, 6, tsrc, headers, mode)

    call set(0, 0, 'subspace_diag module')
    call set_sub(1, LOG_SSDG_ISOLATED_CALC,      'isolated calc.')
    call set_sub(2, LOG_SSDG_ISOLATED_COMM_COLL, 'isolated comm. coll.')
    call set_sub(3, LOG_SSDG_ISOLATED_HPSI,      'isolated hpsi')
    call set_sub(4, LOG_SSDG_PERIODIC_CALC,      'periodic calc.')
    call set_sub(5, LOG_SSDG_PERIODIC_COMM_COLL, 'periodic comm. coll.')
    call set_sub(6, LOG_SSDG_PERIODIC_HPSI,      'periodic hpsi')
    call set_sub(7, LOG_SSDG_PERIODIC_EIGEN,     'periodic eigen')
    call write_loadbalance(fd, 7, tsrc, headers, mode)

    call set(0, 0, 'subspace_diag so module')
    call set_sub(1, LOG_SSDG_SO_ISOLATED_CALC,      'isolated calc.')
    call set_sub(2, LOG_SSDG_SO_ISOLATED_COMM_COLL, 'isolated comm. coll.')
    call set_sub(3, LOG_SSDG_SO_ISOLATED_HPSI,      'isolated hpsi')
    call set_sub(4, LOG_SSDG_SO_PERIODIC_CALC,      'periodic calc.')
    call set_sub(5, LOG_SSDG_SO_PERIODIC_COMM_COLL, 'periodic comm. coll.')
    call set_sub(6, LOG_SSDG_SO_PERIODIC_HPSI,      'periodic hpsi')
    call write_loadbalance(fd, 6, tsrc, headers, mode)

    call set(0, 0, 'force module')
    call set_sub(1, LOG_CALC_ION_FORCE,      'total')
    call set_sub(2, LOG_CALC_FORCE_ION_ION,  'calc ion-ion')
    call set_sub(3, LOG_CALC_FORCE_FOURIER,  'calc fouier')
    call set_sub(4, LOG_CALC_FORCE_ELEC_ION, 'calc electron-ion')
    call set_sub(5, LOG_CALC_FORCE_GTPSI,    '::calc gtpsi')
    call set_sub(6, LOG_CALC_FORCE_DDEN,     '::calc dden')
    call set_sub(7, LOG_CALC_FORCE_NONLOCAL, '::calc nonlocal')
    call set_sub(8, LOG_CALC_FORCE_LOCAL,    'calc local')
    call write_loadbalance(fd, 8, tsrc, headers, mode)

    call set(0, 0, 'init_ps')
    call set_sub(1, LOG_INIT_PS_TOTAL,      'total')
    call set_sub(2, LOG_INIT_PS_CALC_NPS,   '::calc nps')
    call set_sub(3, LOG_INIT_PS_CALC_JXYZ,  '::calc jxyz')
    call set_sub(4, LOG_INIT_PS_LMA_UV,     '::calc lma and uv')
    call set_sub(5, LOG_INIT_PS_CALC_VPSL,  '::calc vpsl')
    call set_sub(6, LOG_INIT_PS_UVPSI,      '::init uvpsi')
    call write_loadbalance(fd, 6, tsrc, headers, mode)


    call set(0, 0, 'hamiltonian module')
    call set_sub(1, LOG_UHPSI_ALL,            'total')
    call set_sub(2, LOG_UHPSI_UPDATE_OVERLAP, 'update overlap (comm.)')
    call set_sub(3, LOG_UHPSI_STENCIL,        'stencil')
    call set_sub(4, LOG_UHPSI_SUBTRACTION,    'subtraction')
    call set_sub(5, LOG_UHPSI_PSEUDO,         'pseudo-pt')
    call set_sub(6, LOG_UHPSI_PSEUDO_COMM,    'pseudo-pt (comm.)')
    call write_loadbalance(fd, 6, tsrc, headers, mode)


    call set(0, 0, 'hpsi-overlapped')
    call set_sub(1, LOG_UHPSI_OVL_PHASE1,      'pack halo')
    call set_sub(2, LOG_UHPSI_OVL_PHASE2,      'comp/comm overlap')
    call set_sub(3, LOG_UHPSI_OVL_PHASE2_COMM, 'halo communication')
    call set_sub(4, LOG_UHPSI_OVL_PHASE3,      'unpack halo')
    call set_sub(5, LOG_UHPSI_OVL_PHASE4,      'halo computation')
    call write_loadbalance(fd, 5, tsrc, headers, mode)

    call set(0, 0, 'checkpoint/restart')
    call set(1, LOG_RESTART_SELF,         'restart self')
    call set(2, LOG_RESTART_SYNC,         'restart sync')
    call set(3, LOG_CHECKPOINT_SELF,      'checkpoint self')
    call set(4, LOG_CHECKPOINT_SYNC,      'checkpoint sync')
    call write_loadbalance(fd, 4, tsrc, headers, mode)

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

    subroutine set_sub(nid, tid, header)
      implicit none
      integer, intent(in)      :: nid, tid
      character(*), intent(in) :: header
      if (nid > 0) then
        tsrc(nid) = timer_get_sub(tid)
      end if
      write (headers(nid),'(a)') header
    end subroutine
  end subroutine write_performance
end module perflog
