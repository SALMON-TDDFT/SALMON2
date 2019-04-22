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

contains
  subroutine write_gs_performance(fd)
    use timer
    implicit none
    integer, intent(in)       :: fd

    write(*,'(a)') "==================== elapsed time ===================="
    call timer_write(fd, 'elapsed time initializing [s]         = ', LOG_INIT_GS)
    call timer_write(fd, 'elapsed time for reading data [s]     = ', LOG_INIT_GS_RESTART)
    call timer_write(fd, 'elapsed time for init. scf iter. [s]  = ', LOG_INIT_GS_ITERATION)
    call timer_write(fd, 'elapsed time for scf iterations [s]   = ', LOG_GS_ITERATION)
    call timer_write(fd, 'elapsed time for deinit. scf iter. [s]= ', LOG_DEINIT_GS_ITERATION)
    call timer_write(fd, 'elapsed time for writing data [s]     = ', LOG_WRITE_RESULTS)
    call timer_write(fd, 'elapsed time for writing LDA data [s] = ', LOG_WRITE_LDA_DATA)
    call timer_write(fd, 'elapsed time for writing infos [s]    = ', LOG_WRITE_INFOS)
    call timer_write(fd, 'total time [s]                        = ', LOG_TOTAL)
    write(*,'(a)') "======================================================"
    write(*,'(a)') "================== in scf iterations ================="
    call timer_write(fd, 'elapsed time for Gram Schmidt [s]     = ', LOG_CALC_GRAM_SCHMIDT)
    call timer_write(fd, 'elapsed time for subspace-diag. [s]   = ', LOG_CALC_SUBSPACE_DIAG)
    call timer_write(fd, 'elapsed time for calculating rho  [s] = ', LOG_CALC_RHO)
    call timer_write(fd, 'elapsed time for Hartree routine  [s] = ', LOG_CALC_HARTREE)
    call timer_write(fd, 'elapsed time for Exc_Cor routine  [s] = ', LOG_CALC_EXC_COR)
    call timer_write(fd, 'elapsed time for calculating Etot [s] = ', LOG_CALC_TOTAL_ENERGY)
    call timer_write(fd, 'elapsed time for writing info. [s]    = ', LOG_WRITE_RESULTS)
    write(*,'(a)') "======================================================"
    write(*,'(a)') "================== in subspace-diag. ================="
    call timer_write(fd, 'elapsed time for initialization [s]   = ', LOG_DIAG_INIT)
    call timer_write(fd, 'elapsed time for Vlocal [s]           = ', LOG_DIAG_VLOCAL)
    call timer_write(fd, 'elapsed time for Amat [s]             = ', LOG_DIAG_AMAT)
    call timer_write(fd, 'elapsed time for allreduce [s]        = ', LOG_DIAG_ALLREDUCE)
    call timer_write(fd, 'elapsed time for eigen [s]            = ', LOG_DIAG_EIGEN)
    call timer_write(fd, 'elapsed time for set orbital [s]      = ', LOG_DIAG_SET_ORBITAL)
    call timer_write(fd, 'elapsed time for set orbital [s]      = ', LOG_DIAG_UPDATE)
    write(*,'(a)') "======================================================"
    write(*,'(a)') "================== in CG. ============================"
    call timer_write(fd, 'total time for CG [s]                 = ', LOG_GSCG_TOTAL)
    call timer_write(fd, 'elapsed time for init. CG [s]         = ', LOG_GSCG_INIT)
    call timer_write(fd, 'elapsed time for init. iter. CG [s]   = ', LOG_GSCG_INIT_ITERATION)
    call timer_write(fd, 'elapsed time for iterations CG [s]    = ', LOG_GSCG_ITERATION)
    call timer_write(fd, 'elapsed time for deinit. CG [s]       = ', LOG_GSCG_DEINIT)
    call timer_write(fd, 'comm. for inner product (5) in CG [s] = ', LOG_ALLREDUCE_INNER_PRODUCT5)
    call timer_write(fd, 'comm. for inner product (7) in CG [s] = ', LOG_ALLREDUCE_INNER_PRODUCT7)
    write(*,'(a)') "======================================================"
  end subroutine

  subroutine write_rt_performance(fd)
    use timer
    implicit none
    integer, intent(in)       :: fd

    write(fd,'(a)') "==================== elapsed time ===================="
    call timer_write(fd, 'elapsed time for rt initialization  = ', LOG_INIT_RT)
    call timer_write(fd, 'elapsed time for reading lda data   = ', LOG_READ_LDA_DATA)
    call timer_write(fd, 'elapsed time for reading rt data    = ', LOG_READ_RT_DATA)
    call timer_write(fd, 'elapsed time for prep. time prop.   = ', LOG_INIT_TIME_PROPAGATION)
    call timer_write(fd, 'elapsed time for rt iterations      = ', LOG_RT_ITERATION)
    call timer_write(fd, 'elapsed time for writing rt data    = ', LOG_WRITE_RT_DATA)
    call timer_write(fd, 'elapsed time aft. writing rt data   = ', LOG_WRITE_RESULTS)
    call timer_write(fd, 'total time                          = ', LOG_TOTAL)

    write(fd,'(a)') "======================================================"
    write(fd,'(a)') "=========== elapsed time for rt iterations ==========="
    call timer_write(fd, 'elapsed time for Vbox               = ', LOG_CALC_VBOX)
    call timer_write(fd, 'elapsed time for time propagation   = ', LOG_CALC_TIME_PROPAGATION)
    call timer_write(fd, 'elapsed time for calculating rho    = ', LOG_CALC_RHO)
    call timer_write(fd, 'elapsed time for Hartree routine    = ', LOG_CALC_HARTREE)
    call timer_write(fd, 'elapsed time for Exc_Cor routine    = ', LOG_CALC_EXC_COR)
    call timer_write(fd, 'elapsed time for Vhxc               = ', LOG_CALC_VLOCAL) ! FIXME: wrong name
    call timer_write(fd, 'elapsed time for calculating Dp     = ', LOG_CALC_DP)
    call timer_write(fd, 'elapsed time for calculating curr   = ', LOG_CALC_CURRENT)
    call timer_write(fd, 'elapsed time for calculating Etot   = ', LOG_CALC_TOTAL_ENERGY)
    call timer_write(fd, 'elapsed time for calc. projection   = ', LOG_CALC_PROJECTION)
    call timer_write(fd, 'elapsed time for calc. quadrupole   = ', LOG_CALC_QUADRUPOLE) ! FIXME: wrong name
    call timer_write(fd, 'elapsed time for writing energies   = ', LOG_WRITE_ENERGIES)
    call timer_write(fd, 'elapsed time for writing info etc.  = ', LOG_WRITE_INFOS)
    write(*,'(a)') "======================================================"
    write(*,'(a)') "======================================================"
    write(*,'(a)') "=========== communication time ======================="
    call timer_write(fd, 'Allreduce in calculating rho        = ', LOG_ALLREDUCE_RHO)
    call timer_write(fd, 'Allreduce in Hartree                = ', LOG_ALLREDUCE_HARTREE)
    call timer_write(fd, 'Allreduce in dipole calc.           = ', LOG_ALLREDUCE_DIPOLE)
    call timer_write(fd, 'Allgatherv                          = ', LOG_ALLGATHERV_TOTAL)
    write(*,'(a)') "======================================================"
    write(*,'(a)') "=========== total_energy_periodic ===================="
    call timer_write(fd, 'sendrecv                            = ', LOG_TEP_SENDRECV)
    call timer_write(fd, 'orbital energy                      = ', LOG_TEP_ORBITAL_ENERGY)
    call timer_write(fd, 'ion-ion                             = ', LOG_TEP_ION_ION)
    call timer_write(fd, 'ion-electron                        = ', LOG_TEP_ION_ELECTRON)
    call timer_write(fd, 'nonlocal 1                          = ', LOG_TEP_NONLOCAL_1)
    call timer_write(fd, 'nonlocal 2                          = ', LOG_TEP_NONLOCAL_2)
    write(*,'(a)') "=========== current =================================="
    call timer_write(fd, 'sendrecv                            = ', LOG_CUR_SENDRECV)
    call timer_write(fd, 'current (except nonlocal)           = ', LOG_CUR_LOCAL)
    call timer_write(fd, 'current nonlocal (1)                = ', LOG_CUR_NONLOCAL1)
    call timer_write(fd, 'Allreduce nonlocal (1)              = ', LOG_CUR_NONLOCAL1_ALLREDUCE)
    call timer_write(fd, 'current nonlocal (2)                = ', LOG_CUR_NONLOCAL2)
    call timer_write(fd, 'Allreduce nonlocal (2)              = ', LOG_CUR_NONLOCAL2_ALLREDUCE)
  end subroutine write_rt_performance
end module write_performance_results
