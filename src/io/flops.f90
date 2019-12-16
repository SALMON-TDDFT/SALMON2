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
module flops
  implicit none

  public get_hamiltonian_flops

  private summation_threads, get_gflops, get_hamiltonian_chunk_size
  private get_stencil_FLOP, get_pseudo_pt_FLOP, get_update_FLOP

contains
  subroutine get_hamiltonian_flops(lgflops,pgflops,mgflops,sgflops)
    use global_variables
    use salmon_parallel
    use communication
    use timer
    implicit none
    real(8),intent(out) :: lgflops(4) ! processor
    real(8),intent(out) :: pgflops(4) ! processor max
    real(8),intent(out) :: mgflops(4) ! macro grid
    real(8),intent(out) :: sgflops(4) ! system total

    type(comm_maxloc_type) :: tin,tout

    call summation_threads(lgflops)

    pgflops  = lgflops
    tin%rank = nproc_id_global
    tin%val  = lgflops(4)
    call comm_get_max(tin, tout, nproc_group_global)
    call comm_bcast(pgflops, nproc_group_global, tout%rank)

    if (calc_mode == calc_mode_ms) then
       call comm_summation(lgflops, mgflops, 4, nproc_group_tdks)
    end if

    call comm_summation(lgflops, sgflops, 4, nproc_group_global)
  end subroutine

  subroutine summation_threads(lgflops)
    use global_variables, only: NUMBER_THREADS, functional, propagator
    use timer
    implicit none
    real(8), intent(out) :: lgflops(4)
    real(8) :: hflop(3), htime(4)
    integer :: i, cnt
    integer :: chunk_size(0:NUMBER_THREADS-1)
    integer :: ncalls_in_loop

    select case(functional)
    case('VS98','TPSS','TBmBJ', 'tbmbj', "BJ_PW", "bj_pw")
        ncalls_in_loop = 3
      case default
        ncalls_in_loop = 2
    end select

    select case(propagator)
      case('middlepoint')
        ncalls_in_loop = ncalls_in_loop - 1
    end select

    call get_hamiltonian_chunk_size(chunk_size)

    lgflops = 0.0d0
    do i=0,NUMBER_THREADS-1
      cnt = chunk_size(i)

      hflop(1) = get_stencil_FLOP(cnt)   * ncalls_in_loop
      hflop(2) = get_pseudo_pt_FLOP(cnt) * ncalls_in_loop
      hflop(3) = get_update_FLOP(cnt)    * ncalls_in_loop

      htime(1) = timer_thread_get(LOG_HPSI_STENCIL, i)
      htime(2) = timer_thread_get(LOG_HPSI_PSEUDO, i)
      htime(3) = timer_thread_get(LOG_HPSI_UPDATE, i)
      htime(4) = timer_thread_get(LOG_HPSI_INIT, i)

      lgflops(1) = lgflops(1) + get_gflops(hflop(1), htime(1))
      lgflops(2) = lgflops(2) + get_gflops(hflop(2), htime(2))
      lgflops(3) = lgflops(3) + get_gflops(hflop(3), htime(3))
      lgflops(4) = lgflops(4) + get_gflops(sum(hflop), sum(htime))
    end do
  end subroutine

  subroutine get_hamiltonian_chunk_size(chunk_size)
    use global_variables, only: NKB,NUMBER_THREADS
    use omp_lib, only: omp_get_thread_num
    implicit none
    integer,intent(out) :: chunk_size(0:NUMBER_THREADS-1)
    integer :: ikb, tid

    tid = 0
    chunk_size(:) = 0
!$omp parallel private(tid) shared(chunk_size)
!$ tid = omp_get_thread_num()
!$omp do private(ikb)
    do ikb=1,NKB
      chunk_size(tid) = chunk_size(tid) + 1
    end do
!$omp end do
!$omp end parallel
  end subroutine

  function get_stencil_FLOP(chunk_size)
    use global_variables, only: NK_s,NK_e,nmacro_s,nmacro_e,NBoccmax,NL,Nt
    integer,intent(in),optional :: chunk_size
    real(8),parameter           :: FLOP = 158

    real(8) :: get_stencil_FLOP
    integer :: nsize

    if(present(chunk_size)) then
      nsize = chunk_size &
            * (nmacro_e - nmacro_s + 1)
    else
      nsize = (NK_e  - NK_s  + 1) * NBoccmax &
            * (nmacro_e - nmacro_s + 1)
    end if
    get_stencil_FLOP = nsize * 4*FLOP*NL * (Nt + 1)
  end function

  function get_pseudo_pt_FLOP(chunk_size)
    use global_variables, only: NK_s,NK_e,nmacro_s,nmacro_e,NBoccmax,Nt,a_tbl,Mps
    implicit none
    integer,intent(in),optional :: chunk_size
    real(8),parameter           :: FLOP_reduction = (2 + 6)     + 2
    real(8),parameter           :: FLOP_scatter   = (2 + 6 + 1) + 2 ! 1 = conjg(z)
    real(8),parameter           :: FLOP_scalar    = 2 + 2

    real(8) :: get_pseudo_pt_FLOP
    real(8) :: FLOP
    integer :: nsize

    FLOP = FLOP_scalar + (FLOP_reduction + FLOP_scatter) * sum(Mps(a_tbl(:)))

    if(present(chunk_size)) then
      nsize = chunk_size &
            * (nmacro_e - nmacro_s + 1)
    else
      nsize = (NK_e  - NK_s  + 1) * NBoccmax &
            * (nmacro_e - nmacro_s + 1)
    endif
    get_pseudo_pt_FLOP = nsize * 4*FLOP * (Nt + 1)
  end function

  function get_update_FLOP(chunk_size)
    use global_variables, only: NK_s,NK_e,nmacro_s,nmacro_e,NBoccmax,NL,Nt
    implicit none
    integer,intent(in),optional :: chunk_size
    real(8),parameter           :: FLOP = 6 + 2

    real(8) :: get_update_FLOP
    integer :: nsize

    if(present(chunk_size)) then
      nsize = chunk_size &
            * (nmacro_e - nmacro_s + 1)
    else
      nsize = (NK_e  - NK_s  + 1) * NBoccmax &
            * (nmacro_e - nmacro_s + 1)
    endif
    get_update_FLOP = nsize * 4*FLOP*NL * (Nt + 1)
  end function

  function get_gflops(FLOP,time) result(x)
    use math_constants
    implicit none
    real(8),intent(in) :: FLOP
    real(8),intent(in) :: time
    real(8) :: x
    if (is_zero(time)) then
      x = 0.d0
    else
      x = FLOP / (time * (10**9))
    end if
  end function
end module
