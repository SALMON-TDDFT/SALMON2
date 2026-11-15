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

!=======================================================================

! Explicit E/H stepping on a supplied field communicator.
! Native E_n clock, A/J at (n+1/2)dt; sample map matches eh_sendrecv('s').
#include "config.h"
module salmonpy_eh_context
  use iso_c_binding
  use structures, only: s_fdtd_system
  use fdtd_eh, only: ls_fdtd_eh,eh_init,eh_calc,eh_finalize
  use lorentz_drude, only: s_ld_state,ld_finalize
  use sendrecv_grid, only: dealloc_cache
  use communication, only: comm_summation
  use parallelization, only: setup_parallel,nproc_group_global
  implicit none
  private
  public :: sp_eh_init,sp_eh_step,sp_eh_close,sp_eh_sample_h
  type :: field_context
    type(s_fdtd_system) :: fs
    type(ls_fdtd_eh) :: fe
    type(s_ld_state) :: ld
    integer,allocatable :: points(:,:)
    real(8),allocatable :: current(:,:,:,:)
    integer :: step=0,capacity
  end type
  type(field_context),allocatable,save :: state
  logical,save :: used=.false.
contains
  integer(c_int) function sp_eh_init(path,npoints,indices) bind(C)
    use inputoutput, only: read_input
    use math_constants, only: set_math_constants
    use salmon_global, only: theory,nt_em,yn_restart
    character(c_char),intent(in) :: path(*)
    integer(c_int),value :: npoints
    integer(c_int),intent(in) :: indices(3,npoints)
    character(4096) :: filename
    integer :: i,p,lo(3),hi(3)
    sp_eh_init=1
    if(used.or.npoints<1) return
    used=.true.;filename=''
    do i=1,len(filename)
      if(path(i)==c_null_char) exit
      filename(i:i)=path(i)
    end do
    if(i>len(filename)) return
    call set_math_constants
    call setup_parallel
    call read_input(trim(filename))
    if(theory/='maxwell'.or.yn_restart/='n') return
    allocate(state)
    call eh_init(state%fs,state%fe,state%ld)
    state%capacity=nt_em
    allocate(state%points(3,npoints))
    do p=1,npoints
      state%points(:,p)=state%fs%lg%is+indices(:,p)
    end do
    lo=state%fs%mg%is_array;hi=state%fs%mg%ie_array
    allocate(state%current(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    state%current=0d0
    sp_eh_init=0
  end function

  integer(c_int) function sp_eh_step(npoints,J,E) bind(C)
    integer(c_int),value :: npoints
    real(c_double),intent(in) :: J(3,npoints)
    real(c_double),intent(out) :: E(3,npoints)
    real(8) :: local_E(3,npoints)
    integer :: p,d,side,index(3),i,jj,k
    sp_eh_step=1
    if(.not.allocated(state)) return
    if(state%step>=state%capacity.or.npoints/=size(state%points,2)) return
    state%current=0d0
    ! S transpose: each macro density contributes half to each adjacent edge.
    ! One macro point per grid cell, equal cell volumes; no extra volume factor.
    do p=1,npoints
      do d=1,3
        do side=0,1
          index=state%points(:,p);index(d)=index(d)-side
          if(all(index>=state%fs%mg%is).and.all(index<=state%fs%mg%ie)) &
            state%current(d,index(1),index(2),index(3))= &
            state%current(d,index(1),index(2),index(3))+.5d0*J(d,p)
        end do
      end do
    end do
    state%step=state%step+1
    state%fe%iter_sta=state%step;state%fe%iter_end=state%step
    call eh_calc(state%fs,state%fe,state%ld,state%current)
    local_E=0d0
    do p=1,npoints
      index=state%points(:,p)
      if(.not.(all(index>=state%fs%mg%is).and.all(index<=state%fs%mg%ie))) cycle
      i=index(1);jj=index(2);k=index(3)
      local_E(1,p)=.5d0*(state%fe%ex_y(i,jj,k)+state%fe%ex_z(i,jj,k) &
                       +state%fe%ex_y(i-1,jj,k)+state%fe%ex_z(i-1,jj,k))
      local_E(2,p)=.5d0*(state%fe%ey_z(i,jj,k)+state%fe%ey_x(i,jj,k) &
                       +state%fe%ey_z(i,jj-1,k)+state%fe%ey_x(i,jj-1,k))
      local_E(3,p)=.5d0*(state%fe%ez_x(i,jj,k)+state%fe%ez_y(i,jj,k) &
                       +state%fe%ez_x(i,jj,k-1)+state%fe%ez_y(i,jj,k-1))
    end do
    call comm_summation(local_E,E,3*npoints,nproc_group_global)
    sp_eh_step=0
  end function

  ! Read raw H at (step+1/2)*dt, centered spatially on each macro cell.
  ! Each rank contributes only owned Yee faces; no diagonal halo assumption.
  integer(c_int) function sp_eh_sample_h(npoints,H) bind(C)
    integer(c_int),value :: npoints
    real(c_double),intent(out) :: H(3,npoints)
    real(8) :: local_H(3,npoints),value
    integer :: p,d,a,b,s1,s2,index(3),i,j,k
    sp_eh_sample_h=1
    if(.not.allocated(state)) return
    if(state%step<1.or.npoints/=size(state%points,2)) return
    local_H=0d0
    do p=1,npoints
      do d=1,3
        a=mod(d,3)+1;b=mod(d+1,3)+1
        do s1=0,1
          do s2=0,1
            index=state%points(:,p);index(a)=index(a)-s1;index(b)=index(b)-s2
            if(.not.(all(index>=state%fs%mg%is).and.all(index<=state%fs%mg%ie))) cycle
            i=index(1);j=index(2);k=index(3)
            select case(d)
            case(1)
              value=state%fe%hx_y(i,j,k)+state%fe%hx_z(i,j,k)
            case(2)
              value=state%fe%hy_z(i,j,k)+state%fe%hy_x(i,j,k)
            case(3)
              value=state%fe%hz_x(i,j,k)+state%fe%hz_y(i,j,k)
            end select
            local_H(d,p)=local_H(d,p)+.25d0*value
          end do
        end do
      end do
    end do
    call comm_summation(local_H,H,3*npoints,nproc_group_global)
    sp_eh_sample_h=0
  end function

  integer(c_int) function sp_eh_close() bind(C)
    sp_eh_close=1
    if(.not.allocated(state)) return
    if(state%step==state%capacity) call eh_finalize(state%fs,state%fe,state%ld)
    call ld_finalize(state%ld)
    call dealloc_cache(state%fs%srg_ng)
    deallocate(state)
    sp_eh_close=0
  end function
end module
