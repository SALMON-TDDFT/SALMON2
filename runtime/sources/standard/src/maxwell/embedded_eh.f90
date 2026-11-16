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
  public :: sp_eh_init,sp_eh_step,sp_eh_close,sp_eh_sample_h,sp_eh_checkpoint
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


  ! Collective on the field communicator; each rank supplies a DISTINCT file path.
  ! Same initialized configuration/topology only. This native binary format is
  ! compiler/build-specific; the caller must verify configuration and file hashes.
  ! restore=0 creates a new file; restore=1 loads without advancing the solver.
  ! A failed restore leaves every rank unchanged. Files from a failed save must
  ! be discarded by the caller (there is no multi-file atomic publication here).
  integer(c_int) function sp_eh_checkpoint(path,restore) bind(C)
    character(c_char),intent(in) :: path(*)
    integer(c_int),value :: restore
    character(4096) :: filename
    type(ls_fdtd_eh) :: fe
    type(s_ld_state) :: ld
    real(8),allocatable :: current(:,:,:,:)
    integer :: i,ios,total,step,unit,close_status
    logical :: loading
    sp_eh_checkpoint=1
    ios=0;filename=''
    if(.not.allocated(state))ios=1
    if(restore/=0.and.restore/=1)ios=1
    do i=1,len(filename)
      if(path(i)==c_null_char)exit
      filename(i:i)=path(i)
    enddo
    if(i>len(filename).or.len_trim(filename)==0)ios=1
    call comm_summation(ios,total,nproc_group_global)
    if(total/=0)return
    loading=restore==1
    ! Work on detached allocatable storage; fs and its communicator/cache remain live.
    fe=state%fe;ld=state%ld;current=state%current;step=state%step
    if(loading)then
      open(newunit=unit,file=trim(filename),access='stream',form='unformatted', &
           status='old',action='read',iostat=ios)
    else
      open(newunit=unit,file=trim(filename),access='stream',form='unformatted', &
           status='new',action='write',iostat=ios)
    endif
    if(ios==0)then
      call eh_checkpoint_stream(unit,loading,fe,ld,current,step,ios)
      close(unit,iostat=close_status)
      if(close_status/=0)ios=close_status
    endif
    if(ios/=0)ios=1
    call comm_summation(ios,total,nproc_group_global)
    if(total/=0)return
    if(loading)then
      state%fe=fe;state%ld=ld;state%current=current;state%step=step
    endif
    sp_eh_checkpoint=0
  end function sp_eh_checkpoint

  subroutine eh_checkpoint_stream(unit,loading,fe,ld,current,step,ios)
    integer,intent(in) :: unit
    logical,intent(in) :: loading
    type(ls_fdtd_eh),intent(inout) :: fe
    type(s_ld_state),intent(inout) :: ld
    real(8),allocatable,intent(inout) :: current(:,:,:,:)
    integer,intent(inout) :: step
    integer,intent(out) :: ios
    character(24),parameter :: signature='SALMONPY_EH_CHECKPOINT_1'
    character(24) :: magic
    integer :: topology(18),expected(18),capacity,tail,ld_meta(3),iteration(3)
    integer,allocatable :: points(:,:)
    real(8) :: maxima(2)
    logical :: obs
    ios=0
    expected=[state%fs%lg%is,state%fs%lg%ie,state%fs%mg%is,state%fs%mg%ie, &
              state%fs%mg%is_array,state%fs%mg%ie_array]
    if(loading)then
      read(unit,iostat=ios)magic,topology,capacity,step,ld_meta,obs
      if(ios/=0)return
      if(magic/=signature.or.any(topology/=expected).or.capacity/=state%capacity.or. &
         step<0.or.step>capacity.or.any(ld_meta/=[merge(1,0,ld%enabled),ld%num_ld,ld%max_pole_num_ld]).or. &
         (obs.neqv.fe%flag_obs))then
        ios=1;return
      endif
    else
      write(unit,iostat=ios)signature,expected,state%capacity,step, &
        [merge(1,0,ld%enabled),ld%num_ld,ld%max_pole_num_ld],fe%flag_obs
      if(ios/=0)return
    endif
    points=state%points
    call checkpoint_i2(unit,loading,points,ios)
    if(ios/=0)return
    if(any(points/=state%points))then
      ios=1;return
    endif
    if(loading)then
      read(unit,iostat=ios)iteration,maxima
      if(ios/=0)return
      fe%iter_sta=iteration(1);fe%iter_end=iteration(2);fe%iter_now=iteration(3)
      if(fe%flag_obs)then
        fe%e_max=maxima(1);fe%h_max=maxima(2)
      endif
    else
      ! iter_now and observation maxima are not initialized on every init path.
      iteration=[fe%iter_sta,fe%iter_end,step];maxima=0d0
      if(fe%flag_obs)maxima=[fe%e_max,fe%h_max]
      write(unit,iostat=ios)iteration,maxima
      if(ios/=0)return
    endif
    call checkpoint_r4(unit,loading,current,ios)
    call checkpoint_r1(unit,loading,fe%c_ebloch_x,ios)
    call checkpoint_r1(unit,loading,fe%c_ebloch_y,ios)
    call checkpoint_r1(unit,loading,fe%c_ebloch_z,ios)
    call checkpoint_r1(unit,loading,fe%c_hbloch_x,ios)
    call checkpoint_r1(unit,loading,fe%c_hbloch_y,ios)
    call checkpoint_r1(unit,loading,fe%c_hbloch_z,ios)
    call checkpoint_r1(unit,loading,fe%rep,ios)
    call checkpoint_r1(unit,loading,fe%rmu,ios)
    call checkpoint_r1(unit,loading,fe%sig,ios)
    call checkpoint_r2(unit,loading,fe%coo,ios)
    call checkpoint_i1(unit,loading,fe%iobs_po_pe,ios)
    call checkpoint_i2(unit,loading,fe%iobs_li_pe,ios)
    call checkpoint_i2(unit,loading,fe%iobs_pl_pe,ios)
    call checkpoint_i2(unit,loading,fe%iobs_po_id,ios)
    call checkpoint_i1(unit,loading,fe%iobs_num_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_ex_xy_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_ex_yz_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_ex_xz_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_ey_xy_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_ey_yz_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_ey_xz_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_ez_xy_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_ez_yz_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_ez_xz_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_hx_xy_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_hx_yz_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_hx_xz_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_hy_xy_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_hy_yz_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_hy_xz_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_hz_xy_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_hz_yz_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_hz_xz_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_jx_xy_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_jx_yz_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_jx_xz_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_jy_xy_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_jy_yz_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_jy_xz_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_jz_xy_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_jz_yz_ene,ios)
    call checkpoint_c5(unit,loading,fe%obs_jz_xz_ene,ios)
    call checkpoint_i1(unit,loading,fe%inc_po_pe,ios)
    call checkpoint_i2(unit,loading,fe%inc_li_pe,ios)
    call checkpoint_i2(unit,loading,fe%inc_pl_pe,ios)
    call checkpoint_i2(unit,loading,fe%inc_po_id,ios)
    call checkpoint_r2(unit,loading,fe%gbeam_width_xy1,ios)
    call checkpoint_r2(unit,loading,fe%gbeam_width_yz1,ios)
    call checkpoint_r2(unit,loading,fe%gbeam_width_xz1,ios)
    call checkpoint_r2(unit,loading,fe%gbeam_width_xy2,ios)
    call checkpoint_r2(unit,loading,fe%gbeam_width_yz2,ios)
    call checkpoint_r2(unit,loading,fe%gbeam_width_xz2,ios)
    call checkpoint_r1(unit,loading,fe%gbeam_width_x1,ios)
    call checkpoint_r1(unit,loading,fe%gbeam_width_y1,ios)
    call checkpoint_r1(unit,loading,fe%gbeam_width_z1,ios)
    call checkpoint_r1(unit,loading,fe%gbeam_width_x2,ios)
    call checkpoint_r1(unit,loading,fe%gbeam_width_y2,ios)
    call checkpoint_r1(unit,loading,fe%gbeam_width_z2,ios)
    call checkpoint_r3(unit,loading,fe%ex_y,ios)
    call checkpoint_r3(unit,loading,fe%c1_ex_y,ios)
    call checkpoint_r3(unit,loading,fe%c2_ex_y,ios)
    call checkpoint_r3(unit,loading,fe%ex_z,ios)
    call checkpoint_r3(unit,loading,fe%c1_ex_z,ios)
    call checkpoint_r3(unit,loading,fe%c2_ex_z,ios)
    call checkpoint_r3(unit,loading,fe%ey_z,ios)
    call checkpoint_r3(unit,loading,fe%c1_ey_z,ios)
    call checkpoint_r3(unit,loading,fe%c2_ey_z,ios)
    call checkpoint_r3(unit,loading,fe%ey_x,ios)
    call checkpoint_r3(unit,loading,fe%c1_ey_x,ios)
    call checkpoint_r3(unit,loading,fe%c2_ey_x,ios)
    call checkpoint_r3(unit,loading,fe%ez_x,ios)
    call checkpoint_r3(unit,loading,fe%c1_ez_x,ios)
    call checkpoint_r3(unit,loading,fe%c2_ez_x,ios)
    call checkpoint_r3(unit,loading,fe%ez_y,ios)
    call checkpoint_r3(unit,loading,fe%c1_ez_y,ios)
    call checkpoint_r3(unit,loading,fe%c2_ez_y,ios)
    call checkpoint_r3(unit,loading,fe%hx_y,ios)
    call checkpoint_r3(unit,loading,fe%c1_hx_y,ios)
    call checkpoint_r3(unit,loading,fe%c2_hx_y,ios)
    call checkpoint_r3(unit,loading,fe%hx_z,ios)
    call checkpoint_r3(unit,loading,fe%c1_hx_z,ios)
    call checkpoint_r3(unit,loading,fe%c2_hx_z,ios)
    call checkpoint_r3(unit,loading,fe%hy_z,ios)
    call checkpoint_r3(unit,loading,fe%c1_hy_z,ios)
    call checkpoint_r3(unit,loading,fe%c2_hy_z,ios)
    call checkpoint_r3(unit,loading,fe%hy_x,ios)
    call checkpoint_r3(unit,loading,fe%c1_hy_x,ios)
    call checkpoint_r3(unit,loading,fe%c2_hy_x,ios)
    call checkpoint_r3(unit,loading,fe%hz_x,ios)
    call checkpoint_r3(unit,loading,fe%c1_hz_x,ios)
    call checkpoint_r3(unit,loading,fe%c2_hz_x,ios)
    call checkpoint_r3(unit,loading,fe%hz_y,ios)
    call checkpoint_r3(unit,loading,fe%c1_hz_y,ios)
    call checkpoint_r3(unit,loading,fe%c2_hz_y,ios)
    call checkpoint_r3(unit,loading,fe%ex_s,ios)
    call checkpoint_r3(unit,loading,fe%ey_s,ios)
    call checkpoint_r3(unit,loading,fe%ez_s,ios)
    call checkpoint_r3(unit,loading,fe%hx_s,ios)
    call checkpoint_r3(unit,loading,fe%hy_s,ios)
    call checkpoint_r3(unit,loading,fe%hz_s,ios)
    call checkpoint_r3(unit,loading,fe%c2_jx,ios)
    call checkpoint_r3(unit,loading,fe%c2_jy,ios)
    call checkpoint_r3(unit,loading,fe%c2_jz,ios)
    call checkpoint_r3(unit,loading,fe%rmedia,ios)
    call checkpoint_r1(unit,loading,fe%time_lr,ios)
    call checkpoint_r2(unit,loading,fe%fr_lr,ios)
    call checkpoint_r2(unit,loading,fe%fi_lr,ios)
    call checkpoint_r2(unit,loading,fe%dip_lr,ios)
    call checkpoint_r2(unit,loading,fe%curr_lr,ios)
    call checkpoint_r2(unit,loading,fe%e_lr,ios)
    call checkpoint_r2(unit,loading,fe%er_lr,ios)
    call checkpoint_r2(unit,loading,fe%ei_lr,ios)
    call checkpoint_r1(unit,loading,fe%ase_ene,ios)
    call checkpoint_i2(unit,loading,fe%iase_mask_xy,ios)
    call checkpoint_i2(unit,loading,fe%iase_mask_yz,ios)
    call checkpoint_i2(unit,loading,fe%iase_mask_xz,ios)
    call checkpoint_c3(unit,loading,fe%ase_ex_inc_pa_ene,ios)
    call checkpoint_c3(unit,loading,fe%ase_ey_inc_pa_ene,ios)
    call checkpoint_c3(unit,loading,fe%ase_ez_inc_pa_ene,ios)
    call checkpoint_c3(unit,loading,fe%ase_hx_inc_pa_ene,ios)
    call checkpoint_c3(unit,loading,fe%ase_hy_inc_pa_ene,ios)
    call checkpoint_c3(unit,loading,fe%ase_hz_inc_pa_ene,ios)
    call checkpoint_c3(unit,loading,fe%ase_ex_inc_bt_ene,ios)
    call checkpoint_c3(unit,loading,fe%ase_ey_inc_bt_ene,ios)
    call checkpoint_c3(unit,loading,fe%ase_ez_inc_bt_ene,ios)
    call checkpoint_c3(unit,loading,fe%ase_hx_inc_bt_ene,ios)
    call checkpoint_c3(unit,loading,fe%ase_hy_inc_bt_ene,ios)
    call checkpoint_c3(unit,loading,fe%ase_hz_inc_bt_ene,ios)
    call checkpoint_c5(unit,loading,fe%ase_ex_xy_sca_ene,ios)
    call checkpoint_c5(unit,loading,fe%ase_ex_xz_sca_ene,ios)
    call checkpoint_c5(unit,loading,fe%ase_ey_xy_sca_ene,ios)
    call checkpoint_c5(unit,loading,fe%ase_ey_yz_sca_ene,ios)
    call checkpoint_c5(unit,loading,fe%ase_ez_yz_sca_ene,ios)
    call checkpoint_c5(unit,loading,fe%ase_ez_xz_sca_ene,ios)
    call checkpoint_c5(unit,loading,fe%ase_hx_xy_sca_ene,ios)
    call checkpoint_c5(unit,loading,fe%ase_hx_xz_sca_ene,ios)
    call checkpoint_c5(unit,loading,fe%ase_hy_xy_sca_ene,ios)
    call checkpoint_c5(unit,loading,fe%ase_hy_yz_sca_ene,ios)
    call checkpoint_c5(unit,loading,fe%ase_hz_yz_sca_ene,ios)
    call checkpoint_c5(unit,loading,fe%ase_hz_xz_sca_ene,ios)
    call checkpoint_r1(unit,loading,fe%art_ene,ios)
    call checkpoint_c3(unit,loading,fe%art_ex_inc_bt_ene,ios)
    call checkpoint_c3(unit,loading,fe%art_ey_inc_bt_ene,ios)
    call checkpoint_c3(unit,loading,fe%art_ez_inc_bt_ene,ios)
    call checkpoint_c3(unit,loading,fe%art_hx_inc_bt_ene,ios)
    call checkpoint_c3(unit,loading,fe%art_hy_inc_bt_ene,ios)
    call checkpoint_c3(unit,loading,fe%art_hz_inc_bt_ene,ios)
    call checkpoint_c5(unit,loading,fe%art_ex_tot_bt_ene,ios)
    call checkpoint_c5(unit,loading,fe%art_ey_tot_bt_ene,ios)
    call checkpoint_c5(unit,loading,fe%art_ez_tot_bt_ene,ios)
    call checkpoint_c5(unit,loading,fe%art_hx_tot_bt_ene,ios)
    call checkpoint_c5(unit,loading,fe%art_hy_tot_bt_ene,ios)
    call checkpoint_c5(unit,loading,fe%art_hz_tot_bt_ene,ios)
    call checkpoint_i1(unit,loading,ld%media_ld,ios)
    call checkpoint_r4(unit,loading,ld%rjx_ld,ios)
    call checkpoint_r4(unit,loading,ld%rjx_old_ld,ios)
    call checkpoint_r4(unit,loading,ld%rjy_ld,ios)
    call checkpoint_r4(unit,loading,ld%rjy_old_ld,ios)
    call checkpoint_r4(unit,loading,ld%rjz_ld,ios)
    call checkpoint_r4(unit,loading,ld%rjz_old_ld,ios)
    call checkpoint_r4(unit,loading,ld%px_ld,ios)
    call checkpoint_r4(unit,loading,ld%py_ld,ios)
    call checkpoint_r4(unit,loading,ld%pz_ld,ios)
    call checkpoint_r4(unit,loading,ld%c1_jx_ld,ios)
    call checkpoint_r4(unit,loading,ld%c2_jx_ld,ios)
    call checkpoint_r4(unit,loading,ld%c3_jx_ld,ios)
    call checkpoint_r4(unit,loading,ld%c1_jy_ld,ios)
    call checkpoint_r4(unit,loading,ld%c2_jy_ld,ios)
    call checkpoint_r4(unit,loading,ld%c3_jy_ld,ios)
    call checkpoint_r4(unit,loading,ld%c1_jz_ld,ios)
    call checkpoint_r4(unit,loading,ld%c2_jz_ld,ios)
    call checkpoint_r4(unit,loading,ld%c3_jz_ld,ios)
    call checkpoint_r3(unit,loading,ld%rjx_fdtd_ld,ios)
    call checkpoint_r3(unit,loading,ld%rjy_fdtd_ld,ios)
    call checkpoint_r3(unit,loading,ld%rjz_fdtd_ld,ios)
    call checkpoint_r3(unit,loading,ld%ex_old_ld,ios)
    call checkpoint_r3(unit,loading,ld%ey_old_ld,ios)
    call checkpoint_r3(unit,loading,ld%ez_old_ld,ios)
    if(ios/=0)return
    if(loading)then
      read(unit,iostat=ios)tail
      if(ios/=0)return
      if(tail/=172903)then
        ios=1;return
      endif
      ! Reject trailing payload rather than silently accepting another format.
      read(unit,iostat=ios)tail
      if(ios==0)then
        ios=1
      else if(is_iostat_end(ios))then
        ios=0
      endif
    else
      write(unit,iostat=ios)172903
    endif
  end subroutine eh_checkpoint_stream

  subroutine checkpoint_c3(unit,loading,array,ios)
    integer,intent(in) :: unit
    logical,intent(in) :: loading
    complex(8),allocatable,intent(inout) :: array(:,:,:)
    integer,intent(inout) :: ios
    logical :: present_array
    integer :: lo(3),hi(3)
    if(ios/=0)return
    if(loading)then
      read(unit,iostat=ios)present_array
      if(ios/=0)return
      if(present_array.neqv.allocated(array))then
        ios=1;return
      endif
      if(.not.present_array)return
      read(unit,iostat=ios)lo,hi
      if(ios/=0)return
      if(any(lo/=lbound(array)).or.any(hi/=ubound(array)))then
        ios=1;return
      endif
      read(unit,iostat=ios)array
    else
      write(unit,iostat=ios)allocated(array)
      if(ios/=0.or..not.allocated(array))return
      write(unit,iostat=ios)lbound(array),ubound(array)
      if(ios/=0)return
      write(unit,iostat=ios)array
    endif
  end subroutine checkpoint_c3

  subroutine checkpoint_c5(unit,loading,array,ios)
    integer,intent(in) :: unit
    logical,intent(in) :: loading
    complex(8),allocatable,intent(inout) :: array(:,:,:,:,:)
    integer,intent(inout) :: ios
    logical :: present_array
    integer :: lo(5),hi(5)
    if(ios/=0)return
    if(loading)then
      read(unit,iostat=ios)present_array
      if(ios/=0)return
      if(present_array.neqv.allocated(array))then
        ios=1;return
      endif
      if(.not.present_array)return
      read(unit,iostat=ios)lo,hi
      if(ios/=0)return
      if(any(lo/=lbound(array)).or.any(hi/=ubound(array)))then
        ios=1;return
      endif
      read(unit,iostat=ios)array
    else
      write(unit,iostat=ios)allocated(array)
      if(ios/=0.or..not.allocated(array))return
      write(unit,iostat=ios)lbound(array),ubound(array)
      if(ios/=0)return
      write(unit,iostat=ios)array
    endif
  end subroutine checkpoint_c5

  subroutine checkpoint_i1(unit,loading,array,ios)
    integer,intent(in) :: unit
    logical,intent(in) :: loading
    integer,allocatable,intent(inout) :: array(:)
    integer,intent(inout) :: ios
    logical :: present_array
    integer :: lo(1),hi(1)
    if(ios/=0)return
    if(loading)then
      read(unit,iostat=ios)present_array
      if(ios/=0)return
      if(present_array.neqv.allocated(array))then
        ios=1;return
      endif
      if(.not.present_array)return
      read(unit,iostat=ios)lo,hi
      if(ios/=0)return
      if(any(lo/=lbound(array)).or.any(hi/=ubound(array)))then
        ios=1;return
      endif
      read(unit,iostat=ios)array
    else
      write(unit,iostat=ios)allocated(array)
      if(ios/=0.or..not.allocated(array))return
      write(unit,iostat=ios)lbound(array),ubound(array)
      if(ios/=0)return
      write(unit,iostat=ios)array
    endif
  end subroutine checkpoint_i1

  subroutine checkpoint_i2(unit,loading,array,ios)
    integer,intent(in) :: unit
    logical,intent(in) :: loading
    integer,allocatable,intent(inout) :: array(:,:)
    integer,intent(inout) :: ios
    logical :: present_array
    integer :: lo(2),hi(2)
    if(ios/=0)return
    if(loading)then
      read(unit,iostat=ios)present_array
      if(ios/=0)return
      if(present_array.neqv.allocated(array))then
        ios=1;return
      endif
      if(.not.present_array)return
      read(unit,iostat=ios)lo,hi
      if(ios/=0)return
      if(any(lo/=lbound(array)).or.any(hi/=ubound(array)))then
        ios=1;return
      endif
      read(unit,iostat=ios)array
    else
      write(unit,iostat=ios)allocated(array)
      if(ios/=0.or..not.allocated(array))return
      write(unit,iostat=ios)lbound(array),ubound(array)
      if(ios/=0)return
      write(unit,iostat=ios)array
    endif
  end subroutine checkpoint_i2

  subroutine checkpoint_r1(unit,loading,array,ios)
    integer,intent(in) :: unit
    logical,intent(in) :: loading
    real(8),allocatable,intent(inout) :: array(:)
    integer,intent(inout) :: ios
    logical :: present_array
    integer :: lo(1),hi(1)
    if(ios/=0)return
    if(loading)then
      read(unit,iostat=ios)present_array
      if(ios/=0)return
      if(present_array.neqv.allocated(array))then
        ios=1;return
      endif
      if(.not.present_array)return
      read(unit,iostat=ios)lo,hi
      if(ios/=0)return
      if(any(lo/=lbound(array)).or.any(hi/=ubound(array)))then
        ios=1;return
      endif
      read(unit,iostat=ios)array
    else
      write(unit,iostat=ios)allocated(array)
      if(ios/=0.or..not.allocated(array))return
      write(unit,iostat=ios)lbound(array),ubound(array)
      if(ios/=0)return
      write(unit,iostat=ios)array
    endif
  end subroutine checkpoint_r1

  subroutine checkpoint_r2(unit,loading,array,ios)
    integer,intent(in) :: unit
    logical,intent(in) :: loading
    real(8),allocatable,intent(inout) :: array(:,:)
    integer,intent(inout) :: ios
    logical :: present_array
    integer :: lo(2),hi(2)
    if(ios/=0)return
    if(loading)then
      read(unit,iostat=ios)present_array
      if(ios/=0)return
      if(present_array.neqv.allocated(array))then
        ios=1;return
      endif
      if(.not.present_array)return
      read(unit,iostat=ios)lo,hi
      if(ios/=0)return
      if(any(lo/=lbound(array)).or.any(hi/=ubound(array)))then
        ios=1;return
      endif
      read(unit,iostat=ios)array
    else
      write(unit,iostat=ios)allocated(array)
      if(ios/=0.or..not.allocated(array))return
      write(unit,iostat=ios)lbound(array),ubound(array)
      if(ios/=0)return
      write(unit,iostat=ios)array
    endif
  end subroutine checkpoint_r2

  subroutine checkpoint_r3(unit,loading,array,ios)
    integer,intent(in) :: unit
    logical,intent(in) :: loading
    real(8),allocatable,intent(inout) :: array(:,:,:)
    integer,intent(inout) :: ios
    logical :: present_array
    integer :: lo(3),hi(3)
    if(ios/=0)return
    if(loading)then
      read(unit,iostat=ios)present_array
      if(ios/=0)return
      if(present_array.neqv.allocated(array))then
        ios=1;return
      endif
      if(.not.present_array)return
      read(unit,iostat=ios)lo,hi
      if(ios/=0)return
      if(any(lo/=lbound(array)).or.any(hi/=ubound(array)))then
        ios=1;return
      endif
      read(unit,iostat=ios)array
    else
      write(unit,iostat=ios)allocated(array)
      if(ios/=0.or..not.allocated(array))return
      write(unit,iostat=ios)lbound(array),ubound(array)
      if(ios/=0)return
      write(unit,iostat=ios)array
    endif
  end subroutine checkpoint_r3

  subroutine checkpoint_r4(unit,loading,array,ios)
    integer,intent(in) :: unit
    logical,intent(in) :: loading
    real(8),allocatable,intent(inout) :: array(:,:,:,:)
    integer,intent(inout) :: ios
    logical :: present_array
    integer :: lo(4),hi(4)
    if(ios/=0)return
    if(loading)then
      read(unit,iostat=ios)present_array
      if(ios/=0)return
      if(present_array.neqv.allocated(array))then
        ios=1;return
      endif
      if(.not.present_array)return
      read(unit,iostat=ios)lo,hi
      if(ios/=0)return
      if(any(lo/=lbound(array)).or.any(hi/=ubound(array)))then
        ios=1;return
      endif
      read(unit,iostat=ios)array
    else
      write(unit,iostat=ios)allocated(array)
      if(ios/=0.or..not.allocated(array))return
      write(unit,iostat=ios)lbound(array),ubound(array)
      if(ios/=0)return
      write(unit,iostat=ios)array
    endif
  end subroutine checkpoint_r4

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
