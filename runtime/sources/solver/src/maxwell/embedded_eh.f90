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
  use filesystem,only:resolve_path
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use iso_c_binding
  use structures, only: s_fdtd_system
  use fdtd_eh, only: ls_fdtd_eh,eh_init,eh_calc,eh_finalize,s_tensor_step_workspace,eh_prepare_tensor_step, &
    eh_tensor_off,s_tensor_off_workspace
  use lorentz_drude, only: s_ld_state,ld_finalize
  use sendrecv_grid, only: dealloc_cache
  use communication, only: comm_summation,comm_bcast,comm_get_max,comm_begin_children, &
    comm_end_children,comm_resume_children,comm_release_children,comm_child_counts,comm_owner_child_count
  use parallelization, only: setup_parallel,nproc_group_global
  implicit none
  private
  public :: sp_eh_step_indices
  public :: sp_eh_init,sp_eh_step,sp_eh_close,sp_eh_sample_h,sp_eh_checkpoint,sp_eh_box_energy_step
  type :: field_context
    type(s_fdtd_system) :: fs
    type(ls_fdtd_eh) :: fe
    type(s_ld_state) :: ld
    integer,allocatable :: points(:,:)
    real(8),allocatable :: current(:,:,:,:)
    integer :: mpi_child_owner=0
    integer :: step=0,capacity
  end type
  type(field_context),allocatable,save :: state
  type :: restore_candidate
    type(ls_fdtd_eh) :: fe
    type(s_ld_state) :: ld
    real(8),allocatable :: current(:,:,:,:)
    integer :: step
  end type
  type(restore_candidate),allocatable,save :: step_saved
  logical,save :: step_restored=.false.
  type(restore_candidate),allocatable,save :: pending
  logical,save :: used=.false.
#ifdef SALMONPY_RECOVERY_TESTING
  integer,save :: fault_stage=0,fault_rank=-1
#endif
  interface accumulate_finite
    module procedure accumulate_finite_1d,accumulate_finite_2d,accumulate_finite_3d,accumulate_finite_4d
  end interface
contains
  subroutine accumulate_finite_1d(values,valid)
    real(8),intent(in)::values(:)
    logical,intent(inout)::valid
    integer::i1
    do i1=1,size(values,1)
      if(.not.ieee_is_finite(values(i1)))valid=.false.
    enddo
  end subroutine
  subroutine accumulate_finite_2d(values,valid)
    real(8),intent(in)::values(:,:)
    logical,intent(inout)::valid
    integer::i1,i2
    do i2=1,size(values,2)
    do i1=1,size(values,1)
      if(.not.ieee_is_finite(values(i1,i2)))valid=.false.
    enddo
    enddo
  end subroutine
  subroutine accumulate_finite_3d(values,valid)
    real(8),intent(in)::values(:,:,:)
    logical,intent(inout)::valid
    integer::i1,i2,i3
    do i3=1,size(values,3)
    do i2=1,size(values,2)
    do i1=1,size(values,1)
      if(.not.ieee_is_finite(values(i1,i2,i3)))valid=.false.
    enddo
    enddo
    enddo
  end subroutine
  subroutine accumulate_finite_4d(values,valid)
    real(8),intent(in)::values(:,:,:,:)
    logical,intent(inout)::valid
    integer::i1,i2,i3,i4
    do i4=1,size(values,4)
    do i3=1,size(values,3)
    do i2=1,size(values,2)
    do i1=1,size(values,1)
      if(.not.ieee_is_finite(values(i1,i2,i3,i4)))valid=.false.
    enddo
    enddo
    enddo
    enddo
  end subroutine
integer(c_int) function sp_eh_step_transaction_protocol() bind(C)
  sp_eh_step_transaction_protocol=1
end function
integer(c_int) function sp_eh_step_transaction(mode) bind(C)
  use communication,only:comm_get_max
  use salmon_global,only:yn_obs_plane_em,yn_obs_plane_integral_em,base_directory
  integer(c_int),value::mode
  integer::reference,status
  status=0;reference=mode
  call comm_bcast(reference,nproc_group_global)
  if(mode/=reference.or.mode<0.or.mode>3)status=11
  if(.not.allocated(state))status=12
  call comm_get_max(status,nproc_group_global)
  if(status/=0)then
    sp_eh_step_transaction=status;return
  endif
  if(mode==0)then
    if(allocated(step_saved).or.allocated(pending))status=14
    if(trim(base_directory)/='./'.and.trim(base_directory)/='.')status=15
    if(state%fe%flag_lr.or.state%fe%flag_ase.or.state%fe%flag_art.or.state%fe%flag_o_restart)status=15
    if(any(yn_obs_plane_em=='y').or.any(yn_obs_plane_integral_em=='y'))status=15
    if(allocated(state%fe%iobs_num_ene))then
      if(any(state%fe%iobs_num_ene>0))status=15
    endif
  else if(mode==1)then
    if(.not.allocated(step_saved).or.step_restored)status=14
  endif
  call comm_get_max(status,nproc_group_global)
  if(status/=0)then
    sp_eh_step_transaction=status;return
  endif
  select case(mode)
  case(0)
    allocate(step_saved)
    step_saved%fe=state%fe;step_saved%ld=state%ld
    step_saved%current=state%current;step_saved%step=state%step
    step_restored=.false.
  case(1)
    state%fe=step_saved%fe;state%ld=step_saved%ld
    call move_alloc(step_saved%current,state%current)
    state%step=step_saved%step;step_restored=.true.
  case(2)
    if(allocated(step_saved))deallocate(step_saved)
    step_restored=.false.
  end select
  sp_eh_step_transaction=0
end function

integer(c_int) function sp_eh_box_energy_step(box,result) bind(C) result(status)
 use phys_constants,only:cspeed_au
 use salmon_global,only:dt_em,media_type
 integer(c_int),intent(in)::box(3,2)
 real(c_double),intent(inout)::result(19)
 real(8)::local(22),total(22),stats(19)
 integer::region(3,2),m,d,a,b,code
 logical::diagonal,cross_background,conducting
 diagonal=.false.;cross_background=.false.;conducting=.false.
 local=0d0;stats=0d0;code=0
 if(.not.allocated(state).or..not.allocated(step_saved))then
  local(22)=1d0
 else if(step_restored.or.state%step/=step_saved%step+1)then
  local(22)=1d0
 else
  region=box+spread(state%fs%lg%is,2,2)
  if(state%fe%tensor_enabled)then
   ! Pole cross currents are already included in rj_old/rj. Only background
   ! cross coupling needs an additional constitutive boundary-energy term.
   diagonal=.true.
   do m=lbound(state%fe%tensor_corr,4),ubound(state%fe%tensor_corr,4)
    do b=1,3
     do a=1,3
      if(a/=b.and.state%fe%tensor_corr(a,b,0,m)/=0d0)cross_background=.true.
     enddo
    enddo
   enddo
  endif
  if(state%fe%subpixel_enabled)cross_background=state%fe%subpixel_cross(0)
  do m=lbound(state%fe%rep,1),ubound(state%fe%rep,1)
   if(state%fe%rep(m)/=1d0.or.state%fe%rmu(m)/=1d0.or.media_type(m)=='ld-diagonal')diagonal=.true.
   if(state%fe%sig(m)/=0d0)conducting=.true.
  enddo
  do d=1,3
   if(region(d,1)<state%fs%lg%is(d)+1.or.region(d,2)>state%fs%lg%ie(d)-1)local(20)=1d0
   if(region(d,1)>region(d,2))local(20)=1d0
   if(state%fs%a_bc(d,1)=='abc'.and.region(d,1)<=state%fs%lg%is(d)+state%fe%ipml_l)local(20)=1d0
   if(state%fs%a_bc(d,2)=='abc'.and.region(d,2)>=state%fs%lg%ie(d)-state%fe%ipml_l)local(20)=1d0
   if(state%fe%inc_num>0)then
    if(state%fe%inc_dist1/='none'.and.index(state%fe%inc_dist1,achar(iachar('x')+d-1))==0)then
     if(state%fe%inc_po_id(1,d)>=region(d,1).and.state%fe%inc_po_id(1,d)<=region(d,2))local(20)=1d0
    endif
    if(state%fe%inc_num>1.and.state%fe%inc_dist2/='none'.and.index(state%fe%inc_dist2,achar(iachar('x')+d-1))==0)then
     if(state%fe%inc_po_id(2,d)>=region(d,1).and.state%fe%inc_po_id(2,d)<=region(d,2))local(20)=1d0
    endif
   endif
  enddo
  if(conducting)diagonal=.true.
  ! Cross interpolation requires matching halo calls. Agree on admission first.
  if(cross_background)then
   code=int(local(20));call comm_get_max(code,nproc_group_global)
   local(20)=real(code,8)
  endif
  if(local(20)==0d0)then
   code=observe(step_saved%fe,state%fe,state%ld,state%current,region,stats)
   if(code==1)then
    local(20)=1d0
   else if(code/=0)then
    local(21)=1d0
   else
    local(1:19)=stats
   endif
  endif
 endif
 ! One owning-field-group reduction; bad values never enter MPI buffers.
 call comm_summation(local,total,22,nproc_group_global)
 status=14
 if(total(22)>0d0)return
 status=1
 if(total(20)>0d0)return
 status=4
 if(total(21)>0d0)return
 do d=1,19
  if(.not.ieee_is_finite(total(d)))return
 enddo
 result=total(1:19);status=0
contains
 integer function observe(before,after,ld,current,region,stats) result(code)
  use salmonpy_fdtd_energy,only:sp_fdtd_box_energy,sp_fdtd_box_energy_diagonal, &
   sp_fdtd_box_energy_conductive
  type(ls_fdtd_eh),target,intent(in)::before,after
  type(s_ld_state),target,intent(in)::ld
  integer,intent(in)::region(3,2)
  real(8),target,intent(in)::current(3,state%fs%mg%is_array(1):state%fs%mg%ie_array(1), &
   state%fs%mg%is_array(2):state%fs%mg%ie_array(2),state%fs%mg%is_array(3):state%fs%mg%ie_array(3))
  real(8),intent(out)::stats(19)
  type(c_ptr)::oldptr(12),newptr(12),jptr(6),coeffptr(12)
  integer(c_int)::lo(3),hi(3),poles
  real(8)::conductive_stats(18)
  integer::cross_code
  lo=lbound(before%ex_y);hi=ubound(before%ex_y)
    oldptr(1)=c_loc(before%ex_y);newptr(1)=c_loc(after%ex_y)
    oldptr(2)=c_loc(before%ex_z);newptr(2)=c_loc(after%ex_z)
    oldptr(3)=c_loc(before%ey_z);newptr(3)=c_loc(after%ey_z)
    oldptr(4)=c_loc(before%ey_x);newptr(4)=c_loc(after%ey_x)
    oldptr(5)=c_loc(before%ez_x);newptr(5)=c_loc(after%ez_x)
    oldptr(6)=c_loc(before%ez_y);newptr(6)=c_loc(after%ez_y)
    oldptr(7)=c_loc(before%hx_y);newptr(7)=c_loc(after%hx_y)
    oldptr(8)=c_loc(before%hx_z);newptr(8)=c_loc(after%hx_z)
    oldptr(9)=c_loc(before%hy_z);newptr(9)=c_loc(after%hy_z)
    oldptr(10)=c_loc(before%hy_x);newptr(10)=c_loc(after%hy_x)
    oldptr(11)=c_loc(before%hz_x);newptr(11)=c_loc(after%hz_x)
    oldptr(12)=c_loc(before%hz_y);newptr(12)=c_loc(after%hz_y)
  poles=0;jptr=c_null_ptr
  if(ld%enabled)then
   poles=ld%max_pole_num_ld
   jptr=[c_loc(ld%rjx_old_ld),c_loc(ld%rjy_old_ld),c_loc(ld%rjz_old_ld), &
         c_loc(ld%rjx_ld),c_loc(ld%rjy_ld),c_loc(ld%rjz_ld)]
  endif
  stats=0d0
  if(diagonal)then
   coeffptr=c_null_ptr
   coeffptr(1:9)=[c_loc(after%c1_ex_y),c_loc(after%c1_ey_z),c_loc(after%c1_ez_x), &
     c_loc(after%c2_jx),c_loc(after%c2_jy),c_loc(after%c2_jz), &
     c_loc(after%c2_hx_y),c_loc(after%c2_hy_z),c_loc(after%c2_hz_x)]
   if(conducting)then
    if(ld%enabled)coeffptr(10:12)=[c_loc(ld%c2_jx_ld),c_loc(ld%c2_jy_ld),c_loc(ld%c2_jz_ld)]
    conductive_stats=0d0
    code=sp_fdtd_box_energy_conductive(lo,hi,state%fs%mg%is,state%fs%mg%ie,region,state%fs%hgs, &
     dt_em,cspeed_au,poles,oldptr,newptr,jptr,c_loc(current),coeffptr,conductive_stats)
    if(code==0)then
     stats(1:17)=conductive_stats(1:17);stats(19)=conductive_stats(18)
    endif
   else
    code=sp_fdtd_box_energy_diagonal(lo,hi,state%fs%mg%is,state%fs%mg%ie,region,state%fs%hgs, &
     dt_em,cspeed_au,poles,oldptr,newptr,jptr,c_loc(current),coeffptr(1:9),stats(1:17))
   endif
  else
  code=sp_fdtd_box_energy(lo,hi,state%fs%mg%is,state%fs%mg%ie,region,state%fs%hgs, &
   dt_em,cspeed_au,poles,oldptr,newptr,jptr,c_loc(current),stats(1:17))
  endif
  if(cross_background)then
   cross_code=add_cross(before,after,region,stats)
   if(cross_code/=0)code=cross_code
  endif
 end function
 integer function add_cross(before,after,region,stats) result(code)
  use salmonpy_fdtd_energy,only:sp_fdtd_cross_energy
  type(ls_fdtd_eh),intent(in)::before,after
  integer,intent(in)::region(3,2)
  real(8),intent(inout)::stats(19)
  real(8),allocatable::eo(:,:,:,:),en(:,:,:,:),ko(:,:,:,:),kn(:,:,:,:)
  type(s_tensor_off_workspace)::scratch
  real(8)::extra(3)
  integer::lo(3),hi(3),ios,d
  lo=state%fs%mg%is_array;hi=state%fs%mg%ie_array
  allocate(eo(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
   en(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
   ko(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
   kn(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
   scratch%edge(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
   scratch%node(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=ios)
  call comm_get_max(ios,nproc_group_global)
  code=4
  if(ios/=0)return
  eo(1,:,:,:)=before%ex_y+before%ex_z;en(1,:,:,:)=after%ex_y+after%ex_z
  eo(2,:,:,:)=before%ey_z+before%ey_x;en(2,:,:,:)=after%ey_z+after%ey_x
  eo(3,:,:,:)=before%ez_x+before%ez_y;en(3,:,:,:)=after%ez_x+after%ez_y
  call eh_tensor_off(state%fs,after,0,eo,ko,scratch)
  call eh_tensor_off(state%fs,after,0,en,kn,scratch)
  extra=0d0
  code=sp_fdtd_cross_energy(lo,hi,state%fs%mg%is,state%fs%mg%ie,region,state%fs%hgs,eo,en,ko,kn,extra)
  if(code/=0)return
  stats(1)=stats(1)+extra(1);stats(2)=stats(2)+extra(2)
  stats(3)=stats(3)+extra(2)-extra(1)
  stats(4)=stats(4)+extra(1);stats(5)=stats(5)+extra(2)
  stats(6)=stats(6)+extra(2)-extra(1)
  stats(18)=extra(3)
  stats(16)=stats(16)+extra(2)-extra(1)+extra(3)
  stats(17)=stats(17)+extra(2)-extra(1)+extra(3)
  do d=1,19
   if(.not.ieee_is_finite(stats(d)))then
    code=4;return
   endif
  enddo
 end function
end function sp_eh_box_energy_step

integer(c_int) function sp_eh_step_indices(values) bind(C)
 integer(c_int),intent(out)::values(2)
 values=-1
 sp_eh_step_indices=1
 if(.not.allocated(state))return
 values(1)=state%step
 if(allocated(pending))values(2)=pending%step
 sp_eh_step_indices=0
end function

integer(c_int) function sp_eh_checkpoint_protocol() bind(C)
  sp_eh_checkpoint_protocol=1
end function

  ! Known layout/shape statuses follow cold-init state/child cleanup before return.
  integer(c_int) function sp_eh_init_recovery_protocol() bind(C)
    sp_eh_init_recovery_protocol=1
  end function

  integer(c_int) function sp_eh_init(path,npoints,indices) bind(C)
    use inputoutput, only: read_input,s_input_parser
    use math_constants, only: set_math_constants
    use salmon_global, only: theory,nt_em,yn_restart
#ifdef SALMONPY_RECOVERY_TESTING
    use parallelization, only: nproc_id_global
#endif
    character(c_char),intent(in) :: path(*)
    integer(c_int),value :: npoints
    integer(c_int),intent(in) :: indices(3,npoints)
    character(4096) :: filename
    integer :: i,p,lo(3),hi(3),child_status,admission
#ifdef SALMONPY_RECOVERY_TESTING
    integer :: allocation_failure_slot
#endif
    logical :: exists,grid_ready
    type(s_input_parser) :: parser
    sp_eh_init=1
    ! All ranks reject preflight together, without consuming the native lifetime.
    call setup_parallel
    admission=0
    if(used.or.npoints<1) admission=1
    filename=''
    do i=1,len(filename)
      if(path(i)==c_null_char) exit
      filename(i:i)=path(i)
    end do
    if(i>len(filename).or.len_trim(filename)==0) admission=1
    exists=.false.
    if(admission==0) inquire(file=resolve_path(trim(filename)),exist=exists)
    if(.not.exists) admission=1
    call comm_get_max(admission,nproc_group_global)
    if(admission/=0) return
    call set_math_constants
    ! Own parser files locally; a parsed-but-unsupported request is retryable.
    call read_input(trim(filename),parser,admission,native_host=.true.)
    call comm_get_max(admission,nproc_group_global)
    if(admission/=0)then
      sp_eh_init=admission
      return
    endif
    admission=0
    if(theory/='maxwell'.or.yn_restart/='n') admission=1
    call comm_get_max(admission,nproc_group_global)
    if(admission/=0) return
    allocate(state,stat=admission)
    call comm_get_max(admission,nproc_group_global)
    if(admission/=0)then
      if(allocated(state))deallocate(state)
      return
    endif
    call comm_begin_children(state%mpi_child_owner,child_status)
    call comm_get_max(child_status,nproc_group_global)
    if(child_status/=0)then
      call comm_end_children(state%mpi_child_owner)
      call comm_release_children(state%mpi_child_owner)
      deallocate(state);return
    endif
    ! Only the documented early parameter failures are retryable. Unhandled
    ! later numerical/I/O/allocation failures still consume this lifetime.
    used=.true.
#ifdef SALMONPY_RECOVERY_TESTING
    allocation_failure_slot=0
    if(fault_stage>=24.and.fault_stage<=39)then
      if(nproc_id_global==fault_rank)then
        select case(fault_stage)
        case(24,25)
          allocation_failure_slot=merge(2,11,fault_stage==24)
        case(26,27)
          allocation_failure_slot=merge(-2,-20,fault_stage==26)
        case(28,29)
          allocation_failure_slot=merge(102,105,fault_stage==28)
        case(30)
          allocation_failure_slot=202
        case(31)
          allocation_failure_slot=204
        case(32)
          allocation_failure_slot=205
        case(33)
          allocation_failure_slot=206
        case(34)
          allocation_failure_slot=209
        case(35)
          allocation_failure_slot=211
        case(36)
          allocation_failure_slot=214
        case(37)
          allocation_failure_slot=216
        case(38)
          allocation_failure_slot=217
        case(39)
          allocation_failure_slot=218
        end select
      endif
      fault_stage=0;fault_rank=-1
    endif
    call eh_init(state%fs,state%fe,state%ld,admission,grid_ready,excluded_response_points=indices, &
                 allocation_failure_slot=allocation_failure_slot)
#else
    call eh_init(state%fs,state%fe,state%ld,admission,grid_ready,excluded_response_points=indices)
#endif
    call comm_get_max(admission,nproc_group_global)
    if(admission/=0)then
      ! Release communication caches before the child communicator they borrow.
      if(grid_ready)call dealloc_cache(state%fs%srg_ng)
      call comm_end_children(state%mpi_child_owner)
      call comm_release_children(state%mpi_child_owner)
      deallocate(state)
      used=.false.
      sp_eh_init=admission
      return
    endif
    state%capacity=nt_em
    allocate(state%points(3,npoints))
    do p=1,npoints
      state%points(:,p)=state%fs%lg%is+indices(:,p)
    end do
    lo=state%fs%mg%is_array;hi=state%fs%mg%ie_array
    allocate(state%current(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    state%current=0d0
    call comm_end_children(state%mpi_child_owner)
    sp_eh_init=0
  end function

  integer(c_int) function sp_eh_recovery_protocol() bind(C)
    sp_eh_recovery_protocol=1
  end function
#ifdef SALMONPY_RECOVERY_TESTING
  integer(c_int) function sp_eh_test_fault(stage,rank) bind(C)
    integer(c_int),value :: stage,rank
    sp_eh_test_fault=1
    if(stage<20.or.stage>42.or.rank<0)return
    fault_stage=stage;fault_rank=rank;sp_eh_test_fault=0
  end function
#endif

  logical function owns_current_point(p) result(owned)
    implicit none
    integer,intent(in) :: p
    integer :: d,side,index(3)
    owned=.false.
    do d=1,3
      do side=0,1
        index=state%points(:,p);index(d)=index(d)-side
        if(all(index>=state%fs%mg%is).and.all(index<=state%fs%mg%ie))owned=.true.
      enddo
    enddo
  end function

  integer(c_int) function sp_eh_current_mask(npoints,mask) bind(C) result(status)
    implicit none
    integer(c_int),value :: npoints
    integer(c_int),intent(out) :: mask(npoints)
    integer :: p
    status=10;mask=0
    if(.not.allocated(state))return
    if(npoints/=size(state%points,2))return
    do p=1,npoints
      if(owns_current_point(p))mask(p)=1
    enddo
    status=0
  end function

  integer(c_int) function sp_eh_step_owned(npoints,nlocal,J,A,E) bind(C) result(status)
    implicit none
    integer(c_int),value :: npoints,nlocal
    real(c_double),intent(in) :: J(3,nlocal),A(3,npoints)
    real(c_double),intent(out) :: E(3,npoints)
    real(8),allocatable :: local_J(:,:)
    integer :: p,k,error,expected
    ! Validate counts collectively before touching input or advancing any state.
    error=0;expected=0;E=0d0
    if(.not.allocated(state))then
      error=12
    else if(npoints/=size(state%points,2).or.nlocal<0)then
      error=10
    else
      do p=1,npoints
        if(owns_current_point(p))expected=expected+1
      enddo
      if(nlocal/=expected)error=10
    endif
    call comm_get_max(error,nproc_group_global)
    status=error
    if(error/=0)return
    allocate(local_J(3,npoints));local_J=0d0;k=0
    do p=1,npoints
      if(.not.owns_current_point(p))cycle
      k=k+1;local_J(:,p)=J(:,k)
    enddo
    status=step_checked(npoints,local_J,E,A,owned=.true.)
  end function

  integer(c_int) function sp_eh_step(npoints,J,E) bind(C)
    integer(c_int),value :: npoints
    real(c_double),intent(in) :: J(3,npoints)
    real(c_double),intent(out) :: E(3,npoints)
    sp_eh_step=step_checked(npoints,J,E)
  end function

  integer(c_int) function sp_eh_step_checked(npoints,J,A,E) bind(C)
    integer(c_int),value :: npoints
    real(c_double),intent(in) :: J(3,npoints),A(3,npoints)
    real(c_double),intent(out) :: E(3,npoints)
    sp_eh_step_checked=step_checked(npoints,J,E,A)
  end function

  ! Rank-local physical currents already located on E components. Only owned
  ! cells are consumed. No macro-to-edge averaging and no current broadcast.
  integer(c_int) function sp_eh_step_edges(npoints,nedges,indices,J,A,E) bind(C)
    integer(c_int),value::npoints,nedges
    integer(c_int),intent(in)::indices(3,nedges)
    real(c_double),intent(in)::J(3,nedges),A(3,npoints)
    real(c_double),intent(out)::E(3,npoints)
    real(8)::zero(3,npoints)
    zero=0d0
    sp_eh_step_edges=step_checked(npoints,zero,E,A,owned=.true.,edge_indices=indices,edge_current=J)
  end function

  ! Read-only geometry and owned raw E values. No MPI or time advance.
  integer(c_int) function sp_eh_edge_geometry(counts,origins,spacing) bind(C) result(status)
    integer(c_int),intent(out)::counts(3)
    real(c_double),intent(out)::origins(3,3),spacing(3)
    integer::d,k
    status=12
    if(.not.allocated(state))return
    if(allocated(pending))return
    counts=state%fs%lg%ie-state%fs%lg%is+1;spacing=state%fs%hgs
    do d=1,3
    do k=1,3
      origins(k,d)=state%fe%coo(state%fs%lg%is(k),k)
      if(k==d)origins(k,d)=origins(k,d)+.5d0*spacing(k)
    enddo
    enddo
    status=0
  end function

  integer(c_int) function sp_eh_edges_read(nedges,indices,E) bind(C) result(status)
    integer(c_int),value::nedges
    integer(c_int),intent(in)::indices(3,nedges)
    real(c_double),intent(out)::E(3,nedges)
    integer::p,index(3),i,j,k
    status=12
    if(.not.allocated(state))return
    if(allocated(pending))return
    status=10
    do p=1,nedges
      if(any(indices(:,p)<0).or.any(indices(:,p)>state%fs%lg%ie-state%fs%lg%is))return
    enddo
    E=0d0
    do p=1,nedges
      index=state%fs%lg%is+indices(:,p)
      if(.not.(all(index>=state%fs%mg%is).and.all(index<=state%fs%mg%ie)))cycle
      i=index(1);j=index(2);k=index(3)
      E(:,p)=[state%fe%ex_y(i,j,k)+state%fe%ex_z(i,j,k), &
              state%fe%ey_z(i,j,k)+state%fe%ey_x(i,j,k), &
              state%fe%ez_x(i,j,k)+state%fe%ez_y(i,j,k)]
    enddo
    status=0
  end function

  ! Explicit unchecked update: caller owns collective order and valid inputs.
  ! No full-state backup, finite scans, consistency votes or retry contract.
  integer(c_int) function sp_eh_step_performance(npoints,nlocal,owned,J,A_old,A_next,E) bind(C) result(status)
    use salmon_global,only:dt_em
    integer(c_int),value::npoints,nlocal,owned
    real(c_double),intent(in)::J(3,nlocal),A_old(3,npoints)
    real(c_double),intent(out)::A_next(3,npoints),E(3,npoints)
    real(8)::local_E(3,npoints)
    real(8),volatile::delta(3,npoints)
    integer::p,d,side,index(3),i,jj,k,row
    logical::entered
    status=12
    if(.not.allocated(state))return
    status=10
    if(npoints/=size(state%points,2).or.nlocal<0.or.(owned/=0.and.owned/=1))return
    if(owned==0.and.nlocal/=npoints)return
    status=13
    if(state%step>=state%capacity)return
    status=14
    if(allocated(pending))return
    call comm_resume_children(state%mpi_child_owner,status)
    entered=status==0
    if(.not.entered)return
    state%current=0d0;row=0
    do p=1,npoints
      if(owned==1)then
        if(.not.owns_current_point(p))cycle
      endif
      row=row+1
      ! Count checks are local ABI bounds protection, never a collective vote.
      if(row>nlocal)then
        status=10;call comm_end_children(state%mpi_child_owner);return
      endif
      do d=1,3
        do side=0,1
          index=state%points(:,p);index(d)=index(d)-side
          if(all(index>=state%fs%mg%is).and.all(index<=state%fs%mg%ie)) &
            state%current(d,index(1),index(2),index(3))= &
            state%current(d,index(1),index(2),index(3))+.5d0*J(d,row)
        enddo
      enddo
    enddo
    state%step=state%step+1
    state%fe%iter_sta=state%step;state%fe%iter_end=state%step
    call eh_calc(state%fs,state%fe,state%ld,state%current)
    local_E=0d0
    do p=1,npoints
      index=state%points(:,p)
      if(.not.(all(index>=state%fs%mg%is).and.all(index<=state%fs%mg%ie)))cycle
      i=index(1);jj=index(2);k=index(3)
      local_E(1,p)=.5d0*(state%fe%ex_y(i,jj,k)+state%fe%ex_z(i,jj,k) &
        +state%fe%ex_y(i-1,jj,k)+state%fe%ex_z(i-1,jj,k))
      local_E(2,p)=.5d0*(state%fe%ey_z(i,jj,k)+state%fe%ey_x(i,jj,k) &
        +state%fe%ey_z(i,jj-1,k)+state%fe%ey_x(i,jj-1,k))
      local_E(3,p)=.5d0*(state%fe%ez_x(i,jj,k)+state%fe%ez_y(i,jj,k) &
        +state%fe%ez_x(i,jj,k-1)+state%fe%ez_y(i,jj,k-1))
    enddo
    ! Physical sampling reduction remains. Volatile intermediate preserves the
    ! existing separate binary64 multiply/subtract for the A integration.
    call comm_summation(local_E,E,3*npoints,nproc_group_global)
    delta=dt_em*E;A_next=A_old-delta
    call comm_end_children(state%mpi_child_owner)
    status=0
  end function sp_eh_step_performance

  integer function step_checked(npoints,J,E,A,owned,edge_indices,edge_current) result(status)
    integer,intent(in)::npoints
    real(8),intent(in)::J(3,npoints)
    real(8),intent(in),optional::A(3,npoints)
    logical,intent(in),optional :: owned
    integer,intent(in),optional::edge_indices(:,:)
    real(8),intent(in),optional::edge_current(:,:)
    real(8),intent(out)::E(3,npoints)
    logical::entered
    status=0;E=0d0
    if(.not.allocated(state))status=12
    call comm_get_max(status,nproc_group_global)
    if(status/=0)return
    call comm_resume_children(state%mpi_child_owner,status)
    entered=status==0
    call comm_get_max(status,nproc_group_global)
    if(status==0)then
      status=step_checked_impl(npoints,J,E,A,owned,edge_indices,edge_current)
    else
      status=16
    endif
    if(entered)call comm_end_children(state%mpi_child_owner)
  end function

  integer function step_checked_impl(npoints,J,E,A,owned,edge_indices,edge_current) result(status)
    use salmon_global,only:dt_em,obs_samp_em
#ifdef SALMONPY_RECOVERY_TESTING
    use parallelization,only:nproc_id_global
#endif
    use communication,only:comm_get_max
    integer,intent(in) :: npoints
    real(8),intent(in) :: J(3,npoints)
    real(8),intent(in),optional :: A(3,npoints)
    logical,intent(in),optional :: owned
    integer,intent(in),optional::edge_indices(:,:)
    real(8),intent(in),optional::edge_current(:,:)
    real(8),intent(out) :: E(3,npoints)
    type(restore_candidate) :: before
    type(s_tensor_step_workspace),target :: tensor_work
#ifdef SALMONPY_RECOVERY_TESTING
    integer::workspace_failure
#endif
    real(8) :: local_E(3,npoints),reference(3,npoints)
    integer :: p,d,side,index(3),i,jj,k,error,mode
    logical :: finite_input,local_current
    status=10;E=0d0;error=0
    if(.not.allocated(state))error=12
    call comm_get_max(error,nproc_group_global)
    if(error/=0)then
      status=error;return
    endif
    local_current=.false.
    if(present(owned))local_current=owned
    mode=merge(1,0,present(A))+merge(2,0,local_current)+merge(4,0,present(edge_indices))
    call comm_bcast(mode,nproc_group_global)
    if(mode/=merge(1,0,present(A))+merge(2,0,local_current)+merge(4,0,present(edge_indices)))error=11
    call comm_get_max(error,nproc_group_global)
    if(error/=0)then
      status=error;return
    endif
    if(npoints/=size(state%points,2))error=10
    if(state%step>=state%capacity)error=13
    if(allocated(pending))error=14
    if(state%fe%flag_lr.or.state%fe%flag_ase.or.state%fe%flag_art)error=15
    if(allocated(state%fe%iobs_num_ene))then
      if(any(state%fe%iobs_num_ene>0))error=15
    endif
    if(present(edge_indices).neqv.present(edge_current))error=10
    if(present(edge_indices).and.present(edge_current))then
      if(size(edge_indices,1)/=3.or.size(edge_current,1)/=3.or. &
         size(edge_indices,2)/=size(edge_current,2))error=10
      do p=1,size(edge_indices,2)
        if(any(edge_indices(:,p)<0).or.any(edge_indices(:,p)>state%fs%lg%ie-state%fs%lg%is))error=10
      enddo
    endif
    finite_input=.true.
    if(present(edge_current))call accumulate_finite(edge_current,finite_input)
    call accumulate_finite(J,finite_input)
    if(present(A))call accumulate_finite(A,finite_input)
    if(.not.finite_input)error=10
    call comm_get_max(error,nproc_group_global)
    if(error/=0)then
      status=error;return
    endif
    if(.not.local_current)then
      reference=J
      call comm_bcast(reference,nproc_group_global)
      if(any(reference/=J))error=11
    endif
    if(present(A))then
      reference=A
      call comm_bcast(reference,nproc_group_global)
      if(any(reference/=A))error=11
    endif
    if(error==0)then
#ifdef SALMONPY_RECOVERY_TESTING
      workspace_failure=0
      if(fault_stage>=40.and.fault_stage<=42)then
        if(nproc_id_global==fault_rank)then
          select case(fault_stage)
          case(40)
            workspace_failure=1
          case(41)
            workspace_failure=7
          case(42)
            workspace_failure=13
          end select
        endif
        fault_stage=0;fault_rank=-1
      endif
      call eh_prepare_tensor_step(state%fs,state%fe,state%ld%max_pole_num_ld,tensor_work,error,workspace_failure)
#else
      call eh_prepare_tensor_step(state%fs,state%fe,state%ld%max_pole_num_ld,tensor_work,error)
#endif
    endif
    call comm_get_max(error,nproc_group_global)
    if(error/=0)then
      status=error;return
    endif
    before%fe=state%fe;before%ld=state%ld
    before%current=state%current;before%step=state%step
    state%current=0d0
    ! Preserve SALMON's S-transpose edge weighting and E -> LD -> H order.
    do p=1,npoints
      do d=1,3
        do side=0,1
          index=state%points(:,p);index(d)=index(d)-side
          if(all(index>=state%fs%mg%is).and.all(index<=state%fs%mg%ie)) &
            state%current(d,index(1),index(2),index(3))= &
            state%current(d,index(1),index(2),index(3))+.5d0*J(d,p)
        enddo
      enddo
    enddo
    if(present(edge_indices).and.present(edge_current))then
      do p=1,size(edge_indices,2)
        index=state%fs%lg%is+edge_indices(:,p)
        if(.not.(all(index>=state%fs%mg%is).and.all(index<=state%fs%mg%ie)))cycle
        state%current(:,index(1),index(2),index(3))= &
          state%current(:,index(1),index(2),index(3))+edge_current(:,p)
      enddo
    endif
    state%step=state%step+1
    state%fe%iter_sta=state%step;state%fe%iter_end=state%step
    call eh_calc(state%fs,state%fe,state%ld,state%current,guard,tensor_workspace=tensor_work)
    if(error/=0)then
      state%fe=before%fe;state%ld=before%ld
      call move_alloc(before%current,state%current)
      state%step=before%step;E=0d0;status=error;return
    endif
    status=0
  contains
    logical function guard(stage)
      use parallelization,only:nproc_id_global
      integer,intent(in) :: stage
      logical :: valid
      real(8) :: H(3,npoints),ev(3),hv(3)
      integer :: h_status,ix,iy,iz
      valid=.true.
      select case(stage)
      case(20)
        if(state%fe%tensor_failed)valid=.false.
        call accumulate_finite(state%fe%ex_y,valid)
        call accumulate_finite(state%fe%ex_z,valid)
        call accumulate_finite(state%fe%ey_x,valid)
        call accumulate_finite(state%fe%ey_z,valid)
        call accumulate_finite(state%fe%ez_x,valid)
        call accumulate_finite(state%fe%ez_y,valid)
      case(21)
        if(state%ld%enabled)then
          call accumulate_finite(state%ld%rjx_ld,valid)
          call accumulate_finite(state%ld%rjy_ld,valid)
          call accumulate_finite(state%ld%rjz_ld,valid)
          call accumulate_finite(state%ld%rjx_old_ld,valid)
          call accumulate_finite(state%ld%rjy_old_ld,valid)
          call accumulate_finite(state%ld%rjz_old_ld,valid)
          call accumulate_finite(state%ld%px_ld,valid)
          call accumulate_finite(state%ld%py_ld,valid)
          call accumulate_finite(state%ld%pz_ld,valid)
          call accumulate_finite(state%ld%rjx_fdtd_ld,valid)
          call accumulate_finite(state%ld%rjy_fdtd_ld,valid)
          call accumulate_finite(state%ld%rjz_fdtd_ld,valid)
          call accumulate_finite(state%ld%ex_old_ld,valid)
          call accumulate_finite(state%ld%ey_old_ld,valid)
          call accumulate_finite(state%ld%ez_old_ld,valid)
        endif
      case(22)
        call accumulate_finite(state%fe%hx_y,valid)
        call accumulate_finite(state%fe%hx_z,valid)
        call accumulate_finite(state%fe%hy_x,valid)
        call accumulate_finite(state%fe%hy_z,valid)
        call accumulate_finite(state%fe%hz_x,valid)
        call accumulate_finite(state%fe%hz_y,valid)
      case(23)
        local_E=0d0
        do p=1,npoints
          index=state%points(:,p)
          if(.not.(all(index>=state%fs%mg%is).and.all(index<=state%fs%mg%ie)))cycle
          i=index(1);jj=index(2);k=index(3)
          local_E(1,p)=.5d0*(state%fe%ex_y(i,jj,k)+state%fe%ex_z(i,jj,k) &
            +state%fe%ex_y(i-1,jj,k)+state%fe%ex_z(i-1,jj,k))
          local_E(2,p)=.5d0*(state%fe%ey_z(i,jj,k)+state%fe%ey_x(i,jj,k) &
            +state%fe%ey_z(i,jj-1,k)+state%fe%ey_x(i,jj-1,k))
          local_E(3,p)=.5d0*(state%fe%ez_x(i,jj,k)+state%fe%ez_y(i,jj,k) &
            +state%fe%ez_x(i,jj,k-1)+state%fe%ez_y(i,jj,k-1))
        enddo
        call comm_summation(local_E,E,3*npoints,nproc_group_global)
        call accumulate_finite(E,valid)
        if(present(A))call accumulate_finite(A-dt_em*E,valid)
        h_status=sp_eh_sample_h(npoints,H)
        valid=valid.and.h_status==0
        call accumulate_finite(H,valid)
        ! Verify the ordinary point-output/maxima arithmetic before any append.
        ! Keep the actual output path and its arithmetic unchanged.
        if(state%fe%flag_obs.and.mod(state%step,obs_samp_em)==0)then
          do iz=state%fs%mg%is(3),state%fs%mg%ie(3)
          do iy=state%fs%mg%is(2),state%fs%mg%ie(2)
          do ix=state%fs%mg%is(1),state%fs%mg%ie(1)
            ev=[state%fe%ex_y(ix,iy,iz)+state%fe%ex_z(ix,iy,iz), &
                state%fe%ey_z(ix,iy,iz)+state%fe%ey_x(ix,iy,iz), &
                state%fe%ez_x(ix,iy,iz)+state%fe%ez_y(ix,iy,iz)]
            hv=[(state%fe%hx_s(ix,iy,iz)+(state%fe%hx_y(ix,iy,iz)+state%fe%hx_z(ix,iy,iz)))/2d0, &
                (state%fe%hy_s(ix,iy,iz)+(state%fe%hy_z(ix,iy,iz)+state%fe%hy_x(ix,iy,iz)))/2d0, &
                (state%fe%hz_s(ix,iy,iz)+(state%fe%hz_x(ix,iy,iz)+state%fe%hz_y(ix,iy,iz)))/2d0]
            call accumulate_finite(ev*state%fe%uVperm_from_au,valid)
            call accumulate_finite(hv*state%fe%uAperm_from_au,valid)
            valid=valid.and.ieee_is_finite(sqrt(sum(ev**2))*state%fe%uVperm_from_au) &
              .and.ieee_is_finite(sqrt(sum(hv**2))*state%fe%uAperm_from_au)
          enddo
          enddo
          enddo
        endif
      end select
#ifdef SALMONPY_RECOVERY_TESTING
      if(fault_stage==stage)then
        if(nproc_id_global==fault_rank)valid=.false.
        fault_stage=0;fault_rank=-1
      endif
#endif
      if(.not.valid)error=stage
      call comm_get_max(error,nproc_group_global)
      guard=error==0
    end function
  end function

  ! Read raw H at (step+1/2)*dt, centered spatially on each macro cell.
  ! Each rank contributes only owned Yee faces; no diagonal halo assumption.
  ! Read-only, rank-local plane tiles. No MPI, files, timers or state changes.
  ! Query with nvalues=0; values then has no referenced elements.
  integer(c_int) function sp_eh_plane_read(axis,offset,nvalues,info,values) bind(C)
    use salmon_global,only:obs_samp_em,nt_em
    use parallelization,only:nproc_id_global
    integer(c_int),value::axis,offset,nvalues
    integer(c_int),intent(out)::info(21)
    real(c_double),intent(out)::values(*)
    integer::d,idx,tangent(2),j,a,b,q,k,pos,ijk(3),i1,i2,i3,n1,n2
    real(8)::v(6)
    sp_eh_plane_read=1;info=0
    if(.not.allocated(state))return
    if(allocated(pending))return
    d=axis
    if(d<1.or.d>3)return
    idx=state%fs%lg%is(d)+offset
    if(idx<=state%fs%lg%is(d).or.idx>state%fs%lg%ie(d))return
    if(.not.state%fe%flag_save)return
    if(mod(state%step,obs_samp_em)/=0)then
      sp_eh_plane_read=2;return
    endif
    info=[16909060,1,d,idx,idx-1,state%fs%lg%is,state%fs%lg%ie, &
          state%fs%mg%is,state%fs%mg%ie,obs_samp_em,nt_em,nproc_id_global,state%step]
    if(nvalues==0)then
      sp_eh_plane_read=0;return
    endif
    j=0
    do a=1,3
      if(a==d)cycle
      j=j+1;tangent(j)=a
    enddo
    n1=state%fs%mg%ie(tangent(1))-state%fs%mg%is(tangent(1))+1
    n2=state%fs%mg%ie(tangent(2))-state%fs%mg%is(tangent(2))+1
    if(nvalues/=12*n1*n2)return
    pos=0
    do q=0,1
      ijk(d)=idx-q
      do b=1,n2
        ijk(tangent(2))=state%fs%mg%is(tangent(2))+b-1
        do a=1,n1
          ijk(tangent(1))=state%fs%mg%is(tangent(1))+a-1
          v=0d0
          if(idx-q>=state%fs%mg%is(d).and.idx-q<=state%fs%mg%ie(d))then
            i1=ijk(1);i2=ijk(2);i3=ijk(3)
            if(state%step==0)then
              v=[state%fe%ex_y(i1,i2,i3)+state%fe%ex_z(i1,i2,i3), &
                 state%fe%ey_z(i1,i2,i3)+state%fe%ey_x(i1,i2,i3), &
                 state%fe%ez_x(i1,i2,i3)+state%fe%ez_y(i1,i2,i3), &
                 state%fe%hx_y(i1,i2,i3)+state%fe%hx_z(i1,i2,i3), &
                 state%fe%hy_z(i1,i2,i3)+state%fe%hy_x(i1,i2,i3), &
                 state%fe%hz_x(i1,i2,i3)+state%fe%hz_y(i1,i2,i3)]
            else
              v=[state%fe%ex_s(i1,i2,i3),state%fe%ey_s(i1,i2,i3),state%fe%ez_s(i1,i2,i3), &
                 state%fe%hx_s(i1,i2,i3),state%fe%hy_s(i1,i2,i3),state%fe%hz_s(i1,i2,i3)]
            endif
          endif
          do k=1,6
            if(.not.ieee_is_finite(v(k)))then
              sp_eh_plane_read=3;return
            endif
            pos=pos+1;values(pos)=v(k)
          enddo
        enddo
      enddo
    enddo
    sp_eh_plane_read=0
  end function sp_eh_plane_read

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
    integer :: i,ios,total,step,unit,close_status,mode
    logical :: loading
    sp_eh_checkpoint=1
    ios=0;filename=''
    if(.not.allocated(state))ios=1
    if(restore<0.or.restore>4)ios=1
    if(allocated(step_saved).and.restore/=0)ios=14
    mode=restore
    call comm_bcast(mode,nproc_group_global)
    if(mode/=restore)ios=1
    if(restore==3.and..not.allocated(pending))ios=1
    call comm_summation(ios,total,nproc_group_global)
    if(total/=0)return
    if(restore==3)then
      if(.not.allocated(pending))return
      state%fe=pending%fe;state%ld=pending%ld
      call move_alloc(pending%current,state%current)
      state%step=pending%step
      deallocate(pending);sp_eh_checkpoint=0;return
    endif
    if(restore==4)then
      if(allocated(pending))deallocate(pending)
      sp_eh_checkpoint=0;return
    endif
    if(restore==2.and.allocated(pending))ios=1
    do i=1,len(filename)
      if(path(i)==c_null_char)exit
      filename(i:i)=path(i)
    enddo
    if(i>len(filename).or.len_trim(filename)==0)ios=1
    call comm_summation(ios,total,nproc_group_global)
    if(total/=0)return
    loading=restore==1.or.restore==2
    ! Work on detached allocatable storage; fs and its communicator/cache remain live.
    fe=state%fe;ld=state%ld;current=state%current;step=state%step
    if(loading)then
      open(newunit=unit,file=resolve_path(trim(filename)),access='stream',form='unformatted', &
           status='old',action='read',iostat=ios)
    else
      open(newunit=unit,file=resolve_path(trim(filename)),access='stream',form='unformatted', &
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
    if(restore==2)then
      allocate(pending)
      pending%fe=fe;pending%ld=ld
      call move_alloc(current,pending%current)
      pending%step=step
    else if(loading)then
      state%fe=fe;state%ld=ld;state%current=current;state%step=step
    endif
    sp_eh_checkpoint=0
  end function sp_eh_checkpoint

  subroutine eh_checkpoint_stream(unit,loading,fe,ld,current,step,ios)
    use salmon_global,only:dt_em
    integer,intent(in) :: unit
    logical,intent(in) :: loading
    type(ls_fdtd_eh),intent(inout) :: fe
    type(s_ld_state),intent(inout) :: ld
    real(8),allocatable,intent(inout) :: current(:,:,:,:)
    integer,intent(inout) :: step
    integer,intent(out) :: ios
    character(24),parameter :: signature='SALMONPY_EH_CHECKPOINT_1'
    character(24) :: magic,expected_magic
    real(8),allocatable::tensor_context(:,:,:,:)
    real(8)::tensor_options(6)
    real(8),allocatable::tensor_rates(:,:,:)
    integer::q,axis
    integer :: topology(18),expected(18),capacity,tail,ld_meta(3),iteration(3)
    integer,allocatable :: points(:,:)
    real(8) :: maxima(2)
    logical :: obs
    ios=0
    expected_magic=signature
    if(fe%tensor_enabled)expected_magic='SALMONPY_EH_TENSOR_CP_1'
    if(fe%subpixel_enabled)expected_magic='SALMONPY_EH_SUBPIXEL_CP1'
    expected=[state%fs%lg%is,state%fs%lg%ie,state%fs%mg%is,state%fs%mg%ie, &
              state%fs%mg%is_array,state%fs%mg%ie_array]
    if(loading)then
      read(unit,iostat=ios)magic,topology,capacity,step,ld_meta,obs
      if(ios/=0)return
      if(magic/=expected_magic.or.any(topology/=expected).or.capacity/=state%capacity.or. &
         step<0.or.step>capacity.or.any(ld_meta/=[merge(1,0,ld%enabled),ld%num_ld,ld%max_pole_num_ld]).or. &
         (obs.neqv.fe%flag_obs))then
        ios=1;return
      endif
    else
      write(unit,iostat=ios)expected_magic,expected,state%capacity,step, &
        [merge(1,0,ld%enabled),ld%num_ld,ld%max_pole_num_ld],fe%flag_obs
      if(ios/=0)return
    endif
    if(fe%tensor_enabled)then
      tensor_context=fe%tensor_raw
      call checkpoint_r4(unit,loading,tensor_context,ios)
      if(ios/=0)return
      if(any(tensor_context/=fe%tensor_raw))then
        ios=1;return
      endif
      tensor_context=real(fe%tensor_owner,8)
      call checkpoint_r4(unit,loading,tensor_context,ios)
      if(ios/=0)return
      if(any(tensor_context/=real(fe%tensor_owner,8)))then
        ios=1;return
      endif
      do q=0,ubound(fe%tensor_weight,5)
        tensor_context=fe%tensor_weight(:,:,:,:,q)
        call checkpoint_r4(unit,loading,tensor_context,ios)
        if(ios/=0)return
        if(any(tensor_context/=fe%tensor_weight(:,:,:,:,q)))then
          ios=1;return
        endif
      enddo
      if(fe%subpixel_enabled)then
        do q=0,ubound(fe%subpixel_corr,6)
          do axis=1,3
            tensor_context=fe%subpixel_corr(:,axis,:,:,:,q)
            call checkpoint_r4(unit,loading,tensor_context,ios)
            if(ios/=0)return
            if(any(tensor_context/=fe%subpixel_corr(:,axis,:,:,:,q)))then
              ios=1;return
            endif
          enddo
        enddo
      endif
      tensor_rates=fe%tensor_rates
      call checkpoint_r3(unit,loading,tensor_rates,ios)
      if(ios/=0)return
      if(any(tensor_rates/=fe%tensor_rates))then
        ios=1;return
      endif
      tensor_options=[fe%tensor_rtol,real(fe%tensor_maxiter,8),dt_em,state%fs%hgs]
      if(loading)then
        read(unit,iostat=ios)tensor_options
        if(ios/=0)return
        if(any(tensor_options/=[fe%tensor_rtol,real(fe%tensor_maxiter,8),dt_em,state%fs%hgs]))then
          ios=1;return
        endif
      else
        write(unit,iostat=ios)tensor_options
        if(ios/=0)return
      endif
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

  integer(c_int) function sp_eh_mpi_resources(values) bind(C)
    integer(c_int),intent(out)::values(3)
    values=0
    call comm_child_counts(values(1),values(2))
    if(allocated(state))call comm_owner_child_count(state%mpi_child_owner,values(3))
    sp_eh_mpi_resources=0
  end function
  integer(c_int) function sp_eh_resources(values) bind(C)
    integer(c_int),intent(out)::values(3)
    values=0
    if(allocated(state))values(1)=1
    if(allocated(step_saved))values(2)=1
    if(allocated(pending))values(3)=1
    sp_eh_resources=0
  end function
  integer(c_int) function sp_eh_close() bind(C)
    integer::status
    logical::entered
    sp_eh_close=1
    if(allocated(state))then
      call comm_resume_children(state%mpi_child_owner,status)
      entered=status==0
      call comm_get_max(status,nproc_group_global)
      if(status/=0)then
        if(entered)call comm_end_children(state%mpi_child_owner)
        return
      endif
    endif
    if(allocated(pending))deallocate(pending)
    if(allocated(step_saved))deallocate(step_saved)
    step_restored=.false.
    if(.not.allocated(state))return
    if(state%step==state%capacity) call eh_finalize(state%fs,state%fe,state%ld)
    call ld_finalize(state%ld)
    call dealloc_cache(state%fs%srg_ng)
    call comm_end_children(state%mpi_child_owner)
    call comm_release_children(state%mpi_child_owner)
    deallocate(state)
    sp_eh_close=0
  end function
end module
