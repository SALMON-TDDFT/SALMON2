! Copyright 2026 SALMONpy contributors. Apache-2.0.
! Stage checks own no numerical state; the active RT context owns its error.
module rt_step_recovery
  use structures
  use iso_c_binding,only:c_int
  use communication,only:comm_get_max,comm_get_globalinfo
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  public::recovery_wave,recovery_scalar,recovery_values,recovery_agree
#ifdef SALMONPY_RECOVERY_TESTING
  integer,save::fault_stage=0,fault_rank=-1
#endif
contains
#ifdef SALMONPY_RECOVERY_TESTING
  integer(c_int) function sp_fast_hse_test_fault(stage,rank) bind(C)
    use iso_c_binding,only:c_int
    integer(c_int),value::stage,rank
    integer::group,id,n
    call comm_get_globalinfo(group,id,n)
    sp_fast_hse_test_fault=1
    if(stage<20.or.stage>33.or.rank<0.or.rank>=n)return
    fault_stage=stage;fault_rank=rank
    sp_fast_hse_test_fault=0
  end function
#endif
  subroutine finish_check(info,stage,bad,status)
    type(s_parallel_info),intent(in)::info
    integer,intent(in)::stage
    integer,intent(inout)::bad,status
#ifdef SALMONPY_RECOVERY_TESTING
    if(fault_stage==stage)then
      if(info%id_rko==fault_rank)bad=stage
      fault_stage=0;fault_rank=-1
    endif
#endif
    call comm_get_max(bad,info%icomm_rko)
    status=max(status,bad)
  end subroutine

  subroutine recovery_agree(info,stage,failed,status)
    type(s_parallel_info),intent(in)::info
    integer,intent(in)::stage
    logical,intent(in)::failed
    integer,intent(inout)::status
    integer::bad
    bad=0
    if(failed)bad=stage
    call finish_check(info,stage,bad,status)
  end subroutine

  subroutine recovery_wave(mg,info,psi,status)
    type(s_rgrid),intent(in)::mg
    type(s_parallel_info),intent(in)::info
    type(s_orbital),intent(in)::psi
    integer,intent(inout)::status
    integer::ix,iy,iz,is,io,ik,im,bad
    bad=0
    if(.not.allocated(psi%zwf))then
      bad=20
    else
      do im=info%im_s,info%im_e
      do ik=info%ik_s,info%ik_e
      do io=info%io_s,info%io_e
      do is=1,size(psi%zwf,4)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        if(.not.ieee_is_finite(real(psi%zwf(ix,iy,iz,is,io,ik,im),8)))bad=20
        if(.not.ieee_is_finite(aimag(psi%zwf(ix,iy,iz,is,io,ik,im))))bad=20
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
    endif
    call finish_check(info,20,bad,status)
  end subroutine

  subroutine recovery_scalar(mg,info,value,stage,status,extra_bad)
    type(s_rgrid),intent(in)::mg
    type(s_parallel_info),intent(in)::info
    type(s_scalar),intent(in)::value
    integer,intent(in)::stage
    integer,intent(inout)::status
    logical,optional,intent(in)::extra_bad
    integer::ix,iy,iz,bad
    bad=0
    if(present(extra_bad))then
      if(extra_bad)bad=stage
    endif
    if(.not.allocated(value%f))then
      bad=stage
    else
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        if(.not.ieee_is_finite(value%f(ix,iy,iz)))bad=stage
      enddo
      enddo
      enddo
    endif
    call finish_check(info,stage,bad,status)
  end subroutine

  subroutine recovery_values(info,value,stage,status)
    type(s_parallel_info),intent(in)::info
    real(8),intent(in)::value(:)
    integer,intent(in)::stage
    integer,intent(inout)::status
    integer::i,bad
    bad=0
    do i=1,size(value)
      if(.not.ieee_is_finite(value(i)))bad=stage
    enddo
    call finish_check(info,stage,bad,status)
  end subroutine
end module
