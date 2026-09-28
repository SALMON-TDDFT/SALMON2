! Shared hybrid settings for conventional and core-partitioned exchange.
module exx_functional
  use salmon_global, only: xc,hse_omega,pbeh_coulomb_radius,rvv10_b,rvv10_c,rvv10_nq,exx_mlwf_radius
  implicit none
  private
  public :: is_hybrid,is_global_hybrid
  public :: exchange_fraction,exchange_screening,lcfo_write_functional,lcfo_check_functional
contains
  pure logical function is_global_hybrid(name)
    implicit none
    character(*),intent(in) :: name
    is_global_hybrid=name=='pbe0'.or.name=='pbeh40'.or.name=='pbeh40_rvv10'
  end function

  pure logical function is_hybrid(name)
    implicit none
    character(*),intent(in) :: name
    is_hybrid=name=='hse06'.or.is_global_hybrid(name)
  end function

  real(8) function exchange_fraction()
    implicit none
    exchange_fraction=.25d0
    if(xc=='pbeh40'.or.xc=='pbeh40_rvv10')exchange_fraction=.4d0
  end function
  real(8) function exchange_screening()
    implicit none
    exchange_screening=hse_omega
    if(is_global_hybrid(xc))exchange_screening=0d0
  end function
  function parameters() result(p)
    implicit none
    real(8) :: p(7)
    p=[exchange_fraction(),exchange_screening(),0d0,0d0,0d0,0d0,exx_mlwf_radius]
    if(is_global_hybrid(xc))p(3)=pbeh_coulomb_radius
    if(xc=='pbeh40_rvv10')p(4:6)=[rvv10_b,rvv10_c,real(rvv10_nq,8)]
  end function
  subroutine lcfo_write_functional(path,run_id,status)
    implicit none
    character(*),intent(in) :: path,run_id
    integer,intent(out) :: status
    integer :: u,ios
    open(newunit=u,file=path,status='replace',action='write',iostat=status)
    if(status/=0)return
    write(u,'(a)',iostat=status)'SLCFO_FUNCTIONAL_V1',trim(run_id),trim(xc)
    if(status==0)write(u,'(7es26.17)',iostat=status)parameters()
    close(u,iostat=ios)
    if(ios/=0)status=ios
  end subroutine
  subroutine lcfo_check_functional(path,run_id,status)
    use ieee_arithmetic, only: ieee_is_finite
    implicit none
    character(*),intent(in) :: path,run_id
    integer,intent(out) :: status
    character(96) :: magic,saved_run,saved_xc
    real(8) :: p(7)
    integer :: u,ios
    logical :: exists
    status=1
    inquire(file=path,exist=exists)
    if(.not.exists)then
      ! Pre-metadata HSE data retain their historical reconstruction route.
      if(.not.is_global_hybrid(xc))status=0
    else
      open(newunit=u,file=path,status='old',action='read',iostat=ios)
      if(ios==0)then
        read(u,'(a)',iostat=ios)magic
        if(ios==0)read(u,'(a)',iostat=ios)saved_run
        if(ios==0)read(u,'(a)',iostat=ios)saved_xc
        if(ios==0)read(u,*,iostat=ios)p
        close(u)
        if(ios==0)then
          if(magic=='SLCFO_FUNCTIONAL_V1'.and.saved_run==run_id.and.saved_xc==xc)then
            if(all(ieee_is_finite(p)))then
              if(all(abs(p-parameters())<=1d-13*max(1d0,abs(parameters()))))status=0
            endif
          endif
        endif
      endif
    endif
    if(status/=0)write(*,'(2a)')'LCFO functional metadata missing or mismatched: ',path
  end subroutine
end module
