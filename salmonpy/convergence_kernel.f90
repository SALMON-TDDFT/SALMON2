! SALMONpy convergence diagnostics. All IEEE inquiries are scalar.
module salmonpy_convergence_kernel
  use iso_c_binding
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
contains
  integer(c_int) function sp_convergence_check(n,values,positive) bind(c)
    integer(c_int),value :: n,positive
    real(c_double),intent(in) :: values(*)
    integer :: i
    sp_convergence_check=4
    if(n<1)return
    do i=1,n
      if(.not.ieee_is_finite(values(i)))return
      if(positive/=0.and.values(i)<=0d0)return
    enddo
    sp_convergence_check=0
  end function

  integer(c_int) function sp_convergence_norm(n,values,reference,atol,rtol,stats,accepted) bind(c)
    integer(c_int),value :: n
    real(c_double),intent(in) :: values(*)
    real(c_double),value :: reference,atol,rtol
    real(c_double),intent(out) :: stats(6)
    integer(c_int),intent(out) :: accepted
    real(c_double) :: loss,gain,limit,delta
    integer :: i
    sp_convergence_norm=4;accepted=0;stats=0d0
    if(.not.ieee_is_finite(reference))return
    if(reference<=0d0)return
    if(sp_convergence_check(n,values,1_c_int)/=0)return
    if(.not.ieee_is_finite(atol).or..not.ieee_is_finite(rtol))return
    if(atol<0d0.or.rtol<0d0)return
    loss=0d0;gain=0d0
    do i=1,n
      delta=reference-values(i)
      if(.not.ieee_is_finite(delta))return
      loss=max(loss,delta);gain=max(gain,-delta)
    enddo
    stats(1)=loss;stats(2)=gain;stats(3)=loss/reference;stats(4)=gain/reference
    stats(5)=max(loss,gain);stats(6)=stats(5)/reference
    do i=1,6
      if(.not.ieee_is_finite(stats(i)))return
    enddo
    limit=atol+rtol*reference
    if(.not.ieee_is_finite(limit))return
    if(stats(5)<=limit)accepted=1
    sp_convergence_norm=0
  end function

  integer(c_int) function sp_convergence_compare(n,coarse,fine,atol,rtol,stats,accepted) bind(c)
    integer(c_int),value :: n
    real(c_double),intent(in) :: coarse(*),fine(*)
    real(c_double),value :: atol,rtol
    real(c_double),intent(out) :: stats(2)
    integer(c_int),intent(out) :: accepted
    integer :: i
    real(c_double) :: delta,scale,limit
    sp_convergence_compare=4;stats=0d0;accepted=1;scale=0d0
    if(sp_convergence_check(n,coarse,0_c_int)/=0)return
    if(sp_convergence_check(n,fine,0_c_int)/=0)return
    if(.not.ieee_is_finite(atol).or..not.ieee_is_finite(rtol))return
    if(atol<0d0.or.rtol<0d0)return
    do i=1,n
      delta=abs(coarse(i)-fine(i));limit=atol+rtol*abs(fine(i))
      if(.not.ieee_is_finite(delta).or..not.ieee_is_finite(limit))return
      stats(1)=max(stats(1),delta);scale=max(scale,abs(fine(i)))
      if(delta>limit)accepted=0
    enddo
    if(scale>0d0)then
      stats(2)=stats(1)/scale
      if(.not.ieee_is_finite(stats(2)))stats(2)=huge(1d0)
    else if(stats(1)>0d0)then
      stats(2)=huge(1d0)
    endif
    sp_convergence_compare=0
  end function
  subroutine scaled_square(value,scale,ssq)
    real(c_double),intent(in) :: value
    real(c_double),intent(inout) :: scale,ssq
    real(c_double) :: magnitude
    magnitude=abs(value)
    if(magnitude==0d0)return
    if(scale<magnitude)then
      ssq=1d0+ssq*(scale/magnitude)**2
      scale=magnitude
    else
      ssq=ssq+(magnitude/scale)**2
    endif
  end subroutine

  integer(c_int) function sp_convergence_waveforms(n,coarse,fine,atol,rtol,stats,accepted) bind(c)
    integer(c_int),value :: n
    real(c_double),intent(in) :: coarse(*),fine(*)
    real(c_double),value :: atol,rtol
    real(c_double),intent(out) :: stats(4)
    integer(c_int),intent(out) :: accepted
    real(c_double) :: difference,es,eq,rs,rq,limit
    integer :: i
    sp_convergence_waveforms=4;stats=0d0;accepted=0
    if(sp_convergence_check(n,coarse,0_c_int)/=0)return
    if(sp_convergence_check(n,fine,0_c_int)/=0)return
    if(.not.ieee_is_finite(atol).or..not.ieee_is_finite(rtol))return
    if(atol<0d0.or.rtol<0d0)return
    es=0d0;eq=1d0;rs=0d0;rq=1d0
    do i=1,n
      difference=coarse(i)-fine(i)
      if(.not.ieee_is_finite(difference))return
      call scaled_square(difference,es,eq)
      call scaled_square(fine(i),rs,rq)
      stats(4)=max(stats(4),abs(difference))
    enddo
    stats(1)=es*sqrt(eq/dble(n))
    stats(3)=rs*sqrt(rq/dble(n))
    if(stats(3)>0d0)then
      stats(2)=stats(1)/stats(3)
      if(.not.ieee_is_finite(stats(2)))stats(2)=huge(1d0)
    endif
    do i=1,4
      if(.not.ieee_is_finite(stats(i)))return
    enddo
    limit=atol+rtol*stats(3)
    ! Inclusive boundary, allowing only machine rounding, not a physical margin.
    limit=limit+8d0*epsilon(1d0)*max(stats(1),stats(3))
    if(.not.ieee_is_finite(limit))return
    if(stats(3)>0d0.or.atol>0d0)then
      if(stats(1)<=limit)accepted=1
    endif
    sp_convergence_waveforms=0
  end function
end module
