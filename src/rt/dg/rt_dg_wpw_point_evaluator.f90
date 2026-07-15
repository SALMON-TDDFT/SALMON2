module rt_dg_wpw_point_evaluator
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: evaluate_windowed_kg_point
  public :: evaluate_wannier_point

contains

  pure subroutine evaluate_windowed_kg_point(chi,grad_chi,gvec,position,omega_cell,value,grad_value,info)
    real(8),intent(in)::chi,grad_chi(3),gvec(3),position(3),omega_cell
    complex(8),intent(out)::value,grad_value(3)
    integer,intent(out)::info
    complex(8)::phase
    complex(8),parameter::zi=(0d0,1d0)
    value=(0d0,0d0); grad_value=(0d0,0d0); info=0
    if(omega_cell<=0d0 .or. .not.ieee_is_finite(omega_cell) .or. .not.ieee_is_finite(chi) .or. &
       any(.not.ieee_is_finite(grad_chi)) .or. any(.not.ieee_is_finite(gvec)) .or. &
       any(.not.ieee_is_finite(position))) then
      info=1
      return
    end if
    phase=exp(zi*sum(gvec*position))/sqrt(omega_cell)
    value=chi*phase
    grad_value=(cmplx(grad_chi,0d0,kind=8)+zi*gvec*chi)*phase
  end subroutine evaluate_windowed_kg_point

  pure subroutine evaluate_wannier_point(coefficients,basis_values,basis_gradients,values,gradients,info)
    complex(8),intent(in)::coefficients(:,:),basis_values(:),basis_gradients(:,:)
    complex(8),intent(out)::values(:),gradients(:,:)
    integer,intent(out)::info
    integer::ib,iw
    values=(0d0,0d0); gradients=(0d0,0d0); info=0
    if(size(coefficients,1)/=size(basis_values) .or. size(basis_gradients,1)/=3 .or. &
       size(basis_gradients,2)/=size(basis_values) .or. size(coefficients,2)/=size(values) .or. &
       size(gradients,1)/=3 .or. size(gradients,2)/=size(values)) then
      info=1
      return
    end if
    if(.not.all_finite_2(coefficients) .or. .not.all_finite_1(basis_values) .or. &
       .not.all_finite_2(basis_gradients)) then
      info=2
      return
    end if
    do iw=1,size(values)
      do ib=1,size(basis_values)
        values(iw)=values(iw)+coefficients(ib,iw)*basis_values(ib)
        gradients(:,iw)=gradients(:,iw)+coefficients(ib,iw)*basis_gradients(:,ib)
      end do
    end do
  end subroutine evaluate_wannier_point

  pure logical function all_finite_1(values) result(ok)
    complex(8),intent(in)::values(:)
    ok=all(ieee_is_finite(real(values,8))).and.all(ieee_is_finite(aimag(values)))
  end function all_finite_1

  pure logical function all_finite_2(values) result(ok)
    complex(8),intent(in)::values(:,:)
    ok=all(ieee_is_finite(real(values,8))).and.all(ieee_is_finite(aimag(values)))
  end function all_finite_2
end module rt_dg_wpw_point_evaluator
