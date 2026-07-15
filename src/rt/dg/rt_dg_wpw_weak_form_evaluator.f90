module rt_dg_wpw_weak_form_evaluator
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  real(8), parameter :: sipg_penalty_factor = 10.0d0
  public :: wpw_volume_weak_pair
  public :: wpw_sipg_face_pair

contains

  pure subroutine wpw_volume_weak_pair(u, grad_u, v, grad_v, potential, weight, overlap_value, h_value, info)
    complex(8), intent(in) :: u, grad_u(3), v, grad_v(3)
    real(8), intent(in) :: potential, weight
    complex(8), intent(out) :: overlap_value, h_value
    integer, intent(out) :: info

    overlap_value = (0.0d0, 0.0d0)
    h_value = (0.0d0, 0.0d0)
    info = 0
    if (weight < 0.0d0 .or. .not. ieee_is_finite(weight) .or. .not. ieee_is_finite(potential) .or. &
        .not. finite_complex(u) .or. .not. finite_complex(v) .or. &
        .not. all_finite_complex(grad_u) .or. .not. all_finite_complex(grad_v)) then
      info = 1
      return
    end if
    overlap_value = weight * conjg(u) * v
    h_value = weight * (0.5d0 * sum(conjg(grad_u) * grad_v) + potential * conjg(u) * v)
  end subroutine wpw_volume_weak_pair

  pure subroutine wpw_sipg_face_pair(u_minus, u_plus, grad_u_minus, grad_u_plus, &
                                      v_minus, v_plus, grad_v_minus, grad_v_plus, &
                                      normal, h_normal, face_weight, h_face, info)
    complex(8), intent(in) :: u_minus, u_plus, grad_u_minus(3), grad_u_plus(3)
    complex(8), intent(in) :: v_minus, v_plus, grad_v_minus(3), grad_v_plus(3)
    real(8), intent(in) :: normal(3), h_normal, face_weight
    complex(8), intent(out) :: h_face
    integer, intent(out) :: info
    complex(8) :: jump_u, jump_v, average_dn_u, average_dn_v
    real(8) :: normal_norm, penalty_over_h

    h_face = (0.0d0, 0.0d0)
    info = 0
    normal_norm = sqrt(sum(normal * normal))
    if (h_normal <= 0.0d0 .or. face_weight < 0.0d0 .or. &
        .not. ieee_is_finite(h_normal) .or. .not. ieee_is_finite(face_weight) .or. &
        .not. ieee_is_finite(normal_norm) .or. abs(normal_norm - 1.0d0) > 1.0d-12 .or. &
        .not. finite_complex(u_minus) .or. .not. finite_complex(u_plus) .or. &
        .not. finite_complex(v_minus) .or. .not. finite_complex(v_plus) .or. &
        .not. all_finite_complex(grad_u_minus) .or. .not. all_finite_complex(grad_u_plus) .or. &
        .not. all_finite_complex(grad_v_minus) .or. .not. all_finite_complex(grad_v_plus)) then
      info = 1
      return
    end if

    penalty_over_h = sipg_penalty_factor / h_normal
    jump_u = u_minus - u_plus
    jump_v = v_minus - v_plus
    average_dn_u = 0.5d0 * (sum(grad_u_minus * normal) + sum(grad_u_plus * normal))
    average_dn_v = 0.5d0 * (sum(grad_v_minus * normal) + sum(grad_v_plus * normal))
    h_face = face_weight * (-0.5d0 * (conjg(average_dn_u) * jump_v + &
      conjg(jump_u) * average_dn_v) + penalty_over_h * conjg(jump_u) * jump_v)
  end subroutine wpw_sipg_face_pair

  pure logical function finite_complex(value) result(ok)
    complex(8), intent(in) :: value
    ok = ieee_is_finite(real(value,8)) .and. ieee_is_finite(aimag(value))
  end function finite_complex

  pure logical function all_finite_complex(values) result(ok)
    complex(8), intent(in) :: values(:)
    integer :: i
    ok = .true.
    do i = 1, size(values)
      if (.not. finite_complex(values(i))) then
        ok = .false.
        return
      end if
    end do
  end function all_finite_complex
end module rt_dg_wpw_weak_form_evaluator
