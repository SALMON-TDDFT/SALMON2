module dg_wpw_lcfo_volume_operator
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  type,public::s_dg_wpw_lcfo_volume_components
    complex(8),allocatable::fixed_kinetic(:,:),fixed_nonlocal(:,:),fixed_face(:,:)
  end type
  public::rebuild_dg_wpw_ww_volume,rebuild_dg_wpw_wp_pp_volume
contains
  subroutine rebuild_dg_wpw_ww_volume(fixed_kinetic,fixed_nonlocal,fixed_face,potential_volume,hamiltonian,info)
    complex(8),intent(in)::fixed_kinetic(:,:),fixed_nonlocal(:,:),fixed_face(:,:),potential_volume(:,:)
    complex(8),intent(out)::hamiltonian(:,:)
    integer,intent(out)::info
    info=1
    if(any(shape(hamiltonian)/=shape(potential_volume)).or.any(shape(hamiltonian)/=shape(fixed_kinetic)).or.&
       any(shape(hamiltonian)/=shape(fixed_nonlocal)).or.any(shape(hamiltonian)/=shape(fixed_face)).or.&
       .not.finite2(fixed_kinetic).or..not.finite2(fixed_nonlocal).or..not.finite2(fixed_face).or.&
       .not.finite2(potential_volume))return
    hamiltonian=fixed_kinetic+fixed_nonlocal+fixed_face+potential_volume;info=0
  end subroutine
  subroutine rebuild_dg_wpw_wp_pp_volume(fixed_kinetic,fixed_nonlocal,potential_volume,hamiltonian,info)
    complex(8),intent(in)::fixed_kinetic(:),fixed_nonlocal(:),potential_volume(:)
    complex(8),intent(out)::hamiltonian(:)
    integer,intent(out)::info
    info=1
    if(size(hamiltonian)/=size(potential_volume).or.size(hamiltonian)/=size(fixed_kinetic).or.&
       size(hamiltonian)/=size(fixed_nonlocal).or..not.finite1(fixed_kinetic).or.&
       .not.finite1(fixed_nonlocal).or..not.finite1(potential_volume))return
    hamiltonian=fixed_kinetic+fixed_nonlocal+potential_volume;info=0
  end subroutine
  logical function finite1(a)result(ok)
    complex(8),intent(in)::a(:);ok=all(ieee_is_finite(real(a,8))).and.all(ieee_is_finite(aimag(a)))
  end function
  logical function finite2(a)result(ok)
    complex(8),intent(in)::a(:,:);ok=all(ieee_is_finite(real(a,8))).and.all(ieee_is_finite(aimag(a)))
  end function
end module
