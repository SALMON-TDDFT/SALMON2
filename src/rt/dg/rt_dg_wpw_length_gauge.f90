module rt_dg_wpw_length_gauge
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: validate_wpw_position_operator, build_wpw_length_gauge_hamiltonian
  public :: evaluate_wpw_polarization, update_wpw_polarization_branch
contains
  subroutine validate_wpw_position_operator(metric,position,tolerance,residual,ok)
    complex(8),intent(in)::metric(:,:),position(:,:)
    real(8),intent(in)::tolerance
    real(8),intent(out)::residual
    logical,intent(out)::ok
    real(8)::scale
    ok=.false.;residual=huge(1d0)
    if(size(metric,1)/=size(metric,2) .or. any(shape(position)/=shape(metric))) return
    if(tolerance<=0d0 .or. .not.all(ieee_is_finite(real(metric))) .or. &
      .not.all(ieee_is_finite(aimag(metric))) .or. .not.all(ieee_is_finite(real(position))) .or. &
      .not.all(ieee_is_finite(aimag(position)))) return
    ! The stored WW/WP/PP volume blocks are covariant matrix elements.  Their
    ! raw Hermiticity residual is also the S-metric residual of S^-1 Z.
    scale=max(1d0,maxval(abs(position)))
    residual=maxval(abs(position-conjg(transpose(position))))/scale
    ok=ieee_is_finite(residual) .and. residual<=tolerance
  end subroutine

  subroutine build_wpw_length_gauge_hamiltonian(h0,position,field,h,ok)
    complex(8),intent(in)::h0(:,:),position(:,:)
    real(8),intent(in)::field
    complex(8),intent(out)::h(:,:)
    logical,intent(out)::ok
    ok=all(shape(h0)==shape(position)) .and. all(shape(h)==shape(h0)) .and. ieee_is_finite(field)
    if(.not.ok)return
    ok=all(ieee_is_finite(real(h0))) .and. all(ieee_is_finite(aimag(h0))) .and. &
      all(ieee_is_finite(real(position))) .and. all(ieee_is_finite(aimag(position)))
    if(ok) h=h0+field*position
  end subroutine

  subroutine evaluate_wpw_polarization(position,coefficients,occupations,pz,ok)
    complex(8),intent(in)::position(:,:),coefficients(:,:)
    real(8),intent(in)::occupations(:)
    real(8),intent(out)::pz
    logical,intent(out)::ok
    integer::state
    complex(8)::value
    pz=0d0;ok=.false.
    if(size(position,1)/=size(position,2) .or. size(coefficients,1)/=size(position,1) .or. &
      size(coefficients,2)/=size(occupations)) return
    if(.not.all(ieee_is_finite(occupations))) return
    do state=1,size(occupations)
      value=dot_product(coefficients(:,state),matmul(position,coefficients(:,state)))
      if(.not.ieee_is_finite(real(value)) .or. .not.ieee_is_finite(aimag(value))) return
      if(abs(aimag(value))>1d-10*max(1d0,abs(real(value)))) return
      pz=pz+occupations(state)*real(value)
    enddo
    ok=ieee_is_finite(pz)
  end subroutine

  subroutine update_wpw_polarization_branch(raw_pz,previous_pz,period,dt,pz,jz,ok)
    real(8),intent(in)::raw_pz,previous_pz,period,dt
    real(8),intent(out)::pz,jz
    logical,intent(out)::ok
    ok=ieee_is_finite(raw_pz).and.ieee_is_finite(previous_pz).and. &
      ieee_is_finite(period).and.ieee_is_finite(dt).and.period>0d0.and.dt>0d0
    if(.not.ok)then;pz=previous_pz;jz=0d0;return;endif
    pz=raw_pz+period*anint((previous_pz-raw_pz)/period)
    jz=(pz-previous_pz)/dt
    ok=ieee_is_finite(pz).and.ieee_is_finite(jz)
  end subroutine
end module
