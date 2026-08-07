module dg_overlapping_wannier_localization
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use,intrinsic::iso_fortran_env,only:real64
  implicit none
  private
  public::evaluate_dg_periodic_localization
contains
  subroutine evaluate_dg_periodic_localization(values,weights,phases,norm,moment,spread,ok,message)
    complex(real64),intent(in)::values(:,:),phases(:,:)
    real(real64),intent(in)::weights(:)
    real(real64),intent(out)::norm(:),spread
    complex(real64),intent(out)::moment(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(real64)::density
    integer::nwannier,npoint,naxis,wannier,point,axis

    ok=.false.;message='';spread=0d0;norm=0d0;moment=(0d0,0d0)
    nwannier=size(values,1);npoint=size(values,2);naxis=size(phases,1)
    if(nwannier<1.or.npoint<1.or.naxis<1.or.size(weights)/=npoint.or.&
        size(phases,2)/=npoint.or.size(norm)/=nwannier.or.&
        any(shape(moment)/=[naxis,nwannier]))then
      message='periodic localization arrays have inconsistent dimensions';return
    end if
    if(.not.all(ieee_is_finite(weights)).or.&
        .not.all(ieee_is_finite(real(values))).or.&
        .not.all(ieee_is_finite(aimag(values))).or.&
        .not.all(ieee_is_finite(real(phases))).or.&
        .not.all(ieee_is_finite(aimag(phases))))then
      message='periodic localization payload is not finite';return
    end if
    if(any(weights<=0d0))then
      message='periodic localization weights must be positive';return
    end if
    if(maxval(abs(abs(phases)-1d0))>64d0*epsilon(1d0))then
      message='periodic localization phase is not unit modulus';return
    end if
    do point=1,npoint;do wannier=1,nwannier
      density=weights(point)*abs(values(wannier,point))**2
      norm(wannier)=norm(wannier)+density
      do axis=1,naxis
        moment(axis,wannier)=moment(axis,wannier)+density*phases(axis,point)
      end do
    end do;end do
    if(any(norm<=tiny(1d0)))then
      message='periodic localization Wannier norm is zero';return
    end if
    do wannier=1,nwannier;do axis=1,naxis
      spread=spread+max(0d0,1d0-abs(moment(axis,wannier)/norm(wannier))**2)
    end do;end do
    if(.not.ieee_is_finite(spread))then
      message='periodic localization spread is not finite';return
    end if
    ok=.true.
  end subroutine evaluate_dg_periodic_localization
end module dg_overlapping_wannier_localization
