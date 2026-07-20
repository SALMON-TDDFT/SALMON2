module dg_wpw_core_w_provider
  use rt_dg_wpw_volume_halo_provider,only:s_dg_wpw_volume_halo_state,read_dg_wpw_volume_halo
  implicit none
  private
  public::evaluate_dg_wpw_core_w_support,reconstruct_dg_wpw_core_w_support
contains
  subroutine evaluate_dg_wpw_core_w_support(owned_ids,owned_values,owned_gradients,support_ids,halos,&
      grid,epoch,values,gradients,info,zero_outside_halo)
    integer,intent(in)::owned_ids(:),support_ids(:),grid(3),epoch
    complex(8),intent(in)::owned_values(:),owned_gradients(:,:)
    type(s_dg_wpw_volume_halo_state),intent(in)::halos(:)
    complex(8),intent(out)::values(:),gradients(:,:)
    integer,intent(out)::info
    logical,intent(in),optional::zero_outside_halo
    complex(8),allocatable::candidate_values(:),candidate_gradients(:,:)
    complex(8)::value,gradient(3)
    integer::i,j,owned_location(1),read_info,nsource
    logical::covered,allow_zero

    values=(0d0,0d0);gradients=(0d0,0d0);info=1
    if(epoch<=0.or..not.strictly_increasing(owned_ids).or..not.strictly_increasing(support_ids).or.&
       size(owned_values)/=size(owned_ids).or.any(shape(owned_gradients)/=[3,size(owned_ids)]).or.&
       size(values)/=size(support_ids).or.any(shape(gradients)/=[3,size(support_ids)]))return
    allocate(candidate_values(size(values)),candidate_gradients(3,size(values)))
    candidate_values=(0d0,0d0);candidate_gradients=(0d0,0d0)
    allow_zero=.false.;if(present(zero_outside_halo))allow_zero=zero_outside_halo
    do i=1,size(support_ids)
      owned_location=findloc(owned_ids,support_ids(i))
      if(owned_location(1)>0)then
        candidate_values(i)=owned_values(owned_location(1))
        candidate_gradients(:,i)=owned_gradients(:,owned_location(1))
      else
        nsource=0;covered=.false.
        do j=1,size(halos)
          if(halos(j)%valid.and.halos(j)%epoch==epoch.and.all(grid>=halos(j)%box_lo).and.&
             all(grid<=halos(j)%box_hi))then
            if(allocated(halos(j)%w_ids))covered=covered.or.any(halos(j)%w_ids==support_ids(i))
          endif
          call read_dg_wpw_volume_halo(halos(j),support_ids(i),grid,epoch,value,gradient,read_info)
          if(read_info==0)then
            nsource=nsource+1;candidate_values(i)=value;candidate_gradients(:,i)=gradient
          endif
        enddo
        if(nsource/=1)then
          if(nsource==0.and.allow_zero.and..not.covered)cycle
          return
        endif
      endif
    enddo
    if(.not.all(finite_complex(candidate_values)).or..not.all(finite_complex(candidate_gradients)))return
    values=candidate_values;gradients=candidate_gradients;info=0
  end subroutine

  subroutine reconstruct_dg_wpw_core_w_support(owned_ids,owned_values,owned_gradients,support_ids,halos,&
      grid,epoch,coefficients,support_values,support_gradients,reconstructed,info,zero_outside_halo)
    integer,intent(in)::owned_ids(:),support_ids(:),grid(3),epoch
    complex(8),intent(in)::owned_values(:),owned_gradients(:,:),coefficients(:,:)
    type(s_dg_wpw_volume_halo_state),intent(in)::halos(:)
    complex(8),intent(inout)::support_values(:),support_gradients(:,:)
    complex(8),intent(out)::reconstructed(:)
    integer,intent(out)::info
    logical,intent(in),optional::zero_outside_halo
    complex(8),allocatable::candidate(:)
    logical::allow_zero

    reconstructed=(0d0,0d0);info=1
    if(size(coefficients,1)/=size(support_ids).or.size(reconstructed)/=size(coefficients,2).or.&
       .not.all(finite_complex(coefficients)))return
    allow_zero=.false.;if(present(zero_outside_halo))allow_zero=zero_outside_halo
    call evaluate_dg_wpw_core_w_support(owned_ids,owned_values,owned_gradients,support_ids,halos,&
      grid,epoch,support_values,support_gradients,info,zero_outside_halo=allow_zero)
    if(info/=0)return
    allocate(candidate(size(reconstructed)))
    candidate=matmul(support_values,coefficients)
    if(.not.all(finite_complex(candidate)))then;info=1;return;endif
    reconstructed=candidate;info=0
  end subroutine reconstruct_dg_wpw_core_w_support

  logical function strictly_increasing(ids)result(ok)
    integer,intent(in)::ids(:);integer::i
    ok=size(ids)>0
    do i=2,size(ids)
      if(ids(i)<=ids(i-1))then;ok=.false.;return;endif
    enddo
  end function

  elemental logical function finite_complex(value)result(ok)
    use,intrinsic::ieee_arithmetic,only:ieee_is_finite
    complex(8),intent(in)::value
    ok=ieee_is_finite(real(value,8)).and.ieee_is_finite(aimag(value))
  end function
end module dg_wpw_core_w_provider
