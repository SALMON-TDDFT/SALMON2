module dg_overlapping_wannier_symmetry
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  integer,parameter::maximum_crystallographic_point_group_order=48
  integer,parameter::maximum_enumerated_subgroups=65536
  public::select_dg_exact_fragment_subgroup
  public::build_dg_fragment_group_representation
  public::project_dg_fragment_covariant_operators
  public::promote_dg_exact_global_subgroup
  public::build_dg_fragment_site_stabilizer
contains
  subroutine build_dg_fragment_site_stabilizer(rotations,translations,fragment_center,allowed, &
      tolerance,selected,product_table,maximum_site_residual,ok,message)
    integer,intent(in)::rotations(:,:,:)
    real(8),intent(in)::translations(:,:),fragment_center(3),tolerance
    logical,intent(in)::allowed(:)
    integer,allocatable,intent(out)::selected(:),product_table(:,:)
    real(8),intent(out)::maximum_site_residual
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,allocatable::work_selected(:)
    integer(8)::product_rotation(3,3)
    real(8)::mapped_center(3),residual,product_translation(3)
    integer::nop,iop,jop,kop,nselected,identity_position,i,j,k
    logical::duplicate,matched
    ok=.false.;message='';maximum_site_residual=0d0;nop=size(rotations,3)
    if(nop<1.or.size(rotations,1)/=3.or.size(rotations,2)/=3.or. &
        size(translations,1)/=3.or.size(translations,2)/=nop.or.size(allowed)/=nop)then
      message='fragment site-stabilizer arrays have inconsistent dimensions';return
    end if
    if(.not.ieee_is_finite(tolerance).or.tolerance<=0d0.or. &
        .not.all(ieee_is_finite(translations)).or..not.all(ieee_is_finite(fragment_center)))then
      message='fragment site-stabilizer coordinates or tolerance are not finite';return
    end if
    allocate(work_selected(min(nop,maximum_crystallographic_point_group_order)))
    nselected=0
    do iop=1,nop
      if(.not.allowed(iop))cycle
      mapped_center=matmul(real(rotations(:,:,iop),8),fragment_center)+translations(:,iop)-fragment_center
      residual=maxval(abs(mapped_center-anint(mapped_center)))
      if(residual>tolerance)cycle
      duplicate=.false.
      do i=1,nselected
        if(all(rotations(:,:,work_selected(i))==rotations(:,:,iop)))then
          duplicate=.true.;exit
        end if
      end do
      if(duplicate)cycle
      if(nselected>=maximum_crystallographic_point_group_order)then
        message='fragment site stabilizer exceeds crystallographic point-group order 48';return
      end if
      nselected=nselected+1;work_selected(nselected)=iop
      maximum_site_residual=max(maximum_site_residual,residual)
    end do
    if(nselected<1)then;message='fragment site stabilizer has no admissible operation';return;end if
    identity_position=0
    do i=1,nselected
      iop=work_selected(i)
      if(all(rotations(:,:,iop)==reshape([1,0,0,0,1,0,0,0,1],[3,3])).and. &
          maxval(abs(translations(:,iop)-anint(translations(:,iop))))<=tolerance)then
        identity_position=i;exit
      end if
    end do
    if(identity_position==0)then;message='fragment site stabilizer is missing identity';return;end if
    if(identity_position/=1)then
      iop=work_selected(1);work_selected(1)=work_selected(identity_position)
      work_selected(identity_position)=iop
    end if
    allocate(selected(nselected),product_table(nselected,nselected));selected=work_selected(1:nselected)
    do i=1,nselected;do j=1,nselected
      iop=selected(i);jop=selected(j)
      product_rotation=matmul(int(rotations(:,:,iop),8),int(rotations(:,:,jop),8))
      product_translation=translations(:,iop)+matmul(real(rotations(:,:,iop),8),translations(:,jop))
      matched=.false.;kop=0
      do k=1,nselected
        if(any(product_rotation/=int(rotations(:,:,selected(k)),8)))cycle
        residual=maxval(abs(product_translation-translations(:,selected(k))- &
          anint(product_translation-translations(:,selected(k)))))
        if(residual<=tolerance)then;matched=.true.;kop=k;exit;end if
      end do
      if(.not.matched)then
        message='fragment site stabilizer affine operations are not closed';return
      end if
      product_table(i,j)=kop
    end do;end do
    ok=.true.
  end subroutine build_dg_fragment_site_stabilizer

  subroutine promote_dg_exact_global_subgroup(product_table,fragment_exact,scalar_block_residual, &
      vector_block_residual,tolerance,subgroup,ok,message)
    integer,intent(in)::product_table(:,:)
    logical,intent(in)::fragment_exact(:,:)
    real(8),intent(in)::scalar_block_residual(:,:),vector_block_residual(:,:),tolerance
    integer,allocatable,intent(out)::subgroup(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(8),allocatable::fragment_residual(:),scalar_residual(:),vector_residual(:),zero_residual(:)
    integer::nop,iop
    nop=size(product_table,1);ok=.false.;message=''
    if(nop<1.or.size(product_table,2)/=nop.or.size(fragment_exact,1)<1.or. &
        size(fragment_exact,2)/=nop.or.size(scalar_block_residual,1)<1.or. &
        size(scalar_block_residual,2)/=nop.or.size(vector_block_residual,1)<1.or. &
        size(vector_block_residual,2)/=nop)then
      message='global symmetry promotion arrays have inconsistent dimensions';return
    end if
    allocate(fragment_residual(nop),scalar_residual(nop),vector_residual(nop),zero_residual(nop))
    zero_residual=0d0
    do iop=1,nop
      fragment_residual(iop)=merge(0d0,2d0*tolerance,all(fragment_exact(:,iop)))
      scalar_residual(iop)=maxval(scalar_block_residual(:,iop))
      vector_residual(iop)=maxval(vector_block_residual(:,iop))
    end do
    call select_dg_exact_fragment_subgroup(product_table,fragment_residual,scalar_residual, &
      vector_residual,zero_residual,tolerance,tolerance,tolerance,tolerance,subgroup,ok,message)
    if(.not.ok)message='global symmetry promotion failed: '//trim(message)
  end subroutine promote_dg_exact_global_subgroup

  subroutine project_dg_fragment_covariant_operators(representation,rotations,scalars,vectors, &
      tolerance,projected_scalars,projected_vectors,pre_projection_defect,post_projection_defect,ok,message)
    complex(8),intent(in)::representation(:,:,:),scalars(:,:,:),vectors(:,:,:,:)
    real(8),intent(in)::rotations(:,:,:),tolerance
    complex(8),allocatable,intent(out)::projected_scalars(:,:,:),projected_vectors(:,:,:,:)
    real(8),intent(out)::pre_projection_defect,post_projection_defect
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(8),allocatable::transformed(:,:),target(:,:)
    real(8)::scale,rotation_defect,determinant,correction_limit
    integer::n,nop,nscalar,nvector,iop,i,j,a,b
    ok=.false.;message='';pre_projection_defect=huge(1d0);post_projection_defect=huge(1d0)
    n=size(representation,1);nop=size(representation,3);nscalar=size(scalars,3);nvector=size(vectors,4)
    if(n<1.or.size(representation,2)/=n.or.nop<1.or.size(rotations,1)/=3.or. &
        size(rotations,2)/=3.or.size(rotations,3)/=nop.or.size(scalars,1)/=n.or. &
        size(scalars,2)/=n.or.nscalar<1.or.size(vectors,1)/=n.or.size(vectors,2)/=n.or. &
        size(vectors,3)/=3.or.nvector<1)then
      message='fragment covariant operator arrays have inconsistent dimensions';return
    end if
    if(.not.ieee_is_finite(tolerance).or.tolerance<=0d0.or. &
        .not.all(ieee_is_finite(real(representation))).or. &
        .not.all(ieee_is_finite(aimag(representation))).or. &
        .not.all(ieee_is_finite(rotations)).or..not.all(ieee_is_finite(real(scalars))).or. &
        .not.all(ieee_is_finite(aimag(scalars))).or..not.all(ieee_is_finite(real(vectors))).or. &
        .not.all(ieee_is_finite(aimag(vectors))))then
      message='fragment covariant operators or tolerance are not finite';return
    end if
    do iop=1,nop
      rotation_defect=maxval(abs(matmul(transpose(rotations(:,:,iop)),rotations(:,:,iop))-identity3()))
      determinant=determinant3(rotations(:,:,iop))
      if(rotation_defect>tolerance.or.abs(abs(determinant)-1d0)>tolerance)then
        message='fragment covariant operator Cartesian rotation is not orthogonal';return
      end if
    end do
    if(maxval(abs(representation(:,:,1)-complex_identity(n)))>tolerance.or. &
        maxval(abs(rotations(:,:,1)-identity3()))>tolerance)then
      message='fragment covariant operator group is not identity-first';return
    end if
    allocate(projected_scalars(n,n,nscalar),projected_vectors(n,n,3,nvector), &
      transformed(n,n),target(n,n))
    if(nop==1)then
      projected_scalars=scalars;projected_vectors=vectors
      pre_projection_defect=0d0;post_projection_defect=0d0;ok=.true.;return
    end if
    scale=max(1d0,max(maxval(abs(scalars)),maxval(abs(vectors))))
    call covariance_defect(representation,rotations,scalars,vectors,scale,pre_projection_defect)
    correction_limit=sqrt(tolerance)
    if(pre_projection_defect>correction_limit)then
      message='fragment operator pre-projection covariance defect exceeds correction limit';return
    end if
    projected_scalars=(0d0,0d0);projected_vectors=(0d0,0d0)
    do iop=1,nop
      do i=1,nscalar
        projected_scalars(:,:,i)=projected_scalars(:,:,i)+matmul(conjg(transpose( &
          representation(:,:,iop))),matmul(scalars(:,:,i),representation(:,:,iop)))
      end do
      do i=1,nvector
        do a=1,3
          target=(0d0,0d0)
          do b=1,3
            transformed=matmul(conjg(transpose(representation(:,:,iop))), &
              matmul(vectors(:,:,b,i),representation(:,:,iop)))
            target=target+rotations(b,a,iop)*transformed
          end do
          projected_vectors(:,:,a,i)=projected_vectors(:,:,a,i)+target
        end do
      end do
    end do
    projected_scalars=projected_scalars/real(nop,8);projected_vectors=projected_vectors/real(nop,8)
    call covariance_defect(representation,rotations,projected_scalars,projected_vectors, &
      scale,post_projection_defect)
    if(post_projection_defect>tolerance)then
      message='fragment operator post-projection covariance exceeds tolerance';return
    end if
    ok=.true.
  end subroutine project_dg_fragment_covariant_operators

  subroutine covariance_defect(representation,rotations,scalars,vectors,scale,defect)
    complex(8),intent(in)::representation(:,:,:),scalars(:,:,:),vectors(:,:,:,:)
    real(8),intent(in)::rotations(:,:,:),scale
    real(8),intent(out)::defect
    complex(8),allocatable::transformed(:,:),target(:,:)
    integer::n,iop,i,a,b
    n=size(representation,1);allocate(transformed(n,n),target(n,n));defect=0d0
    do iop=1,size(representation,3)
      do i=1,size(scalars,3)
        transformed=matmul(conjg(transpose(representation(:,:,iop))), &
          matmul(scalars(:,:,i),representation(:,:,iop)))
        defect=max(defect,maxval(abs(transformed-scalars(:,:,i)))/scale)
      end do
      do i=1,size(vectors,4);do a=1,3
        transformed=matmul(conjg(transpose(representation(:,:,iop))), &
          matmul(vectors(:,:,a,i),representation(:,:,iop)))
        target=(0d0,0d0)
        do b=1,3;target=target+rotations(a,b,iop)*vectors(:,:,b,i);end do
        defect=max(defect,maxval(abs(transformed-target))/scale)
      end do;end do
    end do
  end subroutine covariance_defect

  function complex_identity(n)result(identity)
    integer,intent(in)::n
    complex(8)::identity(n,n)
    integer::i
    identity=(0d0,0d0);do i=1,n;identity(i,i)=1d0;end do
  end function complex_identity

  function identity3()result(identity)
    real(8)::identity(3,3)
    identity=0d0;identity(1,1)=1d0;identity(2,2)=1d0;identity(3,3)=1d0
  end function identity3

  real(8) function determinant3(matrix)
    real(8),intent(in)::matrix(3,3)
    determinant3=matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2))- &
      matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1))+ &
      matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
  end function determinant3

  subroutine build_dg_fragment_group_representation(metric,raw,product_table,tolerance, &
      representation,raw_unitarity_defect,unitarity_defect,closure_defect,ok,message)
    complex(8),intent(in)::metric(:,:),raw(:,:,:)
    integer,intent(in)::product_table(:,:)
    real(8),intent(in)::tolerance
    complex(8),allocatable,intent(out)::representation(:,:,:)
    real(8),intent(out)::raw_unitarity_defect,unitarity_defect,closure_defect
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(8),allocatable::metric_sqrt(:,:),metric_inverse_sqrt(:,:),gram_inverse_sqrt(:,:), &
      transformed(:,:),unitary(:,:),difference(:,:)
    real(8)::metric_scale,representation_scale,correction_limit
    integer::n,nop,iop,jop,kop
    logical::power_ok
    ok=.false.;message='';raw_unitarity_defect=huge(1d0);unitarity_defect=huge(1d0)
    closure_defect=huge(1d0);n=size(metric,1);nop=size(raw,3)
    if(n<1.or.size(metric,2)/=n.or.size(raw,1)/=n.or.size(raw,2)/=n.or.nop<1.or. &
        size(product_table,1)/=nop.or.size(product_table,2)/=nop)then
      message='fragment group representation arrays have inconsistent dimensions';return
    end if
    if(.not.ieee_is_finite(tolerance).or.tolerance<=0d0)then
      message='fragment group representation tolerance must be finite and positive';return
    end if
    if(any(product_table<1).or.any(product_table>nop))then
      message='fragment group representation product table has an invalid index';return
    end if
    metric_scale=max(1d0,maxval(abs(metric)))
    if(maxval(abs(metric-conjg(transpose(metric))))>tolerance*metric_scale)then
      message='fragment group representation metric is not Hermitian';return
    end if
    call hermitian_matrix_power(metric,0.5d0,tolerance,metric_sqrt,power_ok)
    if(.not.power_ok)then;message='fragment group representation metric is not positive definite';return;end if
    call hermitian_matrix_power(metric,-0.5d0,tolerance,metric_inverse_sqrt,power_ok)
    if(.not.power_ok)then;message='fragment group representation metric inverse failed';return;end if
    allocate(representation(n,n,nop),transformed(n,n),unitary(n,n),difference(n,n))
    raw_unitarity_defect=0d0
    do iop=1,nop
      difference=matmul(conjg(transpose(raw(:,:,iop))),matmul(metric,raw(:,:,iop)))-metric
      raw_unitarity_defect=max(raw_unitarity_defect,maxval(abs(difference))/metric_scale)
      transformed=matmul(metric_sqrt,matmul(raw(:,:,iop),metric_inverse_sqrt))
      call hermitian_matrix_power(matmul(conjg(transpose(transformed)),transformed),-0.5d0, &
        tolerance,gram_inverse_sqrt,power_ok)
      if(.not.power_ok)then;message='fragment group representation polar correction is singular';return;end if
      unitary=matmul(transformed,gram_inverse_sqrt)
      representation(:,:,iop)=matmul(metric_inverse_sqrt,matmul(unitary,metric_sqrt))
    end do
    correction_limit=sqrt(tolerance)
    if(raw_unitarity_defect>correction_limit)then
      message='fragment group representation raw unitarity defect exceeds correction limit';return
    end if
    unitarity_defect=0d0
    do iop=1,nop
      difference=matmul(conjg(transpose(representation(:,:,iop))), &
        matmul(metric,representation(:,:,iop)))-metric
      unitarity_defect=max(unitarity_defect,maxval(abs(difference))/metric_scale)
    end do
    representation_scale=max(1d0,maxval(abs(representation)))
    closure_defect=0d0
    do iop=1,nop;do jop=1,nop
      kop=product_table(iop,jop)
      difference=matmul(representation(:,:,iop),representation(:,:,jop))-representation(:,:,kop)
      closure_defect=max(closure_defect,maxval(abs(difference))/representation_scale)
    end do;end do
    if(unitarity_defect>tolerance)then
      message='fragment group representation metric unitarity exceeds tolerance';return
    end if
    if(closure_defect>tolerance)then
      message='fragment group representation closure exceeds tolerance';return
    end if
    ok=.true.
  end subroutine build_dg_fragment_group_representation

  subroutine hermitian_matrix_power(matrix,power,tolerance,result,ok)
    complex(8),intent(in)::matrix(:,:)
    real(8),intent(in)::power,tolerance
    complex(8),allocatable,intent(out)::result(:,:)
    logical,intent(out)::ok
    complex(8),allocatable::vectors(:,:),work(:)
    real(8),allocatable::eigenvalues(:),rwork(:)
    complex(8)::work_query(1)
    integer::n,info,lwork,i
    n=size(matrix,1);ok=.false.
    allocate(vectors(n,n),eigenvalues(n),rwork(max(1,3*n-2)))
    vectors=matrix
    call zheev('V','U',n,vectors,n,eigenvalues,work_query,-1,rwork,info)
    if(info/=0)return
    lwork=max(1,int(real(work_query(1))));allocate(work(lwork));vectors=matrix
    call zheev('V','U',n,vectors,n,eigenvalues,work,lwork,rwork,info)
    if(info/=0.or.minval(eigenvalues)<=tolerance*max(1d0,maxval(abs(eigenvalues))))return
    allocate(result(n,n));result=(0d0,0d0)
    do i=1,n
      result=result+eigenvalues(i)**power*spread(vectors(:,i),2,n)* &
        spread(conjg(vectors(:,i)),1,n)
    end do
    ok=.true.
  end subroutine hermitian_matrix_power

  subroutine select_dg_exact_fragment_subgroup(product_table,atom_residual,boundary_residual, &
      grid_residual,center_residual,atom_tolerance,boundary_tolerance,grid_tolerance, &
      center_tolerance,subgroup,ok,message)
    integer,intent(in)::product_table(:,:)
    real(8),intent(in)::atom_residual(:),boundary_residual(:),grid_residual(:),center_residual(:)
    real(8),intent(in)::atom_tolerance,boundary_tolerance,grid_tolerance,center_tolerance
    integer,allocatable,intent(out)::subgroup(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    logical,allocatable::accepted(:),groups(:,:),candidate(:)
    integer::n,igroup,iop,jop,ngroups,best_group,best_size,candidate_size
    logical::has_inverse

    ok=.false.;message='';n=size(product_table,1)
    if(n<1.or.size(product_table,2)/=n.or.size(atom_residual)/=n.or. &
        size(boundary_residual)/=n.or.size(grid_residual)/=n.or.size(center_residual)/=n)then
      message='exact fragment subgroup arrays have inconsistent dimensions';return
    end if
    if(n>maximum_crystallographic_point_group_order)then
      message='exact fragment subgroup exceeds crystallographic point-group order 48';return
    end if
    if(any(product_table<1).or.any(product_table>n))then
      message='exact fragment subgroup product table has an invalid index';return
    end if
    do iop=1,n
      if(product_table(1,iop)/=iop.or.product_table(iop,1)/=iop)then
        message='exact fragment subgroup operation 1 is not identity';return
      end if
      has_inverse=.false.
      do jop=1,n
        if(product_table(iop,jop)==1.and.product_table(jop,iop)==1)then
          has_inverse=.true.;exit
        end if
      end do
      if(.not.has_inverse)then
        message='exact fragment subgroup product table is not a group';return
      end if
    end do
    if(.not.valid_tolerance(atom_tolerance).or..not.valid_tolerance(boundary_tolerance).or. &
        .not.valid_tolerance(grid_tolerance).or..not.valid_tolerance(center_tolerance))then
      message='exact fragment subgroup tolerances must be finite and positive';return
    end if
    if(.not.valid_residuals(atom_residual).or..not.valid_residuals(boundary_residual).or. &
        .not.valid_residuals(grid_residual).or..not.valid_residuals(center_residual))then
      message='exact fragment subgroup residuals must be finite and nonnegative';return
    end if
    allocate(accepted(n))
    accepted=atom_residual<=atom_tolerance.and.boundary_residual<=boundary_tolerance.and. &
      grid_residual<=grid_tolerance.and.center_residual<=center_tolerance
    if(.not.accepted(1))then
      message='exact fragment subgroup identity operation failed a mapping tolerance';return
    end if

    allocate(groups(n,1));groups=.false.;groups(1,1)=.true.;ngroups=1
    igroup=1
    do while(igroup<=ngroups)
      do iop=2,n
        if(.not.accepted(iop).or.groups(iop,igroup))cycle
        candidate=groups(:,igroup);candidate(iop)=.true.
        call close_generated_subgroup(product_table,candidate)
        if(any(candidate.and..not.accepted))cycle
        if(.not.subgroup_is_known(groups,ngroups,candidate))then
          if(ngroups>=maximum_enumerated_subgroups)then
            message='exact fragment subgroup enumeration limit exceeded';return
          end if
          call append_subgroup(groups,ngroups,candidate)
        end if
      end do
      igroup=igroup+1
    end do

    best_group=1;best_size=1
    do igroup=2,ngroups
      candidate_size=count(groups(:,igroup))
      if(candidate_size>best_size.or.(candidate_size==best_size.and. &
          subgroup_lexically_precedes(groups(:,igroup),groups(:,best_group))))then
        best_group=igroup;best_size=candidate_size
      end if
    end do
    allocate(subgroup(best_size));subgroup=pack([(iop,iop=1,n)],groups(:,best_group))
    ok=.true.
  end subroutine select_dg_exact_fragment_subgroup

  logical function valid_tolerance(value)
    real(8),intent(in)::value
    valid_tolerance=ieee_is_finite(value).and.value>0d0
  end function valid_tolerance

  logical function valid_residuals(values)
    real(8),intent(in)::values(:)
    valid_residuals=all(ieee_is_finite(values)).and.all(values>=0d0)
  end function valid_residuals

  subroutine close_generated_subgroup(product_table,mask)
    integer,intent(in)::product_table(:,:)
    logical,intent(inout)::mask(:)
    integer::i,j
    logical::changed
    changed=.true.
    do while(changed)
      changed=.false.
      do i=1,size(mask)
        if(.not.mask(i))cycle
        do j=1,size(mask)
          if(.not.mask(j))cycle
          if(.not.mask(product_table(i,j)))then
            mask(product_table(i,j))=.true.;changed=.true.
          end if
        end do
      end do
    end do
  end subroutine close_generated_subgroup

  logical function subgroup_is_known(groups,ngroups,candidate)
    logical,intent(in)::groups(:,:),candidate(:)
    integer,intent(in)::ngroups
    integer::i
    subgroup_is_known=.false.
    do i=1,ngroups
      if(all(groups(:,i).eqv.candidate))then
        subgroup_is_known=.true.;return
      end if
    end do
  end function subgroup_is_known

  subroutine append_subgroup(groups,ngroups,candidate)
    logical,allocatable,intent(inout)::groups(:,:)
    integer,intent(inout)::ngroups
    logical,intent(in)::candidate(:)
    logical,allocatable::expanded(:,:)
    allocate(expanded(size(groups,1),ngroups+1))
    expanded(:,1:ngroups)=groups(:,1:ngroups);expanded(:,ngroups+1)=candidate
    call move_alloc(expanded,groups);ngroups=ngroups+1
  end subroutine append_subgroup

  logical function subgroup_lexically_precedes(left,right)
    logical,intent(in)::left(:),right(:)
    integer::i
    subgroup_lexically_precedes=.false.
    do i=1,size(left)
      if(left(i).eqv.right(i))cycle
      subgroup_lexically_precedes=left(i)
      return
    end do
  end function subgroup_lexically_precedes
end module dg_overlapping_wannier_symmetry
