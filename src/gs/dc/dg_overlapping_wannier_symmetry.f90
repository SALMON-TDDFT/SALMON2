module dg_overlapping_wannier_symmetry
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  integer,parameter::maximum_crystallographic_point_group_order=48
  integer,parameter::maximum_enumerated_subgroups=65536
  public::select_dg_exact_fragment_subgroup
contains
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
