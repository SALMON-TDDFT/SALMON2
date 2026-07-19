module dg_wpw_lcfo_operator_adapter
  use dg_wpw_production_context,only:s_dg_wpw_production_context
  use,intrinsic::iso_fortran_env,only:int64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  type,public::s_dg_wpw_lcfo_ww_components
    integer,allocatable::owned_w_ids(:)
    real(8),allocatable::kinetic(:,:),potential(:,:),nonlocal(:,:),face_self(:,:)
    integer,allocatable::cross_face_id(:),cross_row_id(:),cross_col_id(:),cross_axis(:),cross_side(:)
    integer,allocatable::cross_image(:,:)
    real(8),allocatable::cross_value(:)
    character(32)::metric_convention=''
    integer(int64)::provenance_fingerprint=0_int64
    logical::valid=.false.
  end type
  public::import_dg_wpw_lcfo_ww_components
  public::publish_dg_wpw_lcfo_ww_components
  public::validate_dg_wpw_surface_self_policy
contains
  subroutine validate_dg_wpw_surface_self_policy(enabled,info)
    logical,intent(in)::enabled
    integer,intent(out)::info
    info=merge(0,1,enabled)
  end subroutine
  subroutine import_dg_wpw_lcfo_ww_components(local_fragment,n_basis,index_basis,kinetic,potential,nonlocal,&
      face_self,face_id,face_row,face_col,face_axis,face_side,face_image,face_value,metric,components,info)
    integer,intent(in)::local_fragment,n_basis(:),index_basis(:,:)
    real(8),intent(in)::kinetic(:,:),potential(:,:),nonlocal(:,:),face_self(:,:)
    integer,intent(in)::face_id(:),face_row(:),face_col(:),face_axis(:),face_side(:),face_image(:,:)
    real(8),intent(in)::face_value(:)
    character(*),intent(in)::metric
    type(s_dg_wpw_lcfo_ww_components),intent(out)::components
    integer,intent(out)::info
    integer::nlocal,nface,i,j
    info=1
    if(local_fragment<1.or.local_fragment>size(n_basis).or.local_fragment>size(index_basis,2))return
    nlocal=n_basis(local_fragment);nface=size(face_id)
    if(nlocal<0.or.nlocal>size(index_basis,1))return
    if(any(shape(kinetic)/=[nlocal,nlocal]).or.any(shape(potential)/=[nlocal,nlocal]).or.&
       any(shape(nonlocal)/=[nlocal,nlocal]).or.any(shape(face_self)/=[nlocal,nlocal]))return
    if(size(face_row)/=nface.or.size(face_col)/=nface.or.size(face_axis)/=nface.or.&
       size(face_side)/=nface.or.size(face_value)/=nface.or.any(shape(face_image)/=[3,nface]))return
    if(any(index_basis(1:nlocal,local_fragment)<=0).or.any(face_id<=0).or.any(face_row<=0).or.any(face_col<=0))return
    if(any(face_axis<1).or.any(face_axis>3).or.any(abs(face_side)/=1).or.len_trim(metric)==0)return
    do i=1,nface;do j=i+1,nface
      if(face_id(i)==face_id(j))then
        if(face_axis(i)/=face_axis(j).or.face_side(i)/=face_side(j).or.any(face_image(:,i)/=face_image(:,j)))return
        if(face_row(i)==face_row(j).and.face_col(i)==face_col(j))return
      endif
    enddo;enddo
    if(.not.finite_matrix(kinetic).or..not.finite_matrix(potential).or..not.finite_matrix(nonlocal).or.&
       .not.finite_matrix(face_self).or..not.all(ieee_is_finite(face_value)))return
    components%owned_w_ids=index_basis(1:nlocal,local_fragment)
    components%kinetic=kinetic;components%potential=potential;components%nonlocal=nonlocal
    components%face_self=face_self;components%cross_face_id=face_id;components%cross_row_id=face_row
    components%cross_col_id=face_col;components%cross_axis=face_axis;components%cross_side=face_side
    components%cross_image=face_image;components%cross_value=face_value
    components%metric_convention=metric
    components%provenance_fingerprint=fingerprint(components)
    if(components%provenance_fingerprint==0_int64)components%provenance_fingerprint=1_int64
    components%valid=.true.;info=0
  end subroutine

  subroutine publish_dg_wpw_lcfo_ww_components(context,components,info)
    type(s_dg_wpw_production_context),intent(inout)::context
    type(s_dg_wpw_lcfo_ww_components),intent(in)::components
    integer,intent(out)::info
    info=1
    if(.not.components%valid.or.trim(components%metric_convention)/='orthonormal_ww'.or.&
       .not.ids_equal(components%owned_w_ids,context%owned_w_ids))return
    context%pending_ww_owned_w_ids=components%owned_w_ids
    context%pending_ww_kinetic=components%kinetic;context%pending_ww_potential=components%potential
    context%pending_ww_nonlocal=components%nonlocal;context%pending_ww_face_self=components%face_self
    context%pending_ww_cross_face_id=components%cross_face_id
    context%pending_ww_cross_row_id=components%cross_row_id
    context%pending_ww_cross_col_id=components%cross_col_id
    context%pending_ww_cross_axis=components%cross_axis;context%pending_ww_cross_side=components%cross_side
    context%pending_ww_cross_image=components%cross_image;context%pending_ww_cross_value=components%cross_value
    context%pending_ww_metric_convention=components%metric_convention
    context%pending_ww_provenance_fingerprint=components%provenance_fingerprint
    context%pending_ww_valid=.true.;info=0
  end subroutine

  logical function ids_equal(left,right)result(equal)
    integer,intent(in)::left(:),right(:)
    equal=size(left)==size(right)
    if(equal)equal=all(left==right)
  end function

  logical function finite_matrix(a)result(ok)
    real(8),intent(in)::a(:,:);ok=all(ieee_is_finite(a))
  end function
  integer(int64) function fingerprint(c)result(hash)
    type(s_dg_wpw_lcfo_ww_components),intent(in)::c
    integer::i,j
    hash=1469598103934665603_int64
    do i=1,size(c%owned_w_ids);call mix(hash,int(c%owned_w_ids(i),int64));enddo
    do j=1,size(c%kinetic,2);do i=1,size(c%kinetic,1);call mix(hash,transfer(c%kinetic(i,j),hash));enddo;enddo
    do j=1,size(c%potential,2);do i=1,size(c%potential,1);call mix(hash,transfer(c%potential(i,j),hash));enddo;enddo
    do j=1,size(c%nonlocal,2);do i=1,size(c%nonlocal,1);call mix(hash,transfer(c%nonlocal(i,j),hash));enddo;enddo
    do j=1,size(c%face_self,2);do i=1,size(c%face_self,1);call mix(hash,transfer(c%face_self(i,j),hash));enddo;enddo
    do i=1,size(c%cross_face_id)
      call mix(hash,int(c%cross_face_id(i),int64));call mix(hash,int(c%cross_row_id(i),int64))
      call mix(hash,int(c%cross_col_id(i),int64));call mix(hash,int(c%cross_axis(i),int64))
      call mix(hash,int(c%cross_side(i),int64));call mix(hash,int(c%cross_image(1,i),int64))
      call mix(hash,int(c%cross_image(2,i),int64));call mix(hash,int(c%cross_image(3,i),int64))
      call mix(hash,transfer(c%cross_value(i),hash))
    enddo
    do i=1,len_trim(c%metric_convention);call mix(hash,int(iachar(c%metric_convention(i:i)),int64));enddo
  end function
  subroutine mix(hash,value)
    integer(int64),intent(inout)::hash
    integer(int64),intent(in)::value
    hash=ieor(hash,value);hash=ieor(hash,shiftl(hash,13));hash=ieor(hash,shiftr(hash,7))
  end subroutine
end module
