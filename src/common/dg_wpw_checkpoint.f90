module dg_wpw_checkpoint
  use,intrinsic::iso_c_binding,only:c_char,c_int,c_null_char
  use,intrinsic::iso_fortran_env,only:int8,int64,iostat_end
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  character(16),parameter::checkpoint_magic='DG_WPW_CHECKPT_V2'
  character(16),parameter::manifest_magic='DG_WPW_MANIFEST1'
  integer,parameter::checkpoint_version=2,max_dimension=10000000
  type,public::s_dg_wpw_checkpoint_state
    integer::schema_version=checkpoint_version,operator_epoch=-1,n_occ=0
    integer(int64)::layout_fingerprint=0
    character(32)::ownership_kind='',metric_convention='',operator_convention=''
    integer,allocatable::peer_ranks(:),owned_w_ids(:),owned_p_ids(:),support_w_ids(:),support_p_ids(:)
    integer,allocatable::ww_r(:),ww_c(:),wp_w(:),wp_p(:),pp_r(:),pp_c(:)
    real(8),allocatable::eigenvalues(:),occupations(:),potential(:)
    complex(8),allocatable::coeff_w(:,:),coeff_p(:,:)
    complex(8),allocatable::ww_h(:),ww_s(:),wp_h(:),wp_s(:),pp_h(:),pp_s(:)
    complex(8),allocatable::ww_z(:),wp_z(:),pp_z(:)
  end type
  public::write_dg_wpw_checkpoint,read_dg_wpw_checkpoint,dg_wpw_checkpoint_checksum
  public::write_dg_wpw_checkpoint_manifest,read_dg_wpw_checkpoint_manifest
  interface
    integer(c_int) function c_rename(old_path,new_path) bind(C,name='rename')
      import c_char,c_int
      character(c_char),intent(in)::old_path(*),new_path(*)
    end function
  end interface
contains
  subroutine write_dg_wpw_checkpoint(path,state,info)
    character(*),intent(in)::path
    type(s_dg_wpw_checkpoint_state),intent(in)::state
    integer,intent(out)::info
    character(:),allocatable::tmp
    integer::unit,ios,d(18)
    integer(int64)::checksum
    info=1
    if(.not.valid_state(state))return
    call dimensions(state,d);checksum=dg_wpw_checkpoint_checksum(state)
    if(checksum==0)return
    tmp=trim(path)//'.tmp'
    open(newunit=unit,file=tmp,access='stream',form='unformatted',status='replace',action='write',iostat=ios)
    if(ios/=0)return
    write(unit,iostat=ios)checkpoint_magic,state%schema_version,state%operator_epoch,state%n_occ,&
      state%layout_fingerprint,state%ownership_kind,state%metric_convention,state%operator_convention,d,checksum
    if(ios==0)write(unit,iostat=ios)state%peer_ranks,state%owned_w_ids,state%owned_p_ids,state%support_w_ids,state%support_p_ids,&
      state%eigenvalues,state%occupations,state%coeff_w,state%coeff_p,state%potential,&
      state%ww_r,state%ww_c,state%wp_w,state%wp_p,state%pp_r,state%pp_c,&
      state%ww_h,state%ww_s,state%wp_h,state%wp_s,state%pp_h,state%pp_s,&
      state%ww_z,state%wp_z,state%pp_z
    flush(unit);close(unit,iostat=ios)
    if(ios/=0)return
    if(c_rename(trim(tmp)//c_null_char,trim(path)//c_null_char)/=0_c_int)return
    info=0
  end subroutine

  subroutine read_dg_wpw_checkpoint(path,state,info,expected_fingerprint)
    character(*),intent(in)::path
    type(s_dg_wpw_checkpoint_state),intent(out)::state
    integer,intent(out)::info
    integer(int64),intent(in),optional::expected_fingerprint
    character(16)::magic
    integer::unit,ios,d(18)
    integer(int64)::stored_checksum
    integer(int8)::trailing
    info=1
    open(newunit=unit,file=path,access='stream',form='unformatted',status='old',action='read',iostat=ios)
    if(ios/=0)return
    read(unit,iostat=ios)magic,state%schema_version,state%operator_epoch,state%n_occ,&
      state%layout_fingerprint,state%ownership_kind,state%metric_convention,state%operator_convention,d,stored_checksum
    if(ios/=0.or.magic/=checkpoint_magic.or.state%schema_version/=checkpoint_version.or.&
       any(d<0).or.any(d>max_dimension))then;close(unit);return;endif
    if(present(expected_fingerprint))then
      if(state%layout_fingerprint/=expected_fingerprint)then;close(unit);return;endif
    endif
    if(.not.dimension_relations(d,state%n_occ))then;close(unit);return;endif
    allocate(state%peer_ranks(d(18)),state%owned_w_ids(d(1)),state%owned_p_ids(d(2)),&
      state%support_w_ids(d(3)),state%support_p_ids(d(4)))
    allocate(state%eigenvalues(d(5)),state%occupations(d(5)),state%coeff_w(d(1),d(5)),state%coeff_p(d(2),d(5)))
    allocate(state%potential(d(6)))
    allocate(state%ww_r(d(7)),state%ww_c(d(7)),state%ww_h(d(7)),state%ww_s(d(7)))
    allocate(state%wp_w(d(8)),state%wp_p(d(8)),state%wp_h(d(8)),state%wp_s(d(8)))
    allocate(state%pp_r(d(9)),state%pp_c(d(9)),state%pp_h(d(9)),state%pp_s(d(9)))
    allocate(state%ww_z(d(7)),state%wp_z(d(8)),state%pp_z(d(9)))
    read(unit,iostat=ios)state%peer_ranks,state%owned_w_ids,state%owned_p_ids,state%support_w_ids,state%support_p_ids,&
      state%eigenvalues,state%occupations,state%coeff_w,state%coeff_p,state%potential,&
      state%ww_r,state%ww_c,state%wp_w,state%wp_p,state%pp_r,state%pp_c,&
      state%ww_h,state%ww_s,state%wp_h,state%wp_s,state%pp_h,state%pp_s,&
      state%ww_z,state%wp_z,state%pp_z
    if(ios/=0)then;close(unit);return;endif
    read(unit,iostat=ios)trailing
    close(unit)
    if(ios/=iostat_end.or..not.valid_state(state))return
    if(dg_wpw_checkpoint_checksum(state)/=stored_checksum)return
    info=0
  end subroutine

  subroutine write_dg_wpw_checkpoint_manifest(path,checksums,info)
    character(*),intent(in)::path
    integer(int64),intent(in)::checksums(:)
    integer,intent(out)::info
    character(:),allocatable::tmp
    integer::unit,ios
    info=1
    if(size(checksums)<=0.or.size(checksums)>max_dimension.or.any(checksums==0))return
    tmp=trim(path)//'.tmp'
    open(newunit=unit,file=tmp,access='stream',form='unformatted',status='replace',action='write',iostat=ios)
    if(ios/=0)return
    write(unit,iostat=ios)manifest_magic,checkpoint_version,size(checksums),checksums
    flush(unit);close(unit,iostat=ios)
    if(ios/=0)return
    if(c_rename(trim(tmp)//c_null_char,trim(path)//c_null_char)/=0_c_int)return
    info=0
  end subroutine

  subroutine read_dg_wpw_checkpoint_manifest(path,checksums,info)
    character(*),intent(in)::path
    integer(int64),allocatable,intent(out)::checksums(:)
    integer,intent(out)::info
    character(16)::magic
    integer::unit,ios,version,nrank
    integer(int8)::trailing
    info=1
    open(newunit=unit,file=path,access='stream',form='unformatted',status='old',action='read',iostat=ios)
    if(ios/=0)return
    read(unit,iostat=ios)magic,version,nrank
    if(ios/=0.or.magic/=manifest_magic.or.version/=checkpoint_version.or.&
       nrank<=0.or.nrank>max_dimension)then;close(unit);return;endif
    allocate(checksums(nrank));read(unit,iostat=ios)checksums
    if(ios/=0)then;close(unit);return;endif
    read(unit,iostat=ios)trailing;close(unit)
    if(ios/=iostat_end.or.any(checksums==0))return
    info=0
  end subroutine

  subroutine dimensions(s,d)
    type(s_dg_wpw_checkpoint_state),intent(in)::s
    integer,intent(out)::d(18)
    d=[size(s%owned_w_ids),size(s%owned_p_ids),size(s%support_w_ids),size(s%support_p_ids),&
      size(s%eigenvalues),size(s%potential),size(s%ww_r),size(s%wp_w),size(s%pp_r),&
      size(s%coeff_w,1),size(s%coeff_w,2),size(s%coeff_p,1),size(s%coeff_p,2),&
      size(s%occupations),size(s%ww_h),size(s%wp_h),size(s%pp_h),size(s%peer_ranks)]
  end subroutine

  logical function dimension_relations(d,n_occ)result(ok)
    integer,intent(in)::d(18),n_occ
    ok=d(1)>0.and.d(2)>0.and.d(3)>=d(1).and.d(4)>=d(2).and.d(5)>0.and.n_occ>0.and.n_occ<=d(5).and.&
      d(10)==d(1).and.d(11)==d(5).and.d(12)==d(2).and.d(13)==d(5).and.d(14)==d(5).and.&
      d(15)==d(7).and.d(16)==d(8).and.d(17)==d(9)
  end function

  logical function valid_state(s)result(ok)
    type(s_dg_wpw_checkpoint_state),intent(in)::s
    integer::d(18)
    ok=.false.
    if(s%schema_version/=checkpoint_version.or.s%operator_epoch<=0.or.s%layout_fingerprint==0.or.&
       len_trim(s%ownership_kind)==0.or.len_trim(s%metric_convention)==0.or.&
       len_trim(s%operator_convention)==0)return
    if(.not.all_allocated(s))return
    call dimensions(s,d)
    if(any(d<0).or.any(d>max_dimension).or..not.dimension_relations(d,s%n_occ))return
    if(size(s%ww_c)/=d(7).or.size(s%wp_p)/=d(8).or.size(s%pp_c)/=d(9).or.&
       size(s%ww_s)/=d(7).or.size(s%wp_s)/=d(8).or.size(s%pp_s)/=d(9).or.&
       size(s%ww_z)/=d(7).or.size(s%wp_z)/=d(8).or.size(s%pp_z)/=d(9))return
    if(.not.strict_ids(s%owned_w_ids).or..not.strict_ids(s%owned_p_ids).or.&
       .not.strict_ids(s%support_w_ids).or..not.strict_ids(s%support_p_ids))return
    if(.not.increasing_nonnegative(s%peer_ranks))return
    if(.not.finite_real(s%eigenvalues).or..not.finite_real(s%occupations).or..not.finite_real(s%potential).or.&
       .not.finite_complex2(s%coeff_w).or..not.finite_complex2(s%coeff_p).or.&
       .not.finite_complex1(s%ww_h).or..not.finite_complex1(s%ww_s).or.&
       .not.finite_complex1(s%wp_h).or..not.finite_complex1(s%wp_s).or.&
       .not.finite_complex1(s%pp_h).or..not.finite_complex1(s%pp_s).or.&
       .not.finite_complex1(s%ww_z).or..not.finite_complex1(s%wp_z).or..not.finite_complex1(s%pp_z))return
    ok=.true.
  end function

  logical function all_allocated(s)result(ok)
    type(s_dg_wpw_checkpoint_state),intent(in)::s
    ok=allocated(s%peer_ranks).and.allocated(s%owned_w_ids).and.allocated(s%owned_p_ids).and.allocated(s%support_w_ids).and.&
      allocated(s%support_p_ids).and.allocated(s%eigenvalues).and.allocated(s%occupations).and.&
      allocated(s%coeff_w).and.allocated(s%coeff_p).and.allocated(s%potential).and.&
      allocated(s%ww_r).and.allocated(s%ww_c).and.allocated(s%wp_w).and.allocated(s%wp_p).and.&
      allocated(s%pp_r).and.allocated(s%pp_c).and.allocated(s%ww_h).and.allocated(s%ww_s).and.&
      allocated(s%wp_h).and.allocated(s%wp_s).and.allocated(s%pp_h).and.allocated(s%pp_s).and.&
      allocated(s%ww_z).and.allocated(s%wp_z).and.allocated(s%pp_z)
  end function

  integer(int64) function dg_wpw_checkpoint_checksum(s)result(h)
    type(s_dg_wpw_checkpoint_state),intent(in)::s
    integer::i
    h=1469598103934665603_int64
    call mix_i(h,int(s%schema_version,int64));call mix_i(h,int(s%operator_epoch,int64))
    call mix_i(h,int(s%n_occ,int64));call mix_i(h,s%layout_fingerprint)
    do i=1,len(s%ownership_kind);call mix_i(h,int(iachar(s%ownership_kind(i:i)),int64));enddo
    do i=1,len(s%metric_convention);call mix_i(h,int(iachar(s%metric_convention(i:i)),int64));enddo
    do i=1,len(s%operator_convention);call mix_i(h,int(iachar(s%operator_convention(i:i)),int64));enddo
    call mix_int_array(h,s%peer_ranks);call mix_int_array(h,s%owned_w_ids);call mix_int_array(h,s%owned_p_ids)
    call mix_int_array(h,s%support_w_ids);call mix_int_array(h,s%support_p_ids)
    call mix_real_array(h,s%eigenvalues);call mix_real_array(h,s%occupations);call mix_real_array(h,s%potential)
    call mix_complex_array(h,reshape(s%coeff_w,[size(s%coeff_w)]));call mix_complex_array(h,reshape(s%coeff_p,[size(s%coeff_p)]))
    call mix_int_array(h,s%ww_r);call mix_int_array(h,s%ww_c);call mix_int_array(h,s%wp_w);call mix_int_array(h,s%wp_p)
    call mix_int_array(h,s%pp_r);call mix_int_array(h,s%pp_c)
    call mix_complex_array(h,s%ww_h);call mix_complex_array(h,s%ww_s);call mix_complex_array(h,s%wp_h)
    call mix_complex_array(h,s%wp_s);call mix_complex_array(h,s%pp_h);call mix_complex_array(h,s%pp_s)
    call mix_complex_array(h,s%ww_z);call mix_complex_array(h,s%wp_z);call mix_complex_array(h,s%pp_z)
    if(h==0)h=1
  end function
  subroutine mix_i(h,v)
    integer(int64),intent(inout)::h;integer(int64),intent(in)::v
    h=ieor(h,v);h=ieor(h,shiftl(h,13));h=ieor(h,shiftr(h,7));h=ieor(h,shiftl(h,17))
  end subroutine
  subroutine mix_int_array(h,a)
    integer(int64),intent(inout)::h;integer,intent(in)::a(:);integer::i
    call mix_i(h,int(size(a),int64));do i=1,size(a);call mix_i(h,int(a(i),int64));enddo
  end subroutine
  subroutine mix_real_array(h,a)
    integer(int64),intent(inout)::h;real(8),intent(in)::a(:);integer::i
    call mix_i(h,int(size(a),int64));do i=1,size(a);call mix_i(h,transfer(a(i),h));enddo
  end subroutine
  subroutine mix_complex_array(h,a)
    integer(int64),intent(inout)::h;complex(8),intent(in)::a(:);integer::i
    call mix_i(h,int(size(a),int64))
    do i=1,size(a)
      call mix_i(h,transfer(real(a(i),8),h))
      call mix_i(h,transfer(aimag(a(i)),h))
    enddo
  end subroutine
  logical function strict_ids(a)result(ok)
    integer,intent(in)::a(:);integer::i
    ok=size(a)>0.and.all(a>0);do i=2,size(a);if(a(i)<=a(i-1))ok=.false.;enddo
  end function
  logical function increasing_nonnegative(a)result(ok)
    integer,intent(in)::a(:);integer::i
    ok=all(a>=0);do i=2,size(a);if(a(i)<=a(i-1))ok=.false.;enddo
  end function
  logical function finite_real(a)result(ok)
    real(8),intent(in)::a(:);ok=all(ieee_is_finite(a))
  end function
  logical function finite_complex1(a)result(ok)
    complex(8),intent(in)::a(:);ok=all(ieee_is_finite(real(a,8))).and.all(ieee_is_finite(aimag(a)))
  end function
  logical function finite_complex2(a)result(ok)
    complex(8),intent(in)::a(:,:);ok=all(ieee_is_finite(real(a,8))).and.all(ieee_is_finite(aimag(a)))
  end function
end module
